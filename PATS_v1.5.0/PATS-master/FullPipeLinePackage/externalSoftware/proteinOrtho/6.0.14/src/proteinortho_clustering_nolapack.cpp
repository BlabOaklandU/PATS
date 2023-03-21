/*
 *	Clustering algorithm for Proteinortho
 *	Reads edge list and splits connected components
 *	according to algebraic connectivity threshold
 *
 *	Last updated: 2018/03/01		
 *	Author: Marcus Lechner
 */

#ifndef _PROTEINORTHOCLUSTERING
#define _PROTEINORTHOCLUSTERING

//#define DEBUG
//#define timeAnalysis //if set : need -std=c++11 for compiling

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <list> //BFS/DFS
#include <map>
#include <algorithm>
#include <cmath>
#include <vector>
#include <stack>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <memory>
#include <climits> // unsigned int max range

#ifdef _OPENMP
	#include <omp.h>
#endif

using floattype = float;
#define floatprecision_H 1
//the floatprecision_H has to be 1 for float and 2 for double

using namespace std;

#ifdef timeAnalysis
	#include <chrono>
	map<string,floattype> t_master;
#endif



struct wedge {unsigned int edge; unsigned short weight;};
struct protein {vector<wedge> edges; unsigned int species_id; string full_name;};

// Functions
floattype string2floattype(string);
void tokenize(const string& , vector<string>& , const string&);
void parse_file(string);
void remove_edge_index(const unsigned int, const unsigned int);
floattype getConnectivity(vector<unsigned int>&,bool);
void clear_edges(vector<unsigned int>&);
void partition_graph(void);
void print_group(vector<unsigned int>& , floattype );
floattype calc_group(vector<unsigned int>&);
void print_header(void);
void sort_species(void);
void stats(floattype, floattype);
void splitGroups(vector<floattype>&, vector<unsigned int>&, bool);
string getTime(void);

// Parameters
bool param_verbose 		= false;
floattype param_con_threshold 	= 0.1;		// as a reference: a chain a-b-c-d has 0.25
unsigned int debug_level	= 0;
floattype param_sep_purity 	= 1e-7;		// as a reference: a-b-c will give +/-0.707107 and 2.34857e-08 
unsigned int param_max_nodes	= 16777216; // 2^24
floattype param_min_species	= 1;
string param_rmgraph            = "remove.graph";
bool param_useWeights = false;
unsigned int param_minOpenmp = 256;
bool param_useKmereHeuristic = true;
// unsigned int param_maxRam_inKB = 16777216; // 16 GB of memory
bool param_useLapack = false;

// min/max number of alg con iterations
const unsigned int min_iter = 16;			// below this value, results may vary
unsigned int param_max_iter = 8192;			// below this value, results may vary
floattype param_epsilon = 1e-8; // analog to http://people.sc.fsu.edu/~jburkardt/c_src/power_method/power_method_prb.c
const unsigned int kmereHeuristic_minNodes = 1048576; //2^20
const unsigned int kmereHeuristic_protPerSpecies = 3;
const unsigned int kmereHeuristic_minNumberOfGroups = 3;
const unsigned int maxUseWeightsNumNodes = 1048576; //2^20
unsigned int param_lapack_power_threshold_n = 0; //2^8

// Globals
unsigned int species_counter = 0;	// Species
unsigned int protein_counter = 0;	// Proteins
vector<string> species;			// Number -> Name
vector<protein> graph;			// Graph containing all protein data
floattype last_stat = 0;			// For progress stats
unsigned int edges = 0;			// number of edges
vector<shared_ptr<ofstream> > graph_clean;			// File to store graph data
vector<int> reorder_table;		// Tells how proteins/species must be sorted
// unsigned int graph_ram_total_inKB=0;
unsigned int num_cpus=1;

// TMP Globals
map<string,int> species2id;		// Name -> Number
map<string,int> protein2id;		// Name -> Number

//TEST functions
bool test__max_of_diag();
bool test__generate_random_vector();
bool test__get_new_x();
bool test__makeOrthogonal();
bool test__normalize();
bool test__getY();
bool test__splitGroups();

#ifdef DEBUG
	unsigned int total_number_of_iterations_convergence = 0;
	unsigned int total_number_of_kmere_calls = 0;
	void debug__graph_integrity(vector<unsigned int>&);
	void debug__print_edgelist (protein&, const unsigned int, const int);
	void debug__conn_integrity(vector<unsigned int>&, floattype);
	void debug__getConnectivity();
	void debug__print_matrix( int m, int n, floattype* a, int lda );
#endif

///////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////
void printHelp() {
	cerr << "Proteinortho5 - Spectral partitioning algorithm v5.17 (Version without LAPACK!)" << endl;
	cerr << "-----------------------------------------------------" << endl;
	cerr << "This tool is part of Proteinortho" << endl;
	cerr << "" << endl;
	cerr << "Usage:   proteinortho5_clustering [OPTIONS] graph_files..." << endl;
	cerr << "Options: -verbose          report progress" << endl;
	cerr << "         -conn floattype      threshold for connectivity ["<<param_con_threshold<<"]" << endl;
	cerr << "         -purity floattype    threshold for purity ["<<param_sep_purity<<"]" << endl;
	cerr << "         -maxnodes int     max. number of nodes for alg. clustering ["<<param_max_nodes<<"]" << endl;
	cerr << "         -minspecies floattype threshold for species number ["<<param_min_species<<"]" << endl;
	cerr << "         -rmgraph STRING   output file for graph" << endl;
	cerr << "         -seed int         seed value for srand [current unix time]" << endl;
	cerr << "         -epsilon floattype   convergence threshold ["<<param_epsilon<<"]" << endl;
	// cerr << "         -ram int    		maximal used ram threshold for in MB [16384]" << endl;
	cerr << "         -weighted bool    the spectral partition is calculated using the bitscores ["<<param_useWeights<<"]" << endl;
	cerr << "         -cpus int       	the number of threads used for openMP ["<<num_cpus<<"]" << endl;
	cerr << "         -minOpenmp int    the minimum number of nodes for parallel power iteration ["<<param_minOpenmp<<"]" << endl;
	cerr << "         -kmere bool	    use the kmere-split heuristic ["<<param_useKmereHeuristic<<"]" << endl;
	cerr << "         -test 	    	various test-functions are called first [not set]" << endl;
	cerr << "         -maxRunsConvergence int    the maximum number of runs for the calculation of the algebraic connectivity ["<<param_max_iter<<"]" << endl;
}

int main(int argc, char *argv[]) {

	if (argc <= 1) {
		printHelp();
		return EXIT_FAILURE;
	}

	try {
		#ifdef _OPENMP
			omp_set_dynamic(0);     // Explicitly disable dynamic teams
			omp_set_num_threads(num_cpus); 
		#endif

		int rand_seed = 12345; //init randseed 

		// Read parameters
		int paras;
		vector<string> files;
		for (paras = 1; paras < argc; paras++) {
			string parameter = string(argv[paras]);
			if (parameter.substr(0, 1) != "-") {
				files.push_back(parameter);
			}
			else if (parameter == "-verbose") {
				paras++;
				if (string2floattype(string(argv[paras]))) {
					param_verbose = true;
				}
			}
			else if (parameter == "-conn") {
				paras++;
				param_con_threshold = string2floattype(string(argv[paras]));
				if(param_con_threshold<0 || param_con_threshold>1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << endl;throw;}
			}
			else if (parameter == "-purity") {
				paras++;
				param_sep_purity = string2floattype(string(argv[paras]));
				if(param_sep_purity<0){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << endl;throw;}
			}
			else if (parameter == "-ram") {
				paras++;
				// do nothing here is no LAPACK
			}
			else if(parameter == "-kmere"){
				paras++;
				param_useKmereHeuristic = int(string2floattype(string(argv[paras])));
			}
			else if (parameter == "-maxnodes") {
				paras++;
				param_max_nodes = string2floattype(string(argv[paras]));
			}
			else if (parameter == "-minspecies") {
				paras++;
				param_min_species = string2floattype(string(argv[paras]));
				if (param_min_species < 1) {cerr << string("-minspecies must at least be 1. Less than one gene per species is not possible as we only count those that have an entry.").c_str() << endl;throw;}
			}
			else if (parameter == "-debug") {
				paras++;
				debug_level = int(string2floattype(string(argv[paras])));
			}
			else if (parameter == "-epsilon") {
				paras++;
				param_epsilon = string2floattype(string(argv[paras]));
				if(param_epsilon<0||param_epsilon>1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << endl;throw;}
			}
			else if (parameter == "-minOpenmp") {
				paras++;
				param_minOpenmp = int(string2floattype(string(argv[paras])));
			}
			else if (parameter == "-weighted") {
				paras++;
				param_useWeights = int(string2floattype(string(argv[paras])));
			}
			else if (parameter == "-seed") {
				paras++;
				rand_seed = int(string2floattype(string(argv[paras])));
				if(rand_seed<0){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << endl;throw;}
			}
			else if (parameter == "-rmgraph") {
				paras++;
				param_rmgraph = string(argv[paras]);
			}else if(parameter == "-cpus"){
				paras++;
				#ifdef _OPENMP
					if(int(string2floattype(string(argv[paras])))<0){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << endl;throw;}
					omp_set_dynamic(0);     // Explicitly disable dynamic teams
					num_cpus=int(string2floattype(string(argv[paras])));
					omp_set_num_threads(num_cpus); 
				#endif
			}else if(parameter == "-maxRunsConvergence"){
				paras++;
				param_max_iter = int(string2floattype(string(argv[paras])));
			}else if(parameter == "-test"){
				bool test__max_of_diag_result = test__max_of_diag();
				bool test__generate_random_vector_result = test__generate_random_vector();
				bool test__get_new_x_result = test__get_new_x();
				bool test__makeOrthogonal_result = test__makeOrthogonal();
				bool test__normalize_result = test__normalize();
				bool test__getY_result = test__getY();
				bool test__splitGroups_result = test__splitGroups();
				cerr << "- test max_of_diag() : " << test__max_of_diag_result << endl;
				cerr << "- test generate_random_vector() : "<< test__generate_random_vector_result << endl;
				cerr << "- test get_new_x() : " << test__get_new_x_result << endl;
				cerr << "- test makeOrthogonal() : " << test__makeOrthogonal_result << endl;
				cerr << "- test normalize() : " << test__normalize_result << endl;
				cerr << "- test getY() : " << test__getY_result << endl;
				cerr << "- test splitGroups() : " << test__splitGroups_result << endl;
				if( !test__max_of_diag_result || !test__generate_random_vector_result || !test__get_new_x_result || !test__makeOrthogonal_result || !test__normalize_result || !test__getY_result || !test__splitGroups_result){
					cerr << string("Error: tests failed !").c_str() << endl;throw;
				}else{
					cerr << "All test passed." << endl;
					return EXIT_SUCCESS;
				} 
			}
			else {
				printHelp();
				cerr << endl << "Sorry, unknown option '" << string(argv[paras]) << "'!" << endl;
				return EXIT_FAILURE;
			}
		}

		srand(rand_seed);

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Debug level " << debug_level << endl;

		#ifdef DEBUG
			if (debug_level == 43){ return 1;} // ~ 4000 KB Maximum resident set size of the process during its lifetime, in Kbytes.  (memory usage until here)
		#endif

		// Parse files
		for (vector<string>::iterator it=files.begin() ; it != files.end(); it++) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Parsing file " << *it << endl;
			parse_file(*it);
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   I know " << species_counter <<  " species with " << protein_counter << " proteins and " << edges << " edges in sum" << endl;
		}

		// #ifdef DEBUG
		// 	if (debug_level == 42){ cerr << "graph_ram_total_inKB " << graph_ram_total_inKB << endl; return 1;} // 2609492 KB Maximum resident set size and the calculated graph_ram_total_inKB = 2910761.
		// #endif

		// if(graph_ram_total_inKB > param_maxRam_inKB){
		// 	cerr << endl << "ERROR: Memory overflow: allocated ram " << graph_ram_total_inKB/1e+3 << " MB did exceed the maximum ram threshold of "<< param_maxRam_inKB/1e+3 << " MB!"<< endl;
		// 	return EXIT_FAILURE;
		// }

		// Free memory
		files.clear();
		vector<string>().swap(files);
		species2id.clear();
		map<string,int>().swap(species2id);
		protein2id.clear();
		map<string,int>().swap(protein2id);

		// Stats
		if (param_verbose) cerr << species_counter << " species" << endl << protein_counter << " paired proteins" << endl << edges << " bidirectional edges" << endl;

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Maximumum number of nodes for connectivity calculations is " << param_max_nodes << endl;

		// Prepare sort of output
		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Sorting known species" << endl;
		sort_species();

		// Write output header
		print_header();							

		// Open graph-removal file
		string allRMgraphNames="";
		for(unsigned int i = 0 ; i < num_cpus ; i ++){
			stringstream ss;
			ss << i;
			allRMgraphNames+=(param_rmgraph+ss.str())+" ";
			graph_clean.push_back(make_shared<ofstream>((param_rmgraph+ss.str()).c_str()));
		}

		// // Clustering
		// if (debug_level > 0) cerr << getTime() << " [DEBUG]   Clustering" << endl;
		partition_graph();

		for(unsigned int i = 0 ; i < num_cpus ; i ++)
			(*graph_clean[i]).close();

		if(system(("cat "+allRMgraphNames+" >"+param_rmgraph).c_str())!=0 || system(("rm "+allRMgraphNames).c_str())!=0){
			cerr << "[ERROR]   cannot concatenate remove graphs" << endl;
			return EXIT_FAILURE;
		}

		#ifdef timeAnalysis
			for(map<string,floattype>::iterator it = t_master.begin() ; it != t_master.end() ; ++it) cerr << (*it).first << " " << (*it).second << endl;
		#endif
		#ifdef DEBUG
			cout << "conv:" << total_number_of_iterations_convergence << ", kmere_calls:" << total_number_of_kmere_calls<< endl;
		#endif

	}
	catch(string& error) {
		cerr << "[ERROR]   " << error << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}


unsigned int numberOfNodesToMemoryUsageLaplacian_inKB(unsigned int n){
	return (unsigned int)n*(unsigned int)n*sizeof(floattype)/1e+3;
}

typedef vector<unsigned int> ConnectedComponent;

ConnectedComponent DFS(vector<bool> * done, unsigned int cur_node ){

	ConnectedComponent ret; //return vector

	ret.push_back(cur_node);

	for (unsigned int i = 0; i < graph[cur_node].edges.size(); i++) {

		unsigned int adjacency_node = graph[cur_node].edges[i].edge;

		if(adjacency_node > graph.size()){
			cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is reporting an edge/adjacent node, that is not present in the graph.").c_str() << endl;throw;
		}
		if(adjacency_node > (*done).size()){
			cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is not present in done vector.").c_str() << endl;throw;
		}

		if( !(*done)[adjacency_node] ){

			(*done)[adjacency_node] = true;

			ConnectedComponent dfs_down_adjacency_node = DFS(done,adjacency_node);
			ret.insert(ret.end(),dfs_down_adjacency_node.begin(), dfs_down_adjacency_node.end());

		}
	}

	return ret;
}

ConnectedComponent DFS2(vector<bool> * done, unsigned int cur_node ){

	ConnectedComponent ret; //return vector
	vector<unsigned int> q;
	q.push_back(cur_node);

	while(q.size()>0){

		cur_node = q.back();
		q.pop_back();

		ret.push_back(cur_node);

		for (unsigned int i = 0; i < graph[cur_node].edges.size(); i++) {

			unsigned int adjacency_node = graph[cur_node].edges[i].edge;

			if(adjacency_node > graph.size()){
				cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is reporting an edge/adjacent node, that is not present in the graph.").c_str() << endl;throw;
			}
			if(adjacency_node > (*done).size()){
				cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is not present in done vector.").c_str() << endl;throw;
			}

			if( !(*done)[adjacency_node] ){

				(*done)[adjacency_node] = true;
				q.push_back(adjacency_node);
			}
		}
	}
	return ret;
}


ConnectedComponent BFS(vector<bool> * done, unsigned int cur_node ){

	ConnectedComponent ret; //return vector
	list<unsigned int> q;
	q.push_back(cur_node);
	(*done)[cur_node]=true;

	while(q.size()>0){

		list<unsigned int> q_new;

		for(list<unsigned int>::iterator it = q.begin() ; it != q.end() ; ++it){
			cur_node = *it;
			ret.push_back(cur_node);

			for (unsigned int i = 0; i < graph[cur_node].edges.size(); i++) {

				unsigned int adjacency_node = graph[cur_node].edges[i].edge;

				if(adjacency_node > graph.size()){
					cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is reporting an edge/adjacent node, that is not present in the graph.").c_str() << endl;throw;
				}
				if(adjacency_node > (*done).size()){
					cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is not present in done vector.").c_str() << endl;throw;
				}

				if( !(*done)[adjacency_node] ){

					(*done)[adjacency_node] = true;
					q_new.push_back(adjacency_node);
				}
			}
		}

		q=q_new;
	}
	return ret;
}

bool BFS_critical(vector<bool> * done, unsigned int start_node, unsigned int end_node ){

	ConnectedComponent ret; //return vector
	list<unsigned int> q;
	q.push_back(start_node);
	(*done)[start_node]=true;

	while(q.size()>0){

		list<unsigned int> q_new;

		for(list<unsigned int>::iterator it = q.begin() ; it != q.end() ; ++it){
			unsigned int cur_node = *it;
			ret.push_back(cur_node);

			for (unsigned int i = 0; i < graph[cur_node].edges.size(); i++) {

				unsigned int adjacency_node = graph[cur_node].edges[i].edge;

				if(adjacency_node > graph.size()){
					cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is reporting an edge/adjacent node, that is not present in the graph.").c_str() << endl;throw;
				}
				if(adjacency_node > (*done).size()){
					cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is not present in done vector.").c_str() << endl;throw;
				}

				if(cur_node != start_node && adjacency_node == end_node )
					return false;

				if(cur_node==start_node && adjacency_node == end_node)
					continue;

				if( !(*done)[adjacency_node] ){

					(*done)[adjacency_node] = true;
					q_new.push_back(adjacency_node);
				}
			}
		}

		q=q_new;
	}
	return true;
}
bool is_critical(unsigned int node_a, unsigned int node_b /*, unsigned int node_not_a = 0, unsigned int node_not_b = 0*/){

	// a edge (given as a pair of nodes) is critical if:
	// the only way from a to b that is (a,b)
	
	vector<bool> done = vector<bool> (graph.size(), false);	// Keep track on what was done

	list<unsigned int> q;
	q.push_back(node_a);
	done[node_a]=true;

	if( node_a==node_b){
		return false;
	}

	while(q.size()>0){

		unsigned int cur_node = q.back();
		q.pop_back();

		for (unsigned int i = 0; i < graph[cur_node].edges.size(); i++) {

			unsigned int adjacency_node = graph[cur_node].edges[i].edge;
			
			if( adjacency_node==node_b && cur_node==node_a)
				continue;

			if( !done[adjacency_node] ){

				if(adjacency_node==node_b && cur_node != node_a){
					return false;
				}

				done[adjacency_node] = true;
				q.push_back(adjacency_node);
			}
		}
	}
	return true;
}

bool criticalHeuristic(ConnectedComponent *cur_cc){

	bool foundSomeCriticalEdge = false;
	map<pair<unsigned int,unsigned int>, bool> done;
	unsigned int ce_counter=0;

	for (unsigned int i = 0 ; i < (*cur_cc).size() ; i++) {
		unsigned int node_a = (*cur_cc)[i];

		if(graph[node_a].edges.size() < 3)continue;
		for (unsigned int j = 0 ; j < graph[node_a].edges.size() ; j++) {
			unsigned int node_b = graph[node_a].edges[j].edge;

			if(graph[node_b].edges.size() < 3)continue;
			pair<unsigned int,unsigned int> key;
			if(node_a < node_b) key = make_pair(node_a,node_b);
			else key = make_pair(node_b,node_a);

			if(!done.count(key)){
				done[key] = true;

				vector<bool> done_BFS = vector<bool> (graph.size(), false);	// Keep track on what was done (for each node)
				done_BFS[node_a]=true; // mark this node
				
				if(BFS_critical(&done_BFS,node_a,node_b)){

					vector<wedge>::iterator it = graph[node_a].edges.begin();
					for(it = graph[node_a].edges.begin() ; it != graph[node_a].edges.end() ; ++it){
						if((*it).edge == node_b){
							graph[node_a].edges.erase(it);
							break;
						}
					}
					it = graph[node_b].edges.begin();
					for(it = graph[node_b].edges.begin() ; it != graph[node_b].edges.end() ; ++it){
						if((*it).edge == node_a){
							graph[node_b].edges.erase(it);
							break;
						}
					}

					ce_counter++;
					foundSomeCriticalEdge = true;
				}
			}
		}
	}
	cerr << " [WARNING] (critical edge heuristic)  Found and removed "<<ce_counter<<" critical edge(s)." << endl;

	return foundSomeCriticalEdge;
}

struct compare_ConnectedComponents { //sort from large to small
	bool operator() (const ConnectedComponent &a, const ConnectedComponent &b) const {
	    return a.size() > b.size();
	}
};

///////////////////////////////////////////////////////////
// Major partioning algorithm
///////////////////////////////////////////////////////////
void partition_graph() {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	restart_partition_graph:

	vector<bool> done;	// Keep track on what was done
	done = vector<bool> (graph.size(), false);	// Keep track on what was done (for each node)
	bool allNodesAreDone = false;

	while( !allNodesAreDone ){

		vector<ConnectedComponent> CC; // vector of all connected components found

		while( true ){ // CC.size() < num_cpus / gather up to num_cpus connected components

			allNodesAreDone = true;

			for (unsigned int protein_id = 0 ; protein_id < graph.size() ; protein_id++) {
				
				if (done[protein_id]){continue;}// We were here already

				done[protein_id]=true; // mark this node
				ConnectedComponent cur_cc = BFS(&done,protein_id); // get the CC of the current node (protein_id) 

				// Do not report singles
				if (cur_cc.size() < 2) {continue;} // singletons are from no interest

				if (debug_level > 0) cerr << getTime() << " [DEBUG] Found connected component: " << cur_cc.size() << " proteins (ID: " << protein_id << ")" << endl;

				// Skip those that are too large (try critical edge heuristic)
				if (cur_cc.size() > param_max_nodes) {

					cerr << " [WARNING]  Found a very large connected component that contains " << cur_cc.size() << ">" << param_max_nodes << " (maxnodes) elements. This behavior can be adjusted using -maxnodes, still you better consider a higher e-value or coverage thresholds. Now using the Critical-Heuristic: try to identify and remove critical edges." << endl;

					bool foundSomeCriticalEdge = criticalHeuristic(&cur_cc);

					if(!foundSomeCriticalEdge){
						cerr << " [WARNING]  The very large connected component that contains " << cur_cc.size() << " elements, did not have any critical edges." << endl;
						cerr << " [WARNING]  Skipping this connected component. This behavior can be adjusted using -maxnodes, still you better consider a higher e-value or coverage thresholds." << endl;
						clear_edges(cur_cc);
						continue;	
					}else{
						cerr << " [WARNING]  Now restarting the partition algorithm." << endl;
						goto restart_partition_graph;
					}
				}

				CC.push_back(cur_cc);	
				allNodesAreDone=false;	
				break;
			}

			if(allNodesAreDone)break; // no additional CC can be found -> done
		}

		// Connectivity analysis
		sort(CC.begin(), CC.end(), compare_ConnectedComponents());

		vector<floattype> connectivities(CC.size(),0.0); //resulting connectivities

		if (debug_level > 0) cerr << getTime() << " [DEBUG] Found: " << CC.size() << " connected components in total, continue now with the calculation of the algebraic connectivities." << endl;

		vector<bool> done_CC = vector<bool> (CC.size(), false);	// Keep track on what was done

		unsigned int power_calls=0;

		// now do the POWER ITERATION if needed (for large CC)
		// remember: getConnectivity also splits according to the fiedler vector (if the score is low enough)
		for(unsigned int i = 0 ; i < CC.size() ; i++){
			if (debug_level > 0 ) cerr << getTime() << " [DEBUG]   Power iteration is used." << endl;
			connectivities[i] = getConnectivity(CC[i],0);
			done_CC[i]=true;
			if (debug_level > 0)power_calls++;
		}

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   The algebraic connectivity calculation are now done, continue now with the evaluation of the scores." << endl;

		// evaluate the connectivities -> reiterating
		for(unsigned int i = 0 ; i < CC.size() ; i++){

			if (debug_level > 0) cerr << getTime() << " [DEBUG]     Connectivity was " << connectivities[i] << endl;
			// negative values indicate an overwrite. con will be accepted due to good average species coverage
			if (param_con_threshold && connectivities[i] >= 0 && connectivities[i] < param_con_threshold) {
				if (debug_level > 0) cerr << getTime() << " [DEBUG]     Reiterating" << endl;
				// Split groups is done in getConnectivity function
				// Remove flags
				for (unsigned int j = 0; j < CC[i].size(); j++) done[CC[i][j]] = false;
				allNodesAreDone=false;
				continue;
			}
			if (debug_level > 0) cerr << getTime() << " [DEBUG]     Print group:" << endl;

			// Output finial group
			if 	(connectivities[i] < 0) {print_group(CC[i],-connectivities[i]);}				// Connectivity threshold not used
			else if (connectivities[i] >= param_con_threshold) {print_group(CC[i],connectivities[i]);}

			// Print stats
			if (param_verbose) stats(CC[i][0],protein_counter);
			
			// Clean up and go on with the next connected component
			clear_edges(CC[i]);
		}
	}
	if (param_verbose) stats(1,1);
	if (param_verbose) cerr << "\r" << "Done                       " << endl;

}


///////////////////////////////////////////////////////////
// Basic Graph functions
///////////////////////////////////////////////////////////
// Remove all edges from the given list of protein ids 
void clear_edges(vector<unsigned int>& nodes) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	for (unsigned int i = 0; i < nodes.size(); i++) {
		graph[nodes[i]].edges.clear();
		vector<wedge>().swap(graph[nodes[i]].edges);
	}

	#ifdef timeAnalysis
		if(!t_master.count("clear_edges")) t_master["clear_edges"]=0;
		t_master["clear_edges"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

// Remove edges that go beyond the given group groups a, works in a single direction, no overlap in parallel runs if lists are without overlap [TRY]
void removeExternalEdges(map<unsigned int,bool>& a) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif
	// For each node id in Hash
	for (map<unsigned int,bool>::iterator it = a.begin(); it != a.end(); it++) {
		unsigned int node = it->first;
		// For each edge of this node
		for (unsigned int j = 0; j < graph[node].edges.size();) {
			// If target node is not our group, remove
			if (a.find(graph[node].edges[j].edge) == a.end()) {
				remove_edge_index(node,j);				// Works in a single direction, has to be called for outsider group as well
			}
			else {
				j++; // only increment if all is fine
			}
		}
	}
	#ifdef timeAnalysis
		if(!t_master.count("removeExternalEdges")) t_master["removeExternalEdges"]=0;
		t_master["removeExternalEdges"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

// Remove edge in a single direction given an index -- openMP safe
void remove_edge_index(const unsigned int node_id, const unsigned int index) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	protein node = graph[node_id];
	if (index >= graph[node_id].edges.size()) {
		stringstream ss;
		cerr << "Called index out of bounds: " << node.full_name << " sized " << graph[node_id].edges.size() << " called " << index << endl;
		throw;
	}

	unsigned int tid=0;
	#ifdef _OPENMP
		tid=omp_get_thread_num();
	#endif

	if(graph_clean.size()>tid)
		(*graph_clean[tid]) << node.full_name << "\t" << species[node.species_id] << "\t" << graph[node.edges[index].edge].full_name << "\t" << species[graph[node.edges[index].edge].species_id] << endl;
	
	// ofstream os((ss.str()).c_str(),std::ios_base::app);
	// os << node.full_name << "\t" << species[node.species_id] << "\t" << graph[node.edges[index].edge].full_name << "\t" << species[graph[node.edges[index].edge].species_id] << endl;
	// os.close();

	// overwrite this position with the last
	graph[node_id].edges[index] = graph[node_id].edges.back();
	// and remove last (also works if index == last)
	graph[node_id].edges.pop_back();
	
	#ifdef timeAnalysis
		if(!t_master.count("remove_edge_index")) t_master["remove_edge_index"]=0;
		t_master["remove_edge_index"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

///////////////////////////////////////////////////////////
// File parser
///////////////////////////////////////////////////////////
void parse_file(string file) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	// unsigned long graph_ram_total=0;

	if (param_verbose) cerr << "Reading " << file << endl;
	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {
		// For each line
		string file_a = "";	unsigned int file_a_id = 0;
		string file_b = "";	unsigned int file_b_id = 0;

		floattype avg_bitscore_median=255; //for normalization, the bitscore is devided by the median bitscore -> small score (bad for unsigned short cast) -> multiply by factor 255. -> ini 255 -> 1
		//255 is exactly in the middle of ushort range, such that the score can be 255 fold upregulated or down regulated.

		while (!graph_file.eof()) {
			getline(graph_file, line);
			vector<string> fields;
			tokenize(line, fields, "\t");
			// Header line
			if ( (fields.size() == 2 || fields.size() == 6) && fields[0].substr(0, 1) == "#") {
				file_a = fields[0].substr(2, fields[0].size()-2);
				file_b = fields[1];

				if (file_a == "file_a" && file_b == "file_b") continue;	// Header line
				if (file_a == "a" && file_b == "b") continue;	// Header line

				if(fields.size() == 6){
					avg_bitscore_median=(string2floattype(fields[3])+string2floattype(fields[5]))/2;
					if(avg_bitscore_median<100){avg_bitscore_median=100;}
				}

				// Map species a
				if (species2id.find(file_a) == species2id.end())	{
						species.push_back(file_a);
						species2id[file_a] = species_counter++;
				}
				// Map species b
				if (species2id.find(file_b) == species2id.end())	{
						species.push_back(file_b);
						species2id[file_b] = species_counter++;
				}

				file_a_id = species2id[file_a];
				file_b_id = species2id[file_b];
			}
			// Data line
			else if ((fields.size() == 6 || fields.size() == 8) && fields[0].substr(0, 1) != "#") {
				// a b e1 b1 e2 b2 score

				// 5.16 deal with duplicated IDs by adding file ID to protein ID
				string ida = fields[0];
				string idb = fields[1];
				fields[0] += " "; fields[0] += to_string(file_a_id);
				fields[1] += " "; fields[1] += to_string(file_b_id);

				// 5.16 do not point to yourself
				if (!fields[0].compare(fields[1])) {continue;}

				// A new protein
				if (protein2id.find(fields[0]) == protein2id.end())	{
					protein a;
					a.full_name	= ida;
					// graph_ram_total+=a.full_name.size()*sizeof(string);
					a.species_id	= file_a_id;
					// graph_ram_total+=sizeof(unsigned int);
					protein2id[fields[0]] = protein_counter++;
					// graph_ram_total+=fields[0].size()*sizeof(string)+sizeof(unsigned int);

					if( graph.size() >= 1073741824-1 ){
						cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << endl;throw;
					}
					graph.push_back(a);
					// graph_ram_total+=sizeof(protein);
				}
				if (protein2id.find(fields[1]) == protein2id.end())	{
					protein b;
					b.full_name	= idb;
					// graph_ram_total+=b.full_name.size()*sizeof(string);
					b.species_id	= file_b_id;
					// graph_ram_total+=sizeof(unsigned int);
					protein2id[fields[1]] = protein_counter++;
					// graph_ram_total+=fields[1].size()*sizeof(string)+sizeof(unsigned int);

					if( graph.size() >= 1073741824-1 ){
						cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << endl;throw;
					}
					graph.push_back(b);
					// graph_ram_total+=sizeof(protein);
				}

				// Bitscores 
				// check range (in floattype)

//				cerr << string2floattype(fields[3]) << " " << avg_bitscore_median <<" "<< string2floattype(fields[3])/avg_bitscore_median << endl;

				floattype bit_a = 255.0* string2floattype(fields[3])/avg_bitscore_median;
				floattype bit_b = 255.0* string2floattype(fields[5])/avg_bitscore_median;
				if(bit_a<1){bit_a=1;}
				if(bit_b<1){bit_b=1;}

				if(bit_a>USHRT_MAX){
					cerr << " [WARNING] unsigned short overflow " << bit_a <<  ">USHRT_MAX (bitscore of "<< ida<< ") using "<< USHRT_MAX<< " (USHRT_MAX) instead." << endl;
					bit_a=(floattype)USHRT_MAX;
				}
				if(bit_b>USHRT_MAX){
					cerr << " [WARNING] unsigned short overflow " << bit_b <<  ">USHRT_MAX (bitscore of "<< idb<< ") using "<< USHRT_MAX<< " (USHRT_MAX) instead." << endl;
					bit_b=(floattype)USHRT_MAX;
				}

				// assign
				unsigned short bitscore_avg = (bit_a + bit_b)/2;
				if(bitscore_avg<1){bitscore_avg=1;}
		
				// Add link to graph (reciprocal)					
				unsigned int a_id = protein2id[fields[0]];
				unsigned int b_id = protein2id[fields[1]];

				// 5.17, add weight
				wedge w;
				w.edge=b_id;
				// graph_ram_total+=sizeof(unsigned int);
				w.weight=bitscore_avg;
				// graph_ram_total+=sizeof(unsigned int);
				graph[a_id].edges.push_back(w);
				// graph_ram_total+=sizeof(wedge);
				w.edge=a_id;
				// graph_ram_total+=sizeof(unsigned int);
				w.weight=bitscore_avg;
				// graph_ram_total+=sizeof(unsigned int);
				graph[b_id].edges.push_back(w);
				// graph_ram_total+=sizeof(wedge);
				edges++;
			}
		}
		graph_file.close();
	}
	else {
		cerr << string("Could not open file " + file).c_str() << endl;throw;
	}

	// graph_ram_total_inKB += graph_ram_total/1e+3;

	// if (debug_level > 0) cerr << getTime() << " [DEBUG]  Expected Memory Usage of the current input graph: " << graph_ram_total_inKB << " KB = "  << graph_ram_total_inKB/1e+3 << " MB. (current input graph are all the currently loaded files)" << endl;

	if(species_counter==0){species.push_back("0");species_counter++;}

	#ifdef timeAnalysis
		if(!t_master.count("parse_file")) t_master["parse_file"]=0;
		t_master["parse_file"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

///////////////////////////////////////////////////////////
// Output
///////////////////////////////////////////////////////////
// Sort
void sort_species(void) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	reorder_table.reserve(species_counter);
	vector<string> species_sorted (species_counter);
	copy(species.begin(), species.end(), species_sorted.begin());
	sort(species_sorted.begin(), species_sorted.end());
	// find new locations (not efficent but list is small)	
	for (unsigned int i = 0; i < species_counter; i++) {
		for (unsigned int j = 0; j < species_counter; j++) {
			if (species[i] == species_sorted[j]) {reorder_table[j] = i; continue;}
		}
	}

	#ifdef timeAnalysis
		if(!t_master.count("sort_species")) t_master["sort_species"]=0;
		t_master["sort_species"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}


// Progress stats
void stats(floattype i, floattype size) {
	if (!param_verbose) return;
	floattype stat = floattype(i/size*100);
	if (last_stat * 1.01 < stat) {
		last_stat = stat;
		cerr << "\r" << "                          ";
		cerr << "\r" << "Clustering: " << setprecision (2) << fixed << stat << "%";
	}
}

// Header with species names
void print_header() {
	cout << "# Species\tGenes\tAlg.-Conn.";
	for (unsigned int i = 0; i < species_counter; i++) {
		cout << "\t" << species[reorder_table[i]];
	}
	cout << endl;
}

struct protein_degree
{
    inline bool operator() (const pair<string,int> & p1, const pair<string,int>& p2)
    {
        return (p1.second > p2.second);
    }
};

// Group formatting
void print_group(vector<unsigned int>& nodes, floattype connectivity) {

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif
	vector<vector<pair<string,int> > > line(species_counter);	// Output vector

	unsigned int species_number = 0;
	// For each protein in group 	 	
	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int current_protein = nodes[i];
		unsigned int current_species = graph[current_protein].species_id;
		if (line[current_species].size() == 0) 
			species_number++;
		line[current_species].push_back(make_pair( graph[current_protein].full_name , graph[current_protein].edges.size() ));
	}

	cout << species_number << "\t" << nodes.size() << "\t" << setprecision (3) << connectivity;

	// List group data
	for (unsigned int i = 0; i < species_counter; i++) {
		
		string return_line="";

		// sort line 
		if (line[reorder_table[i]].size() > 0) {

			sort( line[reorder_table[i]].begin(), line[reorder_table[i]].end(), protein_degree() );

			return_line = line[reorder_table[i]][0].first;

			for (unsigned int k = 1; k < line[reorder_table[i]].size(); k++) {
				return_line.append(","+line[reorder_table[i]][k].first);
			}
		}
		if(return_line == "")
			return_line = "*";

		// output
		cout << "\t" << return_line;
	}

	cout << endl;

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::print_group")) t_master["partition_graph::print_group"]=0;
		t_master["partition_graph::print_group"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

// Calculate number of species formatting
floattype calc_group(vector<unsigned int>& nodes) {

	map<unsigned int, bool> speciesids;

	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int current_protein = nodes[i];
		speciesids[graph[current_protein].species_id]=1;
	}

	return speciesids.size()==0 ? 99999 : ((floattype)nodes.size())/((floattype)speciesids.size()); // 99999 such that if the species information is missing, then the criterion always fails and the splits are only made based on the alg. connectivity
/*
	vector<unsigned int> line(species_counter,0);	// Species vector
	unsigned int species_number = 0;
	// For each protein in group, count species
	for (unsigned int i = 0; i < nodes.size(); i++) {
		unsigned int current_protein = nodes[i];
		unsigned int current_species = graph[current_protein].species_id;
		if (line[current_species] == 0) {species_number++;}			// new species
		line[current_species]++;						// calc
	}
	line.clear();
	vector<unsigned int>().swap(line);

	unsigned int sum = 0;
	for (unsigned int current_species = 0; current_species < species_counter; current_species++) {
		sum += line[current_species];
	}
	
	floattype avg = floattype(sum)/floattype(species_number);
	return avg;*/
}


///////////////////////////////////////////////////////////
// Misc functions
///////////////////////////////////////////////////////////
// Convert string to floattype
floattype string2floattype(string str) {
	istringstream buffer(str);
	floattype value;
	buffer >> value;
	return value;
}


// Split a string at a certain delim
void tokenize(const string& str, vector<string>& tokens, const string& delimiters = "\t") {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	// Skip delimiters at beginning.
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos = str.find_first_of(delimiters, lastPos);

	while (string::npos != pos || string::npos != lastPos) {
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}

	#ifdef timeAnalysis
		if(!t_master.count("tokenize")) t_master["tokenize"]=0;
		t_master["tokenize"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

///////////////////////////////////////////////////////////
// Algebraic connectivity functions
///////////////////////////////////////////////////////////
// Return maximum degree of given protein_ids -- openMP A ML
unsigned int max_of_diag(vector<unsigned int>& nodes, vector<unsigned int>& diag) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif
	
	unsigned int max = 0;

	bool useOpenMpFlag = (nodes.size() > param_minOpenmp);

	#pragma omp parallel for reduction(max: max) if (useOpenMpFlag)
	for (unsigned int i = 0; i < nodes.size(); i++) {
		if (diag[i] > max) max = diag[i] ;
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::max_of_diag")) t_master["partition_graph::convergence::max_of_diag"]=0;
		t_master["partition_graph::convergence::max_of_diag"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return max;
}

// Generate random vector x of size size
vector<floattype> generate_random_vector(const unsigned int size) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	vector<floattype> x(size);

	x[0] = (floattype)(rand() % 999+1)/1000;	// 0 to 1
	for (unsigned int i = 1; i < size; i++) {
		x[i] = (floattype)(rand() % 999+1)/1000;	// 0 to 1
		if (x[i] == x[i-1]) x[i] /= 3;		// Check: at least one value must be different from the others but still within 0 and 1
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::generate_random_vector")) t_master["partition_graph::convergence::generate_random_vector"]=0;
		t_master["partition_graph::convergence::generate_random_vector"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x;
}



// determine new X, Formula (1) -- openMP B ML
vector<floattype> get_new_x(vector<floattype> x, vector<unsigned int>& nodes, map<unsigned int,unsigned int> &mapping, bool isWeighted) {
	
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	vector<floattype> x_new(x.size(),0);
	
	bool useOpenMpFlag = (nodes.size() > param_minOpenmp);

	if(isWeighted){

		#pragma omp parallel for schedule(dynamic) if (useOpenMpFlag)
		for (unsigned int i = 0; i < nodes.size(); i++) {

			// go through adjacency list of node 
			for (unsigned int j = 0; j < graph[nodes[i]].edges.size(); j++) {
				// y points to z, so take entry z from x
				unsigned int abs_target = graph[nodes[i]].edges[j].edge;
				unsigned int rel_target = mapping[abs_target];

				x_new[i] += x[rel_target]*(floattype)graph[nodes[i]].edges[j].weight;
			}
		}

	}else{

		#pragma omp parallel for schedule(dynamic) if (useOpenMpFlag)
		for (unsigned int i = 0; i < nodes.size(); i++) {

			// go through adjacency list of node 
			for (unsigned int j = 0; j < graph[nodes[i]].edges.size(); j++) {
				// y points to z, so take entry z from x
				unsigned int abs_target = graph[nodes[i]].edges[j].edge;
				unsigned int rel_target = mapping[abs_target];

				x_new[i] += x[rel_target];
			}
		}
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::get_new_x")) t_master["partition_graph::convergence::get_new_x"]=0;
		t_master["partition_graph::convergence::get_new_x"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x_new;
}

// Make vector x orthogonal to 1, Formula (2) -- openMP A ML
vector<floattype> makeOrthogonal(vector<floattype> x) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	floattype sum = 0;

	bool useOpenMpFlag = (x.size() > param_minOpenmp);

	#pragma omp parallel for reduction(+: sum) if (useOpenMpFlag)
	for (unsigned int i = 0; i < x.size(); i++) {sum += x[i];}

	floattype average = sum/x.size();

	#pragma omp parallel for if (useOpenMpFlag)
	for (unsigned int i = 0; i < x.size(); i++) {x[i] -= average;}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::makeOrthogonal")) t_master["partition_graph::convergence::makeOrthogonal"]=0;
		t_master["partition_graph::convergence::makeOrthogonal"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x;
}

// Normalize vector x, Formula (4) -- openMP A ML
vector<floattype> normalize(vector<floattype> x, floattype *length) {

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	floattype sum = 0;

	bool useOpenMpFlag = (x.size() > param_minOpenmp);

	#pragma omp parallel for reduction(+: sum) if (useOpenMpFlag)
	for (unsigned int i = 0; i < x.size(); i++) {sum += x[i]*x[i];}

	*length = sqrt(sum);
	if (*length == 0) {*length = 0.000000001;}// ATTENTION not 0!

	#pragma omp parallel for if (useOpenMpFlag)
	for (unsigned int i = 0; i < x.size(); i++) {x[i] /= *length;}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::normalize")) t_master["partition_graph::convergence::normalize"]=0;
		t_master["partition_graph::convergence::normalize"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x;
}

// Qx, Formula (5) -- openMP A ML
vector<floattype> getY(floattype max_degree, vector<floattype> x_hat, vector<floattype> x_new, vector<unsigned int>& nodes, vector<unsigned int>& diag){

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	bool useOpenMpFlag = (nodes.size() > param_minOpenmp);

	#pragma omp parallel for if (useOpenMpFlag)
	for (unsigned int i = 0; i < x_hat.size(); i++) {

		x_hat[i] *= (2*max_degree - diag[i]);
		x_hat[i] += x_new[i];

	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::getY")) t_master["partition_graph::convergence::getY"]=0;
		t_master["partition_graph::convergence::getY"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x_hat;
}

unsigned int sumOutDegs(const vector<unsigned int>& nodes) {

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	unsigned int sum = 0;
	for (unsigned int a = 0; a < nodes.size(); a++) {
		unsigned int from = nodes[a];
		sum += graph[from].edges.size();
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::sumOutDegs")) t_master["partition_graph::convergence::sumOutDegs"]=0;
		t_master["partition_graph::convergence::sumOutDegs"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return sum/2;
}


floattype getConnectivity(vector<unsigned int>& nodes, bool useLapack) {

	bool useWeights = (param_useWeights && nodes.size() <= maxUseWeightsNumNodes); //param_useWeights = user input whether weights should be used or not. useWeights = the true value, that is true if param_useWeights is true and the maximum number of nodes are not exeeded for the weighted algorithm (maxUseWeightsNumNodes)

	if(!useWeights){ 
		bool isComplete = true;
		for(unsigned int i = 0 ; i < nodes.size() ; i++){
			if(graph[nodes[i]].edges.size() != nodes.size()-1){
				isComplete = false;
				break;
			}
		}
		if(isComplete){
			floattype conn = 1.0;
			if (conn < param_con_threshold) {
				if (param_min_species > 1) {
					floattype avg = calc_group(nodes);
					if (debug_level > 0) cerr << getTime() << " [DEBUG]   Found " << avg << " genes/species on average. User asked for " << param_min_species << endl;
					if (avg <= param_min_species) return -conn;
				}
				clear_edges(nodes);	// destroy group
			}
			return conn;
		}
	}

	unsigned int n=nodes.size();

	if( n>1073741824 ){
		cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << endl;throw;
	}

	floattype maxWeight=-1;
	map<unsigned int,unsigned int> mapping;
	for (unsigned int i = 0; i < (unsigned int)n; i++) {mapping[nodes[i]] = i;}

	floattype connectivity = -1;
	vector<floattype> x_hat(n);

	if(useLapack){
		cerr << " [WARNING] Using dsyevr (lapack) would overflow the given amount of memory. Continue now with the slow standard approach (power iteration)." << endl;
	}

	if(param_useWeights && !useWeights){
		cerr << " [WARNING] The maximum number of nodes for the weighted algorithm is exeeded. Continue now with the faster unweighted algorithm." << endl;
	}
	//diagonal matrix diag : d(u,u)=number of adjacency nodes=deg(u)
	vector<unsigned int> diag(n);

	if(useWeights){
		for (unsigned int i = 0; i < (unsigned int)n; i++) {
			diag[i]=0;
			for (unsigned int j = 0; j < graph[nodes[i]].edges.size(); j++) {
				diag[i] += graph[nodes[i]].edges[j].weight;
				if(useWeights && maxWeight<graph[nodes[i]].edges[j].weight)maxWeight=graph[nodes[i]].edges[j].weight;
			}
		}
	}else{
		for (unsigned int i = 0; i < (unsigned int)n; i++) {
			diag[i]=graph[nodes[i]].edges.size();
			if(useWeights)
				for (unsigned int j = 0; j < graph[nodes[i]].edges.size(); j++) {
					if(maxWeight<graph[nodes[i]].edges[j].weight)maxWeight=graph[nodes[i]].edges[j].weight;
				}
		}
	}

	// Get max degree / sum of weights of nodes
	unsigned int max_d = max_of_diag(nodes,diag);	

	// Init randomized variables. 
	vector<floattype> norm;

	vector<floattype> x = generate_random_vector(n);
	

	// Orthogonalize + normalize vector + get initial lenght
	floattype current_length = 0;
	floattype last_length;

	x_hat = makeOrthogonal(x);
	norm = normalize(x_hat, &last_length);


	// Repeat until difference < param_epsilon
	unsigned int iter = 0;	// catch huge clustering issues by keeping track here

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now();
	#endif	

	while(++iter < param_max_iter) { 

		if (debug_level > 2 && iter%100 == 0) cerr << getTime() << " [DEBUG L2] Step " << iter << " / " << param_max_iter << " error:" <<  abs(current_length-last_length) << endl;
		last_length = current_length;

		// Get a new x
		x = get_new_x(norm, nodes, mapping, useWeights);

		// Get y
		vector<floattype> y = getY(max_d,norm,x,nodes,diag);

		// Orthogonalize
		x_hat = makeOrthogonal(y);

		// Get lenght (lambda) & normalize vector
		norm = normalize(x_hat, &current_length);

		if ( abs(current_length-last_length) < param_epsilon && iter >= min_iter ) break;	// prevent convergence by chance, converge to param_epsilon
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence")) t_master["partition_graph::convergence"]=0;
		t_master["partition_graph::convergence"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	if (debug_level > 0) cerr << getTime() << " [DEBUG]   " << iter << " / " << param_max_iter << " iterations required (error is " << abs(current_length-last_length) << ")" << endl;

	#ifdef DEBUG
		total_number_of_iterations_convergence+=iter;
	#endif

	if(useWeights){
		connectivity = (-current_length+2*max_d)/(maxWeight*(floattype)n);
	}else{
		connectivity = (-current_length+2*max_d)/(n);
	}
	x_hat = normalize(x_hat, &current_length);
	
		
	// 5.17 catch hardly-converging groups
	if (iter >= param_max_iter) {
		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Connectivity score of connected component with " << n << " elements did not converge perfectly in time." << endl;
	}

	if (debug_level > 1){
		cerr << getTime() << " [DEBUG]   Connectivity score:" << connectivity;
		if ( (debug_level > 1 && (unsigned int)n<100 ) || debug_level > 2 ){
			cerr << " eigenvector: (";
			for(unsigned int i = 0 ; i < (unsigned int)n ; i++)
				cerr << x_hat[i] << ","; 
			cerr << ")";
		}
		cerr << endl;
	} 

	// Split groups if connectivity is too low, remove tree like structures that might have arosen
	if (connectivity < param_con_threshold) {
		
		// 5.17 new threshold option overwrites connectivity
		if (param_min_species > 1) {
			floattype avg = calc_group(nodes);
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Found " << avg << " genes/species on average. User asked for " << param_min_species << endl;
			if (avg <= param_min_species) {
				if (debug_level > 0) cerr << getTime() << " [DEBUG]   Group is going to be accepted despite connectivity " << connectivity << endl;
				// just to be safe from infinit loops due to rounding
				if (connectivity == 0) connectivity = 0.001;
				// no splitting despite bad connectivity
				return -connectivity;			
			}
		}

		splitGroups(x_hat, nodes , useLapack);
	}
	
	return connectivity;
}

bool comparator_pairfloattypeUInt ( const pair<floattype,unsigned int>& l, const pair<floattype,unsigned int>& r )
   { return l.first < r.first; }

// Split connected component according to eigenvector -- openMP B ML
void splitGroups(vector<floattype>& y, vector<unsigned int>& nodes , bool useLapack){

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	if(y.size() != nodes.size() || y.size()==0){
		cerr << string("[CRITICAL ERROR]   y.size() != nodes.size() || y.size()==0.").c_str() << endl;throw;
	}

	bool fallback_justdokmerenow=false;

	map<unsigned int,bool> species_of_nodes;
	if(param_useKmereHeuristic && nodes.size() >= kmereHeuristic_minNodes ){
		for(unsigned int i = 0 ; i < nodes.size() ; i++){
			species_of_nodes[graph[nodes[i]].species_id]=1;
		}// now the species_of_nodes.size() = number of different species in the current induced subgraph
	}

	if( param_useKmereHeuristic && ( 
			nodes.size() >= kmereHeuristic_minNodes && 
			(unsigned int)( (floattype)nodes.size()/( (floattype)kmereHeuristic_protPerSpecies * (floattype)species_of_nodes.size() ) ) >= kmereHeuristic_minNumberOfGroups ) ){ // min. 2^20 (~1e+6) nodes total and 4 nodes per bin

		do_kmereAlgorithm:

		if (debug_level > 0) cerr << getTime() << " [DEBUG] (kmere-heuristic) The current connected component is so large that the k-mere heuristic can be used. First: Testing if a normal split would result in a good partition (|.|>20%) of the CC."<< endl;

		map<unsigned int,bool> groupA, groupB, groupZero;
		for (unsigned int i = 0; i < y.size(); i++) {
			if      (abs(y[i]) < param_sep_purity) {
				groupZero[nodes[i]] = true;
			}
			else if (y[i] < 0) {
				groupA[nodes[i]] = true;
			}
			else {
				groupB[nodes[i]] = true;
			}
		}
		if( (groupA.size() == 0 && groupB.size() == 0) ){
			
			if (debug_level > 0) cerr << getTime() << " [DEBUG] (kmere-heuristic) All nodes are below the purty threshold. Continue without purity." << endl;

			for (unsigned int i = 0; i < y.size(); i++) {
				if (y[i] < 0) {
					groupA[nodes[i]] = true;
				}
				else {
					groupB[nodes[i]] = true;
				}
			}
		}
		if((floattype)groupA.size() > 0.2*(floattype)nodes.size() && 
			(floattype)groupB.size() > 0.2*(floattype)nodes.size() && 
			!fallback_justdokmerenow){

			if (debug_level > 0) cerr << getTime() << " [DEBUG] (kmere-heuristic) A normal split would result in a good partition (|.|>20%) of the CC, therefore returning now to the normal algorithm (no k-mere heuristic)." << endl;
			goto do_normalPartitionAlgorithm;
			
		}

		#ifdef DEBUG
			total_number_of_kmere_calls++;
		#endif

		if(species_of_nodes.size() == 0){ //for fallback function of impossible splits (then the species_of_nodes is not set)
			for(unsigned int i = 0 ; i < nodes.size() ; i++){
				species_of_nodes[graph[nodes[i]].species_id]=1;
			}// now the species_of_nodes.size() = number of different species in the current induced subgraph
		}
		
		floattype number_of_groups = (floattype)kmereHeuristic_protPerSpecies * (floattype)species_of_nodes.size() ;

		cerr << " [WARNING] (kmere-heuristic) A normal split would NOT result in a good partition (|.|>20%) of the CC, therefore  the k-mere heuristic is now used. The current connected component will be split in " << number_of_groups << " (= nodes per species <"<< (floattype)kmereHeuristic_protPerSpecies << "> * number of species <"<< (floattype)species_of_nodes.size() <<">) groups greedily accordingly to the fiedler vector."<< endl;

		//get the number of species
		vector<pair<floattype,unsigned int > > d;
		for(unsigned int i = 0 ; i < nodes.size() ; i++)
			d.push_back(make_pair(y[i],nodes[i]));

		sort( d.begin() , d.end() , comparator_pairfloattypeUInt); // sort the vector accordingly to the floattype value

		unsigned int from = 0;
		
		for(unsigned int i = 0 ; from < nodes.size()-1 ; i++){

			unsigned int to = (i+1) * (unsigned int)( (floattype)nodes.size()/number_of_groups ) ;

			map<unsigned int,bool> cur_group;
			for( unsigned int j = from ; j < to ; j++){
				if(j > d.size()-1)
					break;
				cur_group[d[j].second] = true;
			}
			removeExternalEdges(cur_group);

			from = to;
		}

	}else{
		do_normalPartitionAlgorithm:

		if (debug_level > 0) cerr << getTime() << " [DEBUG] perfoming a normal split of " << nodes.size()<< " nodes" << endl;

		// Store data about two groups (Zero cannot be assigned with certainty)
		map<unsigned int,bool> groupA, groupB, groupZero;

		for (unsigned int i = 0; i < y.size(); i++) {
			if      (abs(y[i]) < param_sep_purity) {
				groupZero[nodes[i]] = true;
			}
			else if (y[i] < 0) {
				groupA[nodes[i]] = true;
			}
			else { // y[i] > 0
				groupB[nodes[i]] = true;
			}
		}
		if( (groupA.size() == 0 && groupB.size() == 0) ){
			
			cerr << " [WARNING]  All nodes are below the purty threshold. Continue without purity." << endl;

			for (unsigned int i = 0; i < y.size(); i++) {
				if (y[i] < 0) {
					groupA[nodes[i]] = true;
				}
				else { // y[i] > 0
					groupB[nodes[i]] = true;
				}
			}
		}
		if (debug_level > 0) cerr << getTime() << " [DEBUG] splitting into ("<<groupA.size()<<","<<groupB.size()<<","<<groupZero.size()<<") sized groups!"<< endl;

		// Catch error in laplacien calcs
		if ( (groupA.size() == 0 && groupB.size() == 0) || 
			 ( (groupA.size() == 0 || groupB.size() == 0) && groupZero.size() == 0) ){
			if(useLapack && nodes.size() < 4000){
				cerr << "[WARNING]   Failed to partition subgraph with "<<nodes.size()<<" nodes into ("<<groupA.size()<<","<<groupB.size()<<","<<groupZero.size()<<") sized groups using lapack, now reiterating with power iteration." << endl;
				getConnectivity(nodes, false);
			}else{
				cerr << "[CRITICAL WARNING]   Failed to partition subgraph with "<<nodes.size()<<" nodes into ("<<groupA.size()<<","<<groupB.size()<<","<<groupZero.size()<<") sized groups, now using kmere heuristic as fall-back." << endl;
				fallback_justdokmerenow = true;
				goto do_kmereAlgorithm;
			}

		}else{
			{removeExternalEdges(groupZero);}
			{removeExternalEdges(groupA);}
			{removeExternalEdges(groupB);}
		}
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::splitGroups")) t_master["partition_graph::splitGroups"]=0;
		t_master["partition_graph::splitGroups"] += (floattype)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

}

string getTime(void) {
	time_t t = time(0);   // get time now
	struct tm * now = localtime( & t );
	ostringstream oss;
	if (now->tm_hour < 10) oss << "0";
	oss << now->tm_hour << ":";
	if (now->tm_min < 10) oss << "0";
	oss << now->tm_min << ":";
	if (now->tm_sec < 10) oss << "0";
	oss << now->tm_sec;
	return oss.str();
}


#ifdef DEBUG
	////////////////////// Debug
	void debug__print_edgelist (protein& node, const unsigned int index, const int node_id) {
		cerr << node_id << ": ";	
		for (unsigned int j = 0; j < node.edges.size(); j++) {
			if (j == index) cerr << "*";
			cerr << node.edges[j].edge << " ";
		}
		cerr << endl;
	}

	void debug__conn_integrity(vector<unsigned int>& nodes, floattype conn) {
		if (nodes.size() > 5) return;

		unsigned int sum = 0;
		for (unsigned int a = 0; a < nodes.size(); a++) {
			unsigned int from = nodes[a];
			sum += graph[from].edges.size();
		}

		sum /= 2;

		if (nodes.size() == 3) {
			if (sum == 2 && conn > 0.4) {
				cerr << "gs 3 with 2 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 3a").c_str() << endl;throw;
			}
			if (sum == 3 && conn < 1) {
				cerr << "gs 3 with 3 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 3b").c_str() << endl;throw;
			}
		}
		if (nodes.size() == 4) {
			if (sum == 3 && conn > 0.4) {
				cerr << "gs 4 with 3 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 4a").c_str() << endl;throw;
			}
			if (sum < 6 && conn == 1) {
				cerr << "gs 4 with <6 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 4b").c_str() << endl;throw;
			}
			if (sum == 6 && conn < 1) {
				cerr << "gs 4 with 6 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 4c").c_str() << endl;throw;
			}
		}
		if (nodes.size() == 5) {
			if (sum == 4 && conn > 0.4) {
				cerr << "gs 5 with 4 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 5a").c_str() << endl;throw;
			}
			if (sum < 10 && conn == 1) {
				cerr << "gs 5 with <10 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 5b").c_str() << endl;throw;
			}
			if (sum == 10 && conn < 1) {
				cerr << "gs 5 with 10 edges had conn " << conn << endl;
				cerr << string("integrity issue: connectivity of group size 5c").c_str() << endl;throw;
			}
		}
	}

	void debug__graph_integrity(vector<unsigned int>& nodes) {
		// For each node
		for (unsigned int a = 0; a < nodes.size(); a++) {
			unsigned int from = nodes[a];
	//		cerr << "From: " << from << endl;

			// For each edge in + direction
			for (unsigned int i = 0; i < graph[from].edges.size(); i++) {
				unsigned int to = graph[from].edges[i].edge;
	//			cerr << " To:  " << to << endl;

				if (to == from) {cerr << "ERROR: Edge from " << from << " to " << to << " is selfevident" << endl; cerr << string("integrity issue self hit").c_str() << endl;throw;}

				// Check reverse direction
				// Foreach edge in - direction
				bool found = false;
	//			cerr << " Back:";
				for (unsigned int j = 0; j < graph[to].edges.size(); j++) {
					unsigned int back = graph[to].edges[j].edge;
	//				cerr << " " << back;
					if (back == from) {found = true; continue;}
				}
	//			cerr << endl;
				if (!found) {
					cerr << "ERROR: Edge from " << from << " to " << to << " is unidirectional" << endl; cerr << string("integrity issue direction").c_str() << endl;throw;
				}
			}
		}
	}

	/* Auxiliary routine: printing a matrix */
	void debug__print_matrix( int m, int n, floattype* a, int lda ) {
	    int i, j;
	    for( i = 0; i < m; i++ ) {
	        for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
	        printf( "\n" );
	    }
	}
#endif

bool test__max_of_diag() {
	
	for(unsigned int j = 1 ; j < 100; j++){

		vector<unsigned int> nodes(j);
		for(unsigned int i = 0 ; i < nodes.size() ; i++) nodes[i] = i;
		
		vector<unsigned int> diag(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) diag[i] = i;

		if(max_of_diag(nodes,diag)!=j-1) return false;

	}

	{
		vector<unsigned int> nodes(10);
		for(unsigned int i = 0 ; i < nodes.size() ; i++) nodes[i] = i;
		
		vector<unsigned int> diag(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) diag[i] = 1;
		diag[0] = 2;

		if(max_of_diag(nodes,diag)!=2) return false;	
	}

	return true;
}

bool test__generate_random_vector() {
	
	for(unsigned int j = 1 ; j < 100; j++){
		std::vector<floattype> v = generate_random_vector(j);
		for(unsigned int i = 0 ; i < v.size()-1 ; i++) if(v[i] == v[i+1]) return false;
	}
	return true;
}

bool test__get_new_x() {
	
	for(unsigned int j = 2 ; j < 100; j++){

		graph.clear();
		vector<protein>().swap(graph);
	
		vector<unsigned int> nodes(j);
		map<unsigned int,unsigned int> mapping; 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){ nodes[i] = i; mapping[i]=i; }

		graph.resize(nodes.size());
		// create complete graph K_nodes.size() (all nodes are from nodes.size() different species) with edgeweight 10 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){
			for(unsigned int k = 0 ; k < nodes.size() ; k++){
				wedge w;
				w.edge = i;
				w.weight = 10;
				graph[nodes[k]].edges.push_back(w); 
				graph[nodes[k]].species_id = k; 
				graph[nodes[k]].full_name = "X";
			}
		}

		vector<floattype> x(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) x[i] = i; 

		//test weighted case
		vector<floattype> x_new = get_new_x(x,nodes,mapping,1);
	
		if(x_new.size() != x.size())return false;
		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(x_new[i]!=10*(nodes.size()-1)*(nodes.size())/2) return false;

		//test unweighted case
		x_new = get_new_x(x,nodes,mapping,0);
	
		if(x_new.size() != x.size())return false;
		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(x_new[i]!=(nodes.size()-1)*(nodes.size())/2) return false;
	}

	{
		graph.clear();
		vector<protein>().swap(graph);
		vector<unsigned int> nodes(5);
		graph.resize(nodes.size());
		map<unsigned int,unsigned int> mapping; 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){ 
			nodes[i] = i; mapping[i]=i; 
			graph[nodes[0]].species_id = 0; 
			graph[nodes[0]].full_name = "X";}


		/* 0--1--2
		**    |\
		**    3-4
		** edge weight is between (x,y) is the number xy (if x<y) else yx (not multiply, just concatenate)
		**/

		wedge w;
		
			w.edge = 1;
			w.weight = 10;
			graph[nodes[0]].edges.push_back(w); 
			w.edge = 0;
			graph[nodes[1]].edges.push_back(w); 
			w.edge = 1;
			w.weight = 13;
			graph[nodes[3]].edges.push_back(w); 
			w.edge = 3;
			graph[nodes[1]].edges.push_back(w); 
			w.edge = 1;
			w.weight = 14;
			graph[nodes[4]].edges.push_back(w); 
			w.edge = 4;
			graph[nodes[1]].edges.push_back(w); 
			w.edge = 3;
			w.weight = 34;
			graph[nodes[4]].edges.push_back(w); 
			w.edge = 4;
			graph[nodes[3]].edges.push_back(w); 
			w.edge = 1;
			w.weight = 12;
			graph[nodes[2]].edges.push_back(w); 
			w.edge = 2;
			graph[nodes[1]].edges.push_back(w); 

		vector<floattype> x(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) x[i] = i; 

		//test weighted case
		vector<floattype> x_new = get_new_x(x,nodes,mapping,1);

		if(x_new.size() != x.size())return false;
		if(x_new[0]!=1*10) return false;
		if(x_new[1]!=0*10+2*12+3*13+4*14) return false;
		if(x_new[2]!=1*12) return false;
		if(x_new[3]!=1*13+4*34) return false;
		if(x_new[4]!=1*14+3*34) return false;

		//test unweighted case
		x_new = get_new_x(x,nodes,mapping,0);
	
		if(x_new.size() != x.size())return false;
		if(x_new[0]!=1) return false;
		if(x_new[1]!=0+2+3+4) return false;
		if(x_new[2]!=1) return false;
		if(x_new[3]!=1+4) return false;
		if(x_new[4]!=1+3) return false;
	}

	return true;
}

bool test__makeOrthogonal(){
	for(unsigned int j = 1 ; j < 100; j++){
		std::vector<floattype> x = generate_random_vector(j);
		x=makeOrthogonal(x);
		floattype sum=0;
		for (unsigned int i = 0; i < x.size(); i++) {sum += x[i];}
		//cerr << sum << endl;
		if(sum > 1e-3 || sum < -1e-3)return false;
	}
	return true;
}

bool test__normalize(){
	for(unsigned int j = 1 ; j < 100; j++){
		std::vector<floattype> x = generate_random_vector(j);
		floattype len=0;
		x=normalize(x,&len);
		normalize(x,&len); // second normalize should return a length of 1 
		if(len > 1+1e-3 || len < 1-1e-3)return false;
	}
	return true;
}

bool test__getY(){

	for(unsigned int j = 1 ; j < 100; j++){
		
		vector<unsigned int> nodes(j);
		for(unsigned int i = 0 ; i < nodes.size() ; i++) nodes[i]=i;
		
		vector<unsigned int> diag(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) diag[i]=i;

		floattype max_degree = max_of_diag(nodes,diag);

		vector<floattype> x_hat(nodes.size()),x_new(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++){ x_hat[i]=i; x_new[i] = -(floattype)i*(2*max_degree - diag[i]); }

		vector<floattype> y = getY(max_degree,x_hat,x_new,nodes,diag);
 
 		if(y.size() != x_hat.size())return false;
		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(y[i] > 1e-3 || y[i] < -1e-3)return false;
	}

	{
		vector<unsigned int> nodes(5);
		for(unsigned int i = 0 ; i < nodes.size() ; i++) nodes[i]=i;
		
		vector<unsigned int> diag(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) diag[i]= (floattype)3;
		diag[3]=5;

		floattype max_degree = max_of_diag(nodes,diag); 

		vector<floattype> x_hat(nodes.size()),x_new(nodes.size());
		x_hat[0]= (floattype)2.1; 
		x_hat[1]= (floattype)-3; 
		x_hat[2]= (floattype)10; 
		x_hat[3]= (floattype)9.5; 
		x_hat[4]= (floattype)0; 

		x_new[0] = (floattype)-14.7; //calculated by hand
		x_new[1] = (floattype)21;
		x_new[2] = (floattype)-70;
		x_new[3] = (floattype)-47.5;
		x_new[4] = (floattype)0;

		vector<floattype> y = getY(max_degree,x_hat,x_new,nodes,diag);

		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(y[i] > 1e-3 || y[i] < -1e-3)return false;
	}

	return true;

}


bool test__splitGroups(){
	param_sep_purity = 0.001;

	for(unsigned int j = 6 ; j < 100; j+=2){

		graph.clear();
		vector<protein>().swap(graph);

		vector<unsigned int> nodes(j);
		species.resize(nodes.size());
		map<unsigned int,unsigned int> mapping; 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){ nodes[i] = i; mapping[i]=i; species[i]="X";}

		graph.resize(nodes.size());

		// create complete graph K_nodes.size() (all nodes are from nodes.size() different species) with edgeweight 10 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){
			for(unsigned int k = 0 ; k < nodes.size() ; k++){
				if(i==k)continue;
				wedge w;
				w.edge = i;
				w.weight = 10;
				graph[nodes[k]].edges.push_back(w); 
				graph[nodes[k]].species_id = k; 
				graph[nodes[k]].full_name = "X";
			}
		}

		vector<floattype> y(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size()/3 ; i++) y[i] = -100000; 
		for(unsigned int i = nodes.size()/3 ; i < 2*nodes.size()/3 ; i++) y[i] = 0; 
		for(unsigned int i = 2*nodes.size()/3 ; i < nodes.size() ; i++) y[i] = 100000; 

		splitGroups(y,nodes,true);

		for(unsigned int i = 0 ; i < nodes.size()/3 ; i++)
			for(unsigned int k = 0 ; k < graph[i].edges.size() ; k++) 
				if(graph[i].edges[k].edge >= nodes.size()/3) return false;
		for(unsigned int i = nodes.size()/3 ; i < 2*nodes.size()/3 ; i++)
			for(unsigned int k = 0 ; k < graph[i].edges.size() ; k++) 
				if(graph[i].edges[k].edge < nodes.size()/3 || graph[i].edges[k].edge >= 2*nodes.size()/3) return false; 
		for(unsigned int i = 2*nodes.size()/3 ; i < nodes.size() ; i++)
			for(unsigned int k = 0 ; k < graph[i].edges.size() ; k++) 
				if(graph[i].edges[k].edge < 2*nodes.size()/3) return false;
	}

	for(unsigned int j = 18 ; j < 100; j++){

		graph.clear();
		vector<protein>().swap(graph);
		
		vector<unsigned int> nodes(j);
		map<unsigned int,unsigned int> mapping; 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){ nodes[i] = i; mapping[i]=i; }

		graph.resize(nodes.size());
		// create complete graph K_nodes.size() (all nodes are from nodes.size() different species) with edgeweight 10 
		for(unsigned int i = 0 ; i < nodes.size() ; i++){
			for(unsigned int k = 0 ; k < nodes.size() ; k++){
				if(i==k)continue;
				wedge w;
				w.edge = i;
				w.weight = 10;
				graph[nodes[k]].edges.push_back(w); 
				graph[nodes[k]].species_id = k % 3; 
				graph[nodes[k]].full_name = "X";
			}
		}

		vector<floattype> y(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) y[i] = i % 3;  

		splitGroups(y,nodes,true);

		for(unsigned int i = 0 ; i < nodes.size() ; i++){
			for(unsigned int k = 0 ; k < graph[i].edges.size() ; k++){
				if(i==k)continue;
				if(abs((int)graph[graph[i].edges[k].edge].species_id -(int)graph[i].species_id)>2 ){ return false;}
			}
		}
	}

	return true;
}

#endif /* _PROTEINORTHOCLUSTERING */
