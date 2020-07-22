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


/** @defgroup debug DEBUG
 *  @{
 */

//#define DEBUG
//#define timeAnalysis //if set : need -std=c++11 for compiling

	#ifdef timeAnalysis
		#include <chrono>
		map<string,floattype> t_master;
	#endif

	#ifdef DEBUG
		////////////////////// Debug
		void debug__print_edgelist (protein& node, const unsigned int index, const int node_id);

		void debug__conn_integrity(vector<unsigned int>& nodes, floattype conn);

		void debug__graph_integrity(vector<unsigned int>& nodes);
	#endif

/** @} */ // end of debug



/** @defgroup globalVars Global Variables
 *  @{
 */

	// Parameters (prefix param_*)
	bool param_verbose 		= false; ///< By default no verbose is printed
	floattype param_con_threshold 	= 0.1;		///< as a reference: a chain a-b-c-d has 0.25
	unsigned int debug_level	= 0; ///< Debug stderr level.
	floattype param_sep_purity 	= 1e-7;			///< as a reference: a-b-c will give +/-0.707107 and 2.34857e-08 
	unsigned int param_max_nodes	= 16777216; ///< = 2^24. The absolute maximal number of nodes for a connected component. If exceeded then it is skipped entirely.
	floattype param_min_species	= 1; ///< Minimum number of species per connected component 
	string param_rmgraph            = "remove.graph"; ///< Name of the remove graph, containing all edges that are removed by this clustering algorithm.
	bool param_useWeights = true; ///< If true, then weights are used for the calculation of the algebraic connectivity
	unsigned int param_minOpenmp = 256; ///< the minimum size of a for-loop for openmp to activate (openmp has some initialization costs)
	bool param_useKmereHeuristic = true; ///< If true then kmere-heuristic is used. For ambigous splits (e.g. all but one entrie of the fiedler vector are below 0) then the fiedler vector is clustered in k groups instead. This will lead to k new groups instead of 2. See the variables kmereHeuristic_* for more informations when the kmere-heuristic will be activated.
	unsigned int param_maxRam_inKB = 16777216; ///< = 16 GB of memory as default. Decreasing this may lead to more power iterations and thus to a increase in runtime as a trade-off
	bool param_useLapack = true; ///< If true, then lapack is used accordingly to param_lapack_power_threshold_d

	// min/max number of alg con iterations
	unsigned int critical_min_nodes = 16777216; ///< Depricated
	const unsigned int min_iter = 16;			///< Minimal number of iterations for power iteration
	unsigned int param_max_iter = 8192;			///< Maximal number of iterations for power iteration
	floattype param_epsilon = 1e-8; ///< Epsilon for lapack functions dsyevr/ssyevr as well as the power iteration. Set analog to http://people.sc.fsu.edu/~jburkardt/c_src/power_method/power_method_prb.c
	const unsigned int kmereHeuristic_minNodes = 1048576; ///< = 2^20. The minimum number of nodes of the connected component for the kmere-heuristic.
	const unsigned int kmereHeuristic_protPerSpecies = 1; ///< The mimum number of proteins per species for the kmere-heuristic.
	const unsigned int kmereHeuristic_minNumberOfGroups = 3; ///< The minimum number of species for the kmere-heuristic.
	const unsigned int maxUseWeightsNumNodes = 1048576; ///< = 2^20. The maximum number of nodes for a connected component such that weights are still used for calculations. If the connected component is too large, the unweighted version is used instead to save time.
	floattype param_lapack_power_threshold_d = -1; ///< The minimum graph density for the power iteration method, lapacks (d|s)syevr is used otherwise. If -1 then the linear function is used instead : If d<10^(-5.2)*n, then lapack otherwise power.

	// Globals
	unsigned int species_counter = 0;		///< Number of species (in the input graph)
	unsigned int protein_counter = 0;		///< Number of proteins total (in the input graph)
	vector<string> species;					///< Species ID (number) -> Species name
	vector<protein> graph;					///< Graph containing all protein data (see the class protein for more informations)
	floattype last_stat = 0;				///< For buffering progress stats
	unsigned int edges = 0;					///< Total number of edges
	vector<shared_ptr<ofstream> > graph_clean;			///< File handler to store graph data
	vector<int> reorder_table;				///< Tells how proteins/species must be sorted
	unsigned long graph_ram_total_inKB=0; 	///< The internal size of the input graph in KB (this will be set in parse_file() )
	unsigned int num_cpus=1;				///< By default only one core is used (change this with -cpus)

/** @} */ // end of globalVars





/** @defgroup lapack Interface to the C part of LAPACK
*For more informations : http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_gaeed8a131adf56eaa2a9e5b1e0cce5718.html
*or here http://www.netlib.org/lapack/explore-3.1.1-html/dsyevr.f.html
 *  @{
 */

	extern "C" {
		extern void ssyevr_( char* jobz, char* range, char* uplo, int* n, float* a,
	                int* lda, float* vl, float* vu, int* il, int* iu, float* abstol,
	                int* m, float* w, float* z, int* ldz, int* isuppz, float* work,
	                int* lwork, int* iwork, int* liwork, int* info );
		extern void dsyevr_( char* jobz, char* range, char* uplo, int* n, double* a,
	                int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
	                int* m, double* w, double* z, int* ldz, int* isuppz, double* work,
	                int* lwork, int* iwork, int* liwork, int* info );

	}

	/** 
	*(d|s)syevr LAPACK function
	*(d|s) = floattype (s=float/d=double)
	*sy = symmetric 
	*ev = eigenvalue calculations
	*r = expert modus : more options, i.e. here only the first k=2 eigenvalues are needed
	*For more informations : http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_gaeed8a131adf56eaa2a9e5b1e0cce5718.html
	*or here http://www.netlib.org/lapack/explore-3.1.1-html/dsyevr.f.html
	*/
	template<class T> void dssyevr_( char* jobz, char* range, char* uplo, int* n, T* a,
	            int* lda, T* vl, T* vu, int* il, int* iu, T* abstol,
	            int* m, T* w, T* z, int* ldz, int* isuppz, T* work,
	            int* lwork, int* iwork, int* liwork, int* info ){ // dssyevr_<float>(...) calls the ssyevr_(...) function and dssyevr_<double>(...) calls the dsyevr_(...) function
		#if floatprecision_H == 1
			ssyevr_(jobz, range, uplo, n, a,
		            lda, vl, vu, il, iu, abstol,
		            m, w, z, ldz, isuppz, work,
		            lwork, iwork, liwork, info );
		#elif floatprecision_H == 2
			dsyevr_(jobz, range, uplo, n, a,
		            lda, vl, vu, il, iu, abstol,
		            m, w, z, ldz, isuppz, work,
		            lwork, iwork, liwork, info );
		#else
			cerr << string("Error: invalid floattype (should be either float=1 or double=2)!").c_str() << endl;throw;
		#endif
	}

/** @} */ // end of lapack






/** @defgroup class Auxiliary classes
 *  @{
 */

	/** A graph representation (as vector of indices of the induced subgraph of the variable 'graph') with some graph attributes (graph density, sum of node degrees)
	*
	*/
	class ConnectedComponent
	{
	 public:
		vector<unsigned int> m_content_CC; ///< ids of the induced subgraph of the global variable graph. E.g. m_content_CC=[0,4,6] -> proteins graph[0].full_name,graph[4].full_name,graph[6].full_name
	 	unsigned int d_sum; ///< sum of node degrees 
	 	double density; ///< the graph density calculated in the function BFS (at the end)

	 	/** Constructor: sum of node degrees = 0. density = -1 (will be set in BFS)
	 	*/
	 	ConnectedComponent(){
	 		d_sum=0;
	 		density=-1; 
	 	}
	 	
	 	/** Overloaded [] operator for quick access.
		*/
		unsigned int& operator[](unsigned int i){
			if(i > m_content_CC.size()){cerr << "[CRITICAL ERROR] out of bound in ConnectedComponent[]" << endl; throw;}
			return m_content_CC[i];
		}
		/** Overloaded [] operator for quick access.
		*/
		const unsigned int& operator[](unsigned int i)const{
			if(i > m_content_CC.size()){cerr << "[CRITICAL ERROR] out of bound in ConnectedComponent[]" << endl; throw;} 
			return m_content_CC[i];
		}

		/** Returns the number of nodes of this connected component
		*/
		unsigned int size(){ return m_content_CC.size();}
		unsigned int size()const { return m_content_CC.size();}

		/** Set operator for connected components
		*/
		void operator = (const ConnectedComponent &D ) { 
			m_content_CC = D.m_content_CC;
			d_sum = D.d_sum;
			density = D.density;
		}
		/** Pushes a new node (given as index of variable graph) to this connected component
		*/
		void push_back(unsigned int i) { 
			m_content_CC.push_back(i);
		}
	}; 

	/** Weighted edges
	*/
	struct wedge {
		unsigned int edge; ///< the adjacent node index of variable graph
		unsigned short weight; ///< the weight of this given edge
	};

	/** A protein (the variable 'graph' is just a vector of proteins)
	*/
	struct protein {
		vector<wedge> edges; ///< all the adjacent nodes along with weights -> wedge
		unsigned int species_id; ///< the species id of this protein
		string full_name; ///< the full input name of this protein
	};

	/** For sorting a vector of connected components by density
	*/
	struct compare_ConnectedComponents { //sort from large to small
		bool operator() (const ConnectedComponent &a, const ConnectedComponent &b) const {
		    return a.density < b.density;
		}
	};

	
	/** For sorting by protein degree
	*/
	struct protein_degree
	{
	    inline bool operator() (const pair<string,int> & p1, const pair<string,int>& p2)
	    {
	        return (p1.second > p2.second);
	    }
	};

/** @} */ // end of class


/** @defgroup globalFun Global Functions
 *  @{
 */
	/**
	 * standard breadth-first search (BFS)
	 */
	ConnectedComponent BFS(vector<bool> * done, unsigned int cur_node )

	/**
	 * Returns the peak (maximum so far) resident set size (physical
	 * memory use) measured in bytes, or zero if the value cannot be
	 * determined on this OS.
	 * https://github.com/mapbox/mapbox-gl-native/blob/master/test/src/mbgl/test/getrss.cpp
	 */
	size_t getPeakRSS( );

	/**
	 * Returns the current resident set size (physical memory use) measured
	 * in bytes, or zero if the value cannot be determined on this OS.
	 * https://github.com/mapbox/mapbox-gl-native/blob/master/test/src/mbgl/test/getrss.cpp
	 */
	size_t getCurrentRSS( );

	///////////////////////////////////////////////////////////
	// Main
	///////////////////////////////////////////////////////////
	/**
	 * Auxiliary function that prints the -help
	 */
	void printHelp() ;

	/**
	 * Estimates the memory usage for construction of the quadratic adjacency matrix in KB (the adjacency matrix is used by LAPACK) 
	 */
	unsigned int numberOfNodesToMemoryUsageLaplacian_inKB(unsigned int n);

	/**
	 * depricated
	 */
	unsigned int BFS_not_critical( vector<unsigned int> * done_withBacktrack, unsigned int start_node, unsigned int end_node );

	/**
	 * depricated
	 */
	bool criticalHeuristic(ConnectedComponent *cur_cc);


	///////////////////////////////////////////////////////////
	// Major partioning algorithm
	///////////////////////////////////////////////////////////

	/**
	 * Major partioning algorithm.
	 * First all connected components are identified using the BFS algorithm (see function BFS()).
	 * Then it iterates over all connected components and calculates the algebraic connectivity (either using LAPACK dsyevr or the power iteration) using different multithreading systems.
	 * The algebraic connectivity is calculated in getConnectivity().
	 * If the algebraic connectivity of a connected component is less than the minimum threshold (param_con_threshold) then the connected component is split into 2 (or more). 
	 * Connected component are split in getConnectivity() > splitGroups() > removeExternalEdges()
	 */
	void partition_graph();

	///////////////////////////////////////////////////////////
	// Basic Graph functions
	///////////////////////////////////////////////////////////
	/**
	 * Removes all edges between the given list of protein ids from graph entirely
	 */
	void clear_edges(vector<unsigned int>& nodes);

	/** Remove edges that go beyond the given group groups a, works in a single direction, no overlap in parallel runs if lists are without overlap.
	* This is used for perfoming a cut between groups based on the fiedler vector in partition_graph().
	* map<unsigned int,bool> : maps the node id of a connected component to the split group id (e.g. nodes=[1,5,10]. Splitting 1,5 from 10 : 1=>0, 5=>0, 10=>1)
	*/
	void removeExternalEdges(map<unsigned int,bool>& a);

	/** Remove edge in a single direction given an index -- openMP safe
	*/
	void remove_edge_index(const unsigned int node_id, const unsigned int index);

	///////////////////////////////////////////////////////////
	// File parser
	///////////////////////////////////////////////////////////

	/** Parses the given file and inserts the information into the variable 'graph'
	*/
	void parse_file(string file);

	///////////////////////////////////////////////////////////
	// Output
	///////////////////////////////////////////////////////////
	/** Sort vector of species (variable 'species')
	*/
	void sort_species(void);

	/** Progress stats
	*/
	void stats(floattype i, floattype size);

	/** Print header with species names
	*/
	void print_header();

	/** Group formatting for output
	*/
	void print_group(vector<unsigned int>& nodes, floattype connectivity);

	/** Calculate number of species
	*/
	floattype calc_group(vector<unsigned int>& nodes);

	///////////////////////////////////////////////////////////
	// Misc functions
	///////////////////////////////////////////////////////////
	/** Convert string to floattype (float or double)
	*/
	floattype string2floattype(string str);

	/** Split a string at a certain or multiple delimiters 'delimiters'
	*/
	void tokenize(const string& str, vector<string>& tokens, const string& delimiters = "\t");

	///////////////////////////////////////////////////////////
	// Algebraic connectivity functions
	///////////////////////////////////////////////////////////
	/** Returns maximum degree of given protein_ids -- openMP A ML
	*/
	unsigned int max_of_diag(vector<unsigned int>& nodes, vector<unsigned int>& diag);

	/** Generate random vector x of size size
	*/
	vector<floattype> generate_random_vector(const unsigned int size);

	/** determine new X, Formula (1) -- openMP B ML
	*/
	vector<floattype> get_new_x(vector<floattype> x, vector<unsigned int>& nodes, map<unsigned int,unsigned int> &mapping, bool isWeighted);

	/** Make vector x orthogonal to 1, Formula (2) -- openMP A ML
	*/
	vector<floattype> makeOrthogonal(vector<floattype> x);

	/** Normalize vector x, Formula (4) -- openMP A ML
	*/
	vector<floattype> normalize(vector<floattype> x, floattype *length);

	/** Qx, Formula (5) -- openMP A ML
	*/
	vector<floattype> getY(floattype max_degree, vector<floattype> x_hat, vector<floattype> x_new, vector<unsigned int>& nodes, vector<unsigned int>& diag);

	/** depricated
	*/
	unsigned int sumOutDegs(const vector<unsigned int>& nodes);

	/** The main function for calculating the algebraic connectivity.
	* If useLapack is true, the the lapack function is used (otherwise power iteration is used)
	* If the algebraic connectivity is less than the given threshold then the given connected component is split directly (thus no fiedler vector is returned)
	* The split is done in splitGroups().
	*/
	floattype getConnectivity(vector<unsigned int>& nodes, bool useLapack);

	bool comparator_pairfloattypeUInt ( const pair<floattype,unsigned int>& l, const pair<floattype,unsigned int>& r );

	/** Split connected component according to eigenvector (given in y) -- openMP B ML
	* Removes all edges between positive, negative and around zero entries of the eigenvector y
	* The around zero is defined with param_sep_purity.
	* If the requirements of the kmere-heuristic are satisfied, then the connected component is split in more groups:
	* The number of new groups is defined by = number of nodes / ( kmereHeuristic_protPerSpecies * number of species )
	* The eigenvector y is clustered in k groups by the maginitute of the entries. Then the connected component is split accordingly.
	*/
	void splitGroups(vector<floattype>& y, vector<unsigned int>& nodes , bool useLapack);

	/** gets the current unix time
	*/
	string getTime(void);

	/** Test max_of_diag()
	*/
	bool test__max_of_diag();

	/** Test generate_random_vector()
	*/
	bool test__generate_random_vector();

	/** Test get_new_x()
	*/
	bool test__get_new_x();

	/** Test makeOrthogonal()
	*/
	bool test__makeOrthogonal();

	/** Test normalize()
	*/
	bool test__normalize();

	/** Test getY()
	*/
	bool test__getY();

	/** Test splitGroups()
	*/
	bool test__splitGroups();

/** @} */ // end of debug



#endif /* _PROTEINORTHOCLUSTERING */
