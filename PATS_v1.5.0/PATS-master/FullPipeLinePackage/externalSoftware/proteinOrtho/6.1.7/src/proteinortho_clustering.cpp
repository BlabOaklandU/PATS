/*  
 *	Clustering algorithm for Proteinortho
 *	Reads edge list and splits connected components
 *	according to algebraic connectivity threshold
 *
 *	Last updated: 2022/09/21
 *	Author: Marcus Lechner, Paul Klemm
 */

#ifndef _PROTEINORTHOCLUSTERING
#define _PROTEINORTHOCLUSTERING

//#define DEBUG
//#define timeAnalysis

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
#include <thread>
#include <functional>
#include <cstdint>
#include <cstdio>
#include <queue>
#include <mutex>
#include <condition_variable>

#ifdef _OPENMP
#include <omp.h>
#endif

//the floatprecision_H has to be 1 for float and 2 for double

using namespace std;

#include <chrono>
#ifdef timeAnalysis
//	#include <chrono>
	map<string,float> t_master;
#endif

extern "C" {
	//(d|s)syevr LAPACK function
	//-      float (float/double)
	// --    symmetric 
	//   --- eigenvalue expert (more options, e.g. first k eigenvalues...)
	extern void ssyevr_( char* jobz, char* range, char* uplo, int* n, float* a,
				int* lda, float* vl, float* vu, int* il, int* iu, float* abstol,
				int* m, float* w, float* z, int* ldz, int* isuppz, float* work,
				int* lwork, int* iwork, int* liwork, int* info );
	extern void dsyevr_( char* jobz, char* range, char* uplo, int* n, double* a,
				int* lda, double* vl, double* vu, int* il, int* iu, double* abstol,
				int* m, double* w, double* z, int* ldz, int* isuppz, double* work,
				int* lwork, int* iwork, int* liwork, int* info );
}

struct wedge {unsigned int edge; unsigned short weight=0;};
struct protein {vector<wedge> edges; unsigned int species_id; string full_name=""; bool is_done=false;};

// Functions          
float string2float(string);
void tokenize(const string& , vector<string>& , const string&);
void parse_file(string);
void parse_abc_file(string);
// void remove_edge_index(const unsigned int, const unsigned int);
float getConnectivity_float(vector<unsigned int>*,bool,vector<float>*);
double getConnectivity_double(vector<unsigned int>*,bool,vector<double>*);
void partition_graph(void);
void print_header(void);
void sort_species(void);
void stats(unsigned int,unsigned int,bool);
string getTime(void);
	bool param_verbose 		= true;
	bool param_core 		  = false;
	bool param_abc 		  = false;
	float param_con_threshold 	= 0.1;		// as a reference: a chain a-b-c-d has 0.25
	unsigned int debug_level	= 0;
	float param_sep_purity 	= -1;		// as a reference: a-b-c will give +/-0.707107 and 2.34857e-08 
	unsigned int param_max_nodes	= 16777216; // 2^24
	bool param_max_nodes_was_set=false;
	float param_min_species	= 1;
	string param_rmgraph            = "remove.graph";
	bool param_useWeights = true;
	unsigned int param_minOpenmp = 256; // the minimum size of a for-loop for openmp to activate (openmp has some initialization costs)
	unsigned int param_coreMaxProteinsPerSpecies = 100;
	unsigned int param_coreMinSpecies = 0;
//	bool param_useKmereHeuristic = false; // deprecated
//	unsigned int param_maxRam_inKB = 16777216; // = 16 GB of memory, deprecated
	int param_useLapack = 1;

// min/max number of alg con iterations
// unsigned int critical_min_nodes = 16777216; // replaced
const unsigned int min_iter = 16;			// below this value, results may vary
unsigned int param_max_iter = 8192;			// below this value, results may vary
float param_epsilon = 1e-8; // analog to http://people.sc.fsu.edu/~jburkardt/c_src/power_method/power_method_prb.c
int param_double = 1;
unsigned int param_max_nodes_weight = 1048576; //2^20
float param_lapack_power_threshold_d = -1;//2048; //2^8

// Globals
unsigned int species_counter = 0;	// Species
unsigned int protein_counter = 0;	// Proteins
vector<string> species;			// Number -> Name
vector<protein> graph;			// Graph containing all protein data
unsigned int last_stat_lapack = 0;			// For progress stats
unsigned int last_stat_power = 0;			// For progress stats
bool last_stat_act = false;			// For progress stats
unsigned int edges = 0;			// number of edges
map<size_t,shared_ptr<ofstream> > proteinorthotmp_clean;
map<size_t,shared_ptr<ofstream> > tmp_debug;
map<size_t,shared_ptr<ofstream> > graph_clean;			// File to store graph data
vector<int> reorder_table;		// Tells how proteins/species must be sorted
unsigned long graph_ram_total_inKB=0; // approximately 4MB is the whole program without any input
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

#ifdef DEBUG
	unsigned int total_number_of_iterations_convergence = 0;
	unsigned int total_number_of_kmere_calls = 0;
	void debug__graph_integrity(vector<unsigned int>&);
	void debug__print_edgelist (protein&, const unsigned int, const int);
	void debug__conn_integrity(vector<unsigned int>&, float);
	void debug__print_matrix( int m, int n, float* a, int lda );
#endif

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif

std::string getEnvVar( std::string const & key )
{
	char * val = getenv( key.c_str() );
	return val == NULL ? std::string("") : std::string(val);
}

/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
	/* AIX and Solaris ------------------------------------------ */
	struct psinfo psinfo;
	int fd = -1;
	if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
		return (size_t)0L;      /* Can't open? */
	if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
	{
		close( fd );
		return (size_t)0L;      /* Can't read? */
	}
	close( fd );
	return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
	/* BSD, Linux, and OSX -------------------------------------- */
	struct rusage rusage;
	getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
	return (size_t)rusage.ru_maxrss;
#else
	return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
	/* Unknown OS ----------------------------------------------- */
	return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
	/* Windows -------------------------------------------------- */
	PROCESS_MEMORY_COUNTERS info;
	GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
	return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
	/* OSX ------------------------------------------------------ */
	struct mach_task_basic_info info;
	mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
	if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
		(task_info_t)&info, &infoCount ) != KERN_SUCCESS )
		return (size_t)0L;      /* Can't access? */
	return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
	/* Linux ---------------------------------------------------- */
	long rss = 0L;
	FILE* fp = NULL;
	if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
		return (size_t)0L;      /* Can't open? */
	if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
	{
		fclose( fp );
		return (size_t)0L;      /* Can't read? */
	}
	fclose( fp );
	return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
	/* AIX, BSD, Solaris, and Unknown OS ------------------------ */
	return (size_t)0L;          /* Unsupported. */
#endif 
}

/*void debug_log(string log){
	std::hash<std::thread::id>{}(std::this_thread::get_id());
	if(!log_out.count(tid))
		log_out[tid]=make_shared<ofstream>((param_rmgraph+to_string(tid)).c_str());
}*/

class dispatch_queue{
  // based on: https://github.com/embeddedartistry/embedded-resources/blob/master/examples/cpp/dispatch.cpp

  typedef std::function<void(void)> fp_t;

  public:
  dispatch_queue(std::string name, size_t thread_cnt = 1);
  ~dispatch_queue();

  void start();
  void waitTilDone();

  // dispatch and move
  void dispatch(fp_t&& op);
  void dispatch_allCores(fp_t&& op);

  // Deleted operations
  dispatch_queue(const dispatch_queue& rhs) = delete;
  dispatch_queue& operator=(const dispatch_queue& rhs) = delete;
  dispatch_queue(dispatch_queue&& rhs) = delete;
  dispatch_queue& operator=(dispatch_queue&& rhs) = delete;

  std::string name_;
  std::vector<std::thread> threads_;
  std::queue<fp_t> q_;
  std::queue<fp_t> q_allCores_;
  std::map<size_t,bool> is_computing;
  bool is_doing_allCores;
  bool lock_doing_allCores;
  std::mutex lock_;
  // The mutex class is a synchronization primitive that can be used to protect shared data from being simultaneously accessed by multiple threads.
  // std::unique_lock<std::mutex> lock(lock_); lock.lock() 
  // -> locks the mutex, blocks if the mutex is not available 
  std::condition_variable cv_;
  // cv_.notify_all();

  bool quit_ = false;

  void dispatch_thread_handler(size_t);
};

dispatch_queue::dispatch_queue(std::string name, size_t thread_cnt) : name_{std::move(name)}, threads_(thread_cnt){}

void dispatch_queue::start(){
	is_doing_allCores=false;
	lock_doing_allCores=false;
	if(debug_level>0)
	  std::cerr << "Creating dispatch queue: " << name_ << " threads:" << threads_.size() << " q:" << q_.size() << "\n";
  for(size_t i = 0; i < threads_.size(); i++){
    threads_[i] = std::thread(&dispatch_queue::dispatch_thread_handler,this,i);
  }
}

dispatch_queue::~dispatch_queue(){
	if(debug_level>0) std::cerr << getTime() << " Destructor: Destroying dispatch threads... threads:" << threads_.size() << " q:" << q_.size() << "\n";

  // Signal to dispatch threads that it's time to wrap up
  std::unique_lock<std::mutex> lock(lock_);
  quit_ = true;
  cv_.notify_all();
  lock.unlock();

  // Wait for threads to finish before we exit
  for(size_t i = 0; i < threads_.size(); i++){ 
    if(threads_[i].joinable()){ 
      //printf("Destructor: Joining thread %zu until completion\n", i);
      threads_[i].join(); 
    }
  }
 	if(debug_level>0)
 	  cerr << getTime() << " ~dispatch_queue done" << endl;
}

void dispatch_queue::waitTilDone(){
  while(true){
    bool all_idle=true;
    for (map<size_t,bool>::iterator it=is_computing.begin() ; it != is_computing.end(); it++) {
      if(it->second){all_idle=false;break;}
    }
    if(all_idle && q_.size()==0 && q_allCores_.size()==0){break;}
    sleep(1);
  }
}

void dispatch_queue::dispatch(fp_t&& op){
  std::unique_lock<std::mutex> lock(lock_);
  q_.push(std::move(op));
  cv_.notify_one();
}

void dispatch_queue::dispatch_allCores(fp_t&& op){
  std::unique_lock<std::mutex> lock(lock_);
  q_allCores_.push(std::move(op));
  cv_.notify_all();
}

void dispatch_queue::dispatch_thread_handler(size_t tid){
  std::unique_lock<std::mutex> lock(lock_);
  do{
    // Wait until we have data or a quit signal
    cv_.wait(lock, [this,tid] { 
    	return ( 
    		(
    			//!is_doing_allCores &&
    			( 
    				q_.size()
    				//( (q_.size()>0 && q_allCores_.size()==0) || (q_.size()>num_cpus*10 && q_allCores_.size() > 0) )
    				|| 
    				q_allCores_.size()
    			)
    		)
    		|| quit_ 
    	); 
    });
    // after wait, we own the lock
    is_computing[tid]=false; 
    //if( lock_doing_allCores && (q_.size()>num_cpus*10 || q_allCores_.size()==0) ){lock_doing_allCores=false;}

    if(!quit_ && 
    		q_.size() //&&
    		//( 
    		//	( q_allCores_.size() == 0 ) || 
    		//	( q_allCores_.size() > 0 && !lock_doing_allCores )
    		//)
    	){
    	//lock_doing_allCores=false;
    	unsigned int active_cpus=0;
      for (map<size_t,bool>::iterator it=is_computing.begin() ; it != is_computing.end(); it++) { active_cpus+=it->second; }
    	stats(q_.size()+active_cpus,q_allCores_.size(),1);
      auto op = std::move(q_.front());
      q_.pop();
      is_computing[tid]=true;
      // unlock now that we're done messing with the queue
      lock.unlock();
      op();
      is_computing[tid]=false;
      lock.lock();
    }else if(!quit_ && q_allCores_.size() //&& !is_doing_allCores
    	){
      bool all_idle=true;
      for (map<size_t,bool>::iterator it=is_computing.begin() ; it != is_computing.end(); it++) { if(it->second){all_idle=false;break;} }
      if( all_idle ){
      	//is_doing_allCores=true;
      	//lock_doing_allCores=true;
    		stats(q_.size(),q_allCores_.size()+1,0);
				if (debug_level > 0)
					cerr << tid << " allcores,all_idle:"<< all_idle << endl;
        auto op = std::move(q_allCores_.front());
        q_allCores_.pop();
        is_computing[tid]=true;
        // unlock now that we're done messing with the queue
        lock.unlock();
        op();
        is_computing[tid]=false;
      	//is_doing_allCores=false;
        lock.lock();
      }else{
      	lock.unlock();
      	//cv_.notify_all();
      	lock.lock();
      }
    }
  } while(!quit_);
}

///////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////
void printHelp() {
	cerr << "proteinortho_clustering - Spectral partitioning algorithm (last updated with proteinortho v6.1.6)" << "\n";
	cerr << "-----------------------------------------------------" << "\n";
	cerr << "This tool is part of Proteinortho" << "\n";
	cerr << "" << "\n";
	cerr << "Usage:   proteinortho_clustering [OPTIONS] graph_files..." << "\n";
	cerr << "Options: -verbose          report progress" << "\n";
	cerr << "         -conn float       minimal connectivity: keep dissecting if alg. connectivity is below conn (or if minspecies is satisfied) ["<<param_con_threshold<<"]" << "\n";
	cerr << "         -minspecies float stop clustering if ratio of genes/species of at least minspecies is reached regardless of the connectivity. This overrules the -conn threshold. ["<<param_min_species<<"]" << "\n";
	cerr << "         -core             stop clustering if a split would result in groups that do not span across all species of the inital connected component (unless the connectivity is very low). This overrules the -conn threshold.\n";
	cerr << "         -rmgraph STRING   output file name for the graph" << "\n";
	cerr << "         -seed int         seed value for srand [current unix time]" << "\n";
	cerr << "         -lapack int       use the lapack package for the computation of the algebraic connectivity. 0=no, 1=yes if applicable, 2=always ["<<param_useLapack<<"]" << "\n";
	cerr << "         -cpus int         the number of threads used for openMP ["<<num_cpus<<"]" << "\n";
	cerr << "         -coreMaxProts int  the maximum number of proteins per species for -core ["<<param_coreMaxProteinsPerSpecies<<"]" << "\n";
	cerr << "         -coreMinSpecies int  the minimum number of species for -core ["<<param_coreMinSpecies<<"]" << "\n";
	cerr << "         -abc   flag to indicate the input is a abc formatted graph file instead of a blast-graph (tab-separated). Input is expected to be undirected. c is the similarity score (0-USHRT_MAX) e.g. the blast bitscore" << "\n";
	cerr << "\ntechnical parameters:" << "\n";
	cerr << "         -test             various test-functions are called first [not set]" << "\n";
	cerr << "         -maxnodes int     only consider connected component with up to maxnodes nodes. If exceeded, greedily remove the worst 10 percent of edges (by weight) until satisfied [number of species ** 2]" << "\n";
	cerr << "         -maxweight int    only use the edge weights for connected components with maxweight nodes ["<<param_max_nodes_weight<<"]" << "\n";
	cerr << "         -epsilon float    convergence threshold ["<<param_epsilon<<"]" << "\n";
	cerr << "         -weighted bool    the spectral partition is calculated using the bitscores ["<<param_useWeights<<"]" << "\n";
	cerr << "         -double int      always use double precision. 0=no, 1=yes if applicable, 2=always ["<<param_double<<"]" << "\n";
	cerr << "         -minOpenmp int    the minimum number of nodes for parallel power iteration ["<<param_minOpenmp<<"]" << "\n";
	cerr << "         -powLapD | -power_d float	    the maximal graph density for the power iteration method, lapacks (d|s)syevr is used otherwise [adaptively choose optimal cutoff]" << "\n";
	cerr << "         -maxRunsConvergence int    the maximum number of runs for the calculation of the algebraic connectivity ["<<param_max_iter<<"]" << "\n";
}
// deprecated: 	
// cerr << "         -ram int          maximal used ram threshold for LAPACK and the input graph in MB [16384]" << "\n";
// cerr << "         -kmere bool	    use the kmere-split heuristic ["<<param_useKmereHeuristic<<"]" << "\n";
// cerr << "         -purity float     threshold for purity: treat float values between -purity and purity as 0. -1=adaptively choose ["<<param_sep_purity<<"]" << "\n";

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
				if (string2float(string(argv[paras])) == 0) {
					param_verbose = false; 
				}
			}
			else if (parameter == "-core") {
				param_core = true;
			}
			else if (parameter == "-abc") {
				param_abc = true;
			}
			else if (parameter == "-coreMaxProt" || parameter == "-coreMaxProts") {
				param_core = true;
				paras++;
				param_coreMaxProteinsPerSpecies = string2float(string(argv[paras]));
				if(param_con_threshold<0 || param_con_threshold>1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}
			else if (parameter == "-coreMinSpecies") {
				param_core = true;
				paras++;
				param_coreMinSpecies = string2float(string(argv[paras]));
				if(param_con_threshold<0 || param_con_threshold>1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}
			else if (parameter == "-conn") {
				paras++;
				param_con_threshold = string2float(string(argv[paras]));
				if(param_con_threshold<0 || param_con_threshold>1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}
//			else if (parameter == "-purity") {
//				paras++;
//				param_sep_purity = string2float(string(argv[paras]));
//				if(param_sep_purity<0 && param_sep_purity!=-1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
//			}
//			else if (parameter == "-ram") {
//				paras++;
//				param_maxRam_inKB = string2float(string(argv[paras]))*1e+3;
//			}
			else if (parameter == "-powLapD" || parameter == "-lapPowD" || parameter == "-power_d" || parameter == "-pld") {
				paras++;
				param_lapack_power_threshold_d = (string2float(string(argv[paras])));
			}
//			else if(parameter == "-kmere"){
//				paras++;
//				param_useKmereHeuristic = int(string2float(string(argv[paras])));
//			}
			else if(parameter == "-lapack"){
				paras++;
				param_useLapack = int(string2float(string(argv[paras])));
				if(param_useLapack!=0 && param_useLapack!=1 && param_useLapack!=2){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}
			else if (parameter == "-maxnodes") {
				paras++;
				param_max_nodes = string2float(string(argv[paras]));
				param_max_nodes_was_set=true;
			}
			else if (parameter == "-maxweight") {
				paras++;
				param_max_nodes_weight = string2float(string(argv[paras]));
			}
			else if (parameter == "-minspecies") {
				paras++;
				param_min_species = string2float(string(argv[paras]));
				if (param_min_species < 0) {cerr << string("-minspecies must at least be 0. Less than one gene per species is not possible as we only count those that have an entry.").c_str() << "\n";throw;}
			}
			else if (parameter == "-debug") {
				paras++;
				debug_level = int(string2float(string(argv[paras])));
			}
			else if (parameter == "-epsilon") {
				paras++;
				param_epsilon = string2float(string(argv[paras]));
				if(param_epsilon<0||param_epsilon>1){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}
			else if (parameter == "-double") {
				paras++;
				param_double = string2float(string(argv[paras]));
				if(param_double!=0 && param_double!=1 && param_double!=2){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}
			else if (parameter == "-minOpenmp") {
				paras++;
				param_minOpenmp = int(string2float(string(argv[paras])));
			}
			else if (parameter == "-weighted") {
				paras++;
				param_useWeights = int(string2float(string(argv[paras])));
			}
			else if (parameter == "-seed") {
				paras++;
				rand_seed = int(string2float(string(argv[paras])));
				if(rand_seed<0){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
			}else if (parameter == "-rmgraph") {
				paras++;
				param_rmgraph = string(argv[paras]);
			}else if(parameter == "-cpus"){
				paras++;
				#ifdef _OPENMP
					if(int(string2float(string(argv[paras])))<0){cerr << string("Error: invalid value '"+string(argv[paras])+" for argument "+parameter+"'!").c_str() << "\n";throw;}
					omp_set_dynamic(0);     // Explicitly disable dynamic teams
					num_cpus=int(string2float(string(argv[paras])));
					omp_set_num_threads(num_cpus); 
				#else
					cerr << "Error: missing openMP, please install/activate openMP if you want to use multiple cores.\n";
					throw;
				#endif
			}else if(parameter == "-maxRunsConvergence"){
				paras++;
				param_max_iter = int(string2float(string(argv[paras])));
			}else if(parameter == "-test"){
				bool test__max_of_diag_result = test__max_of_diag();
				bool test__generate_random_vector_result = test__generate_random_vector();
				bool test__get_new_x_result = test__get_new_x();
				bool test__makeOrthogonal_result = test__makeOrthogonal();
				bool test__normalize_result = test__normalize();
				bool test__getY_result = test__getY();
				cerr << "- test max_of_diag() : " << test__max_of_diag_result << "\n";
				cerr << "- test generate_random_vector() : "<< test__generate_random_vector_result << "\n";
				cerr << "- test get_new_x() : " << test__get_new_x_result << "\n";
				cerr << "- test makeOrthogonal() : " << test__makeOrthogonal_result << "\n";
				cerr << "- test normalize() : " << test__normalize_result << "\n";
				cerr << "- test getY() : " << test__getY_result << "\n";
				if( !test__max_of_diag_result || !test__generate_random_vector_result || !test__get_new_x_result || !test__makeOrthogonal_result || !test__normalize_result || !test__getY_result ){ 
					cerr << string("Error: tests failed !").c_str() << "\n";throw;
				}else{
					cerr << "All test passed." << "\n";
					return EXIT_SUCCESS;
				} 
			}
			else {
				printHelp();
				cerr << "\n" << "Sorry, unknown option '" << string(argv[paras]) << "'!" << "\n";
				return EXIT_FAILURE;
			}
		}

		srand(rand_seed);

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Debug level " << debug_level << "\n";

		if(param_core){param_con_threshold=999;} // for -core the -conn parameter is disabled 

		if(getEnvVar("OMP_NUM_THREADS") != "1"){
			// restart with OMP_NUM_THREADS=1
			string cmd="OMP_NUM_THREADS=1";
			for (paras = 0; paras < argc; paras++)
				cmd += " "+string(argv[paras]);
			return system(cmd.c_str());
		}

		// Parse files
		for (vector<string>::iterator it=files.begin() ; it != files.end(); it++) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Parsing file " << *it << "\n";
			if(param_abc){
				parse_abc_file(*it);
			}else{
				parse_file(*it);
			}
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   I know " << species_counter <<  " species with " << protein_counter << " proteins and " << edges << " edges in sum" << "\n";
		}

		graph_ram_total_inKB = getCurrentRSS()/1e+3;

		// #ifdef DEBUG
		if (debug_level == 42){ cerr << "graph_ram_total_inKB " << graph_ram_total_inKB << "\n"; return 1;} // 2609492 KB Maximum resident set size and the calculated graph_ram_total_inKB = 2910761.
		// #endif

//		if(graph_ram_total_inKB >= param_maxRam_inKB){
//			cerr << "\n" << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<"\n" << "WARNING: Putative memory overflow: the given input files ram " << graph_ram_total_inKB/1e+3 << " MB will presumably exceed the maximum ram threshold of "<< param_maxRam_inKB/1e+3 << " MB! You can solve this by giving this proteinortho at least "<< graph_ram_total_inKB/1e+3 << " MB ram with the argument '-ram "<< graph_ram_total_inKB/1e+3 << "' (or more). I will continue anyway, but memory overflow is now a risk."<< "\n" << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" <<"\n";
//			//return EXIT_FAILURE;
//		}

		// Free memory
		files.clear();
		vector<string>().swap(files);
		species2id.clear();
		map<string,int>().swap(species2id);
		protein2id.clear();
		map<string,int>().swap(protein2id);

		// Stats
		if (param_verbose) cerr << species_counter << " species" << "\n" << protein_counter << " paired proteins" << "\n" << edges << " bidirectional edges" << "\n";

		if(param_max_nodes_was_set==false){param_max_nodes = species_counter*species_counter; }

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Maximumum number of nodes for connectivity calculations is " << param_max_nodes << "\n";

		// Prepare sort of output
		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Sorting known species" << "\n";
		sort_species();

		// Write output header
		print_header();							

		// Open graph-removal file
		// string allRMgraphNames="";
		// for(unsigned int i = 0 ; i < num_cpus ; i ++){
		// 	stringstream ss;
		// 	ss << i;
		// 	allRMgraphNames+=(param_rmgraph+ss.str())+" ";
		// 	graph_clean.push_back(make_shared<ofstream>((param_rmgraph+ss.str()).c_str()));
		// }

		// // Clustering
		// if (debug_level > 0) cerr << getTime() << " [DEBUG]   Clustering" << "\n";
		partition_graph();

		if(debug_level>0) cerr << getTime() << "[DEBUG] done with clustering" << endl;

		// concat the remove graph output files
		ofstream OFS((param_rmgraph).c_str());
		for(map<size_t,shared_ptr<ofstream> > ::iterator it = graph_clean.begin() ; it != graph_clean.end() ; ++it){
			it->second->close();
			ifstream IFS((param_rmgraph+to_string(it->first)).c_str());
			if (IFS.is_open()) {
			while (!IFS.eof()) {
					string line;
					getline(IFS, line);
					if(line != "")
						OFS << line << "\n";
				}
			}
			IFS.close();
			// unlink tmp file
			int r = system(("rm '"+param_rmgraph+to_string(it->first)+"'").c_str());
		}
		OFS.close();

		if(debug_level>0) cerr << getTime() << "[DEBUG] done <ofstream>graph_clean" << endl;

		for(map<size_t,shared_ptr<ofstream> > ::iterator it = proteinorthotmp_clean.begin() ; it != proteinorthotmp_clean.end() ; ++it){
			it->second->close(); // close ofstream
			// dump to STDOUT
			ifstream IFS((param_rmgraph+"_proteinortho_tmp_"+to_string(it->first)).c_str());
			if (IFS.is_open()) {
			while (!IFS.eof()) {
					string line;
					getline(IFS, line);
					if(line != "")
						cout << line << "\n";
				}
			}
			IFS.close();
			// unlink tmp file
			int r = system(("rm '"+param_rmgraph+"_proteinortho_tmp_"+to_string(it->first)+"'").c_str());
		}

		if(debug_level>0) cerr << getTime() << "[DEBUG] done <ofstream>proteinorthotmp_clean" << endl;

		for(map<size_t,shared_ptr<ofstream> > ::iterator it = tmp_debug.begin() ; it != tmp_debug.end() ; ++it){
			it->second->close(); // close ofstream
		}
		
		if(debug_level>0) cerr << getTime() << "[DEBUG] done done" << endl;

		// if(system(("cat "+allRMgraphNames+" >"+param_rmgraph).c_str())!=0 || system(("rm "+allRMgraphNames).c_str())!=0){
		// 	cerr << "[ERROR]   cannot concatenate remove graphs" << "\n";
		// 	return EXIT_FAILURE;
		// }

		#ifdef timeAnalysis
			for(map<string,float>::iterator it = t_master.begin() ; it != t_master.end() ; ++it) cerr << (*it).first << " " << (*it).second << "\n";
		#endif
		#ifdef DEBUG
			cout << "conv:" << total_number_of_iterations_convergence << ", kmere_calls:" << total_number_of_kmere_calls<< "\n";
		#endif

	}
	catch(string& error) {
		cerr << "[ERROR]   " << error << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

unsigned int numberOfNodesToMemoryUsageLaplacian_inKB(unsigned int n){
	return (unsigned int)n*(unsigned int)n*sizeof(float)/1e+3;
}

string getCCid(vector<unsigned int> nodes){
	unsigned int n = nodes.size();
	double id=0;
	if(n==0){return "0";}
	map<unsigned int,unsigned int> mapping;
	for (unsigned int i = 0; i < (unsigned int)n; i++) {mapping[nodes[i]] = i;}
	for (unsigned int i = 0 ; i < (unsigned int)n ; i++){
		unsigned int from = nodes[i]; 
		if(!mapping.count(from)){continue;}
		unsigned int sum = 0;			
		for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
			unsigned int to = graph[from].edges[j].edge;
			if(!mapping.count(to)){continue;}
			//id+=to_string(from)+"-"+to_string(graph[from].edges[j].edge)+":"+to_string(graph[from].edges[j].weight)+",";
			id+=(from+1)*(graph[from].edges[j].edge+1)*graph[from].edges[j].weight;
		}
	}
	return to_string((unsigned int)id);
}

class ConnectedComponent // A graph representation (as vector of idx of the induced subgraph of 'graph') with some graph attributes (graph density, sum of node degrees)
{
	public:
	vector<unsigned int> m_content_CC; //ids of the induced subgraph of graph
	unsigned int d_sum; // sum of node degrees 
	double density;
	unsigned int species_num;

	ConnectedComponent(){
		d_sum=0;
		density=-1; // mark for calculation !
		species_num=0;
	}

	void calc_dsum(){
		d_sum=0;
		map<unsigned int,bool> m_content_CC_set;
		for (unsigned int i = 0; i < m_content_CC.size(); i++) {
			m_content_CC_set[m_content_CC[i]]=true;
		}
		for (unsigned int i = 0; i < m_content_CC.size(); i++) {
			for (unsigned int j = 0; j < graph[m_content_CC[i]].edges.size(); j++) {
				if(m_content_CC_set.count(graph[m_content_CC[i]].edges[j].edge)){
					d_sum+=1;
				}
			}
		}
	}
	void calc_density(){
		density=size()>1 ? (((double)d_sum))/(((double)(((double)size()-1.0)*(double)size()))) : 1;
	}
	
	unsigned int& operator[](unsigned int i){
		if(i > m_content_CC.size()){cerr << "[CRITICAL ERROR] out of bound in ConnectedComponent[]" << "\n"; throw;}
		return m_content_CC[i];
	}
	const unsigned int& operator[](unsigned int i)const{
		if(i > m_content_CC.size()){cerr << "[CRITICAL ERROR] out of bound in ConnectedComponent[]" << "\n"; throw;} 
		return m_content_CC[i];
	}

	unsigned int size(){ return m_content_CC.size();}
	unsigned int size()const { return m_content_CC.size();}

	void operator = (const ConnectedComponent &D ) { 
		m_content_CC = D.m_content_CC;
		d_sum = D.d_sum;
		density = D.density;
		species_num = D.species_num;
	}
	void push_back(unsigned int i) { 
		m_content_CC.push_back(i);
	}
};

void partition_CC(ConnectedComponent, dispatch_queue *, bool, unsigned int);
void find_CCs(dispatch_queue*);
void find_CCs_givenNodes(dispatch_queue*,vector<unsigned int>);
void print_group(ConnectedComponent& , float, size_t, bool);
float calc_group(vector<unsigned int>*);

ConnectedComponent BFS( map<unsigned int, bool> * done, unsigned int cur_node , float cut_off ){

	/*
	*
	* Simple BFS implementation
	* done map is used to identify allready visited nodes (can be set before call to indicate forbidden nodes)
	* cut_off = ignore all edges below this cut off
	* 
	*/
	size_t tid=std::hash<std::thread::id>{}(std::this_thread::get_id());
	if(!graph_clean.count(tid))
		graph_clean[tid]=make_shared<ofstream>((param_rmgraph+to_string(tid)).c_str());

	ConnectedComponent ret; //return vector
	list<unsigned int> q;
	q.push_back(cur_node);
	(*done)[cur_node]=true;

	while(q.size()>0){

		list<unsigned int> q_new;

		for(list<unsigned int>::iterator it = q.begin() ; it != q.end() ; ++it){
			cur_node = *it;
			ret.push_back(cur_node);
			ret.d_sum+=graph[cur_node].edges.size();
			(*done)[cur_node] = true;

			if(graph[cur_node].is_done){continue;}

			for (unsigned int j = 0; j < graph[cur_node].edges.size(); j++) {

				if(graph[graph[cur_node].edges[j].edge].is_done){continue;}

				if(graph[cur_node].edges[j].weight < cut_off){
					protein node_i = graph[cur_node];
					protein node_j = graph[node_i.edges[j].edge];
					(*graph_clean[tid]) << node_i.full_name << "\t" << species[node_i.species_id] << "\t" << node_j.full_name << "\t" << species[node_j.species_id] << "\n";
					continue;
				} // ignore

				unsigned int adjacency_node = graph[cur_node].edges[j].edge;

				if(adjacency_node > graph.size()){
					cerr << string("[ERROR] : Input graph is invalid. The node "+graph[cur_node].full_name +" is reporting an edge/adjacent node, that is not present in the graph.").c_str() << "\n";throw;
				}

				if( !done->count(adjacency_node) || !(*done)[adjacency_node] ){
					(*done)[adjacency_node] = true;
					q_new.push_back(adjacency_node);
				}
			}
		}

		q=q_new;
	}
	return ret;
}

struct compare_ConnectedComponents { //sort from large to small
	bool operator() (const ConnectedComponent &a, const ConnectedComponent &b) const {
		return a.density < b.density;
	}
};

void removeLowQualityEdges( ConnectedComponent cur_cc , dispatch_queue *q ){
	/*
	 * remove low quality edges
	 * 
	 * find the range of values of this CC -> define a cut-off = cut_off=min_w+0.1*(max_w-min_w);
	 * remove all edges below that value
	 * redo BFS and start over again for this cluster until -maxnodes are satisfied
	 * 
	 */
	if(debug_level>0) cerr << " [INFO]  removeLowQualityEdges " << cur_cc.size() << " start" << "\n";

	// find the lowest 10% of edge values (min+0.1*(max-min))
	float min_w = -1; 
	float max_w = -1;
	unsigned int original_number_nodes=cur_cc.size();
	for (unsigned int i = 0; i < cur_cc.size(); i++) {
		unsigned int from=cur_cc[i];
		for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
			unsigned int to = graph[from].edges[j].edge;
			unsigned int weight = graph[from].edges[j].weight;
			if(min_w == -1 || weight < min_w){min_w=weight;} 
			if(min_w == -1 || weight > max_w){max_w=weight;}
		}
	}
	float cut_off=min_w+0.1*(max_w-min_w);

	if(debug_level>0) cerr << " [INFO]  removeLowQualityEdges " << cur_cc.size() << " cut_off=" << cut_off << " min_w=" << min_w << " max_w=" << max_w << "\n";

	map<unsigned int, bool> done;	// Keep track on what was done (for each node)
	bool allNodesAreDone = false;

	vector<ConnectedComponent> CC; // vector of all connected components found

	while( true ){ // CC.size() < num_cpus / gather up to num_cpus connected components

		allNodesAreDone = true;

		for (unsigned int id = 0 ; id < cur_cc.size() ; id++) {
		
			unsigned int protein_id = cur_cc[id];

			if (done.count(protein_id) && done[protein_id]){continue;}// We were here already
			if(debug_level>0) cerr << "find_CCs:start at "<< protein_id << "<" << graph.size() << endl;

			done[protein_id]=true; // mark this node
			ConnectedComponent cur_cc = BFS(&done,protein_id,cut_off); // get the CC of the current node (protein_id) 

			if(cur_cc.size()==original_number_nodes){
				// the cutoff is too low -> increase by 10% and restart
				cut_off+=0.1*(max_w-min_w);
				allNodesAreDone=false;
				break;
			}

			// for -core, skip CC with less than param_coreMinSpecies species
			if( param_core && cur_cc.species_num==0 ){
				if(cur_cc.size()<param_coreMinSpecies){allNodesAreDone=false;break;} // if there are less nodes than wanted species -> skip
				map<unsigned int,bool>cc_species;
				for (int i = 0; i < cur_cc.size(); ++i){ cc_species[graph[cur_cc[i]].species_id]=true; }
				cur_cc.species_num=cc_species.size();
				if(cc_species.size()<param_coreMinSpecies){allNodesAreDone=false;break;} // less species as wanted -> skip
			}

			// Skip those that are too large (try heuristic)
			if (cur_cc.size() > param_max_nodes) {

				// reset done vector
				//for (int i = 0; i < cur_cc.size(); ++i){ done[cur_cc[i]]=false; }
				if(debug_level>0) cerr << " [WARNING]  Found a very large connected component that contains " << cur_cc.size() << ">" << param_max_nodes << " (maxnodes) elements. This behavior can be adjusted using -maxnodes. Now using a slow heuristic: try to identify and remove edges." << "\n";
				q->dispatch([cur_cc,q]{ removeLowQualityEdges(cur_cc,q); });
				//protein_id--; 
				continue;
			}

			cur_cc.calc_dsum();
			cur_cc.calc_density();
			if(cur_cc.density > 1){ cerr << "[WARNING] : The input graph has duplicated edges, this lead to an invalid graph density of " << cur_cc.density << " (should be <1). Please clean the .blast-graph with 'proteinortho.pl --cleanblast --step=3 --project=...' or use the cleanupblastgraph tool in src/ to remove the duplicated edges." << "\n"; throw; }
			if (debug_level > 0) cerr << getTime() << " [DEBUG:removeLowQualityEdges] Found connected component: " << cur_cc.size() << " proteins (ID: " << protein_id << "), graph density="<< cur_cc.density << ", sum of degrees="<< cur_cc.d_sum << " ini from " << graph[protein_id].full_name<< "\n";

			q->dispatch([cur_cc,q]{ partition_CC(cur_cc,q,true,false); });
			allNodesAreDone=false;
			break;
		}
		if(allNodesAreDone)break; // no additional CC can be found -> done
	}
	if(debug_level>0) cerr << " [INFO]  removeLowQualityEdges " << cur_cc.size() << " done\n";

	if(debug_level>0) cerr << "cut_off="<<cut_off << " min_w="<<min_w<<" max_w="<<max_w<< endl;
	// for (unsigned int i = 0; i < cur_cc.size(); i++) {
	// 	unsigned int from=cur_cc[i];
	// 	auto it = std::remove_if(
	// 		graph[from].edges.begin(), 
	// 		graph[from].edges.end(), 
	// 		[cut_off,from](wedge cur_el)->bool 
	// 		{
	// 			if(debug_level>0 && cur_el.weight <= cut_off) cerr << "[WARNING] found bad-edge "<< from << "-"<< cur_el.edge << " weight=" << cur_el.weight << endl;
	// 			return cur_el.weight <= cut_off;
	// 		}
	// 	);
	// 	graph[from].edges.erase(it, graph[from].edges.end());
	// }
}

void find_CCs_givenNodes(dispatch_queue *q, vector<unsigned int> todo_work ){

	map<unsigned int, bool> done;	// Keep track on what was done (for each node)
	bool allNodesAreDone = false;

	vector<ConnectedComponent> CC; // vector of all connected components found

	// find all nodes outside todo_work and remove them from the done vector !
	for (unsigned int i = 0 ; i < todo_work.size() ; i++) {
		unsigned int from = todo_work[i]; 
		done[ from ] = 0;
	}
	for (unsigned int i = 0 ; i < todo_work.size() ; i++) {
		unsigned int from = todo_work[i]; 
		for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
			if(!done.count(graph[from].edges[j].edge))
				done[ graph[from].edges[j].edge ]=1;
		}
	}

	while( true ){ // CC.size() < num_cpus / gather up to num_cpus connected components

		allNodesAreDone = true;

		for (unsigned int i = 0 ; i < todo_work.size() ; i++) {

			unsigned int protein_id=todo_work[i];
			if ( graph[protein_id].is_done ){continue;}
			if ( done.count(protein_id) && done[protein_id] ){continue;}// We were here already

			//min_i=protein_id;
			if (debug_level > 0) cerr << getTime() << " [DEBUG:find_CCs_givenNodes] start="<< protein_id << "\n";

			done[protein_id]=true; // mark this node
			ConnectedComponent cur_cc = BFS(&done,protein_id,0); // get the CC of the current node (protein_id) 

			// Do not report singles
			if (cur_cc.size() < 2) {continue;} // singletons are from no interest

			if(debug_level>0) cerr << "ConnectedComponent: @"<< getCCid(todo_work) << "=>" << cur_cc.size() << "@" << getCCid(cur_cc.m_content_CC) << endl;

			// Skip those that are too large (try heuristic)
			if (cur_cc.size() > param_max_nodes) {
				// reset done vector
				// for (int i = 0; i < cur_cc.size(); ++i){ done[cur_cc[i]]=false; }
				if(debug_level>0) cerr << " [WARNING]  Found a very large connected component that contains " << cur_cc.size() << ">" << param_max_nodes << " (maxnodes) elements. This behavior can be adjusted using -maxnodes. Now using a slow heuristic: try to identify and remove edges." << "\n";
				removeLowQualityEdges(cur_cc,q);
				continue;
			}

			if(param_core && cur_cc.species_num==0){
				if(cur_cc.size()<param_coreMinSpecies){allNodesAreDone=false;break;} // if there are less nodes than wanted species -> skip
				map<unsigned int,bool>cc_species;
				for (int i = 0; i < cur_cc.size(); ++i){ cc_species[graph[cur_cc[i]].species_id]=true; }
				cur_cc.species_num=cc_species.size();
				if(cc_species.size()<param_coreMinSpecies){allNodesAreDone=false;break;} // less species as wanted -> skip
			}
			
			cur_cc.calc_dsum();
			cur_cc.calc_density();
			if(cur_cc.density > 1){ cerr << "[WARNING] : The input graph has duplicated edges, this lead to an invalid graph density of " << cur_cc.density << " (should be <1). Please clean the .blast-graph with 'proteinortho.pl --cleanblast --step=3 --project=...' or use the cleanupblastgraph tool in src/ to remove the duplicated edges." << "\n"; throw; }
			if (debug_level > 0) cerr << getTime() << " [DEBUG:find_CCs_givenNodes] Found connected component: " << cur_cc.size() << " proteins (ID: " << protein_id << "), graph density="<< cur_cc.density << ", sum of degrees="<< cur_cc.d_sum << "\n";

			q->dispatch([cur_cc,q]{ partition_CC(cur_cc,q,true,false); });
			
			allNodesAreDone=false;
			//break;
		}
		if(allNodesAreDone)break; // no additional CC can be found -> done
	}
}

void partition_CC(ConnectedComponent cur_cc, dispatch_queue *q, bool do_lapack, unsigned int restarted){
	size_t tid=std::hash<std::thread::id>{}(std::this_thread::get_id());

	vector<float> x_hat;
	float connectivity;

	if( param_double == 2 && restarted==0 ){restarted=1;} // force double
	if( param_double == 0 && restarted>0  ){return;} // force float

	if(cur_cc.size()==0){return;} // sanity check

	if(cur_cc.size()==1){ 
		connectivity=1;return;
	}else{
	
		if(cur_cc.density==-1){
			cur_cc.calc_dsum();
			cur_cc.calc_density();
		}

		if(!do_lapack){
			if (debug_level > 0)
				cerr << "power size:"<< cur_cc.size() << " d:"<< cur_cc.density<<" cut:"<<(double)3.542144e-06*(double)cur_cc.size()-0.01701203<< endl;
			if(debug_level==10){
				cerr << getTime() << " [debug:10] extract cc" << endl;
				ofstream OFS(("cur_cc_n"+to_string(cur_cc.size())+"_d"+to_string(cur_cc.density)+".ids").c_str());
				for (int i = 0; i < cur_cc.size(); ++i){
					OFS << graph[cur_cc[i]].full_name << "\n";
				}
				OFS.close();
			}
			if(restarted>0){
				vector<double> x_hat_dbl;
				connectivity = getConnectivity_double(&cur_cc.m_content_CC,false,&x_hat_dbl);
				x_hat = std::vector<float>(x_hat_dbl.begin(), x_hat_dbl.end());
			}else
				connectivity = getConnectivity_float(&cur_cc.m_content_CC,false,&x_hat);
		}else if(
			do_lapack && (
				cur_cc.size() >= 32768 || 
				( 
					( param_lapack_power_threshold_d> -1 && cur_cc.density < param_lapack_power_threshold_d ) || 
					( param_lapack_power_threshold_d==-1 && cur_cc.density < (double)3.542144e-06*(double)cur_cc.size()-0.01701203 ) 
				)
			)
		){
			if (debug_level > 0)
				cerr << "dispatch all cores size:"<< cur_cc.size() << " d:"<< cur_cc.density<<" cut:"<<(double)3.542144e-06*(double)cur_cc.size()-0.01701203<< endl;
			q->dispatch_allCores([cur_cc,q]{ partition_CC(cur_cc,q,false,false); });
			return;
		}else{
			if (debug_level > 0)
				cerr << "lapack size:"<< cur_cc.size() << " d:"<< cur_cc.density<< endl;
			if(restarted>0){
				vector<double> x_hat_dbl;
				connectivity = getConnectivity_double(&cur_cc.m_content_CC,true,&x_hat_dbl);
				x_hat = std::vector<float>(x_hat_dbl.begin(), x_hat_dbl.end());
			}else
				connectivity = getConnectivity_float(&cur_cc.m_content_CC,true,&x_hat);
		}
	}

	if(param_core && cur_cc.species_num==0){
		map<unsigned int,bool>cc_species;
		for (int i = 0; i < cur_cc.size(); ++i){ cc_species[graph[cur_cc[i]].species_id]=true; }
		cur_cc.species_num=cc_species.size();
		//if(cc_species.size()!=species.size()){allNodesAreDone=false;break;}
	}
	if(param_core && cur_cc.species_num < param_coreMinSpecies){return;}

	if(debug_level>0) cerr << "ConnectedComponent: @"<< getCCid(cur_cc.m_content_CC) << " connectivity=" << connectivity << endl;

	if (connectivity < 0 || connectivity > param_con_threshold) {
		// cerr << " [DEBUG] done "<<connectivity<<"\n";
		print_group(cur_cc,connectivity,tid,false);
		return;
	}
	// 5.17 new threshold option overwrites connectivity
	if (param_min_species >=1 ) {
		if (debug_level > 0) cerr << getTime() << " [DEBUG]  Start the calculation of the average gene/species score " << "\n";
		float avg = calc_group(&cur_cc.m_content_CC);
		if (debug_level > 0) cerr << getTime() << " [DEBUG]   Found " << avg << " genes/species on average. User asked for at least " << param_min_species << "\n";
		if (avg <= param_min_species) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Group is going to be accepted despite connectivity " << connectivity << "\n";
			// just to be safe from infinit loops due to rounding
			if (connectivity == 0) connectivity = 0.001;
			// no splitting despite bad connectivity
			
			if(debug_level > 0)
				cerr << " [DEBUG] param_min_species done "<<connectivity<<"\n";
			print_group(cur_cc,connectivity,tid,false);
			return;
		}
	}

	// Store data about two groups (Zero cannot be assigned with certainty)
	ConnectedComponent groupA, groupB, groupZero;
	map<unsigned int,bool> groupB_set;
	map<unsigned int,bool> groupZero_set;

	if(param_sep_purity==-1){
		// automatic purity detection

		unsigned int n = cur_cc.size();
		map<unsigned int,unsigned int> mapping;
		for(unsigned int i = 0; i < (unsigned int)n; i++){ mapping[cur_cc[i]] = i; }
		
		/*
		 * Automatic prurity detetection (v6.1.2)
		 * purity = threshold of nodes that are neither part of the bisection result of the current split (negative vs positive x_hat entries)
		 *
		 * purity is imporatant if there are multiple clusters that can be split or there is no good split
		 * identified nodes (in purity boundary) are put back to the queue stack (=groupZero)
		 * 
		 * maximal purity = 0.12 -> get all nodes with |x_hat(v)|<0.12
		 * identify all possible purity thresholds 
		 * each purity threshold is evaluated with the total sum of edge weights that are removed with this cut
		 * find most optimal purity cut-off
		 */

		float new_purity=0;
		float best_cut=0; // stores the current best cut (sum of all edges removed by the given purity), first take a purity of 0 (no purity) and set this value (split + and - compartments of x_hat)
		vector<float> x_hat_candidate_purity = x_hat;
		for (unsigned int i = 0; i < x_hat_candidate_purity.size(); i++) { x_hat_candidate_purity[i]=abs(x_hat_candidate_purity[i]); } // purity is a positive threshold
		sort(x_hat_candidate_purity.begin(),x_hat_candidate_purity.end()); // these values are used as purity tests. Later only subset is used of values not too similar to each other
		for (unsigned int i = 0; i < x_hat.size(); i++) {
			if(x_hat[i] < 0) {
				unsigned int from=cur_cc[i];
				for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
					unsigned int to = graph[from].edges[j].edge;
					if(!mapping.count(to)){continue;}
					if(x_hat[mapping[to]] > 0) {
						best_cut += graph[from].edges[j].weight;							
					}
				}
			}
		}
		if(debug_level>0)	cerr << " initial purity=0 best_cut=" << best_cut << endl;
		float last_test_purity = -1; // omit purity test that are very close to previous values, to reduce runtime
		int num_purity_tests=0; // test all of the first 25 purity values anyway
		for(unsigned int pi = 0; pi < (unsigned int)n; pi++){
			//if( x_hat_sort[ pi ] < -0.12 ){continue;} // skip large negative values
			if( x_hat_candidate_purity[ pi ] > 0.12 ){break;} // there is no further valid value (sorted values)
			float test_purity = x_hat_candidate_purity[ pi ];
			if( test_purity == 0 || 
					test_purity == last_test_purity || 
					(num_purity_tests>25 && test_purity-last_test_purity<0.01) ){continue;} // current value is very similar to last purity -> omit
			float cut_value = 0;
			num_purity_tests++;
			last_test_purity = test_purity;
			unsigned int num_nodes = 0; // the number of nodes within the purity bounds
			for (unsigned int i = 0; i < x_hat.size(); i++) {
				if(abs(x_hat[i]) <= test_purity) { // add all the edges between this node (below purity) and nodes above the purity threshold to the cut value
					num_nodes++;
					unsigned int from=cur_cc[i];
					for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
						unsigned int to = graph[from].edges[j].edge;
						if(!mapping.count(to)){continue;} // sanity check if the edge is in the current CC
						if(abs(x_hat[mapping[to]]) > test_purity) {
							cut_value+=graph[from].edges[j].weight;							
						}
					}
				}
			}
			if( cut_value>0 && (best_cut==-1 || best_cut>cut_value) && num_nodes>1){ // update optimum
				new_purity = test_purity;
				best_cut = cut_value;
				if(debug_level>0)	cerr << " initial purity="<<new_purity<<" best_cut=" << best_cut << endl;
			}
		}

		if(new_purity>0.12){new_purity=0;} // if there are no impure nodes -> disable
		unsigned int num_impure_nodes=0;
		for (unsigned int i = 0; i < x_hat.size(); i++) { if(abs(x_hat[i]) < new_purity) { num_impure_nodes++; } }
		if(debug_level>0)	cerr << " num_impure_nodes" << num_impure_nodes << endl;
		if(num_impure_nodes<2 && new_purity > 0.001){new_purity=0; } // disable if there is only one node that is impure with a high purity threshold

		if(debug_level>0)	cerr << " [DEBUG]  detected a purity of "<< (new_purity==0 ? "<disabled>":to_string(new_purity)) << "\n";

		for (unsigned int i = 0; i < x_hat.size(); i++) {
			if(abs(x_hat[i]) <= new_purity) {
				groupZero.m_content_CC.push_back(cur_cc[i]);
				groupZero_set[cur_cc[i]] = true;
			}
			else if (x_hat[i] < 0) {
				groupA.m_content_CC.push_back(cur_cc[i]);
			}
			else { // x_hat[i] > 0
				groupB.m_content_CC.push_back(cur_cc[i]);
				groupB_set[cur_cc[i]] = true;
			}
		}

		if( (groupA.size() == 0 && groupB.size() == 0) ){

			cerr << " [WARNING]  All nodes are below the purity threshold. Continue without purity." << "\n";

			for (unsigned int i = 0; i < groupZero.size(); i++) {
				if (x_hat[ groupZero[i] ] < 0) {
					groupA.m_content_CC.push_back(groupZero[i]);
				}
				else { // x_hat[i] > 0
					groupB.m_content_CC.push_back(groupZero[i]);
				}
			}
		}
	}else{
		for (unsigned int i = 0; i < x_hat.size(); i++) {
			if(abs(x_hat[i]) < param_sep_purity) {
				groupZero.m_content_CC.push_back(cur_cc[i]);
				groupZero_set[cur_cc[i]] = true;
			}
			else if (x_hat[i] < 0) {
				groupA.m_content_CC.push_back(cur_cc[i]);
			}
			else { // x_hat[i] > 0
				groupB.m_content_CC.push_back(cur_cc[i]);
				groupB_set[cur_cc[i]] = true;
			}
		}
		if( (groupA.size() == 0 && groupB.size() == 0) ){
				
			cerr << " [WARNING]  All nodes are below the purity threshold. Continue purity/2." << "\n";

			for (unsigned int i = 0; i < x_hat.size(); i++) {
				if(abs(x_hat[i]) < param_sep_purity/2) {
					groupZero.m_content_CC.push_back(cur_cc[i]);
					groupZero_set[cur_cc[i]] = true;
				}
				else if (x_hat[i] < 0) {
					groupA.m_content_CC.push_back(cur_cc[i]);
				}
				else { // x_hat[i] > 0
					groupB.m_content_CC.push_back(cur_cc[i]);
					groupB_set[cur_cc[i]] = true;
				}
			}

			if( (groupA.size() == 0 && groupB.size() == 0) ){
					
				cerr << " [WARNING]  All nodes are below the purity threshold. Continue purity/10." << "\n";

				for (unsigned int i = 0; i < x_hat.size(); i++) {
					if(abs(x_hat[i]) < param_sep_purity/10) {
						groupZero.m_content_CC.push_back(cur_cc[i]);
						groupZero_set[cur_cc[i]] = true;
					}
					else if (x_hat[i] < 0) {
						groupA.m_content_CC.push_back(cur_cc[i]);
					}
					else { // x_hat[i] > 0
						groupB.m_content_CC.push_back(cur_cc[i]);
						groupB_set[cur_cc[i]] = true;
					}
				}

				if( (groupA.size() == 0 && groupB.size() == 0) ){
						
					cerr << " [WARNING]  All nodes are below the purity threshold. Continue purity/100." << "\n";

					for (unsigned int i = 0; i < x_hat.size(); i++) {
						if(abs(x_hat[i]) < param_sep_purity/100) {
							groupZero.m_content_CC.push_back(cur_cc[i]);
							groupZero_set[cur_cc[i]] = true;
						}
						else if (x_hat[i] < 0) {
							groupA.m_content_CC.push_back(cur_cc[i]);
						}
						else { // x_hat[i] > 0
							groupB.m_content_CC.push_back(cur_cc[i]);
							groupB_set[cur_cc[i]] = true;
						}
					}

					if( (groupA.size() == 0 && groupB.size() == 0) ){
							
						cerr << " [WARNING]  All nodes are below the purity threshold. Continue without purity." << "\n";

						for (unsigned int i = 0; i < groupZero.size(); i++) {
							if (x_hat[ groupZero[i] ] < 0) {
								groupA.m_content_CC.push_back(groupZero[i]);
							}
							else { // x_hat[i] > 0
								groupB.m_content_CC.push_back(groupZero[i]);
							}
						}
					}
				}
			}
		}
	}

	string id;
	string idA;
	string idB;
	string idZ;
	if(debug_level==15){
		unsigned int n = cur_cc.size();
		id=getCCid(cur_cc.m_content_CC);
		idA=getCCid(groupA.m_content_CC);
		idB=getCCid(groupB.m_content_CC);
		idZ=getCCid(groupZero.m_content_CC);
		map<unsigned int,unsigned int> mapping;
		for (unsigned int i = 0; i < (unsigned int)n; i++) {mapping[cur_cc[i]] = i;}
		{
			ofstream OFS(("cluster_n"+to_string(n)+"_id"+id+".edgelist").c_str());
			OFS << "source\ttarget\tweight"<<endl;

			map<pair<unsigned int,unsigned int>,unsigned int> knownedges;
			for(unsigned int i = 0 ; i < (unsigned int)n ; i++){
				unsigned int from = cur_cc[i]; 
				unsigned int sum = 0;
				if(!mapping.count(from)){continue;}
				for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
					unsigned int to = graph[from].edges[j].edge;
					if(!mapping.count(to)){continue;}
					pair<unsigned int,unsigned int> key_cur_edge;
					if(from < to) key_cur_edge = make_pair(from,to);
					else key_cur_edge = make_pair(to,from);
					if(!knownedges.count(key_cur_edge)){
						knownedges[key_cur_edge] = true;
						OFS << graph[from].full_name << "\t" << graph[to].full_name << "\t" << graph[from].edges[j].weight << endl;
					}
				}
			}
			OFS.close();
			string id=getCCid(cur_cc.m_content_CC);
		}
		{
			ofstream OFS(("cluster_n"+to_string(n)+"_id"+id+"_ci0_conn"+to_string(connectivity)+".nodeweight").c_str());
			OFS << "id\teigenvector"<<endl;
			for(unsigned int i = 0 ; i < (unsigned int)n ; i++){
				OFS << graph[cur_cc[i]].full_name << "\t" << x_hat[i] << endl;
			}
			OFS.close();
		}
		{
			ofstream OFS(("cluster_n"+to_string(n)+"_id"+id+"_ci0_conn"+to_string(connectivity)+".split").c_str());
			OFS << "id\tsplit"<<endl;
			for (unsigned int i = 0; i < groupZero.size(); i++) { OFS << graph[groupZero[i]].full_name << "\t0" << endl; }
			for (unsigned int i = 0; i < groupA.size(); i++) { OFS << graph[groupA[i]].full_name << "\t1" << endl; }
			for (unsigned int i = 0; i < groupB.size(); i++) { OFS << graph[groupB[i]].full_name << "\t-1" << endl; }
			OFS.close();
		}
	}

	if (debug_level > 0) cerr << getTime() << " [DEBUG] splitting CC with "<<cur_cc.size()<<" nodes "<<(debug_level==15 ? "@"+id+"" : "")<< " into ("<<groupA.size()<<""<<(debug_level==15 ? " @"+idA+"" : "")<< ","<<groupB.size()<<""<<(debug_level==15 ? " @"+idB+"" : "")<< ","<<groupZero.size()<<""<<(debug_level==15 ? " @"+idZ+"" : "")<< ") sized groups!"<< "\n";

	// Catch error in laplacien calcs
	// all nodes are in one partition -> restart with an increased restart level
	if ( (groupA.size() == 0 && groupB.size() == 0) || 
		 ( (groupA.size() == 0 || groupB.size() == 0) && groupZero.size() == 0) ){
		
		if(restarted==0){
			if(debug_level>0)
				cerr << getTime() << " [WARNING] connected component with "<<cur_cc.size()<<" nodes and density of "<<cur_cc.density<<": the fiedler vector associated to the algebraic connectivity of "<<connectivity<<" would result in invalid split groups ("<<groupA.size()<<","<<groupB.size()<<","<<groupZero.size()<<"). I will retry with double precision, if this persists an error will be thrown."<< endl;
			q->dispatch([cur_cc,q]{ partition_CC(cur_cc,q,true,1); });
		}else if(restarted==1){
			if(debug_level>0)
				cerr << getTime() << " [WARNING] connected component with "<<cur_cc.size()<<" nodes and density of "<<cur_cc.density<<": the fiedler vector associated to the algebraic connectivity of "<<connectivity<<" would result in invalid split groups ("<<groupA.size()<<","<<groupB.size()<<","<<groupZero.size()<<"). I will retry with double precision using the power iteration, if this persists an error will be thrown."<< endl;
			q->dispatch_allCores([cur_cc,q]{ partition_CC(cur_cc,q,false,2); });
		}else{
			cerr << getTime() << " [ERROR] connected component with "<<cur_cc.size()<<" nodes and density of "<<cur_cc.density<<": the fiedler vector associated to the algebraic connectivity of "<<connectivity<<" would result in invalid split groups ("<<groupA.size()<<","<<groupB.size()<<","<<groupZero.size()<<"). Failed even with double precision, the group will contain a star after the algebraic connectivity !"<< endl;
			print_group(cur_cc,connectivity,tid,true);
		}
		
	}else{

		// for the -core algorithm, test if the number of species now decreased by too much or keep on splitting
		if(param_core){
			map<unsigned int,bool>cc_speciesA;
			for (int i = 0; i < groupA.size(); ++i){ cc_speciesA[graph[groupA[i]].species_id]=true; }
			groupA.species_num=cc_speciesA.size();
			map<unsigned int,bool>cc_speciesB;
			for (int i = 0; i < groupB.size(); ++i){ cc_speciesB[graph[groupB[i]].species_id]=true; }
			groupB.species_num=cc_speciesB.size();
			map<unsigned int, bool> done; // for zero additionally test if this group is connected -> if not

			if(groupA.species_num!=cur_cc.species_num && // neither A, B or Zero meet the -core criteria -> test zero and if zero also do not hold -> print the original group cur_cc before the split
					groupB.species_num!=cur_cc.species_num && 
					param_coreMaxProteinsPerSpecies*cur_cc.species_num > cur_cc.size()){

				map<unsigned int,bool>cc_speciesZero;
				for (int i = 0; i < groupZero.size(); ++i){ cc_speciesZero[graph[groupZero[i]].species_id]=true; }
				groupZero.species_num=cc_speciesZero.size();

				if( groupZero.species_num!=cur_cc.species_num ){ // the original CC is still fine -> print
					print_group(cur_cc,connectivity,tid,false);
					return; // done with this group, no need to split into A and B by the definition of -core
				}
			}
		}

		// now add the edges between A and B (and Zero) to the output rmgraph
		// (the rmgraph contains all edges that are removed by the clustering)

		if(!graph_clean.count(tid))
			graph_clean[tid]=make_shared<ofstream>((param_rmgraph+to_string(tid)).c_str());

		if(debug_level>0) cerr << "[INFO] removing between A and B+Zero" << endl;
		// remove edges between A and (B or Zero)
		for(unsigned int i = 0 ; i < groupA.size() ; i++){
			protein node_i = graph[groupA[i]];
			for(unsigned int j = 0 ; j < graph[groupA[i]].edges.size() ; j++){
				protein node_j = graph[graph[groupA[i]].edges[j].edge];
				if(groupB_set.count(graph[groupA[i]].edges[j].edge) || 
					groupZero_set.count(graph[groupA[i]].edges[j].edge) ){
					(*graph_clean[tid]) << node_i.full_name << "\t" << species[node_i.species_id] << "\t" << node_j.full_name << "\t" << species[node_j.species_id] << "\n";
				}
			}
			// auto it = std::remove_if(
			// 	node_i.edges.begin(), 
			// 	node_i.edges.end(), 
			// 	[node_i,groupB_set,groupZero_set](wedge cur_el)->bool 
			// 	{
			// 		if(debug_level>0 && (groupB_set.count(cur_el.edge) || groupZero_set.count(cur_el.edge))) cerr << "[INFO] removing "<< node_i.full_name << "-"<< graph[cur_el.edge].full_name << " weight=" << cur_el.weight << endl;
			// 		return groupB_set.count(cur_el.edge) || groupZero_set.count(cur_el.edge);
			// 	}
			// );
			// node_i.edges.erase(it, node_i.edges.end());
		}

		if(debug_level>0) cerr << "[INFO] removing between B and Zero" << endl;
		// remove edges between B and Zero
		for(unsigned int i = 0 ; i < groupB.size() ; i++){
			protein node_i = graph[groupB[i]];
			for(unsigned int j = 0 ; j < graph[groupB[i]].edges.size() ; j++){
				protein node_j = graph[graph[groupB[i]].edges[j].edge];
				if(groupZero_set.count(graph[groupB[i]].edges[j].edge)){
					(*graph_clean[tid]) << node_i.full_name << "\t" << species[node_i.species_id] << "\t" << node_j.full_name << "\t" << species[node_j.species_id] << "\n";
				}
			}
			// auto it = std::remove_if(
			// 	node_i.edges.begin(), 
			// 	node_i.edges.end(), 
			// 	[node_i,groupZero_set](wedge cur_el)->bool 
			// 	{
			// 		if(debug_level>0 && groupZero_set.count(cur_el.edge)) cerr << "[INFO] removing "<< node_i.full_name << "-"<< graph[cur_el.edge].full_name << " weight=" << cur_el.weight << endl;
			// 		return groupZero_set.count(cur_el.edge);
			// 	}
			// );
			// node_i.edges.erase(it, node_i.edges.end());
		}

		// finally add the new clusters back to the stack of CCs
		// for -core only add those clusters which meet the criteria of number of species

		if(param_core){
			if((groupA.species_num >= param_coreMinSpecies || param_coreMaxProteinsPerSpecies*cur_cc.species_num > cur_cc.size()) && groupA.size()>1)
				q->dispatch([groupA,q]{ find_CCs_givenNodes(q,groupA.m_content_CC); });
			if((groupB.species_num >= param_coreMinSpecies || param_coreMaxProteinsPerSpecies*cur_cc.species_num > cur_cc.size()) && groupB.size()>1)
				q->dispatch([groupB,q]{ find_CCs_givenNodes(q,groupB.m_content_CC); });
			if((groupZero.species_num >= param_coreMinSpecies || param_coreMaxProteinsPerSpecies*cur_cc.species_num > cur_cc.size()) && groupZero.size()>1) // only if zero is connected 
				q->dispatch([groupZero,q]{ find_CCs_givenNodes(q,groupZero.m_content_CC); });
		}else{
			if(groupA.size()>1)
				q->dispatch([groupA,q]{ find_CCs_givenNodes(q,groupA.m_content_CC); });
			if(groupB.size()>1)
				q->dispatch([groupB,q]{ find_CCs_givenNodes(q,groupB.m_content_CC); });
			if(groupZero.size()>1)
				q->dispatch([groupZero,q]{ find_CCs_givenNodes(q,groupZero.m_content_CC); });
		}
	}
}

void find_CCs(dispatch_queue *q){

	map<unsigned int, bool> done;	// Keep track on what was done (for each node)
	bool allNodesAreDone = false;

	vector<ConnectedComponent> CC; // vector of all connected components found
	//unsigned int min_i=0;

	while( true ){ // CC.size() < num_cpus / gather up to num_cpus connected components

		allNodesAreDone = true;

		for (unsigned int protein_id = 0 ; protein_id < graph.size() ; protein_id++) {
			
			if ( graph[protein_id].is_done ){continue;}

			if (done.count(protein_id) && done[protein_id]){continue;}// We were here already

			if(debug_level>0) cerr << getTime() << " [DEBUG:find_CCs] find_CCs:start at "<< protein_id << "<" << graph.size() << endl;

			//min_i=protein_id;

			done[protein_id]=true; // mark this node
			ConnectedComponent cur_cc = BFS(&done,protein_id,0); // get the CC of the current node (protein_id) 

			// Do not report singles
			if (cur_cc.size() < 2) {continue;} // singletons are from no interest

			// for -core, skip CC with less than param_coreMinSpecies species
			if( param_core && cur_cc.species_num==0 ){
				if(cur_cc.size()<param_coreMinSpecies){allNodesAreDone=false;break;} // if there are less nodes than wanted species -> skip
				map<unsigned int,bool>cc_species;
				for (int i = 0; i < cur_cc.size(); ++i){ cc_species[graph[cur_cc[i]].species_id]=true; }
				cur_cc.species_num=cc_species.size();
				if(cc_species.size()<param_coreMinSpecies){allNodesAreDone=false;break;} // less species as wanted -> skip
			}

			// Skip those that are too large (try heuristic)
			if (cur_cc.size() > param_max_nodes) {

				// reset done vector
				//for (int i = 0; i < cur_cc.size(); ++i){ done[cur_cc[i]]=false; }
				if(debug_level>0) cerr << " [WARNING]  Found a very large connected component that contains " << cur_cc.size() << ">" << param_max_nodes << " (maxnodes) elements. This behavior can be adjusted using -maxnodes. Now using a slow heuristic: try to identify and remove edges." << "\n";
				q->dispatch([cur_cc,q]{ removeLowQualityEdges(cur_cc,q); });
				//protein_id--; 
				continue;
			}

			if (debug_level > 0) cerr << getTime() << " [DEBUG:find_CCs] Found connected component: " << cur_cc.size() << " proteins (ID: " << protein_id << "), ini from " << graph[protein_id].full_name<< "\n";
			cur_cc.calc_dsum();
			cur_cc.calc_density();
			if(cur_cc.density > 1){ cerr << "[WARNING] : The input graph has duplicated edges, this lead to an invalid graph density of " << cur_cc.density << " (should be <1). Please clean the .blast-graph with 'proteinortho.pl --cleanblast --step=3 --project=...' or use the cleanupblastgraph tool in src/ to remove the duplicated edges." << "\n"; throw; }
			if (debug_level > 0) cerr << getTime() << " [DEBUG:find_CCs] did calc dsum and density: " << cur_cc.size() << " proteins (ID: " << protein_id << "), graph density="<< cur_cc.density << ", sum of degrees="<< cur_cc.d_sum << " ini from " << graph[protein_id].full_name<< "\n";

			q->dispatch([cur_cc,q]{ partition_CC(cur_cc,q,true,false); });
			allNodesAreDone=false;
			//break;
		}
		if(allNodesAreDone)break; // no additional CC can be found -> done
	}
}

///////////////////////////////////////////////////////////
// Major partioning algorithm
///////////////////////////////////////////////////////////

void partition_graph() {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	dispatch_queue q("proteinortho_clustering_dq", num_cpus);

	q.dispatch([&q]{ find_CCs(&q); });
	q.start();

	if(debug_level>0) cerr << getTime() << " [DEBUG] done with ini dispatch_queue" << endl;

	q.waitTilDone();

	if(debug_level>0) cerr << getTime() << " [DEBUG] quit dispatch_queue" << endl;

	return;
}

///////////////////////////////////////////////////////////
// File parser
///////////////////////////////////////////////////////////
void parse_file(string file) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	// unsigned long graph_ram_total=0;

	if (param_verbose) cerr << "Reading " << file << "\n";
	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {
		// For each line
		string file_a = "";	unsigned int file_a_id = 0;
		string file_b = "";	unsigned int file_b_id = 0;

		float avg_bitscore_median=-1; //for normalization the bitscore is devided by the median bitscore

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

				if(fields.size() == 6){ // either the median scores are directly after the species header OR see <else if>
					avg_bitscore_median=(string2float(fields[3])+string2float(fields[5]))/2;
					if(avg_bitscore_median<1){avg_bitscore_median=-1;}
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
			}else if ( (fields.size() == 4) && fields[0].substr(0, 1) == "#") { // OR the median scores are in an additional line below the species header

				avg_bitscore_median=(string2float(fields[1])+string2float(fields[3]))/2;
				if(avg_bitscore_median<1){avg_bitscore_median=-1;}

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

				//if(debug_level>0) cerr << "parse_file:(" << fields[0]<<","<<fields[1] << ") -> " << (protein2id.find(fields[0]) == protein2id.end()) << endl;

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
						cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << "\n";throw;
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
						cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << "\n";throw;
					}
					graph.push_back(b);
					// graph_ram_total+=sizeof(protein);
				}

				// Bitscores 
				// check range (in float)

				float bit_a,bit_b; 
				if(avg_bitscore_median<0){
					bit_a = string2float(fields[3]);//255 is exactly in the middle of ushort range, such that the score can be 255 fold upregulated or down regulated.
					bit_b = string2float(fields[5]);
				}else{	
					bit_a = 255.0* (string2float(fields[3])/avg_bitscore_median);//255 is exactly in the middle of ushort range, such that the score can be 255 fold upregulated or down regulated.
					bit_b = 255.0* (string2float(fields[5])/avg_bitscore_median);
				}

				if(bit_a<1){bit_a=1;}
				if(bit_b<1){bit_b=1;}

				if(bit_a>USHRT_MAX){
					cerr << " [WARNING] unsigned short overflow " << bit_a <<  ">USHRT_MAX (bitscore of "<< ida<< " adj. to "<< idb<< ") using "<< USHRT_MAX<< " (USHRT_MAX) instead." << "\n";
					bit_a=(float)USHRT_MAX;
				}
				if(bit_b>USHRT_MAX){
					cerr << " [WARNING] unsigned short overflow " << bit_b <<  ">USHRT_MAX (bitscore of "<< idb<< " adj. to "<< ida<< ") using "<< USHRT_MAX<< " (USHRT_MAX) instead." << "\n";
					bit_b=(float)USHRT_MAX;
				}

				// assign
				unsigned short bitscore_avg = (bit_a + bit_b)/2;
				if(bitscore_avg<1){bitscore_avg=1;}
		
				// Add link to graph (reciprocal)					
				unsigned int a_id = protein2id[fields[0]];
				unsigned int b_id = protein2id[fields[1]];

				if(!param_useWeights){bitscore_avg=1;}

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
		cerr << string("Could not open file " + file).c_str() << "\n";throw;
	}

	// graph_ram_total_inKB += graph_ram_total/1e+3;
	// if (debug_level > 0) cerr << getTime() << " [DEBUG]  Expected Memory Usage of the current input graph: " << graph_ram_total_inKB << " KB = "  << graph_ram_total_inKB/1e+3 << " MB. (current input graph are all the currently loaded files)" << "\n";

	if(species_counter==0){species.push_back("0");species_counter++;}

	#ifdef timeAnalysis
		if(!t_master.count("parse_file")) t_master["parse_file"]=0;
		t_master["parse_file"] += (float)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

void parse_abc_file(string file) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif


	species.push_back("input");
	species2id["input"] = species_counter++;

	if (param_verbose) cerr << "Reading " << file << "\n";
	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {
		while (!graph_file.eof()) {
			getline(graph_file, line);
			vector<string> fields;
			tokenize(line, fields, "\t");
			// Header line
			if ((fields.size() == 3) && fields[0].substr(0, 1) != "#") {
				// a b e1 b1 e2 b2 score

				// 5.16 deal with duplicated IDs by adding file ID to protein ID
				string ida = fields[0];
				string idb = fields[1];
				//fields[0] += " "; fields[0] += to_string(file_a_id);
				//fields[1] += " "; fields[1] += to_string(file_b_id);

				// 5.16 do not point to yourself
				if (!fields[0].compare(fields[1])) {continue;}

				//if(debug_level>0) cerr << "parse_file:(" << fields[0]<<","<<fields[1] << ") -> " << (protein2id.find(fields[0]) == protein2id.end()) << endl;

				// A new protein
				if (protein2id.find(fields[0]) == protein2id.end())	{
					protein a;
					a.full_name	= ida;
					// graph_ram_total+=a.full_name.size()*sizeof(string);
					a.species_id = 0;
					// graph_ram_total+=sizeof(unsigned int);
					protein2id[fields[0]] = protein_counter++;
					// graph_ram_total+=fields[0].size()*sizeof(string)+sizeof(unsigned int);

					if( graph.size() >= 1073741824-1 ){
						cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << "\n";throw;
					}
					graph.push_back(a);
					// graph_ram_total+=sizeof(protein);
				}
				if (protein2id.find(fields[1]) == protein2id.end())	{
					protein b;
					b.full_name	= idb;
					// graph_ram_total+=b.full_name.size()*sizeof(string);
					b.species_id	= 0;
					// graph_ram_total+=sizeof(unsigned int);
					protein2id[fields[1]] = protein_counter++;
					// graph_ram_total+=fields[1].size()*sizeof(string)+sizeof(unsigned int);

					if( graph.size() >= 1073741824-1 ){
						cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << "\n";throw;
					}
					graph.push_back(b);
					// graph_ram_total+=sizeof(protein);
				}

				// Bitscores 

				float bit_a = string2float(fields[3]);

				if(bit_a<1){bit_a=1;}

				if(bit_a>USHRT_MAX){
					cerr << " [WARNING] unsigned short overflow " << bit_a <<  ">USHRT_MAX (bitscore of "<< ida<< " adj. to "<< idb<< ") using "<< USHRT_MAX<< " (USHRT_MAX) instead." << "\n";
					bit_a=(float)USHRT_MAX;
				}

				// Add link to graph (reciprocal)					
				unsigned int a_id = protein2id[fields[0]];
				unsigned int b_id = protein2id[fields[1]];

				// 5.17, add weight
				wedge w;
				w.edge=b_id;
				// graph_ram_total+=sizeof(unsigned int);
				w.weight=bit_a;
				// graph_ram_total+=sizeof(unsigned int);
				graph[a_id].edges.push_back(w);
				// graph_ram_total+=sizeof(wedge);
				w.edge=a_id;
				// graph_ram_total+=sizeof(unsigned int);
				w.weight=bit_a;
				// graph_ram_total+=sizeof(unsigned int);
				graph[b_id].edges.push_back(w);
				// graph_ram_total+=sizeof(wedge);
				edges++;
			}
		}
		graph_file.close();
	}
	else {
		cerr << string("Could not open file " + file).c_str() << "\n";throw;
	}

	// graph_ram_total_inKB += graph_ram_total/1e+3;
	// if (debug_level > 0) cerr << getTime() << " [DEBUG]  Expected Memory Usage of the current input graph: " << graph_ram_total_inKB << " KB = "  << graph_ram_total_inKB/1e+3 << " MB. (current input graph are all the currently loaded files)" << "\n";

	if(species_counter==0){species.push_back("0");species_counter++;}

	#ifdef timeAnalysis
		if(!t_master.count("parse_file")) t_master["parse_file"]=0;
		t_master["parse_file"] += (float)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
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
		t_master["sort_species"] += (float)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

// Progress stats
void stats(unsigned int lapack,unsigned int power,bool active_status) {
	if (!param_verbose) return;
	//float stat = float(i/size*100);
	if (last_stat_lapack != lapack || last_stat_power != power || last_stat_act != active_status) {
		last_stat_lapack = lapack; 
		last_stat_power = power;
		last_stat_act = active_status;
		if(debug_level==0) cerr << '\r';
		else cerr << endl;
		cerr << "Clustering: working on (" << lapack << (lapack==1 ? "@BFS/lapack" : "@lapack")<<(active_status?"*":"")<<" + " << power << "@power"<<(!active_status?"*":"")<<") connected component"<< (lapack+power>1?"s":"") <<"          " << std::flush ;
		if(lapack==0 && power==0){ cerr << '\r'<<"Clustering: done                                                       \n"; }
	}
}

// Header with species names
void print_header() {
	cout << "# Species\tGenes\tAlg.-Conn.";
	for (unsigned int i = 0; i < species_counter; i++) {
		cout << "\t" << species[reorder_table[i]];
	}
	cout << "\n";
}

struct protein_degree{
	inline bool operator() (const pair<string,int> & p1, const pair<string,int>& p2){
		return (p1.second > p2.second);
	}
};

// Group formatting
void print_group(ConnectedComponent& cur_cc, float connectivity, size_t tid, bool failed) {

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	for (unsigned int i = 0; i < cur_cc.size(); i++) { graph[cur_cc[i]].is_done=true; }

	if(debug_level>0) cerr << getTime()<< " [DEBUG] print_group start"<< endl;

	if(cur_cc.size()<2){return;}

	map<unsigned int, bool> done;	// Keep track on what was done (for each node)
	for (unsigned int i = 0 ; i < cur_cc.size() ; i++) {
		unsigned int from = cur_cc[i]; 
		done[ from ] = 0;
	}
	for (unsigned int i = 0 ; i < cur_cc.size() ; i++) {
		unsigned int from = cur_cc[i]; 
		for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){
			if(!done.count(graph[from].edges[j].edge))
				done[ graph[from].edges[j].edge ]=1;
		}
	}

	if(!proteinorthotmp_clean.count(tid)){
		proteinorthotmp_clean[tid]=make_shared<ofstream>((param_rmgraph+"_proteinortho_tmp_"+to_string(tid)).c_str());
	}
	//ofstream os("out_tmp_"+to_string(tid),std::ios_base::app);

	vector<vector<pair<string,int> > > line(species_counter);	// Output vector

	unsigned int species_number = 0;
	// For each protein in group 	 	
	for (unsigned int i = 0; i < cur_cc.size(); i++) {
		unsigned int current_protein = cur_cc[i];
		unsigned int current_species = graph[current_protein].species_id;
		if (line[current_species].size() == 0) 
			species_number++;
		line[current_species].push_back(make_pair( graph[current_protein].full_name , graph[current_protein].edges.size() ));
	}

	(*proteinorthotmp_clean[tid]) << species_number << "\t" << cur_cc.size() << "\t" << (debug_level>0 ? "@"+getCCid(cur_cc.m_content_CC)+":":"") << setprecision (3) << (connectivity < 0 ? -connectivity : connectivity) << (failed ? "*": "");
	// cerr << species_number << "\t" << cur_cc.size() << "\t" << setprecision (3) << connectivity;

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
		(*proteinorthotmp_clean[tid]) << "\t" << return_line;
		// cerr << "\t" << return_line;
	}

	(*proteinorthotmp_clean[tid]) << "\n";
	// cerr << "\n";

	if(debug_level>0) cerr << getTime()<< " [DEBUG] print_group @"<<getCCid(cur_cc.m_content_CC)<<" ERROR="<< ((int)cur_cc.size() - (int)cur_cc.size()) << endl;

	//os.close();

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::print_group")) t_master["partition_graph::print_group"]=0;
		t_master["partition_graph::print_group"] += (float)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

// Calculate number of species formatting
float calc_group(vector<unsigned int>* nodes) {

	map<unsigned int, bool> speciesids;

	for (unsigned int i = 0; i < nodes->size(); i++) {
		unsigned int current_protein = (*nodes)[i];
		speciesids[graph[current_protein].species_id]=1;
	}

	return speciesids.size()==0 ? 99999 : ((float)nodes->size())/((float)speciesids.size()); // 99999 such that if the species information is missing, then the criterion always fails and the splits are only made based on the alg. connectivity
	/* SEG FAIL:
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
	
	float avg = (float)sum/(float)species_number;
	return avg;*/
}

///////////////////////////////////////////////////////////
// Misc functions
///////////////////////////////////////////////////////////
// Convert string to float
float string2float(string str) {
	istringstream buffer(str);
	float value;
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
		t_master["tokenize"] += (float)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
}

///////////////////////////////////////////////////////////
// Algebraic connectivity functions
///////////////////////////////////////////////////////////
// Return maximum degree of given protein_ids -- openMP A ML
unsigned int max_of_diag(vector<unsigned int>& nodes, vector<unsigned int>& diag) {

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] max_of_diag start" << endl;

	unsigned int max = 0;

	bool useOpenMpFlag = (nodes.size() > param_minOpenmp);

	#pragma omp parallel for reduction(max: max) if (useOpenMpFlag)
	for (unsigned int i = 0; i < nodes.size(); i++) {
		if (diag[i] > max) max = diag[i] ;
	}

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] max_of_diag end" << endl;
	return max;
}

// Generate random vector x of size size
template<class T> vector<T> generate_random_vector(const unsigned int size) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] generate_random_vector start" << endl;

	vector<T> x(size);

	x[0] = (T)(rand() % 999+1)/1000;	// 0 to 1
	for (unsigned int i = 1; i < size; i++) {
		x[i] = (T)(rand() % 999+1)/1000;	// 0 to 1
		if (x[i] == x[i-1]) x[i] /= 3;		// Check: at least one value must be different from the others but still within 0 and 1
	}

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] generate_random_vector end" << endl;

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::generate_random_vector")) t_master["partition_graph::convergence::generate_random_vector"]=0;
		t_master["partition_graph::convergence::generate_random_vector"] += (T)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x;
}

// determine new X, Formula (1) -- openMP B ML
template<class T> vector<T> get_new_x(vector<T>*x, vector<unsigned int>& nodes, map<unsigned int,unsigned int> &mapping, bool isWeighted) {
	
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] get_new_x start" << endl; 

	vector<T> x_new(x->size(),0);

	bool useOpenMpFlag = (nodes.size() > param_minOpenmp);
	if(debug_level==9) cerr << getTime() << " [DEBUG:9] get_new_x done ini " << useOpenMpFlag << endl;

	if(isWeighted){

		#pragma omp parallel for if (useOpenMpFlag)
		for (unsigned int i = 0; i < nodes.size(); i++) {
			// go through adjacency list of node 
			for (unsigned int j = 0; j < graph[nodes[i]].edges.size(); j++) {
				// y points to z, so take entry z from x
				unsigned int abs_target = graph[nodes[i]].edges[j].edge;
				if(!mapping.count(abs_target)){continue;}
				unsigned int rel_target = mapping[abs_target];

				x_new[i] += (*x)[rel_target]*(T)graph[nodes[i]].edges[j].weight;
			}
		}

	}else{

		#pragma omp parallel for if (useOpenMpFlag)
		for (unsigned int i = 0; i < nodes.size(); i++) {

			// go through adjacency list of node 
			for (unsigned int j = 0; j < graph[nodes[i]].edges.size(); j++) {
				// y points to z, so take entry z from x
				unsigned int abs_target = graph[nodes[i]].edges[j].edge;
				if(!mapping.count(abs_target)){continue;}
				unsigned int rel_target = mapping[abs_target];

				x_new[i] += (*x)[rel_target];
			}
		}
	}

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] get_new_x end" << endl;

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::get_new_x")) t_master["partition_graph::convergence::get_new_x"]=0;
		t_master["partition_graph::convergence::get_new_x"] += (T)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	return x_new;
}

// Make vector x orthogonal to 1, Formula (2) -- openMP A ML
template<class T> vector<T> makeOrthogonal(vector<T> x) {
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] makeOrthogonal start" << endl;

	T sum = 0;

	bool useOpenMpFlag = (x.size() > param_minOpenmp);

	#pragma omp parallel for reduction(+: sum) if (useOpenMpFlag)
	for (unsigned int i = 0; i < x.size(); i++) {sum += x[i];}

	T average = sum/x.size();

	#pragma omp parallel for if (useOpenMpFlag)
	for (unsigned int i = 0; i < x.size(); i++) {x[i] -= average;}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::makeOrthogonal")) t_master["partition_graph::convergence::makeOrthogonal"]=0;
		t_master["partition_graph::convergence::makeOrthogonal"] += (T)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] makeOrthogonal end" << endl;
	return x;
}

// Normalize vector x, Formula (4) -- openMP A ML
template<class T> void normalize(vector<T> *x, T *length) {

	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] normalize start" << endl;
	T sum = 0;

	bool useOpenMpFlag = (x->size() > param_minOpenmp);

	#pragma omp parallel for reduction(+: sum) if (useOpenMpFlag)
	for (unsigned int i = 0; i < x->size(); i++) {sum += (*x)[i]*(*x)[i];}

	*length = (T)sqrt(sum);
	if (*length == 0) {*length = 0.000000001;if(debug_level>0) cerr << "normalize" << "\n";}// ATTENTION not 0!

	#pragma omp parallel for if (useOpenMpFlag)
	for (unsigned int i = 0; i < x->size(); i++) {(*x)[i] /= *length;}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::normalize")) t_master["partition_graph::convergence::normalize"]=0;
		t_master["partition_graph::convergence::normalize"] += (T)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] normalize end" << endl;
	//return x;
}

// Qx, Formula (5) -- openMP A ML
template<class T> vector<T> getY(T max_degree, vector<T> x_hat, vector<T>* x_new, vector<unsigned int>& nodes, vector<unsigned int>& diag){

	if(debug_level==9) cerr << getTime() << " [DEBUG:9] getY start" << endl;
	#ifdef timeAnalysis
		auto t_tmp = std::chrono::steady_clock::now( );
	#endif

	bool useOpenMpFlag = (nodes.size() > param_minOpenmp);

	#pragma omp parallel for if (useOpenMpFlag)
	for (unsigned int i = 0; i < x_hat.size(); i++) {
		x_hat[i] *= ((T)2*max_degree - (T)diag[i]);
		x_hat[i] += (*x_new)[i];
	}

	#ifdef timeAnalysis
		if(!t_master.count("partition_graph::convergence::getY")) t_master["partition_graph::convergence::getY"]=0;
		t_master["partition_graph::convergence::getY"] += (T)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
	#endif
	if(debug_level==9) cerr << getTime() << " [DEBUG:9] getY end" << endl;

	return x_hat;
}

double getConnectivity_double(vector<unsigned int> *nodes, bool useLapack, vector<double> *x_hat) {

	bool useWeights = (param_useWeights && nodes->size() <= param_max_nodes_weight); //param_useWeights = user input whether weights should be used or not. useWeights = the true value, that is true if param_useWeights is true and the maximum number of nodes are not exeeded for the weighted algorithm (param_max_nodes_weight)

	unsigned int n=nodes->size();

	if( n>1073741824 ){
		cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << "\n";throw;
	}

	double maxWeight=-1;
	map<unsigned int,unsigned int> mapping;
	for (unsigned int i = 0; i < (unsigned int)n; i++) {mapping[(*nodes)[i]] = i;}

	double connectivity = -1;
	*x_hat = vector<double>(n);

	if( ( n < 32768 && useLapack && param_useLapack == 1 ) || param_useLapack == 2 ){
		if (debug_level > 0) cerr << getTime() << " [DEBUG]@double using LAPACK" << endl;

		// maximal number of nodes 32768 (2^15 = SHRT_MAX) -> laplace matrix 2^15*2^15=2^30 (max vector size) entries 
		// max vector size = std::vector<int> myvector; cout << myvector.max_size() << "\n"; -> 2^30
		// used ram in MB of lapack = (unsigned int)n*(unsigned int)n*sizeof(double))/1e+6 

		double * laplacian = (double*)calloc( (unsigned int)n*(unsigned int)n,sizeof(double) );

		bool fill_laplacian_return=1; // return value of the fill algorithm -> true : all fine, false : ERROR
		// fill laplacian
		for(unsigned int i = 0 ; i < (unsigned int)n ; i++){

			unsigned int from = (*nodes)[i]; 
			if(!mapping.count(from)){continue;}
			unsigned int sum = 0;

			for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){

				unsigned int to = graph[from].edges[j].edge;

				if(!mapping.count(to)){continue;}
				unsigned int vector_idx = mapping[from] + mapping[to]*n; //transform the 2-d coordinate (i,j) of the nxn matrix to 1-d vector coordinate i+j*n of the 2n vector

				if( vector_idx >= (unsigned int)n*(unsigned int)n){
					fill_laplacian_return = false;
					break;
				}

				if( useWeights
					// && (unsigned int)n <= param_max_nodes_weight
					){
					double w = graph[from].edges[j].weight;
					sum+=w;
					laplacian[vector_idx]=-w;
					if(maxWeight<w)maxWeight=w;
				}else{
					sum++;
					laplacian[vector_idx]=-1.0;
				}
			}

			laplacian[mapping[from]+mapping[from]*n]=sum;
		}
		if(!fill_laplacian_return){
			cerr <<"CRITICAL ERROR : fill_laplacian : out of range" << "\n";
			throw;
		}

		// local variables:
			int il, iu, m = 1, lda = n, ldz = n, info, lwork, liwork, iwkopt;
			double vl, vu;
			double wkopt;
			double* work;
			int* iwork;
			int isuppz[(unsigned int)(2*m)];
			char Vchar='V', Ichar='I', Uchar='U'; // Ichar = for specific range of eigenvalues/vectors
			double eigenvalues[(unsigned int)n]; // need only 1 eigenvalue
			double * eigenvectors = (double*)malloc( (unsigned int)ldz*(unsigned int)m*sizeof(double) ); 
			il = 2; //that is the second one (il=1 -> the first one, il=2 the second one)
			iu = 2; 
			double eps=param_epsilon;

		// Determine optimal workspace 
			lwork = -1;
			liwork = -1;
			int n_int=(int)n;
			
			dsyevr_( &Vchar, &Ichar, &Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu, &eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork, &info );
			
			//dssyevr_<double>( &Vchar, &Ichar, &Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu, &eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork, &info );
			
			lwork = (int)wkopt;
			work = (double*)malloc( lwork*sizeof(double) );
			liwork = iwkopt;
			iwork = (int*)malloc( liwork*sizeof(int) );

		// Solve eigenproblem ...
			
			dsyevr_( &Vchar, &Ichar,&Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu,&eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, work, &lwork, iwork, &liwork , &info );
		
			//dssyevr_<double>( &Vchar, &Ichar,&Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu,&eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, work, &lwork, iwork, &liwork , &info );
			
		// Check for errors in convergence
		if( info > 0 ) {
			cerr << " [ERROR] The algorithm (d|s)syevr failed to compute eigenvalues. Continue now with the slow standard approach (power iteration)." << "\n";
			throw;
			// goto standardComputationOfAlgCon;
		}

		// deallocate
		delete [] laplacian;
		delete [] work;
		delete [] iwork;

		// calculate normalized algebraic connectivity and fill x_hat
		if(useWeights){
			connectivity = eigenvalues[0]/(maxWeight*(double)n);
		}else{
			connectivity = eigenvalues[0]/((double)n);
		}
		for(unsigned int i = 0 ; i < (unsigned int)n ; i++)
			(*x_hat)[i]=-eigenvectors[i];

		// deallocate
		delete [] eigenvectors;	

	}else{
		if (debug_level > 0) cerr << getTime() << " [DEBUG]@double using POWER" << endl;

		if(debug_level==9) cerr << getTime() << " [DEBUG:9] start" << endl;

		if(param_useWeights && !useWeights){
			cerr << " [INFO] The maximum number of nodes for the weighted algorithm is exeeded. Continue now with the faster unweighted algorithm." << "\n";
		}
		//diagonal matrix diag : d(u,u)=number of adjacency nodes=deg(u)
		vector<unsigned int> diag(n);

		if(useWeights){
			for (unsigned int i = 0; i < (unsigned int)n; i++) {
				diag[i]=0;
				unsigned int from = (*nodes)[i]; 
				if(!mapping.count(from)){continue;}

				for (unsigned int j = 0; j < graph[(*nodes)[i]].edges.size(); j++) {

					unsigned int to = graph[from].edges[j].edge;
					if(!mapping.count(to)){continue;}
					diag[i] += graph[(*nodes)[i]].edges[j].weight;
					if(useWeights && maxWeight<graph[(*nodes)[i]].edges[j].weight)maxWeight=graph[(*nodes)[i]].edges[j].weight;
				}
			}
		}
		// else{
		// 	for (unsigned int i = 0; i < (unsigned int)n; i++) {
		// 		diag[i]=graph[(*nodes)[i]].edges.size();
		// 		if(useWeights)
		// 			for (unsigned int j = 0; j < graph[(*nodes)[i]].edges.size(); j++) {
		// 				if(maxWeight<graph[(*nodes)[i]].edges[j].weight)maxWeight=graph[(*nodes)[i]].edges[j].weight;
		// 			}
		// 	}
		// }

		if(debug_level==9) cerr << getTime() << " [DEBUG:9] ini done" << endl;

		// Get max degree / sum of weights of nodes
		unsigned int max_d = max_of_diag((*nodes),diag);	

		// Init randomized variables. 
		vector<double> x = generate_random_vector<double>(n);		

		// Orthogonalize + normalize vector + get initial lenght
		double current_length = 0;
		double last_length;

		(*x_hat) = makeOrthogonal<double>(x);
		normalize(x_hat, &last_length);

		// Repeat until difference < param_epsilon
		unsigned int iter = 0;	// catch huge clustering issues by keeping track here

		#ifdef timeAnalysis
			auto t_tmp = std::chrono::steady_clock::now();
		#endif	

		while(++iter < param_max_iter) { 

			if(debug_level==9) cerr << getTime() << " [DEBUG:9] iter:" << iter << endl;

			last_length = current_length;

			// Get a new x
			x = get_new_x<double>(x_hat, (*nodes), mapping, useWeights);

			// Get y
			vector<double> y = getY<double>(max_d,*x_hat,&x,(*nodes),diag);

			// Orthogonalize
			(*x_hat) = makeOrthogonal<double>(y);

			// Get length (lambda) & normalize vector
			normalize(x_hat, &current_length);

			if ( abs(current_length-last_length) < param_epsilon && iter >= min_iter ) break;	// prevent convergence by chance, converge to param_epsilon
		}
		if(debug_level==9) cerr << getTime() << " [DEBUG:9] post while" << endl;

		#ifdef timeAnalysis
			if(!t_master.count("partition_graph::convergence")) t_master["partition_graph::convergence"]=0;
			t_master["partition_graph::convergence"] += (double)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
		#endif

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   " << iter << " / " << param_max_iter << " iterations required (error is " << abs(current_length-last_length) << ")" << "\n";

		#ifdef DEBUG
			total_number_of_iterations_convergence+=iter;
		#endif

		if(useWeights){
			connectivity = (-current_length+(double)2*max_d)/(maxWeight*(double)n);
		}else{
			connectivity = (-current_length+(double)2*max_d)/((double)n);
		}
		// normalize(x_hat, &current_length);
		//if(debug_level==9) cerr << getTime() << " [DEBUG:9] post last norm" << endl;

		// 5.17 catch hardly-converging groups
		if (iter >= param_max_iter) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Connectivity score of connected component with " << n << " elements did not converge perfectly in time." << "\n";
		}
	}

	if (debug_level > 1){
		cerr << getTime() << " [DEBUG]   Connectivity score:" << connectivity;
		if ( (debug_level > 1 && (unsigned int)n<100 ) || debug_level > 2 ){
			cerr << " eigenvector: (";
			for(unsigned int i = 0 ; i < (unsigned int)n ; i++)
				cerr << (*x_hat)[i] << ","; 
			cerr << ")";
		}
		cerr << "\n";
	}

	// Split groups if connectivity is too low, remove tree like structures that might have arosen
	if (connectivity < param_con_threshold) {
		
		// 5.17 new threshold option overwrites connectivity
		if (param_min_species >=1 ) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]  Start the calculation of the average gene/species score " << "\n";
			double avg = calc_group(nodes);
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Found " << avg << " genes/species on average. User asked for at least " << param_min_species << "\n";
			if (avg <= param_min_species) {
				if (debug_level > 0) cerr << getTime() << " [DEBUG]   Group is going to be accepted despite connectivity " << connectivity << "\n";
				// just to be safe from infinit loops due to rounding
				if (connectivity == 0) connectivity = 0.001;
				// no splitting despite bad connectivity
				return -connectivity;			
			}
		}
	}
	
	return connectivity;
}
float getConnectivity_float(vector<unsigned int> *nodes, bool useLapack, vector<float> *x_hat){

	if (debug_level > 0) cerr << getTime() << " getConnectivity_float"<< endl;

	bool useWeights = (param_useWeights && nodes->size() <= param_max_nodes_weight); //param_useWeights = user input whether weights should be used or not. useWeights = the true value, that is true if param_useWeights is true and the maximum number of nodes are not exeeded for the weighted algorithm (param_max_nodes_weight)

	unsigned int n=nodes->size();

	if( n>1073741824 ){
		cerr << string("[CRITICAL ERROR]   Overflow: number of nodes overflow the maximum number for vector allocation (2^30=1073741824).").c_str() << "\n";throw;
	}

	float maxWeight=-1;
	map<unsigned int,unsigned int> mapping;
	for (unsigned int i = 0; i < (unsigned int)n; i++) {mapping[(*nodes)[i]] = i;}

	float connectivity = -1;
	*x_hat = vector<float>(n);

	if( ( n < 32768 && useLapack && param_useLapack == 1 ) || param_useLapack == 2 ){

		if (debug_level > 0) cerr << getTime() << " [DEBUG:9]@float using LAPACK" << endl;

		// maximal number of nodes 32768 (2^15 = SHRT_MAX) -> laplace matrix 2^15*2^15=2^30 (max vector size) entries 
		// max vector size = std::vector<int> myvector; cout << myvector.max_size() << "\n"; -> 2^30
		// used ram in MB of lapack = (unsigned int)n*(unsigned int)n*sizeof(float))/1e+6 

		float * laplacian = (float*)calloc( (unsigned int)n*(unsigned int)n,sizeof(float) );

		bool fill_laplacian_return=1; // return value of the fill algorithm -> true : all fine, false : ERROR
		// fill laplacian
		for(unsigned int i = 0 ; i < (unsigned int)n ; i++){

			unsigned int from = (*nodes)[i]; 
			if(!mapping.count(from)){continue;}
			unsigned int sum = 0;

			for(unsigned int j = 0 ; j < graph[from].edges.size() ; j++){

				unsigned int to = graph[from].edges[j].edge;
				if(!mapping.count(to)){continue;}
				unsigned int vector_idx = mapping[from] + mapping[to]*n; //transform the 2-d coordinate (i,j) of the nxn matrix to 1-d vector coordinate i+j*n of the 2n vector

				if( vector_idx >= (unsigned int)n*(unsigned int)n){
					fill_laplacian_return = false;
					break;
				}

				if( useWeights
					// && (unsigned int)n <= param_max_nodes_weight
					){
					float w = graph[from].edges[j].weight;
					sum+=w;
					laplacian[vector_idx]=-w;
					if(maxWeight<w)maxWeight=w;
				}else{
					sum++;
					laplacian[vector_idx]=-1.0;
				}
			}

			laplacian[mapping[from]+mapping[from]*n]=sum;
		}
		if(!fill_laplacian_return){
			cerr <<"CRITICAL ERROR : fill_laplacian : out of range" << "\n";
			throw;
		}

		// local variables:
			int il, iu, m = 1, lda = n, ldz = n, info, lwork, liwork, iwkopt;
			float vl, vu;
			float wkopt;
			float* work;
			int* iwork;
			int isuppz[(unsigned int)(2*m)];
			char Vchar='V', Ichar='I', Uchar='U'; // Ichar = for specific range of eigenvalues/vectors
			float eigenvalues[(unsigned int)n]; // need only 1 eigenvalue
			il = 2; //that is the second one (il=1 -> the first one, il=2 the second one)
			iu = 2; 
			//if(debug_level==15){ iu = 5; } // generate more vectors for debug 15
			float * eigenvectors = (float*)malloc( (unsigned int)ldz*(unsigned int)m*(iu-il+1)*sizeof(float) ); 
			float eps=param_epsilon;

		// Determine optimal workspace 
			lwork = -1;
			liwork = -1;
			int n_int=(int)n;
			
			ssyevr_( &Vchar, &Ichar, &Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu, &eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork, &info  );
			
			//dssyevr_<float>( &Vchar, &Ichar, &Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu, &eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork, &info );
			
			lwork = (int)wkopt;
			work = (float*)malloc( lwork*sizeof(float) );
			liwork = iwkopt;
			iwork = (int*)malloc( liwork*sizeof(int) );

		// Solve eigenproblem ...
			ssyevr_( &Vchar, &Ichar,&Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu,&eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, work, &lwork, iwork, &liwork , &info  );
			
			//dssyevr_<float>( &Vchar, &Ichar,&Uchar, &n_int, laplacian, &lda, &vl, &vu, &il, &iu,&eps, &m, eigenvalues, eigenvectors, &ldz, isuppz, work, &lwork, iwork, &liwork , &info );
			
		// Check for errors in convergence
		if( info > 0 ) {
			cerr << " [ERROR] The algorithm (d|s)syevr failed to compute eigenvalues. Continue now with the slow standard approach (power iteration)." << "\n";
			throw;
			// goto standardComputationOfAlgCon;
		}

		// deallocate
		delete [] laplacian;
		delete [] work;
		delete [] iwork;

		// calculate normalized algebraic connectivity and fill x_hat
		if(useWeights){
			connectivity = eigenvalues[0]/(maxWeight*(float)n);
		}else{
			connectivity = eigenvalues[0]/((float)n);
		}
		for(unsigned int i = 0 ; i < (unsigned int)n ; i++)
			(*x_hat)[i]=-eigenvectors[i];

		delete [] eigenvectors;	

	}else{

		if (debug_level > 0) cerr << getTime() << " [DEBUG:9]@float using POWER" << endl;

		if(debug_level==9) cerr << getTime() << " [DEBUG:9]@float start" << endl;

		if(param_useWeights && !useWeights){
			cerr << " [INFO] The maximum number of nodes for the weighted algorithm is exeeded. Continue now with the faster unweighted algorithm." << "\n";
		}
		//diagonal matrix diag : d(u,u)=number of adjacency nodes=deg(u)
		vector<unsigned int> diag(n);

		if(useWeights){
			for (unsigned int i = 0; i < (unsigned int)n; i++) {
				diag[i]=0;
				unsigned int from = (*nodes)[i]; 
				if(!mapping.count(from)){continue;}

				for (unsigned int j = 0; j < graph[(*nodes)[i]].edges.size(); j++) {

					unsigned int to = graph[from].edges[j].edge;
					if(!mapping.count(to)){continue;}
					diag[i] += graph[(*nodes)[i]].edges[j].weight;
					if(useWeights && maxWeight<graph[(*nodes)[i]].edges[j].weight)maxWeight=graph[(*nodes)[i]].edges[j].weight;
				}
			}
		}
		// else{
		// 	#pragma omp parallel for reduction(max: maxWeight)
		// 	for (unsigned int i = 0; i < (unsigned int)n; i++) {
		// 		diag[i]=graph[(*nodes)[i]].edges.size();
		// 		if(useWeights)
		// 			for (unsigned int j = 0; j < graph[(*nodes)[i]].edges.size(); j++) {
		// 				if(maxWeight<graph[(*nodes)[i]].edges[j].weight)maxWeight=graph[(*nodes)[i]].edges[j].weight;
		// 			}
		// 	}
		// }
		if(debug_level==9) cerr << getTime() << " [DEBUG:9]@float ini done" << endl;

		// Get max degree / sum of weights of nodes
		unsigned int max_d = max_of_diag((*nodes),diag);	

		// Init randomized variables. 
		vector<float> x = generate_random_vector<float>(n);		

		// Orthogonalize + normalize vector + get initial lenght
		float current_length = 0;
		float last_length;

		(*x_hat) = makeOrthogonal<float>(x);
		normalize(x_hat, &last_length);

		// Repeat until difference < param_epsilon
		unsigned int iter = 0;	// catch huge clustering issues by keeping track here

		#ifdef timeAnalysis
			auto t_tmp = std::chrono::steady_clock::now();
		#endif	

		while(++iter < param_max_iter) { 

			if(debug_level==9) cerr << getTime() << " [DEBUG:9]@float iter:" << iter << endl;

			last_length = current_length;

			// Get a new x
			x = get_new_x<float>(x_hat, (*nodes), mapping, useWeights);

			// Get y
			vector<float> y = getY<float>(max_d,*x_hat,&x,(*nodes),diag);

			// Orthogonalize
			(*x_hat) = makeOrthogonal<float>(y);

			// Get length (lambda) & normalize vector
			normalize(x_hat, &current_length);

			if ( abs(current_length-last_length) < param_epsilon && iter >= min_iter ) break;	// prevent convergence by chance, converge to param_epsilon
		}
		if(debug_level==9) cerr << getTime() << " [DEBUG:9]@float post while" << endl;

		#ifdef timeAnalysis
			if(!t_master.count("partition_graph::convergence")) t_master["partition_graph::convergence"]=0;
			t_master["partition_graph::convergence"] += (float)std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9;
		#endif

		if (debug_level > 0) cerr << getTime() << " [DEBUG]   " << iter << " / " << param_max_iter << " iterations required (error is " << abs(current_length-last_length) << ")" << "\n";

		#ifdef DEBUG
			total_number_of_iterations_convergence+=iter;
		#endif

		if(useWeights){
			connectivity = (-current_length+(float)2*max_d)/(maxWeight*(float)n);
		}else{
			connectivity = (-current_length+(float)2*max_d)/((float)n);
		}
		//normalize(x_hat, &current_length);
		//if(debug_level==9) cerr << getTime() << " [DEBUG:9]@float post last norm" << endl;
		
		// 5.17 catch hardly-converging groups
		if (iter >= param_max_iter) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Connectivity score of connected component with " << n << " elements did not converge perfectly in time." << "\n";
		}
		if(debug_level>0) cerr << getTime() << " [DEBUG]   length:" << current_length << endl;
	}

	//float cl = 0;
	//normalize(x_hat, &cl);

	if (debug_level > 1){
		cerr << getTime() << " [DEBUG]   Connectivity score:" << connectivity;
		if ( (debug_level > 1 && (unsigned int)n<100 ) || debug_level > 2 ){
			cerr << " eigenvector: (";
			for(unsigned int i = 0 ; i < (unsigned int)n ; i++)
				cerr << (*x_hat)[i] << ","; 
			cerr << ")";
		}
		cerr << "\n";
	}

	// Split groups if connectivity is too low, remove tree like structures that might have arosen
	if (connectivity < param_con_threshold) {
		
		// 5.17 new threshold option overwrites connectivity
		if (param_min_species >=1 ) {
			if (debug_level > 0) cerr << getTime() << " [DEBUG]  Start the calculation of the average gene/species score " << "\n";
			float avg = calc_group(nodes);
			if (debug_level > 0) cerr << getTime() << " [DEBUG]   Found " << avg << " genes/species on average. User asked for at least " << param_min_species << "\n";
			if (avg <= param_min_species) {
				if (debug_level > 0) cerr << getTime() << " [DEBUG]   Group is going to be accepted despite connectivity " << connectivity << "\n";
				// just to be safe from infinit loops due to rounding
				if (connectivity == 0) connectivity = 0.001;
				// no splitting despite bad connectivity
				return -connectivity;			
			}
		}
	}
	
	return connectivity;
}

bool comparator_pairfloatUInt ( const pair<float,unsigned int>& l, const pair<float,unsigned int>& r )
   { return l.first < r.first; }

auto t_tmp = std::chrono::steady_clock::now( );

string getTime(void) {
	//time_t t = time(0);   // get time now
	//struct tm * now = localtime( & t );
	//ostringstream oss;
	//if (now->tm_hour < 10) oss << "0";
	//oss << now->tm_hour << ":";
	//if (now->tm_min < 10) oss << "0";
	//oss << now->tm_min << ":";
	//if (now->tm_sec < 10) oss << "0";
	//oss << now->tm_sec;
	string res = to_string(std::chrono::duration_cast<std::chrono::nanoseconds>( std::chrono::steady_clock::now( ) - t_tmp ).count()/1e+9);
	t_tmp = std::chrono::steady_clock::now( );
	return res;
}

#ifdef DEBUG
	////////////////////// Debug
	void debug__print_edgelist (protein& node, const unsigned int index, const int node_id) {
		cerr << node_id << ": ";	
		for (unsigned int j = 0; j < node.edges.size(); j++) {
			if (j == index) cerr << "*";
			cerr << node.edges[j].edge << " ";
		}
		cerr << "\n";
	}

	void debug__conn_integrity(vector<unsigned int>& nodes, float conn) {
		if (nodes.size() > 5) return;

		unsigned int sum = 0;
		for (unsigned int a = 0; a < nodes.size(); a++) {
			unsigned int from = nodes[a];
			sum += graph[from].edges.size();
		}

		sum /= 2;

		if (nodes.size() == 3) {
			if (sum == 2 && conn > 0.4) {
				cerr << "gs 3 with 2 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 3a").c_str() << "\n";throw;
			}
			if (sum == 3 && conn < 1) {
				cerr << "gs 3 with 3 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 3b").c_str() << "\n";throw;
			}
		}
		if (nodes.size() == 4) {
			if (sum == 3 && conn > 0.4) {
				cerr << "gs 4 with 3 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 4a").c_str() << "\n";throw;
			}
			if (sum < 6 && conn == 1) {
				cerr << "gs 4 with <6 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 4b").c_str() << "\n";throw;
			}
			if (sum == 6 && conn < 1) {
				cerr << "gs 4 with 6 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 4c").c_str() << "\n";throw;
			}
		}
		if (nodes.size() == 5) {
			if (sum == 4 && conn > 0.4) {
				cerr << "gs 5 with 4 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 5a").c_str() << "\n";throw;
			}
			if (sum < 10 && conn == 1) {
				cerr << "gs 5 with <10 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 5b").c_str() << "\n";throw;
			}
			if (sum == 10 && conn < 1) {
				cerr << "gs 5 with 10 edges had conn " << conn << "\n";
				cerr << string("integrity issue: connectivity of group size 5c").c_str() << "\n";throw;
			}
		}
	}

	void debug__graph_integrity(vector<unsigned int>& nodes) {
		// For each node
		for (unsigned int a = 0; a < nodes.size(); a++) {
			unsigned int from = nodes[a];
	//		cerr << "From: " << from << "\n";

			// For each edge in + direction
			for (unsigned int i = 0; i < graph[from].edges.size(); i++) {
				unsigned int to = graph[from].edges[i].edge;
	//			cerr << " To:  " << to << "\n";

				if (to == from) {cerr << "ERROR: Edge from " << from << " to " << to << " is selfevident" << "\n"; cerr << string("integrity issue self hit").c_str() << "\n";throw;}

				// Check reverse direction
				// Foreach edge in - direction
				bool found = false;
	//			cerr << " Back:";
				for (unsigned int j = 0; j < graph[to].edges.size(); j++) {
					unsigned int back = graph[to].edges[j].edge;
	//				cerr << " " << back;
					if (back == from) {found = true; continue;}
				}
	//			cerr << "\n";
				if (!found) {
					cerr << "ERROR: Edge from " << from << " to " << to << " is unidirectional" << "\n"; cerr << string("integrity issue direction").c_str() << "\n";throw;
				}
			}
		}
	}

	/* Auxiliary routine: printing a matrix */
	void debug__print_matrix( int m, int n, float* a, int lda ) {
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
		std::vector<float> v = generate_random_vector<float>(j);
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

		vector<float> x(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) x[i] = i; 

		//test weighted case
		vector<float> x_new = get_new_x<float>(&x,nodes,mapping,1);
	
		if(x_new.size() != x.size())return false;
		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(x_new[i]!=10*(nodes.size()-1)*(nodes.size())/2) return false;

		//test unweighted case
		x_new = get_new_x<float>(&x,nodes,mapping,0);
	
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

		vector<float> x(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) x[i] = i; 

		//test weighted case
		vector<float> x_new = get_new_x<float>(&x,nodes,mapping,1);

		if(x_new.size() != x.size())return false;
		if(x_new[0]!=1*10) return false;
		if(x_new[1]!=0*10+2*12+3*13+4*14) return false;
		if(x_new[2]!=1*12) return false;
		if(x_new[3]!=1*13+4*34) return false;
		if(x_new[4]!=1*14+3*34) return false;

		//test unweighted case
		x_new = get_new_x<float>(&x,nodes,mapping,0);
	
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
		std::vector<float> x = generate_random_vector<float>(j);
		x=makeOrthogonal<float>(x);
		float sum=0;
		for (unsigned int i = 0; i < x.size(); i++) {sum += x[i];}
		//cerr << sum << "\n";
		if(sum > 1e-3 || sum < -1e-3)return false;
	}
	return true;
}

bool test__normalize(){
	for(unsigned int j = 1 ; j < 100; j++){
		std::vector<float> x = generate_random_vector<float>(j);
		float len=0;
		normalize(&x,&len);
		normalize(&x,&len); // second normalize should return a length of 1 
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

		float max_degree = max_of_diag(nodes,diag);

		vector<float> x_hat(nodes.size()),x_new(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++){ x_hat[i]=i; x_new[i] = -(float)i*(2*max_degree - diag[i]); }

		vector<float> y = getY<float>(max_degree,x_hat,&x_new,nodes,diag);
 
		if(y.size() != x_hat.size())return false;
		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(y[i] > 1e-3 || y[i] < -1e-3)return false;
	}

	{
		vector<unsigned int> nodes(5);
		for(unsigned int i = 0 ; i < nodes.size() ; i++) nodes[i]=i;
		
		vector<unsigned int> diag(nodes.size());
		for(unsigned int i = 0 ; i < nodes.size() ; i++) diag[i]= (float)3;
		diag[3]=5;

		float max_degree = max_of_diag(nodes,diag); 

		vector<float> x_hat(nodes.size()),x_new(nodes.size());
		x_hat[0]= (float)2.1; 
		x_hat[1]= (float)-3; 
		x_hat[2]= (float)10; 
		x_hat[3]= (float)9.5; 
		x_hat[4]= (float)0; 

		x_new[0] = (float)-14.7; //calculated by hand
		x_new[1] = (float)21;
		x_new[2] = (float)-70;
		x_new[3] = (float)-47.5;
		x_new[4] = (float)0;

		vector<float> y = getY<float>(max_degree,x_hat,&x_new,nodes,diag);

		for(unsigned int i = 0 ; i < nodes.size() ; i++) if(y[i] > 1e-3 || y[i] < -1e-3)return false;
	}

	return true;

}

#endif /* _PROTEINORTHOCLUSTERING */