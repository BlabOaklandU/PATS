//compile: g++ -std=c++11 -O3 cleanupblastgraph.cpp -o cleanupblastgraph
//Usage: ./cleanupblastgraph BLASTGRAPH (BLASTGRAPH2,BLASTGRAPH3,...)
//-> removes all duplicated edges 
//pk
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <map>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <list>
#include <unordered_set>
#include <cstring>

using namespace std;

struct wedge {unsigned int edge;};
struct protein {vector<wedge> edges;};

// Globals
unsigned int species_counter = 0;	// Species
unsigned int protein_counter = 0;	// Proteins
vector<string> species;			// Number -> Name
vector<protein> graph;			// Graph containing all protein data

// TMP Globals
map<string,unsigned int> species2id;		// Name -> Number
map<string,unsigned int> protein2id;		// Name -> Number

///////////////////////////////////////////////////////////
// Misc functions
///////////////////////////////////////////////////////////
// Convert string to double
double string2double(string str) {
	istringstream buffer(str);
	double value;
	buffer >> value;
	return value;
}
// Convert string to float
float string2float(string str) {
	istringstream buffer(str);
	float value;
	buffer >> value;
	return value;
}


pair<unsigned int, unsigned int> ordPair(unsigned int a,unsigned int b){if(a<b){return make_pair(a,b);}else{return make_pair(b,a);}}

// Split a string at a certain delim
void tokenize(const string& str, vector<string>& tokens, const string& delimiters = "\t") {

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

}


///////////////////////////////////////////////////////////
// File parser
///////////////////////////////////////////////////////////
void parse_file(string file) {

	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {
		// For each line
		string file_a = "";	
		string file_b = "";	
		while (!graph_file.eof()) {
			getline(graph_file, line);
			vector<string> fields;
			tokenize(line, fields, "\t");
			// Header line
			if (fields.size() == 2 && fields[0].substr(0, 1) == "#") {
				file_a = fields[0].substr(2, fields[0].size()-2);
				file_b = fields[1];

				if (file_a == "file_a" && file_b == "file_b") continue;	// Init Header

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

				cout << line<<endl;
			}
			// Data line
			else if ((fields.size() >1) && fields[0].substr(0, 1) != "#") {
				// a b e1 b1 e2 b2 score

				// 5.16 deal with duplicated IDs by adding file ID to protein ID
				string ida = fields[0];
				string idb = fields[1];
				fields[0] += " "; fields[0] += file_a;
				fields[1] += " "; fields[1] += file_b;

				// 5.16 do not point to yourself
				if (!fields[0].compare(fields[1])) {continue;}

				// A new protein 
				map<string,unsigned int>::iterator a_name_it=protein2id.find(fields[0]);
				map<string,unsigned int>::iterator b_name_it=protein2id.find(fields[1]);
				if (a_name_it == protein2id.end())	{
					protein a;
					protein2id[fields[0]] = protein_counter++;
					graph.push_back(a);
					a_name_it=protein2id.find(fields[0]);
				}
				if (b_name_it == protein2id.end())	{
					protein b;
					protein2id[fields[1]] = protein_counter++;
					graph.push_back(b);
					b_name_it=protein2id.find(fields[1]);
				}
				bool isduplicatededge=0;
				if(a_name_it != protein2id.end() && b_name_it != protein2id.end()){
					for(wedge w : graph[protein2id[fields[0]]].edges){
						if(w.edge == protein2id[fields[1]]){
							isduplicatededge=1;
							break;
						}
					}
				}
				if(!isduplicatededge){
					cout << line<<endl;

					// Add link to graph (reciprocal)
					unsigned int a_id = protein2id[fields[0]];
					unsigned int b_id = protein2id[fields[1]];

					// 5.17, add weight
					wedge w;
					w.edge=b_id;
					graph[a_id].edges.push_back(w);
					w.edge=a_id;
					graph[b_id].edges.push_back(w);
					
				}else{
					cerr << "[WARNING] : found duplicated edge : " << line << " (removed)."<< endl;
				}

			}
		}
		graph_file.close();
	}
	else {
		throw string("Could not open file " + file);
	}

}


int main(int argc, char *argv[]) {
	// check for an argument
	
	if(argc < 2 || argv[1] == string("-h") || argv[1] == string("--h") || argv[1] == string("-help") || argv[1] == string("--help") || argv[1] == string("help") || argv[1] == string("h")){
		cerr << "Usage: " << argv[0] << " BLASTGRAPH (BLASTGRAPH2,BLASTGRAPH3,...)" << endl 
		<< "removes all duplicated edges (output goes to STDOUT, duplicated edges are printed to STDERR)" << endl;
		return -1;
	}

	// read in a text file that contains a real matrix stored in column major format
	// but read it into row major format

	string name="";

	for(int argi=1;argi<argc;argi++){
		cerr << "Parsing "<< argv[argi] << endl;
		name=argv[argi];
		parse_file(argv[argi]);
	}
	std::size_t found = string(name).find_last_of("/");
	name=string(name).substr(found+1);

	species2id.clear();
	cerr << "Parsing done ..." << endl;
}
