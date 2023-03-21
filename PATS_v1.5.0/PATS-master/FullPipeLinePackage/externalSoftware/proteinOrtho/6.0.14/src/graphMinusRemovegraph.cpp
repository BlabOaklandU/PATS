/*
 *	Cleanup algorithm for Proteinortho
 *	Reads removegraphfile and graphfile and 
 *	computes the difference-graph
 *
 *	Last updated: 2018/03/01		
 *	Author: Marcus Lechner
 */

#ifndef _GRAPHMINUSREMOVEGRAPH
#define _GRAPHMINUSREMOVEGRAPH

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <cmath>
#include <vector>
#include <map>
#include <stack>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <unordered_set>

using namespace std;

// Functions
double string2double(string);
void tokenize(const string& , vector<string>& , const string&);
void parse_file(string,unordered_set<string>*);
void parse_removegraph(string,unordered_set<string>*);

///////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////


int main(int argc, char *argv[]) {

	try {

		vector<string> files;
		for (int paras = 1; paras < argc; paras++) {
			string parameter = string(argv[paras]);
			if (parameter.substr(0, 1) != "-") {
				files.push_back(parameter);
			}
		}

		if(files.size() < 2){		
			cerr << "USAGE : proteinortho_graphMinusRemovegraph REMOVEGRAPHFILE GRAPHFILE(S)" << endl;
			cerr << "The GRAPHFILE-REMOVEGRAPHFILE will printed out (cout)." << endl;
			cerr << "[ERROR]   invalid number of arguments." << endl;
			return EXIT_FAILURE;
		}

		unordered_set<string> removehash;

		parse_removegraph(files[0],&removehash);

		for(unsigned int i = 1 ; i < files.size() ; i++)
			parse_file(files[i],&removehash);
	}
	catch(string& error) {
		cerr << "USAGE : proteinortho_graphMinusRemovegraph REMOVEGRAPHFILE GRAPHFILE(S)" << endl;
		cerr << "The GRAPHFILE-REMOVEGRAPHFILE will printed out (cout)." << endl;
		cerr << "[ERROR]   " << error << endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

///////////////////////////////////////////////////////////
// File parser
///////////////////////////////////////////////////////////
void parse_file(string file,unordered_set<string>* removehash) {
	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {

		// For each line
		string file_a = "";	string file_a_id = "";
		string file_b = "";	string file_b_id = "";
		while (!graph_file.eof()) {
			getline(graph_file, line);
			vector<string> fields;
			tokenize(line, fields, "\t");
			// Header line
			if (fields.size() == 2 && fields[0].substr(0, 1) == "#") {
				file_a = fields[0].substr(2, fields[0].size()-2);
				file_b = fields[1];

				file_a_id = file_a;//species2id[file_a];
				file_b_id = file_b;//species2id[file_b];

				cout << line<< "\n";
			}
			// Data line
			else if ( (fields.size() == 6 || fields.size() == 8) && fields[0].substr(0, 1) != "#") {
				// a b e1 b1 e2 b2 score

//				cerr << protein_counter << ": " << fields[0] << " <-> " << fields[1] << endl;

				// 5.16 deal with duplicated IDs by adding file ID to protein ID
				string ida = fields[0];
				string idb = fields[1];

				fields[0] += " "; fields[0] += file_a_id;
				fields[1] += " "; fields[1] += file_b_id;

				// 5.16 do not point to yourself
				if (!fields[0].compare(fields[1])) {continue;}

				if(!( fields[0]<fields[1]?(*removehash).count(fields[0]+fields[1]):(*removehash).count(fields[1]+fields[0]) )){
					cout << line<< "\n";
				}
			}else if(line!="")
				cout << line << "\n";
		}
		graph_file.close();
	}
	else {
		throw string("Could not open file " + file);
	}
}


void parse_removegraph(string file,unordered_set<string>* removehash) {

	string line;
	ifstream graph_file(file.c_str());
	if (graph_file.is_open()) {
		// For each line
		while (!graph_file.eof()) {
			getline(graph_file, line);
			vector<string> fields;
			tokenize(line, fields, "\t");
			// Header line
			if (fields.size() == 4 ) {
				// a b e1 b1 e2 b2 score
				// 5.16 deal with duplicated IDs by adding file ID to protein ID
				fields[0] += " "; fields[0] += fields[1];
				fields[2] += " "; fields[2] += fields[3];

				fields[0]<fields[2]?(*removehash).insert(fields[0]+fields[2]):(*removehash).insert(fields[2]+fields[0]);
			}
		}
		graph_file.close();
	}
	else {
		throw string("Could not open file " + file);
	}
}

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
double string2double(string str) {
	istringstream buffer(str);
	double value;
	buffer >> value;
	return value;
}

#endif /* _GRAPHMINUSREMOVEGRAPH */