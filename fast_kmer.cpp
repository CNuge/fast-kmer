#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include "Kmer_seq.h"

using namespace std;


//TODO
//currently have the data being spat to cout, see if instead you can specify a file for 
//cout, maybe a global option or something of that sort?

/*
compile:
g++ fast_kmer.cpp Kmer_seq.cpp -o fast_kmer.o

use:
./fast_kmer.o -i data/small_unittest.fasta -k 6 -k 3
./fast_kmer.o -i data/small_unittest.fasta -o output_kmers.csv -k 5  --gc
./fast_kmer.o -i data/small_unittest.fasta -o output_kmers.csv -k 3 --gc
*/

void process_line(string name, string content,
					vector<int> k_vals, bool calc_gc,char out_delim){

	//construct the Kmer_seq record
	Kmer_seq s(name, content);
	
	//write the name to the outfile
	cout << s.name();

	//for each kmer size, write the values to the file
	for(int k : k_vals){
		//set the k for the instance and given loop
		s.k_size = k;
		s.build_kmap();
		s.count();

		//write kmer vals to the file
		print_vals(cout, s, out_delim);

		//clear the kmer map for the next loop
		s.kmer_map.clear();
	}

	if(calc_gc == true){
		cout << out_delim << s.gc_content();
	}
	s.kmer_map.clear();
	//endline to complete the entry
	cout << endl;

}

////////////
// reading and processing of sequences to kmer df - line-by-line low overhead
////////////
int main(int argc, char **argv){


	string help_info = 
	"fast_kmer - a command line tool for generating dataframes of kmer counts from fasta files."
	"The program takes a fasta file filled with DNA sequences of interest, and generates a" 
	"data frame with alignment free information of interest for each sequence in the file.\n"
	"Arguments:\n"
	"--in -i  : The input fasta file.\n"
	"--out -o : The name of file to output data to (either .csv or .tsv).\n" 
	"           Default output name is 'output_kmers.tsv'.\n"
	"-k       : The size of the kmers to count. Multiple k values can be passed, each with their own -k flag.\n"
	"           Default k is 6.\n"
	"--gc     : Passing this flag with append a column with the gc content to the output dataframe.\n"
	"           Default is false.\n"
	"--help -h: Display the help menu.\n";

	//initate the variables to store the data we are trying to obtain
	string in_filename;
	string out_filename = "output_kmers.tsv";
	vector<int> k_vals;
	bool calc_gc = false;
	char out_delim = '\t';

	//iterate through the array - parse arguments based on flags
	for (int pos = 0; pos < argc; pos++){
		if(string(argv[pos]) == "-i" || string(argv[pos]) == "--in"){
			pos++;
			in_filename = string(argv[pos]);
		}else if(string(argv[pos]) == "-o" || string(argv[pos]) == "--out"){
			pos++;
			out_filename = string(argv[pos]);
		}else if(string(argv[pos]) == "-k"){
			pos++;
			k_vals.push_back(stoi(string(argv[pos])));
		}else if(string(argv[pos]) == "--gc"){
			calc_gc = true;
		}else if(string(argv[pos]) == "-h" || string(argv[pos]) == "--help"){
			cout << help_info << endl;
			return 0;
		}
	}

	//this sets the default k if no value was passed by user
	if(k_vals.empty()){
		k_vals.push_back(6);
	}

	// report the arguments that have been passed
	cout << "input filename is: " << in_filename << endl; 
	cout << "output filename is: " << out_filename << endl;
	cout << "kvals are:";
	for(int k : k_vals){
		cout << " " << k ;
	}
	cout << endl; 
	cout << "gc calculation: " << calc_gc << endl;

	//change delimiter if csv outfile indicated
	if(out_filename.substr((out_filename.size()-4),4) == ".csv"){
		out_delim = ',';
	}

	std::ifstream input(in_filename);
	if(!input.good()){
		std::cerr << "Error opening file: " << in_filename << endl;		
	}

	// after opening arg parsing and report, set cout to the output file
	std::ofstream out(out_filename);
	std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
	std::cout.rdbuf(out.rdbuf()); //redirect std::cout to out.txt!
	//note at the end there is a reset to normal buffer, then a print to say done

	// build an empty map, use it to write the header line
	Kmer_seq header("sequence", "");

	//write the kmer keys to file
	cout << header.name();
	for(int k : k_vals){
		header.k_size = k;
		header.build_kmap();
		for(string s : header.keys()){
			cout << out_delim << s;
		}
		//write the keys to file
		header.kmer_map.clear();
	}
	//if gc, write final field to file
	if(calc_gc == true){
		cout << out_delim << "gc_content";
	}
	cout << endl;

	//TODO - see if you can make this into a cin loop instead of the clunky readline
	//
	std::string line, name, content;
	// iterate through the file, getting lines y
	while( std::getline(input, line).good()){
		// check if currently no entry, or if newline has the start flag
		if(line.empty() || line[0] == '>'){
			// if existing entry, send it to the needed location
			if( !name.empty()){
				//if record, then process
				process_line(name, content, k_vals, calc_gc, out_delim);
				name.clear();
				content.clear();				
			}
			// set the name to the current line, less the leading >
			if( !line.empty()){
				name = line.substr(1); // after building record, take newline and override the field
			}
		} else if (!name.empty()){ // append sequence, only if there is a name it goes with
			if (line.find(' ') != std::string::npos){// halt if spaces in sequence
				name.clear();
				content.clear();
			}else {
				content += line;
			}	
		} 
	}
	//this block ensures the final line in processed following the loop
	if ( !name.empty()){
		process_line(name, content, k_vals, calc_gc, out_delim);
		name.clear();
		content.clear();
	}
	std::cout.rdbuf(coutbuf);
	cout << "done!" << endl;
}

