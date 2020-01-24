#include <iostream>
#include <fstream>


#include <string>
#include <vector>
#include <map>

#include "Kmer_seq.h"

using namespace std;



////////////
// friend functions
////////////

// This recrusive function builds a vector containing
// all the dna string kmers up to the specified k size
vector<string> kmer_vec_build(int k , vector<string> dna_vec){

	vector<string> nts = {"A", "T", "G", "C"};

	int new_k = k-1;

	if(new_k == 0){
		return dna_vec;
	}else{
		//new vector of kmers to be constructed
		vector<string> new_dna_vec;
		// for new kmer construction
		string new_kmer;

		for(string kmer : dna_vec){
			for (string n : nts){
				 new_kmer = kmer + n ;
				 new_dna_vec.push_back(new_kmer);
			}
		}
		return kmer_vec_build(new_k, new_dna_vec);
	}
}


////////////
// input and output
////////////

// called by the istram constructor, reads in the istream and stores the header then the sequence
// for now won't strip the leading > off of the id, can do this outside or add in here later 
istream& read(istream &is, Kmer_seq &item){
	is >> item.id >> item.id; 
	return is;
}

//called to write to file, returns an id and then a tab delimited list
//of the Kmer count values
ostream& print_vals(ostream &os, const Kmer_seq &item, char sep){
	for (int k: item.values()){
		os << sep << k;		
	}
	return os;
}


////////////
// member functions
////////////


//build a hash map for kmers of the specified size
Kmer_seq &Kmer_seq::build_kmap(){

	//get a vector of all the kmers
	vector<string> kmer_seed = {"A", "C", "G", "T"};
	vector<string> kmers;

	kmers = kmer_vec_build(k_size, kmer_seed);//initiate with the basic nts

	for(string k_str : kmers){
		//k_map.insert(pair<string, int>(k_str, 0));
		kmer_map[k_str] = 0;
	}

	return *this;
}

//TODO - make sure this is robust to Ns and weird chars 
// 	   - i.e. they're skipped and not added to the map
// iterate through the string and count the kmer occurances
Kmer_seq &Kmer_seq::count(){
	string k_str; 

	for(int i = 0; i < (dna.size() - k_size); ++i){
		k_str = dna.substr(i, k_size);
		kmer_map[k_str]++;
	}

	return *this;
}


//return the kmer_map keys - can be used to obtain the header line
vector<string> Kmer_seq::keys() const {

	vector<string> k_keys;

	auto it = kmer_map.begin();
	while(it != kmer_map.end()){
		k_keys.push_back(it->first);
		it++;
	}

	return k_keys;	
}

//return the kmer_map values - can be used to obtain the column values
vector<int> Kmer_seq::values() const {
	vector<int> k_vals;

	auto it = kmer_map.begin();
	while(it != kmer_map.end()){
		k_vals.push_back(it->second);
		it++;
	}
	return k_vals;	
}


// take the string of DNA in the Kmer_seq instances and
// calculate the gc content
double Kmer_seq::gc_content(){
	int gc_count = 0; //note if not initalized to 0 weird behaviour ensues!
	int at_count = 0;
	int total_bp = 0;
	char val;
	double gc_perc;

	//for (int i = 0; i < dna.size() ; i++){
	for (char v : dna){
		val = toupper(v);
		if ((val == 'C') || (val == 'G')) {
			++gc_count;
		}else if ((val == 'A') || (val == 'T')){
			++at_count;
		}
	}

	total_bp = gc_count + at_count; //this omits Ns from the total
	gc_perc = gc_count/(double)total_bp; //convert one part to double and all change

	return gc_perc;
}





