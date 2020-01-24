#ifndef KMER_SEQ_H
#define KMER_SEQ_H


#include <iostream>
#include <string>
#include <vector>
#include <map>

//ideas 
//- possibly do away with the istream and o stream, seems like a square peg in a round hole

class Kmer_seq{
friend std::istream &read(std::istream&, Kmer_seq&);
friend std::ostream &print_vals(std::ostream&, const Kmer_seq&, char sep);
friend std::vector<std::string> kmer_vec_build(int k, std::vector<std::string> dna_vec);
public:
	//constructors
	Kmer_seq() = default; //builds empty instance
	Kmer_seq(const std::string &h, const std::string &s):
				id(h), dna(s) {} //fill with passed header and dna sequence
	Kmer_seq(const std::string &h, const std::string &s, const int &k):
				id(h), dna(s), k_size(k) {} //fill with passed header and dna sequence
	Kmer_seq(std::istream &); // initalize with istream, calls read friend func

	// k_size is public, can change it if not provided by constructor
	int k_size = 0;
	// map is public so we can clear it and run multiple k counts 
	std::map<std::string, int> kmer_map; //a map to store the kmer count data

	//declarations of member operations
	Kmer_seq &build_kmap();//reference means in place
	Kmer_seq &count();// iterate through the string and count the kmer occurances
	//these let the relevant parts of the class be accessed
	inline std::string name() const { return id; }; 
	inline std::string seq() const{ return dna; };
	double gc_content();
	std::vector<std::string> keys() const;
	std::vector<int> values() const;
private:
	//private data components
	std::string id; //the header id for the fasta entry
	std::string dna;//the dna string

};


//friend functions defined in Kmer_seq.cpp
std::istream &read(std::istream&, Kmer_seq&); //reading Kmer_seq from ostream
std::ostream &print_vals(std::ostream&, const Kmer_seq&, char sep);//writing Kmer_seq to ostream

#endif
