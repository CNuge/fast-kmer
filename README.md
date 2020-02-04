# fast-kmer
## A command line function for generating dataframes of kmer counts for a series of sequences

As a side project I have recently been learning how to program in C++, and was looking for a way I could apply the skills I've acquired thus far. I decided to [reinvent the wheel](https://www.cbcb.umd.edu/software/jellyfish/) and write a program to generate kmer count data from 

## How to compile
In the future I will add a make file, for now you can compile fast_kmer like so.

	g++ fast_kmer.cpp Kmer_seq.cpp -o fast_kmer.o

## How it works

Here is the help menu which shows the arguments that can be passed to the program

	fast_kmer - a command line tool for generating dataframes of kmer counts from fasta files.
	The program takes a fasta file filled with DNA sequences of interest, and generates a
	data frame with alignment free information of interest for each sequence in the file.
	Arguments:
	--in -i  : The input fasta file.
	--out -o : The name of file to output data to (either .csv or .tsv).
	           Default output name is 'output_kmers.tsv'.
	-k       : The size of the kmers to count. Multiple k values can be passed, each with their own -k flag.
	           Default k is 6.
	--gc     : Passing this flag with append a column with the gc content to the output dataframe.
	           Default is false.
	--help -h: Display the help menu.

## Examples

Count 3mers and 6mers for input file

	./fast_kmer.o -i data/small_unittest.fasta -k 6 -k 3

Count 5mers, also calculate gc content, and write to csv instead of default tsv.

	./fast_kmer.o -i data/small_unittest.fasta -o output_kmers.csv -k 5  --gc
