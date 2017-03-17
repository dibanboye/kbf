#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <boost/math/special_functions/prime.hpp>
#include "cpp-src/KBFUtil.hpp"
#include "hash/generate_hash.h"

using namespace std;

int main(int argc, char* argv[]) {

   // fasta filename
   //string filename = argv[1];

   // get reads from file
   //vector<string> reads = parseFasta(filename); 
   // get kmers from file
   //vector<string> kmers = vector<string>();
   
   vector<string> reads;
   reads.push_back("ATTCG");
   unordered_set<kmer_t> kmers = getKmers(reads, 2);

   cout << endl;
   generate_hash f( kmers, kmers.size(), 2 );
   for (unordered_set<kmer_t>::iterator it1 = kmers.begin();
	it1 != kmers.end();
	++it1) {
     print_kmer( *it1, 2, cout );
     cout << ' ';
     cout << f.get_hash_value( *it1 ) << endl;
   }
   
}

