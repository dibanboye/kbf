#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "cpp-src/KBFUtil.hpp"
#include "FDBG.cpp"

using namespace std;

int main(int argc, char* argv[]) {

  if (argc < 3) {
        cerr << "\tMissing required arguments." << endl;
        cerr << "\tUsage:" << endl;
        cerr << "\tkbf <reads.fa> <k>" << endl;
        exit(1);
  }
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);
   // get reads from file
   vector<string> reads = parseFasta(filename); 
   unordered_set<kmer_t> kmers = getKmers(reads, k);

   cout << endl;
   
    // generate_hash f( kmers, kmers.size(), k );
    // for (unordered_set<kmer_t>::iterator it1 = kmers.begin();
    // 	it1 != kmers.end();
    // 	++it1) {
    //   print_kmer( *it1, k, cout );
    //   cout << ' ';
    //   cout << f( *it1 ) << endl;
    // }

   FDBG Graph( reads, kmers, kmers.size(), k, true );

   return 0;
}

