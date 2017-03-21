#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#define BOOST_LOG_DYN_LINK 1
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include "cpp-src/KBFUtil.hpp"
#include "FDBG.cpp"


using namespace std;

/**
 * Run like "./a.out kmers.fasta 8"
 */
int main(int argc, char* argv[]) {

   // Set debug level
   boost::log::core::get()->set_filter(boost::log::trivial::severity
      >= boost::log::trivial::trace);

   BOOST_LOG_TRIVIAL(info) << "Beginning to build data structure ...";   

   // Check if the user put in the correct command line arguments
   if (argc < 3) {
         BOOST_LOG_TRIVIAL(fatal) << "Missing required arguments. Usage:"
                                     "kbf <reads.fa> <k>";
         exit(1);
   }
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);

   // get reads from file
   vector<string> reads = parseFasta(filename); 
   unordered_set<kmer_t> kmers = getKmers(reads, k);
   BOOST_LOG_TRIVIAL(info) << "Read in " << kmers.size() << " kmers of size " << k;
 
   // print out the kmers 
   unordered_set<kmer_t>::iterator i; 
   for (i = kmers.begin(); i != kmers.end(); ++i) {
      BOOST_LOG_TRIVIAL(debug) << get_kmer_str(*i, k);
   }
 
   BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
   FDBG Graph( reads, kmers, kmers.size(), k, false );
   return 0;

   print_kmer( 0, k, cout );
   if (Graph.detect_membership( 0 ))
     cout << " member" << endl;
    else
     cout << " not member" << endl;

   unordered_set<kmer_t> kmers2 = getKmers(reads, k);

   for (unordered_set<kmer_t>::iterator it1 = kmers2.begin();
     	it1 != kmers2.end();
     	++it1) {
      print_kmer( *it1, k, cout );

      if (Graph.detect_membership( *it1 ))
	cout << " member" << endl;
      else
	cout << " not member" << endl;
    }


   return 0;
}

