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
#include <chrono>
#include "proc_size.cpp"

#include "FDBG.cpp"
#include "formatutil.cpp"
#include "KBFUtil.hpp"


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
				  << argv[0] << " <reads.fa> <k>";
         exit(1);
   }
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);

   // get k-mers and edgemers from file
   unordered_set<kmer_t> kmers;
   unordered_set<kmer_t> edgemers;

   auto start = std::chrono::system_clock::now();
   handle_mers( filename, k, kmers, edgemers );
   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
      
   BOOST_LOG_TRIVIAL(info) << "Getting " << kmers.size() + edgemers.size() << " mers took " << elapsed_seconds.count() << " s";

   // print out the kmers 
   //   unordered_set<kmer_t>::iterator it1; 
   //   for (it1 = kmers.begin(); it1 != kmers.end(); ++it1) {
   //      BOOST_LOG_TRIVIAL(debug) << get_kmer_str(*it1, k);
   //   }

   BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
   start = std::chrono::system_clock::now();
   FDBG Graph( kmers, edgemers, kmers.size(), k, false );
   end = std::chrono::system_clock::now();
   elapsed_seconds = end-start;
   
   BOOST_LOG_TRIVIAL(info) << "Data structure built in " << elapsed_seconds.count() << " s";
   BOOST_LOG_TRIVIAL(info) << "Size(bytes):" << getCurrentRSS();
   BOOST_LOG_TRIVIAL(info) << "Size(Mb):" << getCurrentRSS() / 1024.0;

   BOOST_LOG_TRIVIAL(info) << "Size(bits):" << Graph.getBitSize();
   BOOST_LOG_TRIVIAL(info) << "Bits per element:" << Graph.getBitSize() / static_cast<double>( Graph.n );
   BOOST_LOG_TRIVIAL(info) << "Size(Mb):" << Graph.getBitSize() / (8.0 * 1024);

   BOOST_LOG_TRIVIAL(info) << "Membership check...";
   unordered_set<kmer_t>::iterator i; 
   for (i = kmers.begin(); i != kmers.end(); ++i) {
      //      BOOST_LOG_TRIVIAL(debug) << "Testing k-mer: " << *i << ' ' << get_kmer_str( *i, k );
      //      BOOST_LOG_TRIVIAL(debug) << Graph.inefficient_detect_membership( *i ) << ' ' << Graph.detect_membership( *i );
            if (!Graph.inefficient_detect_membership( *i )) {
            	 BOOST_LOG_TRIVIAL(fatal) << "Membership test failed.";
            	 exit(1);
            }
   }


   BOOST_LOG_TRIVIAL(info) << "Membership tests passed!";
   
   return 0;
}

