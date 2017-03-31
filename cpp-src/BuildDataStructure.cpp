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
#include "KBFUtil.hpp"

using namespace std;

/**
 * Run like "./a.out kmers.fasta 8"
 */
int main(int argc, char* argv[]) {

   // Set debug level
   boost::log::core::get()->set_filter(boost::log::trivial::severity
      >= boost::log::trivial::trace);

   cout << sizeof( uintmax_t ) << endl;

   uint128_t aa = 0;
   cout << aa << ' ' << sizeof( aa ) <<  endl;
   
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

   // get reads from file
   vector<string> reads = parseFasta(filename); 
   unordered_set<kmer_t> kmers = getKmers(reads, k);
   unordered_set<kmer_t> edgemers = getKmers(reads, k + 1);
   BOOST_LOG_TRIVIAL(info) << "Read in " << kmers.size() << " kmers of size " << k;
 
   // print out the kmers 
   //   unordered_set<kmer_t>::iterator i; 
   //   for (i = kmers.begin(); i != kmers.end(); ++i) {
   //      BOOST_LOG_TRIVIAL(debug) << get_kmer_str(*i, k);
   //   }
 
   BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
   auto start = std::chrono::system_clock::now();
   FDBG Graph( kmers, edgemers, kmers.size(), k, false );
   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
   
   BOOST_LOG_TRIVIAL(info) << "Data structure built in " << elapsed_seconds.count() << " s";
   reads.clear();
   BOOST_LOG_TRIVIAL(info) << "Size(bytes):" << getCurrentRSS();
   BOOST_LOG_TRIVIAL(info) << "Size(Mb):" << getCurrentRSS() / 1024.0;

   BOOST_LOG_TRIVIAL(info) << "Size(est. bytes):" << Graph.get_size();
   BOOST_LOG_TRIVIAL(info) << "Size(est. Mb):" << Graph.get_size() / 1024.0;

   BOOST_LOG_TRIVIAL(info) << "Size(bits):" << Graph.getBitSize();
   BOOST_LOG_TRIVIAL(info) << "Bits per element:" << Graph.getBitSize() / static_cast<double>( Graph.n );
   BOOST_LOG_TRIVIAL(info) << "Size(Mb):" << Graph.getBitSize() / (8.0 * 1024);
   
   return 0;
}

