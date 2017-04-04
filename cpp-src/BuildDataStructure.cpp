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
   //   BOOST_LOG_TRIVIAL(info) << "Size of uint128_t ..." << sizeof( uint128_t );   

   // Check if the user put in the correct command line arguments
   if (argc < 3) {
         BOOST_LOG_TRIVIAL(fatal) << "Missing required arguments. Usage:"
				  << argv[0] << " <reads.fa> <k>";
         exit(1);
   }
  
   // fasta filename
   string filename = argv[1];
   int k = stoi(argv[2]);
   FDBG Graph;
   
   //Have we constructed this dataset before on these parameters?
   string dsfile = filename.substr( 0, filename.find_last_of( '.' ) ) + "fdbg" + to_string( k ) + ".bin";

   // get k-mers and edgemers from file
   unordered_set<kmer_t> kmers;
   unordered_set<kmer_t> edgemers;
   auto start = std::chrono::system_clock::now();
   handle_mers( filename, k, kmers, edgemers );
   auto end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
   
   if (file_exists( dsfile )) {
     //Yes, so avoid reconstructing
     BOOST_LOG_TRIVIAL(info) << "Loading from " << dsfile;
     ifstream ifile_ds( dsfile.c_str(), ios::in | ios::binary );
     Graph.load( ifile_ds );
     ifile_ds.close();
   } else {
   
     BOOST_LOG_TRIVIAL(info) << "Getting " << kmers.size() + edgemers.size() << " mers took " << elapsed_seconds.count() << " s";

     BOOST_LOG_TRIVIAL(info) << "Building De Bruijn Graph ...";
     start = std::chrono::system_clock::now();
     Graph.build( kmers, edgemers, kmers.size(), k);
     end = std::chrono::system_clock::now();
     elapsed_seconds = end-start;

     BOOST_LOG_TRIVIAL(info) << "Data structure built in " << elapsed_seconds.count() << " s";
     BOOST_LOG_TRIVIAL(info) << "Writing data structure to file " << dsfile;
     ofstream ofile( dsfile.c_str(), ios::out | ios::binary );
     Graph.save( ofile );
     ofile.close();
   }

   //   BOOST_LOG_TRIVIAL(info) << "Estimated size(bits) (Mb):" << Graph.estimateBitSize()/(8.0 * 1024 * 1024);
   //   BOOST_LOG_TRIVIAL(info) << "Size(Mb):" << Graph.bitSize() / (8.0 * 1024 * 1024);
   //   BOOST_LOG_TRIVIAL(info) << "Bits per element:" << Graph.bitSize() / static_cast<double>( Graph.n );

   BOOST_LOG_TRIVIAL(info) << "Membership check...";
   unordered_set<kmer_t>::iterator i; 
   for (i = kmers.begin(); i != kmers.end(); ++i) {
      //      BOOST_LOG_TRIVIAL(debug) << "Testing k-mer: " << *i << ' ' << get_kmer_str( *i, k );
      //      BOOST_LOG_TRIVIAL(debug) << Graph.inefficient_detect_membership( *i ) << ' ' << Graph.detect_membership( *i );
            if (!Graph.detect_membership( *i )) {
            	 BOOST_LOG_TRIVIAL(fatal) << "Membership test failed.";
            	 exit(1);
            }
   }


   BOOST_LOG_TRIVIAL(info) << "Membership tests passed!";
   
   return 0;
}

