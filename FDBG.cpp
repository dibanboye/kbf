#ifndef FDBG_class
#define FDBG_class

//#include "hash/HashUtil.cpp"
#include "hash/generate_hash.h"
#include <vector>
#include <string>

using namespace std;

class FDBG {
public:

  vector< vector< bool > > IN; //size n x sigma, in edges
  vector< vector< bool > > OUT; //size n x sigma, out edges
  unsigned sigma; //alphabet-size. For now, only 4 is supported
  unsigned n; //number of nodes in graph
  unsigned k; //length of each mer (string in alphabet)
  generate_hash f;

  FDBG( vector< string >& reads, 
	unordered_set<kmer_t>& kmers,
	unsigned n, //number of kmers
	unsigned k, //mer size
	bool b_verify = false, //if true, print summary
	ostream& os = cout
	) 
  {
    sigma = 4;
    
    this->n = n;
    this->k = k;
    
    //construct hash function f
    f.construct_hash_function( kmers, n, k );

    //initialize IN, OUT to zero (false)
    vector< bool > vzero( sigma, false );
    IN.assign( n, vzero );
    OUT.assign( n, vzero );

    add_edges( reads, b_verify, os );
    if (b_verify) {
      for (unordered_set<kmer_t>::iterator it1 = kmers.begin();
	   it1 != kmers.end(); ++it1) {
	print_kmer( *it1, k, cout );
	os << ' ';
	os << f( *it1 ) << endl;
      }
    }
  }

  void add_edges( vector< string >& reads, bool b_verify = false, ostream& os = cout ) {
    string read;
    string kplusone;
    for (unsigned i = 0; i < reads.size(); ++i) {
      read = reads[i];
      unsigned index1 = 0;

      while (index1 + k < read.size()) {
	//we can get a k + 1 - mer
	kplusone = read.substr( index1, k + 1 );
	add_edge( kplusone );
	++index1;
      }
    }

    if (b_verify) {
      for (unsigned i = 0; i < reads.size(); ++i) {
	os << reads[i] << endl;
      }
      os << "IN: " << endl;
      print_matrix( IN, os );

      os << "OUT: " << endl;
      print_matrix( OUT, os );
      
    }
  }

  void add_edge( string& edge ) {
    kmer_t u,v;
    split_edge( edge, u, v );
    unsigned first, last;
    first = access_kmer( u, k, 0 );
    last = access_kmer( v, k, k - 1 );
    // OUT[ f(u), last ] = 1
    OUT[ f(u) ][ last ] = true;
    // IN[ f(v), first ] = 1
    IN[ f(v) ][ first ] = true;
  }

  void split_edge( string& edge, kmer_t& u, kmer_t& v ) {
    unsigned index = 1;
    while (index <= k) {
      set_kmer( u, k, index - 1, edge[ index - 1 ] );
      set_kmer( v, k, index - 1, edge[ index ] );

      ++index;
    }
  }

  void print_matrix( vector< vector< bool > >& mat, ostream& os = cout ) {
    for (unsigned i = 0; i < mat.size(); ++i) {
      for (unsigned j = 0; j < mat[i].size(); ++j) {
	os << mat[i][j] << ' ';
      }
      os << endl;
    }
  }

};


#endif

