#ifndef FDBG_class
#define FDBG_class

//#include "hash/HashUtil.cpp"
#include "hash/generate_hash.h"
#include <vector>
#include <string>
#include <unordered_set>
#include <map>

using namespace std;

class nei_bits {
public:
  bool bits[3]; // xyz , where z specifies IN = 0, or OUT = 1, xy specifies A,C,G, or T

  nei_bits() {
    for (unsigned i = 0; i < 3; ++i) {
      bits[i] = 0;
    }
  }

  nei_bits( const nei_bits& rhs ) {
    for (unsigned i = 0; i < 3; ++i) {
      this->bits[i] = rhs.bits[i];
    }
  }
  
  nei_bits& operator=( const nei_bits& rhs ) {
    for (unsigned i = 0; i < 3; ++i) {
      this->bits[i] = rhs.bits[i];
    }
    return *this;
  }

  bool operator==( const nei_bits& rhs ) {
    for (unsigned i = 0; i < 3; ++i) {
      if (this->bits[i] != rhs.bits[i])
	return false;
    }

    return true;
  }

  friend ostream& operator<<(ostream& os, const nei_bits& rhs );
  
};

  ostream& operator<<(ostream& os, const nei_bits& rhs ) {
    for (unsigned i = 0; i < 3; ++i) {
      os << rhs.bits[i];
    }

    return os;
  }

class tree_info {
public:
  nei_bits parent;
  bool is_stored;

  tree_info() {
    is_stored = false;
  }
	       
};

class forest {
public:
  map< u_int64_t, kmer_t > stored_mers;
  vector< tree_info > parents;
};


/**
 * This class represents the fully dynamic De Bruijn graph
 */
class FDBG {

public:

  vector< vector< bool > > IN; //size n x sigma, in edges
  vector< vector< bool > > OUT; //size n x sigma, out edges
  unsigned sigma; //alphabet-size. For now, only 4 is supported
  unsigned n; //number of nodes in graph
  unsigned k; //length of each mer (string in alphabet)
  generate_hash f; //hash function that takes each kmer to 1..n
  forest myForest; //the forest (described in paper)
  
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

    //debugging
    //print each value and mapping
    unordered_set<kmer_t>::iterator i;
    for (i = kmers.begin(); i!= kmers.end(); ++i) {
       BOOST_LOG_TRIVIAL(debug) << get_kmer_str(*i, k) << " maps to "
       << f(*i);
    }

    //initialize IN, OUT to zero (false)
    BOOST_LOG_TRIVIAL(info) << "Initializing IN and OUT.";
    vector< bool > vzero( sigma, false );
    IN.assign( n, vzero );
    OUT.assign( n, vzero );

    //add edges to IN and OUT using reads
    BOOST_LOG_TRIVIAL(info) << "Adding edges to IN and OUT ...";
    add_edges( reads, b_verify, os );

    //Perform the forest construction
    construct_forest( kmers, 3*k*2 ); //alpha = 3klg(sigma)
  }

  /**
   * Add edges to IN and OUT using reads (not kmers)
   */
  void add_edges( vector< string >& reads, bool b_verify = false, ostream& os = cout ) {
    string read;
    string kplusone;
 
   //for each read, get k+1 length pieces that we can
   //figure out an edge between two kmers
   for (unsigned i = 0; i < reads.size(); ++i) {
      read = reads[i];
      unsigned index1 = 0;

      while (index1 + k < read.size()) {
	//get a k + 1 - mer
	kplusone = read.substr( index1, k + 1 );
        BOOST_LOG_TRIVIAL(debug) << "Read in k+1-mer " << kplusone;
	add_edge( kplusone );
	++index1;
      }
    }


    BOOST_LOG_TRIVIAL(debug) << "Printing IN ...";
    for (int i = 0; i < IN.size(); i++) {
        if (IN[i].size() == 4) {
           BOOST_LOG_TRIVIAL(debug) << " Node " << i << ": " << IN[i][0] << ", "
              << IN[i][1] << ", "
              << IN[i][2] << ", " << IN[i][3];
        }
        else {
           BOOST_LOG_TRIVIAL(error) << "IN does not have the correct dimensions.";
        }
    }

    BOOST_LOG_TRIVIAL(debug) << "Printing OUT ...";
    for (int i = 0; i < OUT.size(); i++) {
        if (OUT[i].size() == 4) {
           BOOST_LOG_TRIVIAL(debug) << " Node " << i << ": " << OUT[i][0] << ", "
              << OUT[i][1] << ", "
              << OUT[i][2] << ", " << OUT[i][3];
        }
        else {
           BOOST_LOG_TRIVIAL(error) << "OUT does not have the correct dimensions.";
        }
    }

  }

  //add edge implied by the k+1-mer (between the first k and last k characters).
  void add_edge( string& edge ) {

    // figure out the two kmers
    kmer_t u,v;
    split_edge( edge, u, v );
    BOOST_LOG_TRIVIAL(debug) << "Adding an edge from " << get_kmer_str(u, k)
        << " to " << get_kmer_str(v, k);

    // get which column in sigma to put in (corresponds to which letter is the first/last)
    // number 0..3 represent each alphabet letter
    unsigned first, last;
    first = access_kmer( u, k, 0 );
    last = access_kmer( v, k, k - 1 );

    // set edge in IN/OUT 
    OUT[ f(u) ][ last ] = true;
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

  void initialize_forest() {

    BOOST_LOG_TRIVIAL(info) << "Initializing the forest ...";
    tree_info empty_tree;
    empty_tree.is_stored = false;
    myForest.parents.assign( n, empty_tree );
  }

  bool detect_membership( kmer_t m ) {
    u_int64_t f_m = f(m);
    if (f_m > n)
      return false;

    vector< kmer_t > neis;
    vector< nei_bits > nbits;
    while ( !(myForest.parents[ f_m ].is_stored) ) {
      //we need to travel through the tree
      get_neighbors( m, neis, nbits );

      for (unsigned ii = 0; ii < nbits.size(); ++ii) {
	//	cout << ' ' << f_m << ' ' <<  nbits[ii] << ' ' << myForest.parents[ f_m ].parent << endl;
	if (nbits[ii] == myForest.parents[ f_m ].parent) {
	  m = neis[ii]; //move to parent in the tree
	  f_m = f(m);
	  if (f_m > n)
	    return false;
	  break;
	}
      }

    }

    //if we make it here, we have made it to the root
    kmer_t root = myForest.stored_mers[ f_m ];
    if (root == m)
      return true;
    else
      return false;
  }

  /*
   * Initially constructs the forest
   * Requires IN, OUT, f to be already constructed
   * alpha is the tree height parameter
   * Requires non-used bits of kmers to be zero
   */
  void construct_forest( unordered_set< kmer_t >& kmers, int alpha ) {
  
    initialize_forest();

    unordered_set< kmer_t > visited_mers;
    vector< int > h( n, -1 );  // height for each node below its tree root
    vector< kmer_t > p1( n, 0 ); // p1,p2 needed to tell when to store a tree root
    vector< kmer_t > p2( n, 0 ); //
    vector< kmer_t > p( n, 0 );  // parent in BFS
    
    while (visited_mers.size() != n ) { //there are still connected components to traverse
      kmer_t root = *kmers.begin(); //this is why we need to delete elements from the kmers set
      //update visited_mers, kmers with root
      update_kmer_sets( visited_mers, kmers, root );
      store( root );
      u_int64_t r = f( root );
      p1[ r ] = root;
      p2[ r ] = root;
      h[ r ] = 0;

      //BFS queue 
      queue< kmer_t > Q;
      Q.push( root );

      //conduct BFS
      vector< kmer_t > neis;
      vector< nei_bits > neis_bits;
      while (!Q.empty()) {
	kmer_t c = Q.front();
	Q.pop();
	get_neighbors( c, neis, neis_bits );

	// Debugging code
	// cout << "Neighbors of: ";
	// print_kmer(c, k, cout);
	// cout << endl;
	// for (unsigned ii = 0; ii < neis.size(); ++ii) {
	//   print_kmer( neis[ii], k, cout);
	//   cout << endl;
	// }
	
	for (unsigned ii = 0; ii < neis.size(); ++ii) {
	  kmer_t m = neis[ii]; //this is 'n' in the pseudocode
	  if ( visited_mers.find( m ) == visited_mers.end() ) {
	    //haven't visited m yet
	    Q.push(m);
	    update_kmer_sets( visited_mers, kmers, m );
	    u_int64_t f_m = f(m); //save these values so we don't have to recompute all the time
	    u_int64_t f_c = f(c);
	    p[ f_m ] = c;

	    //set the parent bit field correctly; want c to be parent of m
	    unsigned letter;
	    if (neis_bits[ii].bits[2]) {
	      //m is an OUT neighbor of c
	      //therefore C is in an IN-neighbor of m
	      myForest.parents[ f_m ].parent.bits[2] = 0;
	      //the character is determined by c's first letter
	      letter = access_kmer( c, k, 0 );
	    } else {
	      //m is an IN neighbor of c
	      //therefore c is in an OUT-neighbor of m
	      myForest.parents[ f_m ].parent.bits[2] = 1;
	      //the character is determined by c's last letter
	      letter = access_kmer( c, k, k - 1 );
	    }

	    //map this value to bits
	    switch (letter) {
	    case 0:
	      myForest.parents[ f_m ].parent.bits[0] = 0;
	      myForest.parents[ f_m ].parent.bits[1] = 0;
	      break;
	    case 1:
	      myForest.parents[ f_m ].parent.bits[0] = 0;
	      myForest.parents[ f_m ].parent.bits[1] = 1;
	      break;
	    case 2:
	      myForest.parents[ f_m ].parent.bits[0] = 1;
	      myForest.parents[ f_m ].parent.bits[1] = 0;
	      break;
	    case 3:
	      myForest.parents[ f_m ].parent.bits[0] = 1;
	      myForest.parents[ f_m ].parent.bits[1] = 1;
	      break;
	    }
	    //	    cout << f_m << ' ' <<  myForest.parents[ f_m ].parent << ' ' << neis_bits[ii] << endl;
	    
	    h[ f_m ] = h[ f_c ] + 1;
	    int height_m = h[ f_m ];
	    if (height_m <= alpha) {
	      p1[ f_m ] = p1[ f_c ];
	      p2[ f_m ] = p2[ f_c ];	      
	    }
	    if ( (alpha < height_m) && (height_m <= 2*alpha)) {
	      store( p1[ f_c ] );
	      p1[ f_m ] = p1[ f_c ];
	      p2[ f_m ] = p1[ f_c ];	      
	    }
	    if ( height_m == (2 * alpha + 1) ) {
	      h[ f_m ] = 0;
	      p1[ f_m ] = m;
	      p2[ f_m ] = p1[ f_c ];
	    }
	  }
	}
      }
    }
  }

  /*
   * Finds all (undirected) neighbors of c, using IN, OUT and kmer representation
   * Returns neighbors in neis
   * and how to get to each corresponding neighbor in nbits
   * Assumes kmers have 0's in unused bits
   */
  
  void get_neighbors( kmer_t c,
		      vector< kmer_t >& neis,
		      vector< nei_bits >& nbits ) {
    neis.clear(); nbits.clear();
    u_int64_t fc = f(c);
    kmer_t d; //for bit operations
    //////IN neighbors
    //copy c for bit operations, shift right by two
    d = c >> 2; 

    //means that in the 0th spot, we now have "00"
    if ( IN[ fc ][ 0 ] ) {
      //have in-neighbor with 'A'
      kmer_t e = d;
      set_kmer( e, k, 0, 'A' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 0;
      nb.bits[1] = 0;
      nb.bits[2] = 0;
      nbits.push_back(nb);
    }
    if ( IN[ fc ][ 1 ] ) {
      //have in-neighbor with 'C'
      kmer_t e = d;
      set_kmer( e, k, 0, 'C' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 0;
      nb.bits[1] = 1;
      nb.bits[2] = 0;
      nbits.push_back(nb);
    }
    if ( IN[ fc ][ 2 ] ) {
      //have in-neighbor with 'G'
      kmer_t e = d;
      set_kmer( e, k, 0, 'G' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 1;
      nb.bits[1] = 0;
      nb.bits[2] = 0;
      nbits.push_back(nb);
    }
    if ( IN[ fc ][ 3 ] ) {
      //have in-neighbor with 'T'
      kmer_t e = d;
      set_kmer( e, k, 0, 'T' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 1;
      nb.bits[1] = 1;
      nb.bits[2] = 0;
      nbits.push_back(nb);
    }
    //OUT
    //similar procedure except need to clear unused bits
    //copy c for bit operations, shift left by two
    d = c << 2; 

    //means that in the kth spot, we now have "00"...
    //but we need to zero the -1th spot
    //clear -1th position
    kmer_t op = 3; //0...011
    op = op << 2*(k);  //11 in i'th spot, zeros elsewhere, no issues even if k=32
    op = ~op;      //00 in i'th spot, ones elsewhere
    d = d & op; //i'th position of mer cleared.
    
    if ( OUT[ fc ][ 0 ] ) {
      //have out-neighbor with 'A'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'A' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 0;
      nb.bits[1] = 0;
      nb.bits[2] = 1;
      nbits.push_back(nb);
    }
    if ( OUT[ fc ][ 1 ] ) {
      //have out-neighbor with 'C'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'C' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 0;
      nb.bits[1] = 1;
      nb.bits[2] = 1;
      nbits.push_back(nb);
    }
    if ( OUT[ fc ][ 2 ] ) {
      //have out-neighbor with 'G'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'G' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 1;
      nb.bits[1] = 0;
      nb.bits[2] = 1;
      nbits.push_back(nb);
    }
    if ( OUT[ fc ][ 3 ] ) {
      //have out-neighbor with 'T'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'T' );
      neis.push_back( e );
      nei_bits nb;
      nb.bits[0] = 1;
      nb.bits[1] = 1;
      nb.bits[2] = 1;
      nbits.push_back(nb);
    }
  }
  
  void store( kmer_t mer ) {
    u_int64_t val = f( mer );
    myForest.stored_mers[ val ] = mer;
    myForest.parents[ val ].is_stored = true;
  }
  
  void update_kmer_sets( unordered_set< kmer_t >& visited,
			 unordered_set< kmer_t >& kmers,
			 kmer_t mer ) {
    kmers.erase( mer );
    visited.insert( mer );
  }

};


#endif

