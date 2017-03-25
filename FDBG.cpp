#ifndef FDBG_class
#define FDBG_class

//#include "hash/HashUtil.cpp"
#include "hash/generate_hash.h"
#include <vector>
#include <string>
#include <unordered_set>
#include <map>

using namespace std;

void set_kmer( kmer_t& mer, unsigned k, unsigned i, Letter c );
kmer_t pushOnFront(kmer_t& orig, Letter& letter);
kmer_t pushOnBack(kmer_t& orig, Letter& letter);

// Representation of A, C, G, and T in bits 00, 01, 10, 11
class Letter {

public:
  bool[2] bits;

  Letter() {
     bits[0] = 0;
     bits[1] = 0;
  }

  Letter(bool first, bool second) {
     bits[0] = first;
     bits[1] = second;
  }

  // Return number 0..3 for this letter
  unsigned getNum() {
    return 2*((int)bits[0]) + ((int)bits[1]);
  }     

  Letter& operator=( const Letter& rhs ) {

    this->bits[0] = rhs->bits[0];
    this->bits[1] = rhs->bits[1];
    return *this;
  }

}

// Class that represents a node in the forest
// The only information encoded is whether it was reached during
// the forest creation via IN or OUT (which determines how to recover
// the kmer) and by what letter it is different (from the node that reached it)
class ForestNode {
public:
  // Get to parent by IN (true) or OUT (false)
  bool INorOUT;

  // What letter must be added to get to parent
  Letter letter;

  bool is_stored; // whether the kmer for this is stored or not

  ForestNode() {
    this->is_stored = false;
  }

  ForestNode(bool INorOUT, Letter letter) {

    this->INorOUT = INorOUT;
    this->letter = letter;

    this->is_stored = false;
  }

  ForestNode( const ForestNode& rhs ) {
    this->letter = rhs->letter;
    this->INorOUT = rhs->INorOUT;

    this.is_stored = rhs.is_stored;
  }
  
  // Given this node's kmer string, figure out parent's kmer string
  kmer_t getNext(kmer_t& mer) {

    if (this->INorOUT) {
      return pushOnBack(mer, letter);
    }
    else {
      return pushOnFront(mer, letter);
    }
  }

  ForestNode& operator=( const ForestNode& rhs ) {

    this->letter = rhs->letter;
    this->is_stored = rhs->is_stored;
    this->INorOUT = rhs->INorOUT;
    return *this;
  }

};


class forest {
public:
  map< u_int64_t, kmer_t > stored_mers;
  vector< ForestNode > nodes;

  // Create a forest that has num_nodes nodes
  forest(int num_nodes) {
    nodes = vector(num_nodes, ForestNode()); 
  }
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
    IN.assign( this.n, vzero );
    OUT.assign( this.n, vzero );

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

  // Take a k+1-mer and split into beginning and end k-mers
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

  // Given a kmer, decide if it is one in our graph
  bool detect_membership( kmer_t m ) {

    BOOST_LOG_TRIVIAL(debug) << "Detecting membership of " << get_kmer_str(m, this.k);

    // The hash value of our kmer
    u_int64_t hash = this.f(m);

    BOOST_LOG_TRIVIAL(debug) << "It has hash value " << hash;

    // If it is a real kmer value, it must map to 0..1-n
    if (hash >= n)
      return false;

    // get forest node
    ForestNode fn = myForest.nodes[hash];

    // fn is where we are in the tree. Keep going until we know its kmer value.
    while ( !(fn.is_stored) ) {

      // deduce the parent's kmer
      m = fn.getNext(m); 

      // get the parent's hash
      hash = this.f(m);

      BOOST_LOG_TRIVIAL(debug) << "... which goes to hash value " << hash;

      // hash must be in 0...n-1
      if (hash >= n) return false;

      fn = this.myForest.nodes[hash];
    }

    BOOST_LOG_TRIVIAL(debug) << "Root " << get_kmer_str(m, this.k) << " is deduced";

    // now we have a forest node that we have the kmer of stored
    // So we just have to test if it is accurate or not
    if (m == myForest.stored_mers[hash])
      return true;
    else
      return false;
  }

  /*
   * Initially constructs the forest
   * Requires IN, OUT, f to be already constructed
   * alpha is the tree height parameter
   * Requires non-used bits of kmers to be zero
   * WARNING: kmers will be modified
   */
  void construct_forest( unordered_set< kmer_t >& kmers, int alpha ) {
 
    // construct a forest with n empty nodes 
    this.myForest = forest(this.n);

    // kmers that we have looked at for the forest construction
    unordered_set< kmer_t > visited_mers;

    vector< int > h( n, -1 );  // height for each node below its tree root
    vector< kmer_t > p1( n, 0 ); // p1,p2 needed to tell when to store a tree root
    vector< kmer_t > p2( n, 0 ); //
    vector< kmer_t > p( n, 0 );  // parent in BFS

    // Do a BFS through each of the UNDIRECTED graph components
    // Keep going until we have looked at all kmers
    // This will only have one loop unless the graph is not connected
    while (visited_mers.size() != n ) {

      // Pick initial root, will be stored in our forest
      kmer_t root = *kmers.begin();
      store( root );

      BOOST_LOG_TRIVIAL(debug) << "Building forest from root " + get_kmer_str(root, this.k);

      // We have visited root
      move_kmer( kmers, visited_mers, root );

      // The hash value of the root
      u_int64_t r = f( root );

      // It is the root, both parents are itself
      p1[ r ] = root;
      p2[ r ] = root;

      // Root has height 0
      h[ r ] = 0;

      // Now do a BFS around the root in order to put all nodes into a tree
      queue< kmer_t > Q;
      Q.push( root );

      // Kmers of this nodes's neighbors
      vector< kmer_t > neighbor_kmers;

      // The neighbors in the form we will need for our forest
      vector< ForestNode > neighbor_nodes;

      // BFS 
      while (!Q.empty()) {
        // Next k-mer to visit neighbors of
	kmer_t c = Q.front();
	Q.pop();
        BOOST_LOG_TRIVIAL(debug) << "Visiting neighbors of node " + get_kmer_str(c, this.k);

        // Neighbors of c
	get_neighbors( c, neighbor_kmers, neighbor_nodes );
	
	for (unsigned ii = 0; ii < neis.size(); ++ii) {
	  kmer_t m = neis[ii]; //this is 'n' in the pseudocode
	  if ( visited_mers.find( m ) == visited_mers.end() ) {
	    //haven't visited m yet
	    Q.push(m);
	    move_kmer( kmers, visited_mers, m );
	    u_int64_t f_m = f(m); //save these values so we don't have to recompute all the time
	    u_int64_t f_c = f(c);
	    p[ f_m ] = c;

	    //set the parent bit field correctly; want c to be parent of m
	    unsigned letter;
	    if (neis_bits[ii].bits[2]) {
	      //m is an OUT neighbor of c
	      //therefore C is in an IN-neighbor of m
	      myForest.nodes[ f_m ].parent.bits[2] = 0;
	      //the character is determined by c's first letter
	      letter = access_kmer( c, k, 0 );
	    } else {
	      //m is an IN neighbor of c
	      //therefore c is in an OUT-neighbor of m
	      myForest.nodes[ f_m ].parent.bits[2] = 1;
	      //the character is determined by c's last letter
	      letter = access_kmer( c, k, k - 1 );
	    }

	    //map this value to bits
	    switch (letter) {
	    case 0:
	      myForest.nodes[ f_m ].parent.bits[0] = 0;
	      myForest.nodes[ f_m ].parent.bits[1] = 0;
	      break;
	    case 1:
	      myForest.nodes[ f_m ].parent.bits[0] = 0;
	      myForest.nodes[ f_m ].parent.bits[1] = 1;
	      break;
	    case 2:
	      myForest.nodes[ f_m ].parent.bits[0] = 1;
	      myForest.nodes[ f_m ].parent.bits[1] = 0;
	      break;
	    case 3:
	      myForest.nodes[ f_m ].parent.bits[0] = 1;
	      myForest.nodes[ f_m ].parent.bits[1] = 1;
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

  /**
   * Given a kmer c, get all neighbor kmers (neighbor_kmers) and their ForestNode
   * representation (where c is the parent). 
   */
  void get_neighbors( kmer_t& c,
		      vector< kmer_t >& neighbor_kmers,
		      vector< ForestNode >& neighbor_nodes ) {

    neis.clear();
    nbits.clear();

    BOOST_LOG_TRIVIAL(debug) << "Finding neighbors of " + get_kmer_str(c, this.k);

    // Need hash value of our kmer to look into IN and OUT
    u_int64_t fc = f(c);

    // First, get neighbors from the IN matrix

    // Neighbors coming in do not have the last letter of the kmer (c)

    if ( IN[ fc ][ 0 ] ) {
      //have in-neighbor with 'A'
      Letter letter (0, 0);

      // Deduce that kmer by tacking 'A' at the beginning of d
      kmer_t e = pushOnFront(c, letter);
      neighbor_kmers.push_back( e );

      // Add to neighbors
      ForestNode nb (0, letter);
      neighbor_nodes.push_back(nb);
    }
    if ( IN[ fc ][ 1 ] ) {
      //have in-neighbor with 'C'
      kmer_t e = d;
      set_kmer( e, k, 0, 'C' );
      neis.push_back( e );
      ForestNode nb (0, 1, 0);
      nbits.push_back(nb);
    }
    if ( IN[ fc ][ 2 ] ) {
      //have in-neighbor with 'G'
      kmer_t e = d;
      set_kmer( e, k, 0, 'G' );
      neis.push_back( e );
      ForestNode nb (1, 0, 0);
      nbits.push_back(nb);
    }
    if ( IN[ fc ][ 3 ] ) {
      //have in-neighbor with 'T'
      kmer_t e = d;
      set_kmer( e, k, 0, 'T' );
      neis.push_back( e );
      ForestNode nb (1, 1, 0);
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
      ForestNode nb (0, 0, 1);
      nbits.push_back(nb);
    }
    if ( OUT[ fc ][ 1 ] ) {
      //have out-neighbor with 'C'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'C' );
      neis.push_back( e );
      ForestNode nb (0, 1, 1);
      nbits.push_back(nb);
    }
    if ( OUT[ fc ][ 2 ] ) {
      //have out-neighbor with 'G'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'G' );
      neis.push_back( e );
      ForestNode nb (1, 0, 1);
      nbits.push_back(nb);
    }
    if ( OUT[ fc ][ 3 ] ) {
      //have out-neighbor with 'T'
      kmer_t e = d;
      set_kmer( e, k, k - 1, 'T' );
      neis.push_back( e );
      ForestNode nb (1, 1, 1);
      nbits.push_back(nb);
    }
  }
  
  void store( kmer_t mer ) {
    u_int64_t val = f( mer );
    myForest.stored_mers[ val ] = mer;
    myForest.Nodes[ val ].is_stored = true;
  }
  
  // Move mer from kmers and into visited
  void move_kmer(unordered_set< kmer_t >& kmers, unordered_set< kmer_t >& visited,
			 kmer_t mer ) {
    kmers.erase( mer );
    visited.insert( mer );
  }

};

/*
 * Sets the i'th position of a mer of length k as indicated by character c
 * c \in {A,C,G,T}
 */
void set_kmer( kmer_t& mer, unsigned k, unsigned i, Letter& c ) {
  //clear i-th position
  kmer_t op = 3; //0...011
  op = op << 2*(k - i - 1);  //11 in i'th spot, zeros elsewhere
  op = ~op;      //00 in i'th spot, ones elsewhere
  mer = mer & op; //i'th position of mer cleared.

  //set i'th position
  kmer_t val = c->getNum();

  val = val << 2*(k - i - 1);  //correct bits in i'th spot, zeros elsewhere
  mer = mer | val;
}

// Push letter onto front of kmer and return kmer
// Going backwards along an edge (the relationship comes from orig's IN).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnFront(kmer_t& orig, Letter& letter) {

  // Get the kmer with back pushed off orig
  kmer_t new_kmer = orig >> 2;

  set_kmer( new_kmer, k, 0, letter);

} 

// Push letter onto back of kmer and return kmer
// Going forward along an edge (the relationship comes from orig's OUT).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnBack(kmer_t& orig, Letter& letter) {

  // Get the kmer with back pushed off orig
  kmer_t new_kmer = orig << 2;

  set_kmer( new_kmer, k - 1, 0, letter);

} 



#endif

