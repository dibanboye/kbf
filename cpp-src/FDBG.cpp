#ifndef FDBG_class
#define FDBG_class

//#include "hash/HashUtil.cpp"
#include "hash/generate_hash.h"
#include <vector>
#include <string>
#include <unordered_set>
#include <map>

using namespace std;

class Letter; //needed for set_kmer declaration

void set_kmer( kmer_t& mer, unsigned k, unsigned i, Letter& c );
void set_kmer( kmer_t& mer, unsigned k, unsigned i, char c );
kmer_t pushOnFront(kmer_t& orig, Letter& letter, unsigned);
kmer_t pushOnBack(kmer_t& orig, Letter& letter, unsigned);

// Representation of A, C, G, and T in bits 00, 01, 10, 11
class Letter {
  bool bits[2];
  void set_bits( char letter ) {
    switch ( letter ) {
    case 'A':
      bits[0] = 0;
      bits[1] = 0;
      break;
    case 'C':
      bits[0] = 0;
      bits[1] = 1;
      break;
    case 'G':
      bits[0] = 1;
      bits[1] = 0;
      break;
    case 'T':
      bits[0] = 1;
      bits[1] = 1;
      break;
    }

  }
  
public:

  void set( unsigned ii ) {

  }
  

  Letter()  {
    set_bits( 'A' );
  }

  //  Letter(bool first, bool second) {
  //     bits[0] = first;
  //     bits[1] = second;
  //  }

  Letter( unsigned letter ) {
    switch (letter) {
    case 0:
      set_bits( 'A' );
      break;
    case 1:
      set_bits( 'C' );
      break;
    case 2:
      set_bits( 'G' );
      break;
    case 3:
      set_bits( 'T' );
      break;
    }
  }
  
  Letter( char letter ) {
    set_bits( letter );
  }

  
  // Return number 0..3 for this letter
  unsigned getNum() {
    return 2*((int)bits[0]) + ((int)bits[1]);
  }     

  Letter& operator=( const Letter& rhs ) {

    this->bits[0] = rhs.bits[0];
    this->bits[1] = rhs.bits[1];
    return *this;
  }

};


  
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

  ForestNode(bool INorOUT, unsigned lett) : letter( lett ) {

    this->INorOUT = INorOUT;

    this->is_stored = false;
  }

  ForestNode( const ForestNode& rhs ) {
    this->letter = rhs.letter;
    this->INorOUT = rhs.INorOUT;

    this->is_stored = rhs.is_stored;
  }
  
  // Given this node's kmer string, figure out parent's kmer string
  kmer_t getNext(kmer_t& mer, unsigned k) {

    if (this->INorOUT) {
      return pushOnFront(mer, letter, k);
    }
    else {
      return pushOnBack(mer, letter, k);
    }
  }

  ForestNode& operator=( const ForestNode& rhs ) {

    this->letter = rhs.letter;
    this->is_stored = rhs.is_stored;
    this->INorOUT = rhs.INorOUT;
    return *this;
  }

};


class forest {
public:
  map< u_int64_t, kmer_t > stored_mers;
  vector< ForestNode > nodes;

  // Create a forest that has num_nodes nodes
  forest(int num_nodes) : nodes( num_nodes ) {

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
  unsigned alpha;   //each tree in forest is guaranteed to be of height alpha to 3alpha

  double get_size() { //in bites
    u_int64_t ssize = 0;

    //Compute size of IN, OUT
    
    double sizeofrow = sizeof( IN[0] ) + 4 / 8.0 ; //each row takes half a byte + overhead
    double sizeofIN = sizeofrow * n;
    double sizeofINOUT = 2*sizeofIN;

    double forestsize = n * 4 / 8.0 + myForest.stored_mers.size() * 2*k / 8.0; //not too accurate

    double hashsize = f.bphf->totalBitSize() / 8.0;  //size in bytes

    return forestsize + hashsize + sizeofINOUT;
    
  }

  u_int64_t getBitSize() {
    u_int64_t res = 8 * n; //IN,OUT
    res += 4 * n + myForest.stored_mers.size() * 2*k; //not too accurate, adding forest bits
    res += f.bphf->totalBitSize();
    return res;
  }
  
  FDBG( vector< string >& reads, 
	unordered_set<kmer_t>& kmers,
	unsigned n, //number of kmers
	unsigned k, //mer size
	bool b_verify = false, //if true, print summary
	ostream& os = cout
	) : myForest( n )
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
      //       BOOST_LOG_TRIVIAL(debug) << get_kmer_str(*i, k) << " maps to "
      //       << f(*i);
    }

    //initialize IN, OUT to zero (false)
    BOOST_LOG_TRIVIAL(info) << "Initializing IN and OUT.";
    vector< bool > vzero( sigma, false );
    IN.assign( this->n, vzero );
    OUT.assign( this->n, vzero );

    //add edges to IN and OUT using reads
    BOOST_LOG_TRIVIAL(info) << "Adding edges to IN and OUT ...";
    add_edges( reads, b_verify, os );

    BOOST_LOG_TRIVIAL(info) << "Shrinking IN and OUT ...";
    //    for (unsigned i = 0; i < n; ++i) {
    //      for (unsigned j = 0; j < 4; ++j) {
    //	IN[i][j].shrink_to_fit();
    //	OUT[i][j].shrink_to_fit();
    //      }
    //    }
    IN.shrink_to_fit();
    OUT.shrink_to_fit();
    
    //Perform the forest construction
    construct_forest( kmers, k*2 ); //alpha = k * lg(sigma)

    myForest.nodes.shrink_to_fit();
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
	//        BOOST_LOG_TRIVIAL(debug) << "Read in k+1-mer " << kplusone;
	add_edge( kplusone );
	++index1;
      }
    }


   //    BOOST_LOG_TRIVIAL(debug) << "Printing IN ...";
    for (int i = 0; i < IN.size(); i++) {
        if (IN[i].size() == 4) {
	  //           BOOST_LOG_TRIVIAL(debug) << " Node " << i << ": " << IN[i][0] << ", "
	  //              << IN[i][1] << ", "
	  //              << IN[i][2] << ", " << IN[i][3];
        }
        else {
	  BOOST_LOG_TRIVIAL(error) << "IN does not have the correct dimensions.";
        }
    }

    //    BOOST_LOG_TRIVIAL(debug) << "Printing OUT ...";
    for (int i = 0; i < OUT.size(); i++) {
        if (OUT[i].size() == 4) {
	  //           BOOST_LOG_TRIVIAL(debug) << " Node " << i << ": " << OUT[i][0] << ", "
	  //              << OUT[i][1] << ", "
	  //              << OUT[i][2] << ", " << OUT[i][3];
        }
        else {
           BOOST_LOG_TRIVIAL(error) << "OUT does not have the correct dimensions.";
        }
    }

  }

  /*
   * Add an edge dynamically to the data structure.
   * From u to v
   * Updates the forest dynamically.
   * Nothing happens if the k-mers aren't compatible.
   */
  void dynamicAddEdge( const kmer_t& u, const kmer_t& v ) {
    //check if u, v are compatible
    unsigned ui, vi;
    for (unsigned i = 0; i < (k - 1); ++i) {
      ui = access_kmer( u, k, i + 1 );
      vi = access_kmer( v, k, i );
      if (ui != vi)
	return;
    }

    //making it this far means that an edge can be added between them
    //Add the edge to OUT[ f(u) ] and to IN[ f(v) ]
    //if edge was already present, quit
    u_int64_t hashU = f( u );
    u_int64_t hashV = f( v );
    unsigned outIndex = access_kmer( v, k, k - 1 );
    unsigned inIndex = access_kmer( u, k, 0 );
    if ( OUT[ hashU ][ outIndex ] )
      return; // edge already exists

    //making it here means that edge is compatible and edge is not already in graph
    //So: begin logic for adding edge
    OUT[ hashU ][ outIndex ] = 1;
    IN[ hashV ][ inIndex ] = 1;

    //Now, need to update the forest
    //TODO

    return;
  }
  
  //add edge implied by the k+1-mer (between the first k and last k characters).
  void add_edge( string& edge ) {

    // figure out the two kmers
    kmer_t u = 0;
    kmer_t v = 0;
    split_edge( edge, u, v );
    //    BOOST_LOG_TRIVIAL(debug) << "Adding an edge from " << get_kmer_str(u, k)
    //        << " to " << get_kmer_str(v, k);

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
    //    BOOST_LOG_TRIVIAL(debug) << "Detecting membership of " << get_kmer_str(m, this->k);

    // The hash value of our kmer
    // Need to keep track of KRval, so it can be updated
    u_int64_t KR_val = f.generate_KRHash_val( m, k, f.Prime ) ;
    u_int64_t hash = f.perfect_from_KR( KR_val );

    //    BOOST_LOG_TRIVIAL(debug) << "It has hash value " << hash;

    // If it is a real kmer value, it must map to 0..1-n
    if (hash >= n)
      return false;

    // get forest node
    ForestNode fn = myForest.nodes[hash];

    //number of times we have traveled in the tree
    unsigned hopcount = 0;
    bool in;    //true = IN, false = OUT
    unsigned letter; //needed to confirm edge from parent's side
    
    // fn is where we are in the tree. Keep going until we know its kmer value.
    while ( !(fn.is_stored) ) {
      ++hopcount;
      if (hopcount > 3*alpha + 1)
	return false; //we have encountered a cycle in our attempt to travel to a root
	
      //do we progress with IN or OUT
      in = fn.INorOUT;
      if (in) {
	//the parent is an IN-neighbor of m
	//so we need m's last character
	letter = access_kmer( m, k, k - 1 );

      } else {
	//the parent is an OUT-neighbor of m
	//so we need m's first character
	letter = access_kmer( m, k, 0 );
      }

      // deduce the parent's kmer
      m = fn.getNext(m, k);
      
      // get the parent's hash
      if (in) {
	unsigned letter_front_parent = access_kmer( m, k, 0 );
	KR_val = f.update_KRHash_val_IN( KR_val, letter_front_parent, letter );
	hash = f.perfect_from_KR( KR_val );
      } else {
	unsigned letter_back_parent = access_kmer( m, k, k - 1 );
	KR_val = f.update_KRHash_val_OUT( KR_val, letter, letter_back_parent );
	hash = f.perfect_from_KR( KR_val );
      }
      //  	hash = f(m); <----left over from old inefficient hashing
	

      // hash must be in 0...n-1
      if (hash >= n) return false;

      //confirm the edge from parent side
      if (in) {
	if (!(OUT[ hash ][ letter ]))
	  return false;
      } else {
	if (!(IN[ hash ][ letter ]))
	  return false;
      }
      
      fn = myForest.nodes[hash];
    }

    //    BOOST_LOG_TRIVIAL(debug) << "Root " << get_kmer_str(m, k) << " is deduced";

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
    BOOST_LOG_TRIVIAL(info) << "Beginning the forest construction, with alpha = " << alpha;
    this->alpha = alpha;
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

      //      BOOST_LOG_TRIVIAL(debug) << "Building forest from root " + get_kmer_str(root, k);
      //      BOOST_LOG_TRIVIAL(debug) << "Forest size, n " << myForest.nodes.size() << ' ' << n;

      // We have visited root
      move_kmer( kmers, visited_mers, root );

      // The hash value of the root
      u_int64_t r = f( root );

      // See pseudocode for p1,p2
      p1[ r ] = root;
      p2[ r ] = root;

      // Root has height 0
      h[ r ] = 0;

      // Now do a BFS around the root in order to put all nodes into a tree
      queue< kmer_t > Q;
      Q.push( root );

      // Kmers of this nodes's neighbors, and if they are in or out neighbors
      vector< kmer_t > neis;
      vector< bool > v_inorout;
      
      // BFS 
      while (!Q.empty()) {
        // Next k-mer to visit neighbors of
	kmer_t c = Q.front();
	Q.pop();
	//        BOOST_LOG_TRIVIAL(debug) << "Visiting neighbors of node " + get_kmer_str(c, k);

        // Neighbors of c
	get_neighbors( c, neis, v_inorout );

       	ForestNode INc( false, access_kmer( c, k, k - 1 ) ); //the ForestNode for an IN-neighbor of c
	ForestNode OUTc( true, access_kmer( c, k, 0 ) ); //the ForestNode for an OUT-neighbor of c
	
	for (unsigned ii = 0; ii < neis.size(); ++ii) {
	  kmer_t m = neis[ii]; //this is 'n' in the pseudocode
	  //	  BOOST_LOG_TRIVIAL(debug) << "Checking neighbor m: " + get_kmer_str(m, k);
	  if ( visited_mers.find( m ) == visited_mers.end() ) {
	    //haven't visited m yet
	    Q.push(m);
	    move_kmer( kmers, visited_mers, m );
	    //	    BOOST_LOG_TRIVIAL(debug) << "Number of kmers visited: " << visited_mers.size();
	    u_int64_t f_m = f(m); //save these values so we don't have to recompute all the time
	    u_int64_t f_c = f(c);
 	    p[ f_m ] = c;

	    if (v_inorout[ ii ]){
	      //m is IN neighbor of c
	      myForest.nodes[ f_m ] = INc;
	    } else {
	      myForest.nodes[ f_m ] = OUTc;
	    }
	    
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
   * Given a kmer c, get all neighbor kmers (neighbor_kmers) 
   * for neighbor_kmers[ i ], v_inorout[ i ] is true if in-neighbor, false otherwise
   */
  void get_neighbors( kmer_t& c,
		      vector< kmer_t >& neighbor_kmers,
		      vector< bool >& v_inorout ) {

    neighbor_kmers.clear();
    v_inorout.clear();

    //    BOOST_LOG_TRIVIAL(debug) << "Finding neighbors of " + get_kmer_str(c, k);

    // Need hash value of our kmer to look into IN and OUT
    u_int64_t fc = f(c);
    //    BOOST_LOG_TRIVIAL(trace) << "c, f(c): " << get_kmer_str(c, k  ) << ' ' << fc;
    //    BOOST_LOG_TRIVIAL(trace) << "Size of IN,OUT: " << IN[ fc ].size() << ' ' << OUT[ fc ].size();
    for ( unsigned ii = 0; ii < 4; ++ii ) {
      if ( IN[ fc ][ ii ] ) {
	//have in-neighbor with letter
	Letter letter ( ii );
	//	BOOST_LOG_TRIVIAL(trace) << "letter: " << ii << ' ' << letter.getNum();
	// Deduce that kmer by tacking letter at the beginning of c
	kmer_t e = pushOnFront(c, letter, k);
	neighbor_kmers.push_back( e );

	//	BOOST_LOG_TRIVIAL(trace) << "neighbor: " + get_kmer_str(e, k);
	
	v_inorout.push_back( true ); //true=IN
      }
      if ( OUT[ fc ][ ii ] ) {
	//have out-neighbor with letter
	Letter letter ( ii );

	// Deduce that kmer by tacking letter at the end of c
	kmer_t e = pushOnBack(c, letter, k);
	neighbor_kmers.push_back( e );

	v_inorout.push_back( false ); //false=OUT

      }
    }

    //    BOOST_LOG_TRIVIAL(trace) << "Neighbors found: ";
    for (unsigned i = 0; i < neighbor_kmers.size(); ++i) {
      //      BOOST_LOG_TRIVIAL(trace) << i << ' ' << get_kmer_str( neighbor_kmers[i], k ) << " is an "<<	v_inorout[ i ] << " neighbor.";
      
    }
  }
  
  void store( kmer_t mer ) {
    u_int64_t val = f( mer );
    myForest.stored_mers[ val ] = mer;
    myForest.nodes[ val ].is_stored = true;
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
  kmer_t val = c.getNum();

  val = val << 2*(k - i - 1);  //correct bits in i'th spot, zeros elsewhere
  mer = mer | val;
}

/*
 * Accesses the i'th position of a mer of length k 
 * Stores result in c 
 */
void access_kmer( kmer_t& mer, unsigned k, unsigned i, Letter& c ) {
  mer = mer >> 2*(k - i - 1);
  kmer_t mask = static_cast<kmer_t>(3);
  mer = mask & mer;
  c.set( static_cast<unsigned>(mer) );
}

// Push letter onto front of kmer and return kmer
// Going backwards along an edge (the relationship comes from orig's IN).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnFront(kmer_t& orig, Letter& letter, unsigned k) {
  //  BOOST_LOG_TRIVIAL(trace) << "pushOnFront " + get_kmer_str(orig, k) << ' ' << letter.getNum();
  // Get the kmer with back pushed off orig
  kmer_t new_kmer = orig >> 2;

  set_kmer( new_kmer, k, 0, letter);

  return new_kmer;
} 

// Push letter onto back of kmer and return kmer
// Going forward along an edge (the relationship comes from orig's OUT).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnBack(kmer_t& orig, Letter& letter, unsigned k) {
  //  BOOST_LOG_TRIVIAL(trace) << "pushOnBack " + get_kmer_str(orig, k) << ' ' << letter.getNum();
  // Get the kmer with front pushed off orig
  kmer_t new_kmer = orig << 2;
  set_kmer( new_kmer, k, k - 1, letter);
  //but we need to zero the -1th spot
  //clear -1th position
  kmer_t op = 3; //0...011
  op = op << 2*(k);  //11 in -1'st spot, zeros elsewhere, no issues even if k=32
  op = ~op;      //00 in i'th spot, ones elsewhere
  new_kmer = new_kmer & op; //i'th position of mer cleared.

  return new_kmer;
} 

/*
 * Sets the i'th position of a mer of length k as indicated by character c
 * c \in {A,C,G,T}
 */
void set_kmer( kmer_t& mer, unsigned k, unsigned i, char c ) {
  //clear i-th position
  kmer_t op = 3; //0...011
  op = op << 2*(k - i - 1);  //11 in i'th spot, zeros elsewhere
  op = ~op;      //00 in i'th spot, ones elsewhere
  mer = mer & op; //i'th position of mer cleared.

  //set i'th position
  kmer_t val;
  switch( c ) {
    case 'A':
      val = 0;
      break;
    case 'C':
      val = 1;
      break;
    case 'G':
      val = 2;
      break;
    case 'T':
      val = 3;
      break;
  }

  val = val << 2*(k - i - 1);  //correct bits in i'th spot, zeros elsewhere
  mer = mer | val;
}

#endif

