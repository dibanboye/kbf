#ifndef FDBG_class
#define FDBG_class

//#include "hash/HashUtil.cpp"
#include "hash/generate_hash.h"
#include <vector>
#include <string>
#include <unordered_set>
#include <map>
#include <chrono>
#include <queue>
#include "BitArray.cpp"
#include "formatutil.cpp"

using namespace std;

class Letter; //needed for set_kmer declaration

kmer_t getMaxVal(unsigned);
void push_last_letter(const kmer_t&, kmer_t&);
void remove_front_letter(const kmer_t&, kmer_t&, const unsigned&);
void set_kmer( kmer_t& mer, unsigned k, unsigned i, Letter& c );
kmer_t pushOnFront(const kmer_t& orig, Letter& letter, unsigned);
kmer_t pushOnBack(const kmer_t& orig, Letter& letter, unsigned);
static vector< uint64_t > DEFAULT_VECTOR;

// Representation of A, C, G, and T in bits 00, 01, 10, 11
class Letter {

public:

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
  
  void set( unsigned ii ) {

  }
  

  Letter()  {
    set_bits( 'A' );
  }

  Letter(bool a, bool b)  {

     this->bits[0] = a;
     this->bits[1] = b;
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


/**
 * Class that represents the forest
 * each group of 4 entries in the vector is a single node
 */
class Forest {

  public:
    /**
     * Each four pieces of data describe a node in the forest
     * One spot for whether it is stored as a root, one for whether
     * its parent can be reached by IN or not, and two for letter to
     * reach its parent by
     * The ith node is the kmer that hashes to i
     */
    BitArray bitarray;

    // The number of nodes
    u_int64_t n;

    // The kmers of our roots
    map<u_int64_t, kmer_t> roots;

    /**Constructor*/
    Forest(u_int64_t n): bitarray(4*n) {

      this->n = n;
    }

      /**Constructor*/
  Forest() {
    //default constructor
  }

  void allocate(u_int64_t n) {
    this->n = n;
    bitarray.allocate( 4*n );
  }
  
    /**
     * Set value of the ith node
     * Not stored
     */
    void setNode(const u_int64_t& i, bool IN, const Letter& l) {

       u_int64_t index = this->nodeIndex(i);

       this->bitarray.set(index, false);
       this->bitarray.set(index + 1, IN); 

       this->bitarray.set(index + 2, l.bits[0]);
       this->bitarray.set(index + 3, l.bits[1]);

    }

    /**
     * Set the letter of the ith node
     */
    void setLetter(const u_int64_t& i, const Letter& l) {

       u_int64_t index = this->nodeIndex(i);

       this->bitarray.set(index + 2, l.bits[0]);
       this->bitarray.set(index + 3, l.bits[1]);

    }

    /**
     * Get the letter to the parent of the ith node in the form of an unsigned
     */
    unsigned getLetter(const u_int64_t& i) {

       unsigned l;

       u_int64_t index = this->nodeIndex(i);

       l += 2*this->bitarray.get(index + 2);
       l += this->bitarray.get(index + 3);

       return l;
    }

    // Return index in data of the beginning of the ith node
    u_int64_t nodeIndex(const u_int64_t& i) {
      return 4*i;
    }

    /**
     * Get the data for the ith node, put in data vector
     */
    void getNodeData(const u_int64_t& i, vector<bool>& data) {

       u_int64_t index = nodeIndex(i);

       data.push_back(this->bitarray.get(index));
       data.push_back(this->bitarray.get(index + 1));
       data.push_back(this->bitarray.get(index + 2));
       data.push_back(this->bitarray.get(index + 3));

    }

    /**
     * Store the kmer of the ith node
     */
    void storeNode(const u_int64_t& i, const kmer_t& str) {

       this->roots[i] = str;

       this->bitarray.set(nodeIndex(i), true);
    }

    /**
     * Unstore the kmer of the ith node
     */
    void unstoreNode(const u_int64_t& i) {

       this->roots.erase(i);

       this->bitarray.set(nodeIndex(i), false);
    }

    // Returns whether the ith node is stored
    bool isStored(u_int64_t i) {
       return this->bitarray.get(nodeIndex(i));
    }

    // Whether the ith node has IN
    bool parent_in_IN(u_int64_t i) {
       return this->bitarray.get(nodeIndex(i) + 1);
    }

    // Set whether the ith node has IN
    bool set_parent_in_IN(const u_int64_t& i, bool in) {
       this->bitarray.set(nodeIndex(i) + 1, in);
    }

    // Deduce the parent of the ith node's kmer given the ith node's
    // kmer mer and the length of the kmers k
    kmer_t getNext(u_int64_t i, const kmer_t& mer, unsigned k) {

      u_int64_t index = nodeIndex(i);
      Letter l (this->bitarray.get(index + 2), this->bitarray.get(index + 3));

      if (this->bitarray.get(index + 1)) {
        // reach via IN
        return pushOnFront(mer, l, k);
      }
      else {
        return pushOnBack(mer, l, k);
      }
    }

    u_int64_t getBitSize() {
        return this->bitarray.total_bit_size() + this->roots.size()*8*sizeof(kmer_t);
    }

  void save( ostream& of ) {
    of.write ( (char*) &n, sizeof( u_int64_t ) );
    unsigned number_stored = roots.size();
    of.write ( (char*) &number_stored, sizeof( unsigned ) );
    for (auto i = roots.begin(); i != roots.end(); ++i) {
      u_int64_t hash = i -> first;
      kmer_t label = i -> second;
      of.write ( (char*) &hash, sizeof( u_int64_t ) );
      of.write ( (char*) &label, sizeof( kmer_t ) );
    }
    bitarray.save( of );
  }

  void load( istream& of ) {
    of.read ( (char*) &n, sizeof( u_int64_t ) );
    unsigned number_stored;
    of.read ( (char*) &number_stored, sizeof( unsigned ) );
    roots.clear();
    for (unsigned i = 0; i < number_stored; ++i) {
      u_int64_t hash;
      kmer_t label;
      of.read( (char*) &hash, sizeof( u_int64_t ) );
      of.read( (char*) &label, sizeof( kmer_t ) );
      roots[ hash ] = label;
    }
    bitarray.load( of );
  }

  
};

/**
 * Class meant to represent IN or OUT
 * Has n rows, each row has 4 columns
 */
class INorOUT {

  private:
    BitArray bitarray;
    u_int64_t n; // number of rows

  public:

  INorOUT(u_int64_t n) : bitarray (4*n) {
    this->n = n;
  }

  INorOUT()  {
    //default constructor
  }

  void  allocate(u_int64_t n) {
    this->n = n;
    bitarray.allocate( 4*n );
  }
  
    // Set row, col value
    void set(unsigned row, unsigned col, bool val) {

        unsigned index = 4*row + col;
        this->bitarray.set(index, val);

    }

    // Get row, col value
    bool get(unsigned row, unsigned col) {

        unsigned index = 4*row + col;
        return this->bitarray.get(index);

    }

  size_t getBitSize() {
        return this->bitarray.total_bit_size() + 8*sizeof(u_int64_t);
    }

  void save( ostream& of ) {
    of.write ( (char*) &n, sizeof( u_int64_t ) );
    bitarray.save( of );
  }

  void load( istream& of ) {
    of.read ( (char*) &n, sizeof( u_int64_t ) );
    //    cerr << "INorOUT n:" << n << endl;
    bitarray.load( of );
  }
};


/**
 * This class represents the fully dynamic De Bruijn graph
 */
class FDBG {

public:

  INorOUT IN;
  INorOUT OUT;
  unsigned sigma; //alphabet-size. For now, only 4 is supported
  u_int64_t n; //number of nodes in graph
  unsigned k; //length of each mer (string in alphabet)
  generate_hash f; //hash function that takes each kmer to 1..n
  Forest fo; // the forest NEW
  unsigned alpha;   //each tree in forest is guaranteed to be of height alpha to 3alpha
  double construction_time;

  void save( ostream& of ) {
    BOOST_LOG_TRIVIAL(debug) << "Saving small variables to file...";
    of.write ( (char*) &n, sizeof( u_int64_t ) );
    of.write ( (char*) &sigma, sizeof( unsigned ) );
    of.write ( (char*) &k, sizeof( unsigned ) );
    of.write ( (char*) &alpha, sizeof( unsigned ) );
    of.write ( (char*) &construction_time, sizeof( double ) );

    BOOST_LOG_TRIVIAL(debug) << "Saving IN to file...";
    IN.save( of );
    BOOST_LOG_TRIVIAL(debug) << "Saving OUT to file...";
    OUT.save( of );
    BOOST_LOG_TRIVIAL(debug) << "Saving Forest to file...";
    fo.save( of );
    BOOST_LOG_TRIVIAL(debug) << "Saving hash to file...";
    f.save( of );
  }

  void load( istream& of ) {
    BOOST_LOG_TRIVIAL(debug) << "Loading variables from file...";
    
    of.read ( (char*) &n, sizeof( u_int64_t ) );
    of.read ( (char*) &sigma, sizeof( unsigned ) );
    of.read ( (char*) &k, sizeof( unsigned ) );
    of.read ( (char*) &alpha, sizeof( unsigned ) );
    of.read ( (char*) &construction_time, sizeof( double ) );
    //    BOOST_LOG_TRIVIAL(debug) << n << ' ' << sigma << ' ' << k << ' ' << alpha;
    
    BOOST_LOG_TRIVIAL(debug) << "Loading IN from file...";
    IN.load( of );
    BOOST_LOG_TRIVIAL(debug) << "Loading OUT from file...";
    OUT.load( of );
    BOOST_LOG_TRIVIAL(debug) << "Loading fo from file...";
    fo.load( of );

    //    printForest();
      
    BOOST_LOG_TRIVIAL(debug) << "Loading hash from file...";
    f.load( of );
  }

  
  /**
   * The number of bits that our data should be using
   */
  u_int64_t estimateBitSize() {
    u_int64_t res = 8 * n; //IN,OUT
    res += 4 * n + fo.roots.size() * 2*k; //not too accurate, adding forest bits
    res += f.bphf->totalBitSize();
    return res;
  }

  size_t bitSize() {

    // IN and OUT
    u_int64_t res = this->IN.getBitSize() + this->OUT.getBitSize();

    // forest
    res += this->fo.getBitSize();

    // mphf
    res += f.bphf->totalBitSize();

    return res;
  }

  FDBG() {
    //default constructor
  }

  void build(
	     unordered_set<kmer_t>& kmers,
	     unordered_set<kmer_t>& edgemers,
	     u_int64_t n, //number of kmers
	     unsigned k ) { //mer size
    auto start = std::chrono::system_clock::now();
    sigma = 4;
    
    this->n = n;
    this->k = k;

    IN.allocate( n );
    OUT.allocate( n );
    fo.allocate( n );
    
    //construct hash function f
    f.construct_hash_function( kmers, n, k );

    //printHashFunction(kmers);

    //initialize IN, OUT to zero (false)
    BOOST_LOG_TRIVIAL(info) << "Initializing IN and OUT.";

    //add edges to IN and OUT
    BOOST_LOG_TRIVIAL(info) << "Adding edges to IN and OUT ...";

    add_edges(edgemers);

    //For debugging, print IN and OUT given kmers
    //printINandOUT(kmers);
    
    //Perform the forest construction
    construct_forest( kmers, k*2 ); //alpha = k * lg(sigma)

    //cout << "SIZE: " << this->fo.roots.size() << endl;

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    construction_time = elapsed_seconds.count();
  }

  
  FDBG( unordered_set<kmer_t>& kmers,
        unordered_set<kmer_t>& edgemers,
	u_int64_t n, //number of kmers
	unsigned k, //mer size
	bool b_verify = false, //if true, print summary
	ostream& os = cout
	) : IN (n), OUT (n), fo (n)
  {
    build( kmers, edgemers, n, k );
  }

  /**
   * Add edges to IN and OUT using K+1-mers
   */
  void add_edges(unordered_set<kmer_t>& edges) {

    unordered_set<kmer_t>::iterator i;
    for (i = edges.begin(); i!= edges.end(); ++i) {
        add_edge(*i);
    }

  }

  void add_edge(const kmer_t& edge ) {

    // figure out the two kmers
    kmer_t u = 0;
    kmer_t v = 0;

    split_edge( edge, u, v );

    // get which column in sigma to put in (corresponds to which letter is the first/last)
    // number 0..3 represent each alphabet letter
    unsigned first, last;
    first = access_kmer( u, k, 0 );
    last = access_kmer( v, k, k - 1 );

    u_int64_t hash_u = f(u);
    u_int64_t hash_v = f(v);

    if (hash_u >= this->n) {
      BOOST_LOG_TRIVIAL(error) << "ERROR: Edge " << get_kmer_str(edge, this->k + 1)
        <<  " has produced kmer "
        << get_kmer_str(u, this->k)
        << "("<< u << ")"
        << " with invalid hash function value.";
    }

    if (hash_v >= this->n) {
      BOOST_LOG_TRIVIAL(error) << "Edge " << get_kmer_str(edge, this->k + 1)
        <<  " has produced kmer "
        << get_kmer_str(v, this->k)
        << "("<< v << ")"
        << " with invalid hash function value.";
    }

    this->OUT.set(hash_u, last, true);
    this->IN.set(hash_v, first, true);
 
  }

  // Take a k+1-mer and split into beginning and end k-mers
  // u is the beginning, v is the end
  void split_edge( const kmer_t& edge, kmer_t& u, kmer_t& v ) {

     push_last_letter(edge, u);

     remove_front_letter(edge, v, this->k);

  }

  /**
   * Computes the average tree height and how many are above/below max/min
   */
  void getTreeData(unsigned& num_trees, double& avg_height,
     unsigned& num_above, unsigned& num_below) {

     num_trees = this->fo.roots.size();

     avg_height = 0;
     num_above = 0;
     num_below = 0;
     double weight = 1.0/(this->fo.roots.size());

     //BOOST_LOG_TRIVIAL(debug) << "There are " << this->fo.roots.size()
     //   << " trees and so the weight is " << weight;

     //BOOST_LOG_TRIVIAL(debug) << "Trees should have heights between "
     //   << this->alpha << " and " << 3*this->alpha;

     unsigned height;
     kmer_t root;

     map<kmer_t, unsigned> heights;
     vector<kmer_t> sorted_kmers;

     // go through each tree

     map<u_int64_t, kmer_t>::iterator iter;
     for (iter = this->fo.roots.begin(); iter != this->fo.roots.end(); ++iter) {

        root = iter->second;

        //BOOST_LOG_TRIVIAL(debug) << "Looking at tree with root "
        //   << get_kmer_str(root, this->k);

        height = getTreeHeightRoot(root, heights, sorted_kmers);

        avg_height += (height*weight);

        if (height < this->alpha) {
           num_below++;
        }
        else if (height > 3*this->alpha + 1) {
           num_above++;
        }

     }

  }


  /**
   * Given a node, get the other nodes for all edges
   * that don't exist. For example, if we have ATTC and there
   * does not exist an edge ATTCA, we return TTCA. The nodes may or
   * may not exist in the graph.
   * This is for generating random edges.
   */
  void getNonExistentEdges(const kmer_t& node, vector<kmer_t>& not_neighbors_in,
     vector<kmer_t>& not_neighbors_out) {

     not_neighbors_in.clear();
     not_neighbors_out.clear();

     u_int64_t hash = this->f(node);

     Letter l;
     kmer_t neighbor;
     for (unsigned i = 0; i < 4; i++) {

        if (this->IN.get(hash,i) == 0) {
           l = Letter (i);
           neighbor = pushOnFront(node, l, this->k);
           //BOOST_LOG_TRIVIAL(debug) << "Non existent edge from "
           //   << get_kmer_str(neighbor, this->k);
           not_neighbors_in.push_back(neighbor);
        }

        if (this->OUT.get(hash,i) == 0) {
           l = Letter (i);
           neighbor = pushOnBack(node, l, this->k);
           //BOOST_LOG_TRIVIAL(debug) << "Non existent edge to "
           //   << get_kmer_str(neighbor, this->k);
           not_neighbors_out.push_back(neighbor);
        }

     }

  }

  /*
   * Add an edge dynamically to the data structure.
   * From u to v
   * Updates the forest dynamically.
   * Nothing happens if the k-mers aren't compatible.
   * Returns bool of whether an edge is actually added or not
   */
  bool dynamicAddEdge( const kmer_t& u, const kmer_t& v ) {

    //BOOST_LOG_TRIVIAL(debug) << "Adding an edge from " << get_kmer_str(u, this->k)
    //   << " to " << get_kmer_str(v, this->k) << "...";

    //check if u, v are compatible
    //unsigned ui, vi;
    //for (unsigned i = 0; i < (this->k - 1); ++i) {
    //  ui = access_kmer( u, this->k, i + 1 );
    //  vi = access_kmer( v, this->k, i );
      //BOOST_LOG_TRIVIAL(debug) << "Checking the " << (i+1) << "/" << i << " spots";
    //  if (ui != vi) {
     //   BOOST_LOG_TRIVIAL(debug) << "Cannot add an edge because there is not "
     //      << " k-1 length overlap.";
//	return false;
  //    }
   // }

    // For testing
    //if (!detect_membership(u)) {
       //BOOST_LOG_TRIVIAL(debug) << "The node " << get_kmer_str(u, this->k)
       //   << " is not in the graph.";
    //   return false;
    //}

    //if (!detect_membership(v)) {
       //BOOST_LOG_TRIVIAL(debug) << "The node " << get_kmer_str(v, this->k)
       //   << " is not in the graph.";
    //   return false;
    //}
    

    //making it this far means that an edge can be added between them
    //Add the edge to OUT[ f(u) ] and to IN[ f(v) ]
    //if edge was already present, quit
    u_int64_t hashU = f( u );
    u_int64_t hashV = f( v );
    unsigned outIndex = access_kmer( v, k, k - 1 );
    unsigned inIndex = access_kmer( u, k, 0 );
    if ( OUT.get(hashU, outIndex) ) {
      BOOST_LOG_TRIVIAL(debug) << "This edge already exists";
      return false; // edge already exists
    }

    //making it here means that edge is compatible and edge is not already in graph
    //So: begin logic for adding edge
    OUT.set(hashU, outIndex, true);
    IN.set(hashV, inIndex, true);

    // Now, need to update the forest
    // Basically, there is the potential to merge too small trees

    //BOOST_LOG_TRIVIAL(debug) << "Finding the trees that these two edges are in ...";

    // Find the heights of the trees these two are in, and their roots
    kmer_t root_u;
    u_int64_t root_u_hash;
    unsigned height_u = getRoot(u, root_u, root_u_hash);
    kmer_t root_v;
    u_int64_t root_v_hash;
    unsigned height_v = getRoot(v, root_v, root_v_hash);

    //BOOST_LOG_TRIVIAL(debug) << "root " << get_kmer_str(root_u, this->k)
    //   << " and root " << get_kmer_str(root_v, this->k);

    map<kmer_t, unsigned> u_heights;
    vector<kmer_t> u_sorted_kmers;
    map<kmer_t, unsigned> v_heights;
    vector<kmer_t> v_sorted_kmers;
    unsigned treeheight_u = getTreeHeightRoot(root_u, u_heights, u_sorted_kmers);
    unsigned treeheight_v = getTreeHeightRoot(root_v, v_heights, v_sorted_kmers);

    //map<kmer_t, unsigned>::iterator iter;
    //for (iter = u_heights.begin(); iter != u_heights.end(); ++iter) {
    //   BOOST_LOG_TRIVIAL(debug) << "The height of kmer " << get_kmer_str(iter->first, this->k)
    //      << " is " << iter->second;
    //}

    //BOOST_LOG_TRIVIAL(debug) << "One is in a tree of height " << height_u
    //   << " with root " << get_kmer_str(root_u, this->k)
    //   << " and the other a tree of height " << height_v
    //   << " with root " << get_kmer_str(root_v, this->k);

    if (root_u == root_v) {
       // No trees to merge
       //BOOST_LOG_TRIVIAL(debug) << "Both trees already have the same root. No merging.";
       return true;
    }

    // The letter u has, v doesn't
    Letter u_letter (inIndex);
    // similar
    Letter v_letter (outIndex);

    // One of the trees is too small, try and make big enough ones
    if ((treeheight_u < this->alpha) || (treeheight_v < this->alpha)) {

       mergeTrees(u, v, u_heights, v_heights, u_sorted_kmers, v_sorted_kmers,
          root_u_hash, root_v_hash, hashU, hashV, u_letter, v_letter);

    }

    return true;
  }

   /*
    * Utility function for dynamicRemoveEdge
    * tree_mers is all the mers from a
    * tree of height < alpha
    * Looks for adjacent tree to merge with
    */
   void removalFixTree( vector< kmer_t >& tree_mers,
			vector< uint64_t >& tree_mers_hash,
			map< kmer_t, unsigned >& tree_heights ) {
      //make unordered_sets from the tree hash values
      unordered_set< uint64_t > uo_hash ( tree_mers_hash.begin(), tree_mers_hash.end() );
      for (unsigned i = 0; i < tree_mers_hash.size(); ++i ) {
	 uint64_t tree_hash = tree_mers_hash[i];
	 kmer_t tree_mer = tree_mers[i];
	 //get all edges incident with tree_mer
	 for (unsigned j = 0; j < 4; ++j) {
	    if (IN.get( tree_hash, j )) {
	       Letter L( j );
	       kmer_t neighbor = pushOnFront( tree_mer, L, this->k );

	       if ( uo_hash.find( f(neighbor) ) == uo_hash.end()  ) { //TODO:could use hash update here
		  //we have found a neighbor not in this tree
		  //need to get the tree data of neighbor
		  vector< kmer_t > neighborSortedKmers;
		  map< kmer_t, unsigned > neighborHeights;

		  getTreeHeight( neighbor, neighborHeights, neighborSortedKmers );
		  //mergeTrees parameters will change
		  kmer_t neighbor_root;
		  uint64_t neighbor_root_hash;
		  getRoot( neighbor, neighbor_root, neighbor_root_hash );
		  kmer_t tree_root;
		  uint64_t tree_root_hash;
		  getRoot( tree_mer, tree_root, tree_root_hash );
		  Letter nei_letter = access_kmer( neighbor, k, 0 );
		  Letter tree_letter = access_kmer( tree_mer, k, k - 1 );
		  //		  BOOST_LOG_TRIVIAL(debug) << "Merging trees after removal";
		  //		  cerr << "here";
		  mergeTrees( neighbor, tree_mer, 
			      neighborHeights, tree_heights,
			      neighborSortedKmers, tree_mers,
			      neighbor_root_hash, tree_root_hash,
			      f(neighbor), tree_hash,
			      nei_letter, tree_letter
			      );

		  return;
	       }
	    }
	    if (OUT.get( tree_hash, j )) {
	       Letter L( j );
	       kmer_t neighbor = pushOnBack( tree_mer, L, this->k );
	       if ( uo_hash.find( f(neighbor) ) == uo_hash.end()  ) { //TODO: inefficient hash comp.
		  //we have found a neighbor not in this tree
		  //need to get the tree data of neighbor
		  vector< kmer_t > neighborSortedKmers;
		  map< kmer_t, unsigned > neighborHeights;
		  getTreeHeight( neighbor, neighborHeights, neighborSortedKmers );
		  kmer_t neighbor_root;
		  uint64_t neighbor_root_hash;
		  getRoot( neighbor, neighbor_root, neighbor_root_hash );
		  kmer_t tree_root;
		  uint64_t tree_root_hash;
		  getRoot( tree_mer, tree_root, tree_root_hash );
		  Letter nei_letter = access_kmer( neighbor, k, k - 1 );
		  Letter tree_letter = access_kmer( tree_mer, k, 0 );
		  //		  BOOST_LOG_TRIVIAL(debug) << "Merging trees after removal";
		  mergeTrees( tree_mer,  neighbor,
			      tree_heights, neighborHeights,
			      tree_mers, neighborSortedKmers,
			      tree_root_hash, neighbor_root_hash,
			      tree_hash, f( neighbor ),
			      tree_letter, nei_letter );

		  return;
	       }
	    }
	 }
      }
   }
   
   /*
    * removalUpdateForest
    * An edge has been removed between a child and parent
    * in one of the forest's trees 
    * This function performs the necessary changes to the forest
    * in an attempt to maintain tree height bounds
    *
    * This is a utility function for dynamicRemoveEdge
    */
   void removalUpdateForest( const kmer_t& child, const uint64_t& child_hash,
			     const kmer_t& parent, const uint64_t& parent_hash ) {

      //BOOST_LOG_TRIVIAL(debug) << "Updating the forest after edge between"
      //   << get_kmer_str(child, this->k) << " and " << get_kmer_str(parent, this->k)
      //   << " was deleted";

      //the child's subtree is without a root
      //temporarily, make it the root of its subtree
      //TODO: replace this by sampling best possible root
      fo.storeNode( child_hash, child );

      //get the heights of the two new trees
      //First, look at parent's tree
      vector< kmer_t > parentSortedKmers;
      vector< uint64_t > parentSortedKmers_hash;
      map< kmer_t, unsigned > parentHeights;
      unsigned parentTreeHeight = getTreeHeight( parent, parentHeights, parentSortedKmers, parentSortedKmers_hash );
      if (alpha > parentTreeHeight) {
	 //Fix this tree if possible
         //BOOST_LOG_TRIVIAL(debug) << "The tree of " << get_kmer_str(parent, this->k) << " must be fixed";
	 removalFixTree( parentSortedKmers, parentSortedKmers_hash, parentHeights );
      }
      //Second, look at child's tree
      vector< kmer_t > childSortedKmers;
      vector< uint64_t > childSortedKmers_hash;
      map< kmer_t, unsigned > childHeights;
      unsigned childTreeHeight = getTreeHeight( child, childHeights, childSortedKmers, childSortedKmers_hash );
      if (alpha > childTreeHeight) {
	 //Fix this tree if possible
         //BOOST_LOG_TRIVIAL(debug) << "The tree of " << get_kmer_str(child, this->k) << " must be fixed";
	 removalFixTree( childSortedKmers, childSortedKmers_hash, childHeights );
      }
   }
   

   
  /**
   * Either merge two trees into one at the edge u,v
   * Or break off one tree into the other
   * It is assumed that the De Bruijn graph has an edge from u to v
   * New root is sampled
   * u_heights and v_heights give the heights of all kmers in each tree
   * Returns whether the two trees were successfully merged
   */
  // TODO
     bool mergeTrees(const kmer_t& u, const kmer_t& v,
		     const map<kmer_t, unsigned>& u_heights, const map<kmer_t, unsigned>& v_heights,
		     const vector<kmer_t>& u_sorted_kmers, const vector<kmer_t>& v_sorted_kmers,
		     const u_int64_t& root_u_hash, const u_int64_t& root_v_hash,
		     const u_int64_t& u_hash, const u_int64_t& v_hash,
		     const Letter& u_letter, const Letter& v_letter) {

     //BOOST_LOG_TRIVIAL(debug) << "Merging the tree of " << get_kmer_str(u, this->k)
     //   << " with the tree of " << get_kmer_str(v, this->k);

     // Heights of u and v in their trees
     unsigned height_u = u_heights.at(u);
     unsigned height_v = v_heights.at(v);
     
     // Heights of trees
     unsigned treeheight_u = u_heights.at(u_sorted_kmers[u_sorted_kmers.size() - 1]);
     unsigned treeheight_v = v_heights.at(v_sorted_kmers[v_sorted_kmers.size() - 1]);

     if (treeheight_u > 3*this->alpha + 1) {
        BOOST_LOG_TRIVIAL(warning) << "Tree too big before merging.";
     }

     if (treeheight_v > 3*this->alpha + 1) {
        BOOST_LOG_TRIVIAL(warning) << "Tree too big before merging.";
     }

     // roots
     kmer_t root_u = u_sorted_kmers[0];
     kmer_t root_v = v_sorted_kmers[0];

     if (root_u == root_v) {
        BOOST_LOG_TRIVIAL(warning) << "You are attempting to merge a tree with itself.";
       return false;
     }

     unsigned merge_case = -1;

     // If one of the heights is greater than 2*alpha we will break off one into the other
     if ((height_u > 2*this->alpha) && (treeheight_v < this->alpha)) {
        // break off some of u's tree into v's

        merge_case = 1;

        // Find the node that is up alpha in the tree
        kmer_t breakpoint;
        u_int64_t breakpoint_hash;
        travelUp(u, this->alpha, breakpoint, breakpoint_hash);

        // TORI TEST
        kmer_t above_breakpoint;
        u_int64_t above_breakpoint_hash;
        travelUp(u, this->alpha + 1, above_breakpoint, above_breakpoint_hash);

        // temporarily store this as a node to effectively break it off
        this->fo.storeNode(breakpoint_hash, breakpoint);

        // TORI
        kmer_t root_test;
        kmer_t root_test_hash;
        getRoot(u, root_test, root_test_hash);
        if (root_test != breakpoint) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 1. We did not break at the breakpoint";
        }

        map<kmer_t, unsigned> test_broken;
        vector<kmer_t> testbrokenkmers;
        unsigned testheight = getTreeHeight(breakpoint, test_broken, testbrokenkmers);
        if (testheight > 2*this->alpha) {
           BOOST_LOG_TRIVIAL(warning) << "The broken off tree is too big,";
        }


        // reverse edges from breakpoint down to v
        reverseEdgesToRoot(u);

        // Now unstore the breakpoint        
        this->fo.unstoreNode(breakpoint_hash);

        // forest edge from u to v
        this->fo.setNode(u_hash, false, v_letter);

        // TORI
        root_test;
        root_test_hash;
        getRoot(u, root_test, root_test_hash);
        if (root_test != root_v) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 1. We did not merge u into v's tree.";
        }
       
        // TORI
        root_test;
        root_test_hash;
        getRoot(breakpoint, root_test, root_test_hash);
        if (root_test != root_v) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 1. We did not merge breakpoint into v's tree.";
        }

        kmer_t root_test2;
        u_int64_t root_test2_hash;
        getRoot(above_breakpoint, root_test2, root_test2_hash);
        if (root_test2 != root_u) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 1. u's tree is not still intact.";
        }


       if (root_test == root_test2) {
          BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 1. The two trees are now 1!";
       }
     
     }
     else if ((height_v > 2*this->alpha) && (treeheight_u < this->alpha)) {

        merge_case = 2;

        // break off some of v's tree into u's

        // Find the node that is up alpha in the tree
        kmer_t breakpoint;
        u_int64_t breakpoint_hash;
        travelUp(v, this->alpha, breakpoint, breakpoint_hash);

        // TORI TEST
        kmer_t above_breakpoint;
        u_int64_t above_breakpoint_hash;
        travelUp(v, this->alpha + 1, above_breakpoint, above_breakpoint_hash);

        // temporarily store this as a node to effectively break it off
        this->fo.storeNode(breakpoint_hash, breakpoint);

        // TORI
        kmer_t root_test;
        kmer_t root_test_hash;
        getRoot(v, root_test, root_test_hash);
        if (root_test != breakpoint) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 2. We did not break at the breakpoint";
        }

        map<kmer_t, unsigned> test_broken;
        vector<kmer_t> testbrokenkmers;
        unsigned testheight = getTreeHeight(breakpoint, test_broken, testbrokenkmers);
        if (testheight > 2*this->alpha) {
           BOOST_LOG_TRIVIAL(warning) << "The broken off tree is too big,";
        }


        // reverse edges from breakpoint down to v
        reverseEdgesToRoot(v);

        // Now unstore the breakpoint        
        this->fo.unstoreNode(breakpoint_hash);

        // forest edge from v to u
        this->fo.setNode(v_hash, true, u_letter);

        // TORI
        root_test;
        root_test_hash;
        getRoot(v, root_test, root_test_hash);
        if (root_test != root_u) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 2. We did not merge u into v's tree.";
        }
       
        // TORI
        root_test;
        root_test_hash;
        getRoot(breakpoint, root_test, root_test_hash);
        if (root_test != root_u) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 2. We did not merge breakpoint into v's tree.";
        }

        kmer_t root_test2;
        u_int64_t root_test2_hash;
        getRoot(above_breakpoint, root_test2, root_test2_hash);
        if (root_test2 != root_v) {
           BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 2. u's tree is not still intact.";
        }


       if (root_test == root_test2) {
          BOOST_LOG_TRIVIAL(warning) << "Something is wrong in merge case 2. The two trees are now 1!";
       }


     }
     else if ((treeheight_u < this->alpha) || (treeheight_v < this->alpha)) {

         merge_case = 3;

	 // At least one of our trees is too small, and the other is of height less than or equal to 2*alpha

	 // The distance from the furthest leaf in u along the path from that leaf to root_u
	 // then down to u, over to v, up to v's root, and then down v's furthest path that
	 // the root should be
	 unsigned root_hops = 0.5*(treeheight_u + height_u + 1 + height_v + treeheight_v);

	 kmer_t new_root;
	 u_int64_t new_root_hash; 

	 // The root is somewhere in u's tree
	 if (root_hops <= treeheight_u + height_u) {
	   //BOOST_LOG_TRIVIAL(debug) << "The root should be in " << get_kmer_str(u, this->k)
	   //   << "'s tree";
		
		// How much up from u the new root is
		unsigned hops = treeheight_u + height_u - root_hops;
		//BOOST_LOG_TRIVIAL(debug) << "It is " << hops << " hops up";

		if (hops > height_u ) {
		   // For now, we don't go past the root
		   //BOOST_LOG_TRIVIAL(warning) << "New root is past an old root. Will just use old root.";
		   hops = height_u;
		}

		// Make sure that root can reach v's nodes in less than or equal to 3*alpha
		if ((hops + height_v + treeheight_v + 1) > 3*alpha) {
		   BOOST_LOG_TRIVIAL(warning) << "Trees will not be merged because we might end up with one too big.";
		   return false;
		}

		// Go up to new root
		travelUp(u, hops, new_root, new_root_hash);

		// Make this our new root
		reverseEdgesToRoot(new_root);
		reverseEdgesToRoot(v);
		
		if (new_root != root_u) {
		   //BOOST_LOG_TRIVIAL(debug) << "Unstoring root " << root_u_hash;
		   this->fo.unstoreNode(root_u_hash);
		}

		this->fo.unstoreNode(root_v_hash);

		this->fo.storeNode(new_root_hash, new_root);

		// forest edge from v to u
		this->fo.setNode(v_hash, true, u_letter);
	     }
	     // The root is somewhere in v's tree
	     else {
		//BOOST_LOG_TRIVIAL(debug) << "The root should be in " << get_kmer_str(v, this->k)
		//   << "'s tree";

		// How much up from v the new root is
		unsigned hops = root_hops - treeheight_u - height_u - 1;
		//BOOST_LOG_TRIVIAL(debug) << "It is " << hops << " hops up";

		if (hops > height_v ) {
		   // For now, we don't go past the root
		   //BOOST_LOG_TRIVIAL(warning) << "New root is past an old root. Will just use old root.";
		   hops = height_v;
		}

		// Make sure that root can reach u's nodes in less than or equal to 3*alpha
		if ((treeheight_u + height_u + 1 + hops) > 3*alpha) {
		   BOOST_LOG_TRIVIAL(warning) << "Trees will not be merged because we might end up with one too big.";
		   return false;
		}

		// Go up to new root
		travelUp(v, hops, new_root, new_root_hash);

		reverseEdgesToRoot(new_root);
		reverseEdgesToRoot(u);

		if (new_root != root_v) {
		   //BOOST_LOG_TRIVIAL(debug) << "Unstoring root " << root_v_hash;
		   this->fo.unstoreNode(root_v_hash);
		}

		this->fo.unstoreNode(root_u_hash);

		this->fo.storeNode(new_root_hash, new_root);

		// forest edge from u to v
		this->fo.setNode(u_hash, false, v_letter);

	     }
     } 
     else {
        BOOST_LOG_TRIVIAL(warning) << "You have attempted to merge two trees that should not be merged.";
     }
    
     map<kmer_t, unsigned> testheights;
     vector<kmer_t> testkmers;
     unsigned testheight = getTreeHeight(u, testheights, testkmers);

     if (testheight > 3*this->alpha + 1) {
        BOOST_LOG_TRIVIAL(warning) << "u's tree too big after merging! u's tree was found to have height "
           << testheight << " and the maximum is " << 3*this->alpha << ". Merge case " << merge_case
           << ". It began with height " << height_u;
     }

     testheight = getTreeHeight(v, testheights, testkmers);

     if (testheight > 3*this->alpha + 1) {
        BOOST_LOG_TRIVIAL(warning) << "v's tree too big after merging! A tree was found to have height "
           << testheight << " and the maximum is " << 3*this->alpha << ". Merge case " << merge_case
           << ". It began with height " << height_v;
     }

     return true;

  }


   /**
    * Find the longest leaf in the tree that is not on the same path to the root
    * as u
    * heights is a map of all kmers in the tree to their heights
    * sorted_kmers is the kmers in the tree sorted from smallest to largest height
    * Stores results in leaf, leaf_height
    */
   void longestLeaf(const kmer_t& u, const map<kmer_t, unsigned>& heights,
      const vector<kmer_t>& sorted_kmers, kmer_t& leaf, unsigned& leaf_height) {

     // Find longest leaf that is not on the same path as u to the root
     for (int i = (sorted_kmers.size() - 1); i >= 0 ; --i) {

        //BOOST_LOG_TRIVIAL(debug) << "Checking kmer " << get_kmer_str(sorted_kmers[i], this->k);

        if (i == 0) {
           // We've gotten to the root. So no need to do tests.
            leaf = sorted_kmers[i];
            leaf_height = 0;
            break;
        }

        if ((!isAncestor(sorted_kmers[i], u)) && (!isAncestor(u, sorted_kmers[i]))) {
           leaf = sorted_kmers[i];
           leaf_height = heights.at(sorted_kmers[i]);
           break;
        }

     }
   }


  /*
   * Remove an edge from the data structure.
   * From u to v
   * Updates the forest.
   * Returns bool of whether an edge is actually deleted or not
   */
  bool dynamicRemoveEdge( const kmer_t& u, const kmer_t& v ) {

    //BOOST_LOG_TRIVIAL(debug) << "Deleting an edge from " << get_kmer_str(u, this->k)
    //   << " to " << get_kmer_str(v, this->k) << "...";

    //check if u, v are compatible
    unsigned ui, vi;
    for (unsigned i = 0; i < (this->k - 1); ++i) {
      ui = access_kmer( u, this->k, i + 1 );
      vi = access_kmer( v, this->k, i );
      //BOOST_LOG_TRIVIAL(debug) << "Checking the " << (i+1) << "/" << i << " spots";
      if (ui != vi) {
        //BOOST_LOG_TRIVIAL(debug) << "Cannot add an edge because there is not "
        //   << " k-1 length overlap.";
	return false;
      }
    }

    // For testing
    //if (!detect_membership(u)) {
    //   BOOST_LOG_TRIVIAL(debug) << "The node " << get_kmer_str(u, this->k)
    //      << " is not in the graph.";
    //   return false;
    //}

    //if (!detect_membership(v)) {
    //   BOOST_LOG_TRIVIAL(debug) << "The node " << get_kmer_str(v, this->k)
    //      << " is not in the graph.";
    //   return false;
    //}
    
    //making it this far means that an edge can be removed between them
    u_int64_t hashU = f( u );
    u_int64_t hashV = f( v );
    unsigned outIndex = access_kmer( v, k, k - 1 );
    unsigned inIndex = access_kmer( u, k, 0 );

    if ( !OUT.get(hashU, outIndex) ) {
      BOOST_LOG_TRIVIAL(debug) << "This edge doesn't exist.";
      return false; // edge doesn't exist
    }

    // Remove this edge from IN and OUT
    OUT.set(hashU, outIndex, false);
    IN.set(hashV, inIndex, false);

    // Now, need to update the forest if it includes this edge
    if ((!this->fo.isStored(hashU)) && (this->fo.getNext(hashU, u, this->k) == v)) {
       // u's parent is v
       removalUpdateForest( u, hashU, v, hashV );
    }
    else if ((!this->fo.isStored(hashV)) && (this->fo.getNext(hashV, v, this->k) == u)) {
       // v's parent is u
       removalUpdateForest( v, hashV, u, hashU );
    }
    else {
       // this edge is not in the forest
    } 

    return true;
  }

  

  /**
   * Reverse all forest edges along the path from node to its root
   * Used for combining two trees
   */
  void reverseEdgesToRoot(kmer_t node) {

     //BOOST_LOG_TRIVIAL(debug) << "Reversing edges to root of node "
     //   << get_kmer_str(node, this->k);

     u_int64_t node_kr = this->f.generate_KRHash_val_mod(node, this->k);
     u_int64_t node_hash = f.perfect_from_KR_mod(node_kr);

     // we are already at the root, nothing to store
     if (this->fo.isStored(node_hash)) {
        //BOOST_LOG_TRIVIAL(debug) << "That node is a root. Nothing to reverse.";
        return;
     }

     // Get all parent data
     kmer_t parent;
     u_int64_t parent_kr, parent_hash;
     getParentInfo(node, node_kr, node_hash, parent, parent_kr, parent_hash);
     // Whether forest edges are via IN or OUT
     bool in = this->fo.parent_in_IN(node_hash);
     bool parent_in;

     // Need to save this to move forward after we've changed link to parent
     kmer_t grandparent;
     u_int64_t grandparent_kr, grandparent_hash;

     // Keep reversing edges until we get to the last one
     while (!this->fo.isStored(parent_hash)) {

        //BOOST_LOG_TRIVIAL(debug) << "Looking at edge between " <<
        //   get_kmer_str(node, this->k) << " and " << get_kmer_str(parent, this->k);

        // Save links to the next place on the tree
        getParentInfo(parent, parent_kr, parent_hash, grandparent,
           grandparent_kr, grandparent_hash);

        // Need to save how to get from parent to grandparent, since forest
        // edge will be removed
        parent_in = this->fo.parent_in_IN(parent_hash);

        //BOOST_LOG_TRIVIAL(debug) << grandparent_hash;
        // Reverse the edge between node and its parent
        flipEdge(parent_hash, node_hash, node, in);

        // move forward
        node = parent;
        node_kr = parent_kr;
        node_hash = parent_hash;
        parent = grandparent;
        parent_kr = grandparent_kr;
        parent_hash = grandparent_hash;
        in = parent_in;

     }

     // Last edge going to root
     flipEdge(parent_hash, node_hash, node, in);

  }


  /**
   * Add a forest edge from node to new_parent.
   * "in" gives whether you can get from new_parent
   * to node via IN or not.
   * So this can be used to "flip" a forest edge, although it doesn't
   * assume that that edge is currently in the forest (since "in" is
   * passed in, where
   * the child was once new_parent, and the parent was once node.
   */
  void flipEdge(const u_int64_t& node_hash, const u_int64_t& new_parent_hash,
     const kmer_t& new_parent_kmer, bool in) {

     //BOOST_LOG_TRIVIAL(debug) << "Switching parent of "
     //   << node_hash
     //   << " to parent " << get_kmer_str(new_parent_kmer, this->k);

     // If the child went to parent via OUT, then the parent will get to
     // the child from IN, etc.
     this->fo.set_parent_in_IN(node_hash, !in);

     //BOOST_LOG_TRIVIAL(debug) << "Node " << node_hash
     //   << " reaches its parent via IN? " << this->fo.parent_in_IN(node_hash);

     Letter l;
     if (in) {
        // Need the last letter of the new_parent
        l = access_kmer(new_parent_kmer, this->k, this->k - 1);
        this->fo.setLetter(node_hash, l);
     }
     else {
        // Need the first letter of the new_parent
        l = access_kmer(new_parent_kmer, this->k, 0);
        this->fo.setLetter(node_hash, l);

     }

     //BOOST_LOG_TRIVIAL(debug) << "By what letter? " << this->fo.getLetter(node_hash);


   }


  /**
   * Get parent info from child info
   */
  void getParentInfo(const kmer_t& node, const u_int64_t& node_kr, const u_int64_t& node_hash,
                     kmer_t& parent, u_int64_t& parent_kr, u_int64_t& parent_hash) {

      //BOOST_LOG_TRIVIAL(debug) << "Getting info for the parent of "
      //   << get_kmer_str(node, this->k);

      parent = this->fo.getNext(node_hash, node, this->k);

      //BOOST_LOG_TRIVIAL(debug) << "That's parent " << get_kmer_str(parent, this->k);

      // The front and back letters of the edge between node and parent
      // used to compute hash values of parent
      unsigned front, back;

      // Get data depending on how we access parent
      if (this->fo.parent_in_IN(node_hash)) {

         // We got it in IN, so the edge goes from parent to node
        unsigned back = access_kmer(node, this->k, this->k - 1);
	unsigned front = access_kmer( parent, this->k, 0 );

	parent_kr = f.update_KRHash_val_IN_mod(node_kr, front, back);

	parent_hash = f.perfect_from_KR_mod(parent_kr);

        //BOOST_LOG_TRIVIAL(debug) << "This edge starts with letter "
        //   << front << " and ends with letter " << back;
        //BOOST_LOG_TRIVIAL(debug) << "It has hash value " << parent_hash;

      } else {
         // We got it in OUT, so the edge goes from node to parent

        unsigned front = access_kmer(node, this->k, 0);
	unsigned back = access_kmer(parent, this->k, this->k - 1 );

        //BOOST_LOG_TRIVIAL(debug) << "This edge starts with letter "
        //   << front << " and ends with letter " << back;

	parent_kr = f.update_KRHash_val_OUT_mod(node_kr, front, back);
	parent_hash = f.perfect_from_KR_mod(parent_kr);
        //BOOST_LOG_TRIVIAL(debug) << "It has hash value " << parent_hash;
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

  void printHashFunction(unordered_set<kmer_t>& kmers) {

    unordered_set<kmer_t>::iterator i;
    for (i = kmers.begin(); i != kmers.end(); ++i) {
        cout << setw(10) << get_kmer_str(*i, this->k);
        cout << setw(10) << this->f(*i); 
        cout << endl;
    } 

  }

  // Given a kmer, decide if it is one in our graph
  bool detect_membership( kmer_t m ) {
    //    BOOST_LOG_TRIVIAL(debug) << "Detecting membership of " << get_kmer_str(m, this->k);

    // The hash value of our kmer
    // Need to keep track of KRval, so it can be updated
    // largeUnsigned type is provided by 'generate_hash.h'
     //    largeUnsigned KR_val = f.generate_KRHash_raw( m, k ) ;
     u_int64_t KR_val = f.generate_KRHash_val_mod( m, k ) ;
     u_int64_t hash = f.perfect_from_KR_mod( KR_val );
    
     //     BOOST_LOG_TRIVIAL(debug) << "Correct hash, computed: " << f(m) << ' ' << hash;
    
    //    BOOST_LOG_TRIVIAL(debug) << "It has hash value " << hash;

    // If it is a real kmer value, it must map to 0..1-n
    if (hash >= n) {
       //       BOOST_LOG_TRIVIAL(debug)  << "Returning false because of hash function blowup" << endl;
      return false;
    }

    //number of times we have traveled in the tree
    unsigned hopcount = 0;
    bool in;    //true = IN, false = OUT
    unsigned letter; //needed to confirm edge from parent's side
    
    // Keep going in the tree until we reach a hash that is stored as a root
    while ( !(this->fo.isStored(hash)) ) {
      ++hopcount;
      if (hopcount > 3*alpha + 10) {
	 //BOOST_LOG_TRIVIAL( debug ) << "Returning false because of hopcount" << endl;
	 return false; //we have encountered a cycle in our attempt to travel to a root
      }
	
      //do we progress with IN or OUT
      in = this->fo.parent_in_IN(hash);
      if (in) {
	//the parent is an IN-neighbor of m
	//so we need m's last character
	letter = access_kmer( m, k, k - 1 );

      } else {
	//the parent is an OUT-neighbor of m
	//so we need m's first character
	letter = access_kmer( m, k, 0 );
      }

      //      BOOST_LOG_TRIVIAL(debug) << "Child k-mer: " << get_kmer_str( m, k );
      
      // deduce the parent's kmer
      m = this->fo.getNext(hash, m, k);

      //      BOOST_LOG_TRIVIAL(debug) << "Parent k-mer: " << get_kmer_str( m, k );
      
      // get the parent's hash
      if (in) {
	unsigned letter_front_parent = access_kmer( m, k, 0 );
	KR_val = f.update_KRHash_val_IN_mod( KR_val, letter_front_parent, letter );
	hash = f.perfect_from_KR_mod( KR_val );
      } else {
	unsigned letter_back_parent = access_kmer( m, k, k - 1 );
	KR_val = f.update_KRHash_val_OUT_mod( KR_val, letter, letter_back_parent );
	hash = f.perfect_from_KR_mod( KR_val );
      }
      //hash = f(m); 
      //      BOOST_LOG_TRIVIAL(debug) << "Correct hash, computed: " << f(m) << ' ' << hash;
      //            BOOST_LOG_TRIVIAL(debug) << "Correct, computed KR_val: "
      //            			       << f.generate_KRHash_val_mod( m, k ) << ' ' << KR_val;

      // hash must be in 0...n-1
      if (hash >= n) {
	 //	 BOOST_LOG_TRIVIAL(debug) << "Returning false because of hash function blowup" << endl;
	return false;
      }
      //confirm the edge from parent side
      if (in) {
	if (!(OUT.get(hash, letter))) {
	   //BOOST_LOG_TRIVIAL(debug) << "Returning false because IN,OUT verification failed" << endl;
	   return false;
	}
      } else {
	if (!(IN.get(hash, letter))) {
	   //BOOST_LOG_TRIVIAL(debug) << "Returning false because IN,OUT verification failed" << endl;
	  return false;
	}
      }
      
    }

    //    BOOST_LOG_TRIVIAL(debug) << "Root " << get_kmer_str(m, k) << " is deduced";

    // now we have a forest node that we have the kmer of stored
    // So we just have to test if it is accurate or not
    if (m == fo.roots[hash])
      return true;
    else
      return false;
  }

  // Set the parent's kmer given the child's kmer
  kmer_t getParent(const kmer_t& m) {

     kmer_t parent_kr = 0;
     u_int64_t m_kr = f.generate_KRHash_val_mod( m, this->k );
     kmer_t parent_kmer;

     getParent(m, m_kr, parent_kmer, parent_kr);

     return parent_kmer;

     //BOOST_LOG_TRIVIAL(debug) << "The parent's hash value is " << f.perfect_from_KR_mod(parent_kr);
  }


  // Set the parent's kmer, and karp-rabin value 
  // given the child's kmer and krval
  void getParent(const kmer_t& m, const u_int64_t& kr_val,
     kmer_t& parent_kmer, u_int64_t& parent_kr) {

     //BOOST_LOG_TRIVIAL(debug) << "Finding parent data for kmer "
     //   << get_kmer_str(m, this->k);

     // Get the hash value of m
     u_int64_t hash = f.perfect_from_KR_mod( kr_val );
     //BOOST_LOG_TRIVIAL(debug) << "The hash value was found to be " << hash; 

     // whether we get to the parent via IN
     bool in  = this->fo.parent_in_IN(hash);

     unsigned letter; // what letter m has but its parent doesn't

     if (in) {	//the parent is an IN-neighbor of m
	//the parent is an IN-neighbor of m
	//so we need m's last character
	letter = access_kmer( m, this->k, this->k - 1 );
     } else {
	//the parent is an OUT-neighbor of m
	//so we need m's first character
	letter = access_kmer( m, this->k, 0 );
     }

     //BOOST_LOG_TRIVIAL(debug) << "Reached via IN? " << in << ". What letter? " << letter;

     // deduce the parent's kmer
     parent_kmer = this->fo.getNext(hash, m, k);


     //BOOST_LOG_TRIVIAL(debug) << "The parent's kmer is " << get_kmer_str(parent_kmer, this->k);

     if (in) {
	unsigned letter_front_parent = access_kmer( parent_kmer, k, 0 );
	parent_kr = f.update_KRHash_val_IN_mod( kr_val, letter_front_parent, letter );
     } else {
	unsigned letter_back_parent = access_kmer( parent_kmer, k, k - 1 );
	parent_kr = f.update_KRHash_val_OUT_mod( kr_val, letter, letter_back_parent );
     }
  }


  /**
   * Checks whether ancestor is an ancestor in the tree of node
   */
  bool isAncestor(kmer_t node, const kmer_t& ancestor) {

    // First check whether these are the same before moving up
    if (node == ancestor) {
       //BOOST_LOG_TRIVIAL(debug) << "Node an ancestor are equal. So ancestor is an ancestor.";
       return true;
    }

    // ancestor's hash
    u_int64_t ancestor_hash = this->f(ancestor);

    // node's hash and other stuff to travel up efficiently
    u_int64_t kr = f.generate_KRHash_val_mod(node, this->k);
    u_int64_t hash = f.perfect_from_KR_mod(kr);
    u_int64_t parent;
    u_int64_t parent_kr;

    //BOOST_LOG_TRIVIAL(debug) << "Finding whether " << get_kmer_str(ancestor, this->k)
    //   << " is an ancestor of " << get_kmer_str(node, this->k);

    // node could be a root
    if (this->fo.isStored(hash)) {
       //BOOST_LOG_TRIVIAL(debug) << "The node is a root and not equal to ancestor";
       // We already checked whether node and ancestor were the same
       return false;
    } 

    // Keep getting parent until we've found a root or ancestor
    while (!this->fo.isStored(hash)) {

       getParent(node, kr, parent, parent_kr);

       // is this the ancestor
       if (parent == ancestor) {
          //BOOST_LOG_TRIVIAL(debug) << "This is ancestor. So it is true.";
          return true;
       }

       // move up
       hash = f.perfect_from_KR_mod(parent_kr);
       node = parent;
       kr = parent_kr;
       //BOOST_LOG_TRIVIAL(debug) << "Moved to parent " << get_kmer_str(node, this->k);
    }

    //BOOST_LOG_TRIVIAL(debug) << "Ancestor was not found. Not an ancestor.";
    return false;
  }



  /**
   * Get the root in the forest that this node goes to
   */
  unsigned getRoot(kmer_t node, kmer_t& root, u_int64_t& root_hash) {

    unsigned node_height = 0;

    u_int64_t kr = f.generate_KRHash_val_mod(node, this->k);
    u_int64_t hash = f.perfect_from_KR_mod(kr);
    u_int64_t parent;
    u_int64_t parent_kr;

    //BOOST_LOG_TRIVIAL(debug) << "Finding the root of " << node;

    if (this->fo.isStored(hash)) {
      //BOOST_LOG_TRIVIAL(debug) << "Root is " << node;
      root = node;
      root_hash = hash;
      return node_height;
    } 

    // Keep getting parent until we've found a root
    while (!this->fo.isStored(hash)) {

       getParent(node, kr, parent, parent_kr);
       //BOOST_LOG_TRIVIAL(debug) << "Moved to parent " << parent;
       //BOOST_LOG_TRIVIAL(debug) << "With kr val " << parent_kr;
       hash = f.perfect_from_KR_mod(parent_kr);
       //BOOST_LOG_TRIVIAL(debug) << "Hash value " << hash;
       node = parent;
       //BOOST_LOG_TRIVIAL(debug) << "Now our node is " << node;
       kr = parent_kr;
       node_height++;
    }

    //BOOST_LOG_TRIVIAL(debug) << "Root is " << node;
    root = node;
    root_hash = hash;

    return node_height;
  }

  /**
   * Travel up k hops from the node
   * If gets to root before finished hopping just returns the root
   */
  // TODO
  // WORK IN PROGRESS
  void travelUp(kmer_t node, const unsigned& hops, kmer_t& ancestor, u_int64_t& ancestor_hash) {

    u_int64_t kr = f.generate_KRHash_val_mod(node, this->k);
    u_int64_t hash = f.perfect_from_KR_mod(kr);
    u_int64_t parent;
    u_int64_t parent_kr;

    //BOOST_LOG_TRIVIAL(debug) << "Travelling up " << hops << " hops from "
    //   << get_kmer_str(node, this->k);

    if (this->fo.isStored(hash)) {
      //BOOST_LOG_TRIVIAL(debug) << "Travelled 0 hops";
      ancestor = node;
      ancestor_hash = hash;
      return;
    } 

    unsigned hop_count = 0;
    // Keep getting parent until we've travelled enough hops or found a root
    while ((hop_count < hops) && !this->fo.isStored(hash)) {

       getParent(node, kr, parent, parent_kr);
       //BOOST_LOG_TRIVIAL(debug) << "Moved to parent " << get_kmer_str(parent, this->k);
       hash = f.perfect_from_KR_mod(parent_kr);
       node = parent;
       kr = parent_kr;
       hop_count++;
    }

    //BOOST_LOG_TRIVIAL(debug) << "Finished travelling. Found " << get_kmer_str(node, this->k);
    ancestor = node;
    ancestor_hash = hash;

    return;
  }


  /**
   * Given a node, get all of its children in the forest
   * children_hash contains hash values of children in the
   * the same as order as children
   */
   void getChildren(const kmer_t& node,
		    //		    uint64_t& node_kr,   //The input node's karp-rabin hash value
		    vector<kmer_t>& children,
		    vector< uint64_t >& children_hash) {
     children.clear();
     children_hash.clear();
     
     // get KR value to make getting neighbor hash values easier
     u_int64_t node_kr = f.generate_KRHash_val_mod(node, this->k);

     vector<kmer_t> neighbors;
     vector<bool> inorout;
     //BOOST_LOG_TRIVIAL(debug) << "Getting children of "
     //  << get_kmer_str(node, this->k);

     // get all neighbors
     get_neighbors(node, neighbors, inorout);

     //BOOST_LOG_TRIVIAL(debug) << "Found neighbors: ";
     //for (int i = 0; i < neighbors.size(); i++) {
     //  BOOST_LOG_TRIVIAL(debug) << " " << neighbors[i];
     //}
     
     // for each edge, the first and last characters
     unsigned first, last;
     unsigned node_first = access_kmer(node, this->k, 0);
     unsigned node_last = access_kmer(node, this->k, k-1);
     u_int64_t neighbor_hash;
     kmer_t parent; // holds the parent of a neighbor for comparison

     // get the hash values for each neighbor, check if node is their parent in forest
     for (int i = 0; i < neighbors.size(); i++) {

       if (inorout[i] == 1) {
         // this neighbor has an edge going towards node
         // so we need its first letter, node's last
         first = access_kmer(neighbors[i], this->k, 0);
         neighbor_hash = f.perfect_from_KR_mod(f.update_KRHash_val_IN_mod(node_kr,
           first, node_last));
         //BOOST_LOG_TRIVIAL(debug) << "The hash for neighbor "
         //  << get_kmer_str(neighbors[i], this->k)
         //  << " is " << child_hash;

       }
       else {
         // this neighbor has an edge from node towards it
         // so we need its last letter, node's first
         last = access_kmer(neighbors[i], this->k, k-1);
         neighbor_hash = f.perfect_from_KR_mod(f.update_KRHash_val_OUT_mod(node_kr,
           node_first, last));
        // BOOST_LOG_TRIVIAL(debug) << "The hash for neighbor "
        //   << get_kmer_str(neighbors[i], this->k)
        //   << " is " << child_hash;
       }

       //BOOST_LOG_TRIVIAL(debug) << "That child's parent is " <<
       //  get_kmer_str(this->fo.getNext(child_hash, neighbors[i], this->k), this->k);

       // Get the parent of this neighbor
       parent = this->fo.getNext(neighbor_hash, neighbors[i], this->k);

       // this node is its parent, and it is not a root
       if ((node == parent) && !(this->fo.isStored(neighbor_hash))) {
         // this is a child in the forest
         children.push_back(neighbors[i]);
	 children_hash.push_back( neighbor_hash );
         //BOOST_LOG_TRIVIAL(debug) << get_kmer_str(neighbors[i], this->k)
         //  << " is a child";
	 //print_kmer( neighbors[i], k, cerr );
	 //cerr << ' ' << f(neighbors[i]) << ' ' << neighbor_hash << endl;
       }

       
     }
  

  }

  /**
   * Given a kmer, find the height of the tree that node is in in the forest
   */
  unsigned getTreeHeight(const kmer_t& node,
			 map<kmer_t, unsigned>& heights,
			 vector<kmer_t>& sorted_kmers,
			 vector<uint64_t>& sorted_kmers_hash = DEFAULT_VECTOR ) {
     //BOOST_LOG_TRIVIAL(debug) << "Getting the height of tree with node "
     //   << get_kmer_str(node, this->k);

     // get the root of this tree
     kmer_t root;
     u_int64_t root_hash;
     getRoot(node, root, root_hash);

     unsigned res = getTreeHeightRoot(root, heights, sorted_kmers, sorted_kmers_hash );

     return res;
  }

  
   
  /**
   * Same as above, but it is assumed that the input node is the root of the tree
   * Keeps track of each kmer in tree and its height
   * Also makes vector of kmers sorted from smallest to biggest height
   * sorted_kmers_hash is the hash value of the sorted_kmers in the same order
   * 
   */
  // TODO: Make array copying more efficient
   //TODO: is not fully using hash function update. Needs to be fixed
  unsigned getTreeHeightRoot(const kmer_t& root,
			     map<kmer_t, unsigned>& heights,
			     vector<kmer_t>& sorted_kmers,
			     vector< uint64_t>& sorted_kmers_hash = DEFAULT_VECTOR ) {

     //BOOST_LOG_TRIVIAL(debug) << "Getting the height of tree with root "
     //  << get_kmer_str(root, this->k);

     // keep descending down the tree by getting children until
     // there doesn't exist anymore. Count how many times this
     // must happen

     heights.clear();
     sorted_kmers.clear();
     sorted_kmers_hash.clear();

     // add root to heights
     heights[root] = 0;
     sorted_kmers.push_back(root);
     sorted_kmers_hash.push_back( f(root) );
				 
     vector<kmer_t> children;
     vector<uint64_t> children_hash;
     getChildren(root, children, children_hash);

     vector<kmer_t> children_children;        // the children of the children
     vector<uint64_t> children_children_hash; // the children of the children
     vector<kmer_t> children_node; // holds the children of a single node
     vector<uint64_t> children_node_hash; // holds the children of a single node

     unsigned height = 0;

     unsigned warnings = 0;

     // Until we have reached the most far node in the tree
     while (children.size() != 0) {
        height++;

        // Get children of each child
        for (int i = 0; i < children.size(); i++) {
           // add child with height
           heights[children[i]] = height;
           sorted_kmers.push_back(children[i]);
	   sorted_kmers_hash.push_back( children_hash[i] );

           //BOOST_LOG_TRIVIAL(debug) << "Getting the children of node "
           //  << get_kmer_str(children[i], this->k);

           // get the ith child's children
           getChildren(children[i], children_node, children_node_hash );

           // add each one to children_children
           for (int j = 0; j < children_node.size(); j++) {
              //BOOST_LOG_TRIVIAL(debug) << "Found child "
              //  << get_kmer_str(children_node[j], this->k);
              children_children.push_back(children_node[j]);
	      children_children_hash.push_back(children_node_hash[j]);
           }
        }

        children.swap(children_children);
        children_children.clear();
	children_hash.swap( children_children_hash );
	children_children_hash.clear();
        //BOOST_LOG_TRIVIAL(debug) << "There are now " << children.size() << " children";
     }

      //BOOST_LOG_TRIVIAL(debug) << "Height is determined to be " << height;

      return height;

  }


   // Given a kmer, decide if it is one in our graph
   // recomputes hash value at each step, so for testing purposes only
  bool inefficient_detect_membership( kmer_t m ) {
    //    BOOST_LOG_TRIVIAL(debug) << "Detecting membership of " << get_kmer_str(m, this->k);

    // The hash value of our kmer
    // Need to keep track of KRval, so it can be updated
     u_int64_t hash = f(m);

    //    BOOST_LOG_TRIVIAL(debug) << "It has hash value " << hash;

    // If it is a real kmer value, it must map to 0..1-n
    if (hash >= n) {

      return false;
    }

    //number of times we have traveled in the tree
    unsigned hopcount = 0;
    bool in;    //true = IN, false = OUT
    unsigned letter; //needed to confirm edge from parent's side
    
    while ( !(this->fo.isStored(hash)) ) {
      ++hopcount;
      if (hopcount > 3*alpha + 10) {

	return false; //we have encountered a cycle in our attempt to travel to a root
      }
	
      //do we progress with IN or OUT
      in = this->fo.parent_in_IN(hash);
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
      m = this->fo.getNext(hash, m, k);
      hash = f(m); 

      // hash must be in 0...n-1
      if (hash >= n) {

	return false;
      }
      //confirm the edge from parent side
      if (in) {
	if (!(OUT.get(hash, letter))) {

	  return false;
	}
      } else {
	if (!(IN.get(hash, letter))) {

	  return false;
	}
      }
      
    }

    //    BOOST_LOG_TRIVIAL(debug) << "Root " << get_kmer_str(m, k) << " is deduced";

    // now we have a forest node that we have the kmer of stored
    // So we just have to test if it is accurate or not
    if (m == fo.roots[hash])
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

        // If a neighbor can be reached via OUT from c, this is the letter to reach it
        Letter first_c = access_kmer(c, k, 0);
        // Now for IN
        Letter last_c = access_kmer(c, k, k - 1);

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
              /**
               * The iith neighbor of c was reached by c via IN. Therefore to reach
               * c from the iith neighbor, we must go OUT. We will need the last
               * character of c
               */
              fo.setNode(f_m, false, last_c);
	    } else {
              // similar logic to the first case
              fo.setNode(f_m, true, first_c);
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

    kmers.swap( visited_mers );
  }

  /**
   * Given a kmer c, get all neighbor kmers (neighbor_kmers) 
   * for neighbor_kmers[ i ], v_inorout[ i ] is true if in-neighbor, false otherwise
   */
  void get_neighbors(const kmer_t& c,
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
      if ( IN.get(fc, ii) ) {
	//have in-neighbor with letter
	Letter letter ( ii );
	//	BOOST_LOG_TRIVIAL(trace) << "letter: " << ii << ' ' << letter.getNum();
	// Deduce that kmer by tacking letter at the beginning of c
	kmer_t e = pushOnFront(c, letter, k);
	neighbor_kmers.push_back( e );

	//	BOOST_LOG_TRIVIAL(trace) << "neighbor: " + get_kmer_str(e, k);
	
	v_inorout.push_back( true ); //true=IN
      }
      if ( OUT.get(fc, ii) ) {
	//have out-neighbor with letter
	Letter letter ( ii );

	// Deduce that kmer by tacking letter at the end of c
	kmer_t e = pushOnBack(c, letter, k);
	neighbor_kmers.push_back( e );

	v_inorout.push_back( false ); //false=OUT

      }
    }

    //    BOOST_LOG_TRIVIAL(trace) << "Neighbors found: ";
    //for (unsigned i = 0; i < neighbor_kmers.size(); ++i) {
      //      BOOST_LOG_TRIVIAL(trace) << i << ' ' << get_kmer_str( neighbor_kmers[i], k ) << " is an "<<	v_inorout[ i ] << " neighbor.";
      
    //}
  }
 
  void store( kmer_t mer ) {
    u_int64_t val = f( mer );
    this->fo.storeNode(val, mer);
  }
  
  // Move mer from kmers and into visited
  void move_kmer(unordered_set< kmer_t >& kmers, unordered_set< kmer_t >& visited,
			 kmer_t mer ) {
    kmers.erase( mer );
    visited.insert( mer );
  }




  // Print IN and OUT with the Kmer strings down the side, and the characters on the top
  void printINandOUT(unordered_set<kmer_t>& kmers) {
 
    cout << setw(30) << "===== IN =====" << endl << endl;
 
    cout << setw(10) << ""; 
    cout << setw(5) << "A";
    cout << setw(5) << "C";
    cout << setw(5) << "G";
    cout << setw(5) << "T" << endl;

    unordered_set<kmer_t>::iterator i;

    for (i = kmers.begin(); i != kmers.end(); ++i) {
        cout << setw(10) << get_kmer_str(*i, this->k);

        u_int64_t hash = this->f(*i); 
        cout << setw(5) << this->IN.get(hash, 0);
        cout << setw(5) << this->IN.get(hash, 1);
        cout << setw(5) << this->IN.get(hash, 2);
        cout << setw(5) << this->IN.get(hash, 3);
        cout << endl;
    } 

    cout << endl;



    cout << setw(30) << "===== OUT =====" << endl << endl;
 
    cout << setw(10) << ""; 
    cout << setw(5) << "A";
    cout << setw(5) << "C";
    cout << setw(5) << "G";
    cout << setw(5) << "T" << endl;

    for (i = kmers.begin(); i != kmers.end(); ++i) {
        cout << setw(10) << get_kmer_str(*i, this->k);

        u_int64_t hash = this->f(*i); 
        cout << setw(5) << this->OUT.get(hash, 0);
        cout << setw(5) << this->OUT.get(hash, 1);
        cout << setw(5) << this->OUT.get(hash, 2);
        cout << setw(5) << this->OUT.get(hash, 3);
        cout << endl;
    } 

    cout << endl;

  }

  // Print the forest with the Kmer strings down the side
  void printForest(unordered_set<kmer_t>& kmers) {
 
    cout << setw(30) << "===== FOREST =====" << endl << endl;

 
    cout << setw(10) << ""; 
    cout << setw(10) << "Parent";
    cout << setw(10) << "Is root?" << endl;;

    unordered_set<kmer_t>::iterator i;

    for (i = kmers.begin(); i != kmers.end(); ++i) {
        cout << setw(10) << get_kmer_str(*i, this->k);

        u_int64_t hash = this->f(*i);

        if (!this->fo.isStored(hash)) {
            // This is not a root
           cout << setw(10) << get_kmer_str(this->fo.getNext(hash, *i, this->k), this->k);
           cout << setw(10) << "No" << endl;
        }
        else {
            // This is a root
           cout << setw(10) << "None";
           cout << setw(10) << "Yes" << endl;
        }
 
    } 

    cout << endl;

  }



// END CLASS
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

// Push letter onto front of kmer and return kmer
// Going backwards along an edge (the relationship comes from orig's IN).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnFront(const kmer_t& orig, Letter& letter, unsigned k) {
  //  BOOST_LOG_TRIVIAL(trace) << "pushOnFront " + get_kmer_str(orig, k) << ' ' << letter.getNum();
  // Get the kmer with back pushed off orig
  kmer_t new_kmer = orig >> 2;

  set_kmer( new_kmer, k, 0, letter);

  return new_kmer;
} 

// Push letter onto back of kmer and return kmer
// Going forward along an edge (the relationship comes from orig's OUT).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnBack(const kmer_t& orig, Letter& letter, unsigned k) {
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

/**
 * Given a kmer length k, returns the maximum value that kmer
 * is allowed to be
 */
kmer_t getMaxVal(unsigned k) {

   kmer_t max = 0;

   for (int i = 0; i < k; ++i) {
       max <<= 2;
       max |= 3;
   }

   return max;
}

// Push off the last letter
void push_last_letter(const kmer_t& edge, kmer_t& u) {

   u = edge >> 2;

}

// Remove front letter from kmer of length k+1
// To make it a kmer of length k
void remove_front_letter(const kmer_t& edge, kmer_t& v, const unsigned& k) {

  //v = edge & ~(3 << 2*(k));

  kmer_t op = 3; 
  op = op << 2*(k);
  op = ~op;
  v = edge & op;

  //BOOST_LOG_TRIVIAL(debug) << "Removed front letter of " << edge << " to get " << v;
}

#endif

