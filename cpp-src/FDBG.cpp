#ifndef FDBG_class
#define FDBG_class

//#include "hash/HashUtil.cpp"
#include "hash/generate_hash.h"
#include <vector>
#include <string>
#include <unordered_set>
#include <map>
#include <chrono>
#include "BitArray.cpp"

using namespace std;

class Letter; //needed for set_kmer declaration

kmer_t getMaxVal(unsigned);
void push_last_letter(const kmer_t&, kmer_t&);
void remove_front_letter(const kmer_t&, kmer_t&, const unsigned&);
void set_kmer( kmer_t& mer, unsigned k, unsigned i, Letter& c );
void set_kmer( kmer_t& mer, unsigned k, unsigned i, char c );
kmer_t pushOnFront(const kmer_t& orig, Letter& letter, unsigned);
kmer_t pushOnBack(const kmer_t& orig, Letter& letter, unsigned);

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
    void setNode(u_int64_t& i, bool IN, Letter& l) {

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

     // go through each tree

     map<u_int64_t, kmer_t>::iterator iter;
     for (iter = this->fo.roots.begin(); iter != this->fo.roots.end(); ++iter) {

        root = iter->second;

        //BOOST_LOG_TRIVIAL(debug) << "Looking at tree with root "
        //   << get_kmer_str(root, this->k);

        height = getTreeHeightRoot(root);

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
    getRoot(u, root_u, root_u_hash);
    kmer_t root_v;
    u_int64_t root_v_hash;
    getRoot(v, root_v, root_v_hash);

    //BOOST_LOG_TRIVIAL(debug) << "root " << get_kmer_str(root_u, this->k)
    //   << " and root " << get_kmer_str(root_v, this->k);



    unsigned height_u = getTreeHeightRoot(root_u);
    unsigned height_v = getTreeHeightRoot(root_v);

    //BOOST_LOG_TRIVIAL(debug) << "One is in a tree of height " << height_u
    //   << " with root " << get_kmer_str(root_u, this->k)
    //   << " and the other a tree of height " << height_v
    //   << " with root " << get_kmer_str(root_v, this->k);

    if (root_u == root_v) {
       // No trees to merge
       //BOOST_LOG_TRIVIAL(debug) << "Both trees already have the same root. No merging.";
       return true;
    }

    // Both trees are too small. Merge them.
    if ((height_u < this->alpha) && (height_v < this->alpha)) {

       //BOOST_LOG_TRIVIAL(debug) << "Both trees are below the min height. Merging.";
       // the TO kmer's (v) tree is left alone
       // the FROM kmer's tree is added to the TO kmer's tree
       // by reversing all edges from the FROM kmer to its root
       reverseEdgesToRoot(u);

       // the root of the FROM kmer is unstored
       this->fo.unstoreNode(root_u_hash);

       // Add forest edge from u to v
       // u's parent is switched
       Letter l (outIndex);
       this->fo.setNode(hashU, false, l);

       //BOOST_LOG_TRIVIAL(debug) << "Merged trees.";
    }
    else {
        //BOOST_LOG_TRIVIAL(debug) << "Both trees were at least the minimum height. No merging.";
    }

    return true;
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
    if (!detect_membership(u)) {
       BOOST_LOG_TRIVIAL(debug) << "The node " << get_kmer_str(u, this->k)
          << " is not in the graph.";
       return false;
    }

    if (!detect_membership(v)) {
       BOOST_LOG_TRIVIAL(debug) << "The node " << get_kmer_str(v, this->k)
          << " is not in the graph.";
       return false;
    }
    
    //making it this far means that an edge can be removed between them
    //Add the edge to OUT[ f(u) ] and to IN[ f(v) ]
    //if edge was already present, quit
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
       // u is now made root of its subtree
       //BOOST_LOG_TRIVIAL(debug) << get_kmer_str(u, this->k) << " is now root of its subtree.";
       this->fo.storeNode(hashU, u);
    }
    else if ((!this->fo.isStored(hashV)) && (this->fo.getNext(hashV, v, this->k) == u)) {
       // v's parent is u
       // v is now made root of its subtree
       //BOOST_LOG_TRIVIAL(debug) << get_kmer_str(v, this->k) << " is now root of its subtree.";
       this->fo.storeNode(hashV, v);
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
   * Get the root in the forest that this node goes to
   */
  void getRoot(kmer_t node, kmer_t& root, u_int64_t& root_hash) {

    u_int64_t kr = f.generate_KRHash_val_mod(node, this->k);
    u_int64_t hash = f.perfect_from_KR_mod(kr);
    u_int64_t parent;
    u_int64_t parent_kr;

    //BOOST_LOG_TRIVIAL(debug) << "Finding the root of " << node;

    if (this->fo.isStored(hash)) {
      //BOOST_LOG_TRIVIAL(debug) << "Root is " << node;
      root = node;
      root_hash = hash;
      return;
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
    }

    //BOOST_LOG_TRIVIAL(debug) << "Root is " << node;
    root = node;
    root_hash = hash;

  }

  /**
   * Given a node, get all of its children in the forest
   */
  void getChildren(const kmer_t& node, vector<kmer_t>& children) {

     children.clear();

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
         //BOOST_LOG_TRIVIAL(debug) << get_kmer_str(neighbors[i], this->k)
         //  << " is a child";
       }

       
     }
  

  }

  /**
   * Given a kmer, find the height of the tree that node is in in the forest
   */
  unsigned getTreeHeight(const kmer_t& node) {

     //BOOST_LOG_TRIVIAL(debug) << "Getting the height of tree with node "
     //   << get_kmer_str(node, this->k);

     // get the root of this tree
     kmer_t root;
     u_int64_t root_hash;
     getRoot(node, root, root_hash);

     return getTreeHeightRoot(root);
  }

  /**
   * Same as above, but it is assumed that the input node is the root of the tree
   */
  unsigned getTreeHeightRoot(const kmer_t& root) {

     //BOOST_LOG_TRIVIAL(debug) << "The root of this tree is node "
     //  << get_kmer_str(root, this->k);

     // keep descending down the tree by getting children until
     // there doesn't exist anymore. Count how many times this
     // must happen

     vector<kmer_t> children;
     getChildren(root, children);

     vector<kmer_t> children_children; // the children of the children
     vector<kmer_t> children_node; // holds the children of a single node

     unsigned height = 0;

     unsigned warnings = 0;

     // Until we have reached the most far node in the tree
     while (children.size() != 0) {

        // Get children of each child
        for (int i = 0; i < children.size(); i++) {

           //BOOST_LOG_TRIVIAL(debug) << "Getting the children of node "
           //  << get_kmer_str(children[i], this->k);

           // get the ith child's children
           getChildren(children[i], children_node);

           // add each one to children_children
           for (int j = 0; j < children_node.size(); j++) {

              //BOOST_LOG_TRIVIAL(debug) << "Found child "
              //  << get_kmer_str(children_node[j], this->k);
              children_children.push_back(children_node[j]);
           }
        }

        children = children_children;
        children_children.clear();

        height++;

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

