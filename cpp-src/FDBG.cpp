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

kmer_t getMaxVal(unsigned);
void push_last_letter(const kmer_t &, kmer_t &);
void remove_front_letter(const kmer_t &, kmer_t &, const unsigned &);
void set_kmer(kmer_t &mer, unsigned k, unsigned i, Letter &c);
void set_kmer(kmer_t &mer, unsigned k, unsigned i, char c);
kmer_t pushOnFront(kmer_t &orig, Letter &letter, unsigned);
kmer_t pushOnBack(kmer_t &orig, Letter &letter, unsigned);

// Representation of A, C, G, and T in bits 00, 01, 10, 11
class Letter
{

public:
  bool bits[2];
  void set_bits(char letter)
  {
    switch (letter)
    {
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

  void set(unsigned ii)
  {
  }

  Letter()
  {
    set_bits('A');
  }

  Letter(bool a, bool b)
  {

    this->bits[0] = a;
    this->bits[1] = b;
  }

  //  Letter(bool first, bool second) {
  //     bits[0] = first;
  //     bits[1] = second;
  //  }

  Letter(unsigned letter)
  {
    switch (letter)
    {
    case 0:
      set_bits('A');
      break;
    case 1:
      set_bits('C');
      break;
    case 2:
      set_bits('G');
      break;
    case 3:
      set_bits('T');
      break;
    }
  }

  Letter(char letter)
  {
    set_bits(letter);
  }

  // Return number 0..3 for this letter
  unsigned getNum()
  {
    return 2 * ((int)bits[0]) + ((int)bits[1]);
  }

  Letter &operator=(const Letter &rhs)
  {

    this->bits[0] = rhs.bits[0];
    this->bits[1] = rhs.bits[1];
    return *this;
  }
};

/**
 * Class that represents the forest
 * each group of 4 entries in the vector is a single node
 */
class Forest
{

public:
  /**
     * Each four pieces of data describe a node in the forest
     * One spot for whether it is stored as a root, one for whether
     * its parent can be reached by IN or not, and two for letter to
     * reach its parent by
     * The ith node is the kmer that hashes to i
     */
  vector<bool> data;

  // The number of nodes
  u_int64_t n;

  // The kmers of our roots
  map<u_int64_t, kmer_t> roots;

  /**Constructor*/
  Forest(u_int64_t n) : data(4 * n)
  {

    this->n = n;
  }

  /**
     * Set value of the ith node
     * Not stored
     */
  void setNode(u_int64_t &i, bool IN, Letter &l)
  {

    u_int64_t index = this->nodeIndex(i);

    this->data[index] = false;
    this->data[index + 1] = IN;

    this->data[index + 2] = l.bits[0];
    this->data[index + 3] = l.bits[1];
  }

  // Return index in data of the beginning of the ith node
  u_int64_t nodeIndex(u_int64_t i)
  {
    return 4 * i;
  }

  /**
     * Store the kmer of the ith node
     */
  void storeNode(u_int64_t i, kmer_t str)
  {

    this->roots[i] = str;

    this->data[nodeIndex(i)] = 1;
  }

  // Returns whether the ith node is stored
  bool isStored(u_int64_t i)
  {
    return this->data[nodeIndex(i)];
  }

  // Whether the ith node has IN
  bool parent_in_IN(u_int64_t i)
  {
    return this->data[nodeIndex(i) + 1];
  }

  // Deduce the parent of the ith node's kmer given the ith node's
  // kmer mer and the length of the kmers k
  kmer_t getNext(u_int64_t i, kmer_t &mer, unsigned k)
  {

    u_int64_t index = nodeIndex(i);
    Letter l(this->data[index + 2], this->data[index + 3]);

    if (this->data[index + 1])
    {
      // reach via IN
      return pushOnFront(mer, l, k);
    }
    else
    {
      return pushOnBack(mer, l, k);
    }
  }

  /**Print out a plain array*/
  void print()
  {

    int index = 0;

    for (int i = 0; i < this->n; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        cout << this->data[index];
        index++;
      }
      cout << endl;
    }
  }
};

/**
 * Class meant to represent IN or OUT
 * Has n rows, each row has 4 columns
 */
class INorOUT
{

private:
  vector<bool> data; // Holds entire matrix in one vector
  u_int64_t n;       // number of rows

public:
  INorOUT(u_int64_t n) : data(4 * n)
  {
    this->n = n;
  }

  vector<bool>::reference operator()(u_int64_t row, unsigned col)
  {

    u_int64_t index = 4 * row + col;

    return data[index];
  }

  void print()
  {

    int index = 0;

    for (int i = 0; i < this->n; i++)
    {
      for (int j = 0; j < 4; j++)
      {
        cout << this->data[index];
        index++;
      }
      cout << endl;
    }
  }
};

/**
 * This class represents the fully dynamic De Bruijn graph
 */
class FDBG
{

public:
  INorOUT IN;
  INorOUT OUT;
  unsigned sigma;  //alphabet-size. For now, only 4 is supported
  u_int64_t n;     //number of nodes in graph
  unsigned k;      //length of each mer (string in alphabet)
  generate_hash f; //hash function that takes each kmer to 1..n
  Forest fo;       // the forest NEW
  unsigned alpha;  //each tree in forest is guaranteed to be of height alpha to 3alpha

  /**
   * The number of bits that our data should be using
   */
  u_int64_t getBitSize()
  {
    u_int64_t res = 8 * n;                  //IN,OUT
    res += 4 * n + fo.roots.size() * 2 * k; //not too accurate, adding forest bits
    res += f.bphf->totalBitSize();
    return res;
  }

  FDBG(unordered_set<kmer_t> &kmers,
       unordered_set<kmer_t> &edgemers,
       u_int64_t n,           //number of kmers
       unsigned k,            //mer size
       bool b_verify = false, //if true, print summary
       ostream &os = cout) : IN(n), OUT(n), fo(n)
  {
    sigma = 4;

    this->n = n;
    this->k = k;

    //construct hash function f
    f.construct_hash_function(kmers, n, k);

    //printHashFunction(kmers);

    //initialize IN, OUT to zero (false)
    BOOST_LOG_TRIVIAL(info) << "Initializing IN and OUT.";

    //add edges to IN and OUT
    BOOST_LOG_TRIVIAL(info) << "Adding edges to IN and OUT ...";

    add_edges(edgemers);

    //For debugging, print IN and OUT given kmers
    //printINandOUT(kmers);

    //Perform the forest construction
    construct_forest(kmers, k * 2); //alpha = k * lg(sigma)

    //cout << "SIZE: " << this->fo.roots.size() << endl;

    //printHashFunction(kmers);
    //printForest();
  }

  /**
   * Add edges to IN and OUT using K+1-mers
   */
  void add_edges(unordered_set<kmer_t> &edges)
  {

    unordered_set<kmer_t>::iterator i;
    for (i = edges.begin(); i != edges.end(); ++i)
    {
      add_edge(*i);
    }
  }

  void add_edge(const kmer_t &edge)
  {

    // figure out the two kmers
    kmer_t u = 0;
    kmer_t v = 0;

    split_edge(edge, u, v);

    // get which column in sigma to put in (corresponds to which letter is the first/last)
    // number 0..3 represent each alphabet letter
    unsigned first, last;
    first = access_kmer(u, k, 0);
    last = access_kmer(v, k, k - 1);

    u_int64_t hash_u = f(u);
    u_int64_t hash_v = f(v);

    if (hash_u >= this->n)
    {
      BOOST_LOG_TRIVIAL(error) << "Edge " << get_kmer_str(edge, this->k + 1)
                               << " has produced kmer "
                               << get_kmer_str(u, this->k)
                               << "(" << u << ")"
                               << " with invalid hash function value "
                               << hash_u;
    }

    if (hash_v >= this->n)
    {
      BOOST_LOG_TRIVIAL(error) << "Edge " << get_kmer_str(edge, this->k + 1)
                               << " has produced kmer "
                               << get_kmer_str(v, this->k)
                               << "(" << v << ")"
                               << " with invalid hash function value "
                               << hash_v;
    }

    this->OUT(hash_u, last) = true;
    this->IN(hash_v, first) = true;
  }

  // Take a k+1-mer and split into beginning and end k-mers
  // u is the beginning, v is the end
  void split_edge(const kmer_t &edge, kmer_t &u, kmer_t &v)
  {

    push_last_letter(edge, u);

    remove_front_letter(edge, v, this->k);
  }

  /*
   * Add an edge dynamically to the data structure.
   * From u to v
   * Updates the forest dynamically.
   * Nothing happens if the k-mers aren't compatible.
   */
  void dynamicAddEdge(const kmer_t &u, const kmer_t &v)
  {
    //check if u, v are compatible
    unsigned ui, vi;
    for (unsigned i = 0; i < (k - 1); ++i)
    {
      ui = access_kmer(u, k, i + 1);
      vi = access_kmer(v, k, i);
      if (ui != vi)
        return;
    }

    //making it this far means that an edge can be added between them
    //Add the edge to OUT[ f(u) ] and to IN[ f(v) ]
    //if edge was already present, quit
    u_int64_t hashU = f(u);
    u_int64_t hashV = f(v);
    unsigned outIndex = access_kmer(v, k, k - 1);
    unsigned inIndex = access_kmer(u, k, 0);
    if (OUT(hashU, outIndex))
      return; // edge already exists

    //making it here means that edge is compatible and edge is not already in graph
    //So: begin logic for adding edge
    OUT(hashU, outIndex) = 1;
    IN(hashV, inIndex) = 1;

    //Now, need to update the forest
    //TODO

    return;
  }

  void print_matrix(vector<vector<bool>> &mat, ostream &os = cout)
  {
    for (unsigned i = 0; i < mat.size(); ++i)
    {
      for (unsigned j = 0; j < mat[i].size(); ++j)
      {
        os << mat[i][j] << ' ';
      }
      os << endl;
    }
  }

  void printHashFunction(unordered_set<kmer_t> &kmers)
  {

    unordered_set<kmer_t>::iterator i;
    for (i = kmers.begin(); i != kmers.end(); ++i)
    {
      cout << setw(10) << get_kmer_str(*i, this->k);
      cout << setw(10) << this->f(*i);
      cout << endl;
    }
  }

  // Given a kmer, decide if it is one in our graph
  bool detect_membership(kmer_t m)
  {
    //    BOOST_LOG_TRIVIAL(debug) << "Detecting membership of " << get_kmer_str(m, this->k);

    // The hash value of our kmer
    // Need to keep track of KRval, so it can be updated
    // largeUnsigned type is provided by 'generate_hash.h'
    largeUnsigned KR_val = f.generate_KRHash_raw(m, k);
    // u_int64_t KR_val = f.generate_KRHash_val_mod( m, k ) ;
    u_int64_t hash = f.perfect_from_KR(KR_val);

    //     BOOST_LOG_TRIVIAL(debug) << "Correct hash, computed: " << f(m) << ' ' << hash;

    //    BOOST_LOG_TRIVIAL(debug) << "It has hash value " << hash;

    // If it is a real kmer value, it must map to 0..1-n
    if (hash >= n)
    {
      //       BOOST_LOG_TRIVIAL(debug)  << "Returning false because of hash function blowup" << endl;
      return false;
    }

    //number of times we have traveled in the tree
    unsigned hopcount = 0;
    bool in;         //true = IN, false = OUT
    unsigned letter; //needed to confirm edge from parent's side

    // Keep going in the tree until we reach a hash that is stored as a root
    while (!(this->fo.isStored(hash)))
    {
      ++hopcount;
      if (hopcount > 3 * alpha + 10)
      {
        //BOOST_LOG_TRIVIAL( debug ) << "Returning false because of hopcount" << endl;
        return false; //we have encountered a cycle in our attempt to travel to a root
      }

      //do we progress with IN or OUT
      in = this->fo.parent_in_IN(hash);
      if (in)
      {
        //the parent is an IN-neighbor of m
        //so we need m's last character
        letter = access_kmer(m, k, k - 1);
      }
      else
      {
        //the parent is an OUT-neighbor of m
        //so we need m's first character
        letter = access_kmer(m, k, 0);
      }

      //      BOOST_LOG_TRIVIAL(debug) << "Child k-mer: " << get_kmer_str( m, k );

      // deduce the parent's kmer
      m = this->fo.getNext(hash, m, k);

      //      BOOST_LOG_TRIVIAL(debug) << "Parent k-mer: " << get_kmer_str( m, k );

      // get the parent's hash
      if (in)
      {
        unsigned letter_front_parent = access_kmer(m, k, 0);
        f.update_KRHash_val_IN(KR_val, letter_front_parent, letter);
        hash = f.perfect_from_KR(KR_val);
      }
      else
      {
        unsigned letter_back_parent = access_kmer(m, k, k - 1);
        f.update_KRHash_val_OUT(KR_val, letter, letter_back_parent);
        hash = f.perfect_from_KR(KR_val);
      }
      //hash = f(m);
      //      BOOST_LOG_TRIVIAL(debug) << "Correct hash, computed: " << f(m) << ' ' << hash;
      //            BOOST_LOG_TRIVIAL(debug) << "Correct, computed KR_val: "
      //            			       << f.generate_KRHash_val_mod( m, k ) << ' ' << KR_val;

      // hash must be in 0...n-1
      if (hash >= n)
      {
        //	 BOOST_LOG_TRIVIAL(debug) << "Returning false because of hash function blowup" << endl;
        return false;
      }
      //confirm the edge from parent side
      if (in)
      {
        if (!(OUT(hash, letter)))
        {
          //BOOST_LOG_TRIVIAL(debug) << "Returning false because IN,OUT verification failed" << endl;
          return false;
        }
      }
      else
      {
        if (!(IN(hash, letter)))
        {
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

  // Given a kmer, decide if it is one in our graph
  // recomputes hash value at each step, so for testing purposes only
  bool inefficient_detect_membership(kmer_t m)
  {
    //    BOOST_LOG_TRIVIAL(debug) << "Detecting membership of " << get_kmer_str(m, this->k);

    // The hash value of our kmer
    // Need to keep track of KRval, so it can be updated
    u_int64_t hash = f(m);

    //    BOOST_LOG_TRIVIAL(debug) << "It has hash value " << hash;

    // If it is a real kmer value, it must map to 0..1-n
    if (hash >= n)
    {

      return false;
    }

    //number of times we have traveled in the tree
    unsigned hopcount = 0;
    bool in;         //true = IN, false = OUT
    unsigned letter; //needed to confirm edge from parent's side

    while (!(this->fo.isStored(hash)))
    {
      ++hopcount;
      if (hopcount > 3 * alpha + 10)
      {

        return false; //we have encountered a cycle in our attempt to travel to a root
      }

      //do we progress with IN or OUT
      in = this->fo.parent_in_IN(hash);
      if (in)
      {
        //the parent is an IN-neighbor of m
        //so we need m's last character
        letter = access_kmer(m, k, k - 1);
      }
      else
      {
        //the parent is an OUT-neighbor of m
        //so we need m's first character
        letter = access_kmer(m, k, 0);
      }

      // deduce the parent's kmer
      m = this->fo.getNext(hash, m, k);
      hash = f(m);

      // hash must be in 0...n-1
      if (hash >= n)
      {

        return false;
      }
      //confirm the edge from parent side
      if (in)
      {
        if (!(OUT(hash, letter)))
        {

          return false;
        }
      }
      else
      {
        if (!(IN(hash, letter)))
        {

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
  void construct_forest(unordered_set<kmer_t> &kmers, int alpha)
  {
    BOOST_LOG_TRIVIAL(info) << "Beginning the forest construction, with alpha = " << alpha;
    this->alpha = alpha;
    // kmers that we have looked at for the forest construction
    unordered_set<kmer_t> visited_mers;

    vector<int> h(n, -1);    // height for each node below its tree root
    vector<kmer_t> p1(n, 0); // p1,p2 needed to tell when to store a tree root
    vector<kmer_t> p2(n, 0); //
    vector<kmer_t> p(n, 0);  // parent in BFS

    // Do a BFS through each of the UNDIRECTED graph components
    // Keep going until we have looked at all kmers
    // This will only have one loop unless the graph is not connected
    while (visited_mers.size() != n)
    {

      // Pick initial root, will be stored in our forest
      kmer_t root = *kmers.begin();
      store(root);

      //      BOOST_LOG_TRIVIAL(debug) << "Building forest from root " + get_kmer_str(root, k);

      // We have visited root
      move_kmer(kmers, visited_mers, root);

      // The hash value of the root
      u_int64_t r = f(root);

      // See pseudocode for p1,p2
      p1[r] = root;
      p2[r] = root;

      // Root has height 0
      h[r] = 0;

      // Now do a BFS around the root in order to put all nodes into a tree
      queue<kmer_t> Q;
      Q.push(root);

      // Kmers of this nodes's neighbors, and if they are in or out neighbors
      vector<kmer_t> neis;
      vector<bool> v_inorout;

      // BFS
      while (!Q.empty())
      {
        // Next k-mer to visit neighbors of
        kmer_t c = Q.front();
        Q.pop();
        //        BOOST_LOG_TRIVIAL(debug) << "Visiting neighbors of node " + get_kmer_str(c, k);

        // Neighbors of c
        get_neighbors(c, neis, v_inorout);

        // If a neighbor can be reached via OUT from c, this is the letter to reach it
        Letter first_c = access_kmer(c, k, 0);
        // Now for IN
        Letter last_c = access_kmer(c, k, k - 1);

        for (unsigned ii = 0; ii < neis.size(); ++ii)
        {
          kmer_t m = neis[ii]; //this is 'n' in the pseudocode
          //	  BOOST_LOG_TRIVIAL(debug) << "Checking neighbor m: " + get_kmer_str(m, k);
          if (visited_mers.find(m) == visited_mers.end())
          {
            //haven't visited m yet
            Q.push(m);
            move_kmer(kmers, visited_mers, m);
            //	    BOOST_LOG_TRIVIAL(debug) << "Number of kmers visited: " << visited_mers.size();
            u_int64_t f_m = f(m); //save these values so we don't have to recompute all the time
            u_int64_t f_c = f(c);
            p[f_m] = c;

            if (v_inorout[ii])
            {
              /**
               * The iith neighbor of c was reached by c via IN. Therefore to reach
               * c from the iith neighbor, we must go OUT. We will need the last
               * character of c
               */
              fo.setNode(f_m, false, last_c);
            }
            else
            {
              // similar logic to the first case
              fo.setNode(f_m, true, first_c);
            }

            h[f_m] = h[f_c] + 1;
            int height_m = h[f_m];
            if (height_m <= alpha)
            {
              p1[f_m] = p1[f_c];
              p2[f_m] = p2[f_c];
            }
            if ((alpha < height_m) && (height_m <= 2 * alpha))
            {
              store(p1[f_c]);
              p1[f_m] = p1[f_c];
              p2[f_m] = p1[f_c];
            }
            if (height_m == (2 * alpha + 1))
            {
              h[f_m] = 0;
              p1[f_m] = m;
              p2[f_m] = p1[f_c];
            }
          }
        }
      }
    }

    kmers.swap(visited_mers);
  }

  /**
   * Given a kmer c, get all neighbor kmers (neighbor_kmers) 
   * for neighbor_kmers[ i ], v_inorout[ i ] is true if in-neighbor, false otherwise
   */
  void get_neighbors(kmer_t &c,
                     vector<kmer_t> &neighbor_kmers,
                     vector<bool> &v_inorout)
  {

    neighbor_kmers.clear();
    v_inorout.clear();

    //    BOOST_LOG_TRIVIAL(debug) << "Finding neighbors of " + get_kmer_str(c, k);

    // Need hash value of our kmer to look into IN and OUT
    u_int64_t fc = f(c);
    //    BOOST_LOG_TRIVIAL(trace) << "c, f(c): " << get_kmer_str(c, k  ) << ' ' << fc;
    //    BOOST_LOG_TRIVIAL(trace) << "Size of IN,OUT: " << IN[ fc ].size() << ' ' << OUT[ fc ].size();
    for (unsigned ii = 0; ii < 4; ++ii)
    {
      if (IN(fc, ii))
      {
        //have in-neighbor with letter
        Letter letter(ii);
        //	BOOST_LOG_TRIVIAL(trace) << "letter: " << ii << ' ' << letter.getNum();
        // Deduce that kmer by tacking letter at the beginning of c
        kmer_t e = pushOnFront(c, letter, k);
        neighbor_kmers.push_back(e);

        //	BOOST_LOG_TRIVIAL(trace) << "neighbor: " + get_kmer_str(e, k);

        v_inorout.push_back(true); //true=IN
      }
      if (OUT(fc, ii))
      {
        //have out-neighbor with letter
        Letter letter(ii);

        // Deduce that kmer by tacking letter at the end of c
        kmer_t e = pushOnBack(c, letter, k);
        neighbor_kmers.push_back(e);

        v_inorout.push_back(false); //false=OUT
      }
    }

    //    BOOST_LOG_TRIVIAL(trace) << "Neighbors found: ";
    for (unsigned i = 0; i < neighbor_kmers.size(); ++i)
    {
      //      BOOST_LOG_TRIVIAL(trace) << i << ' ' << get_kmer_str( neighbor_kmers[i], k ) << " is an "<<	v_inorout[ i ] << " neighbor.";
    }
  }

  void store(kmer_t mer)
  {
    u_int64_t val = f(mer);
    this->fo.storeNode(val, mer);
  }

  // Move mer from kmers and into visited
  void move_kmer(unordered_set<kmer_t> &kmers, unordered_set<kmer_t> &visited,
                 kmer_t mer)
  {
    kmers.erase(mer);
    visited.insert(mer);
  }

  // Print IN and OUT with the Kmer strings down the side, and the characters on the top
  void printINandOUT(unordered_set<kmer_t> &kmers)
  {

    cout << setw(30) << "===== IN =====" << endl
         << endl;

    cout << setw(10) << "";
    cout << setw(5) << "A";
    cout << setw(5) << "C";
    cout << setw(5) << "G";
    cout << setw(5) << "T" << endl;

    unordered_set<kmer_t>::iterator i;
    for (i = kmers.begin(); i != kmers.end(); ++i)
    {
      cout << setw(10) << get_kmer_str(*i, this->k);

      u_int64_t hash = this->f(*i);
      cout << setw(5) << this->IN(hash, 0);
      cout << setw(5) << this->IN(hash, 1);
      cout << setw(5) << this->IN(hash, 2);
      cout << setw(5) << this->IN(hash, 3);
      cout << endl;
    }

    cout << endl;
    cout << setw(30) << "===== OUT =====" << endl
         << endl;

    cout << setw(10) << "";
    cout << setw(5) << "A";
    cout << setw(5) << "C";
    cout << setw(5) << "G";
    cout << setw(5) << "T" << endl;

    for (i = kmers.begin(); i != kmers.end(); ++i)
    {
      cout << setw(10) << get_kmer_str(*i, this->k);

      u_int64_t hash = this->f(*i);
      cout << setw(5) << this->OUT(hash, 0);
      cout << setw(5) << this->OUT(hash, 1);
      cout << setw(5) << this->OUT(hash, 2);
      cout << setw(5) << this->OUT(hash, 3);
      cout << endl;
    }
  }

  // Print the forest with the Kmer strings down the side
  void printForest()
  {

    cout << setw(30) << "===== FOREST =====" << endl
         << endl;

    for (int i = 0; i < this->n; ++i)
    {
      cout << setw(10) << i;

      u_int64_t index = this->fo.nodeIndex(i);

      cout << setw(5) << this->fo.data[index];
      cout << setw(5) << this->fo.data[index + 1];
      cout << setw(5) << this->fo.data[index + 2];
      cout << setw(5) << this->fo.data[index + 3];
      cout << endl;
    }
  }

  // END CLASS
};

/*
 * Sets the i'th position of a mer of length k as indicated by character c
 * c \in {A,C,G,T}
 */
void set_kmer(kmer_t &mer, unsigned k, unsigned i, Letter &c)
{
  //clear i-th position
  kmer_t op = 3;              //0...011
  op = op << 2 * (k - i - 1); //11 in i'th spot, zeros elsewhere
  op = ~op;                   //00 in i'th spot, ones elsewhere
  mer = mer & op;             //i'th position of mer cleared.

  //set i'th position
  kmer_t val = c.getNum();

  val = val << 2 * (k - i - 1); //correct bits in i'th spot, zeros elsewhere
  mer = mer | val;
}

// Push letter onto front of kmer and return kmer
// Going backwards along an edge (the relationship comes from orig's IN).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnFront(kmer_t &orig, Letter &letter, unsigned k)
{
  //  BOOST_LOG_TRIVIAL(trace) << "pushOnFront " + get_kmer_str(orig, k) << ' ' << letter.getNum();
  // Get the kmer with back pushed off orig
  kmer_t new_kmer = orig >> 2;

  set_kmer(new_kmer, k, 0, letter);

  return new_kmer;
}

// Push letter onto back of kmer and return kmer
// Going forward along an edge (the relationship comes from orig's OUT).
// For example, orig = AGCT, then it returns GAGC
kmer_t pushOnBack(kmer_t &orig, Letter &letter, unsigned k)
{
  //  BOOST_LOG_TRIVIAL(trace) << "pushOnBack " + get_kmer_str(orig, k) << ' ' << letter.getNum();
  // Get the kmer with front pushed off orig
  kmer_t new_kmer = orig << 2;
  set_kmer(new_kmer, k, k - 1, letter);
  //but we need to zero the -1th spot
  //clear -1th position
  kmer_t op = 3;            //0...011
  op = op << 2 * (k);       //11 in -1'st spot, zeros elsewhere, no issues even if k=32
  op = ~op;                 //00 in i'th spot, ones elsewhere
  new_kmer = new_kmer & op; //i'th position of mer cleared.

  return new_kmer;
}

/*
 * Sets the i'th position of a mer of length k as indicated by character c
 * c \in {A,C,G,T}
 */
void set_kmer(kmer_t &mer, unsigned k, unsigned i, char c)
{
  //clear i-th position
  kmer_t op = 3;              //0...011
  op = op << 2 * (k - i - 1); //11 in i'th spot, zeros elsewhere
  op = ~op;                   //00 in i'th spot, ones elsewhere
  mer = mer & op;             //i'th position of mer cleared.

  //set i'th position
  kmer_t val;
  switch (c)
  {
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

  val = val << 2 * (k - i - 1); //correct bits in i'th spot, zeros elsewhere
  mer = mer | val;
}

/**
 * Given a kmer length k, returns the maximum value that kmer
 * is allowed to be
 */
kmer_t getMaxVal(unsigned k)
{

  kmer_t max = 0;

  for (int i = 0; i < k; ++i)
  {
    max <<= 2;
    max |= 3;
  }

  return max;
}

// Push off the last letter
void push_last_letter(const kmer_t &edge, kmer_t &u)
{

  u = edge >> 2;
}

// Remove front letter from kmer of length k+1
// To make it a kmer of length k
void remove_front_letter(const kmer_t &edge, kmer_t &v, const unsigned &k)
{

  //v = edge & ~(3 << 2*(k));

  kmer_t op = 3;
  op = op << 2 * (k);
  op = ~op;
  v = edge & op;

  //BOOST_LOG_TRIVIAL(debug) << "Removed front letter of " << edge << " to get " << v;
}

#endif
