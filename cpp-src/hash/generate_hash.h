#ifndef GEN_HASH
#define GEN_HASH

#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <unordered_set>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;
#include "HashUtil.cpp"

using namespace std;

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

/**
 * Take in a set of k-mers, generate a hash function
 */
class generate_hash {

    public:
        u_int64_t n_kmer; //the number of k-mers 
        unsigned k_kmer; //the lengths of the k-mers (max 32)

        vector<string> kmer_data; // pointer to kmer_data TODO: should read from file

        std::unordered_set<u_int64_t> KRHash; //image of k-mers through Karp-Rabin hash function
        //std::vector<u_int64_t> KRHash_vec; //vector form

        //the image of our k-mers under our Karp-Rabin hash function
        //vector<u_int64_t> KR_hash_val;

	/* 
	 * Will store the precomputed powers of 
	 * Needs to be a vector will have k distinct powers
	 * Using 128 bit type to prevent overflow
	 * Stores powers in order r^1, r^2, ..., r^k
	 */
	vector< uint128_t > powersOfR; 

        boophf_t* bphf; //MPHF we will generate

        u_int64_t r; // the base for our Karp-Rabin Hash function
        u_int64_t Prime; // the prime for our Karp-Rabin Hash function

        const static short sigma = 4; // alphabet size

        /**
         * Create hash function out of n k-mers of length k
         */
        generate_hash( unordered_set< kmer_t >& kmers , u_int64_t n , unsigned k ) {
	  std::srand(std::time(NULL));
	  
	  construct_hash_function( kmers,  n,  k );
        }

	generate_hash() {
	  //default constructor
	  std::srand(std::time(NULL));
	}

	void construct_hash_function( unordered_set< kmer_t >& kmers , u_int64_t n , unsigned k ) {
	  this->n_kmer = n; // number of k-mers
	  this->k_kmer = k; // length of each k-mer
	  
          BOOST_LOG_TRIVIAL(info) << "Constructing the hash function ...";
	  build_KRHash(kmers); // build KR hash function
	  build_minimalPerfectHash(); // build minimal perfect hash function
	  
	}

	/*
	 * Once we know k (k_kmer) and r
	 * we can precompute the powers of r
	 */
	void precomputePowers() {
	  uint128_t ri;
	  powersOfR.clear();
	  //NEED 1 to k. Not 0 to (k - 1)
	  for (unsigned i = 1; i <= k_kmer; ++i) {
	    ri = mypower( r, i );
	    powersOfR.push_back( ri );
	  }
	}
	

        /**
         * Find the hash value of a k-mer
         */
        u_int64_t get_hash_value(const kmer_t& seq)
        {
            u_int64_t krv = generate_KRHash_val(seq, k_kmer, Prime);
            u_int64_t res = this->bphf->lookup(krv); // still need only 64 bits for kmer_t
            return res;
        }

	/**
         * Find the hash value of a k-mer
	 * Allows f( v ) notation
         */
	u_int64_t operator()(const kmer_t& seq) {
	  return get_hash_value( seq );
	}
					       

        // Task4: generate_KRHash_val
        // data is a k-mer
        // k is the length of the k-mer
        // r is the base 
        // P is the prime
        void build_KRHash( unordered_set< kmer_t >& kmers ){

            BOOST_LOG_TRIVIAL(info) << "Constructing Karp-Rabin hash function ...";

            u_int64_t v; // holder for KRH value

            // prime we will mod out by
	    const u_int64_t tau = 1;
	    BOOST_LOG_TRIVIAL(debug) << "Minimum prime: " << tau*k_kmer*n_kmer*n_kmer;
            u_int64_t P = getPrime(max((u_int64_t)this->sigma, (u_int64_t)tau*k_kmer*n_kmer*n_kmer));

            // Find satisfied base and prime combination
            //memset(KR_hash_val, 0, n_kmer * sizeof(u_int64_t));

            // keep generating new base until we find one that is injective over our k-mers
	    bool f_injective;
	    do
	      {
	      f_injective = true; //assume f is injective until evidence otherwise
	      this->r = randomNumber((u_int64_t) 1, P-1);
	      //Once we have a candidate base r
	      //we should avoid recomputing its powers all the time
	      precomputePowers();
	      
	      for ( unordered_set< kmer_t >::iterator
		      it1 = kmers.begin(); it1 != kmers.end();
		    ++it1 ) {
		v = generate_KRHash_val( *it1, k_kmer, P);
		//		BOOST_LOG_TRIVIAL(trace) << "hash of kmer: " << v;
		if (this->KRHash.find(v) == this->KRHash.end())
		  {
		    // this is a new value   
		    this->KRHash.insert(v);
		  }
		else // not injective
		  {
                    BOOST_LOG_TRIVIAL(trace) << "Base " <<this->r << " with prime "
                       << P << " failed injectivity.";
		    this->KRHash.clear(); // clear it out and start over
		    f_injective = false;
		    break;
		  }
		
	      }
	    } while (!f_injective);

	    
            BOOST_LOG_TRIVIAL(info) << "Base " << this->r << " with prime " << P << " is injective.";
	    this->Prime = P;

	    return;
	}
        

	/*
	 * Computes powers with u_int64_t and integer exponents
	 */
	uint128_t mypower( const u_int64_t& base, unsigned exponent ) {
	  uint128_t rvalue( 1 );
	  while (exponent > 0) {
	    rvalue *= static_cast< uint128_t >(base);
	    --exponent;
	  }

	  return rvalue;
	}

	
        /**
         * Given a kmer, find out its KRH using base r and prime P
         */
        u_int64_t generate_KRHash_val(const kmer_t& kmer,
				      const unsigned& k,
				      const u_int64_t& P){
	  //	  BOOST_LOG_TRIVIAL(trace) << "Generating KRHash val";
	  //use 128 bits to prevent overflow
	  uint128_t val = 0; // what will be the KRH value

            // go through each bp and add value
	  for (unsigned i = 0;
	       i < k;
	       ++i) {
	      // val += baseNum(kmer.at(i)) * pow(r, i);
	      val +=
		static_cast< uint128_t > ( access_kmer( kmer, k, static_cast<unsigned>(i)) ) *
		powersOfR[i]; //powersOfR[i] = r^{i + 1}
	  }

	  val = val % P;

	  return static_cast<u_int64_t>(val);
        }

	/*
	 * This function takes as input a Karp-Rabin value (KR_val)
	 * Then updates it by subtracting the value from 'first' character source kmer
	 * Then dividing by r (at this point, it has shifted last k-1 characters up)
	 * Finally adding the last term corresponding to the 'last' character
	 *
	 * target k-mer is OUT neighbor of source k-mer
	 *
	 */
	u_int64_t update_KRHash_val_OUT
	  ( u_int64_t& KR_val_in,       //KR hash of source kmer
	    const unsigned& first,   //character at front of source k-mer
	    const unsigned& last ) { //last character in target k-mer
	  uint128_t KR_val = KR_val_in;
	  KR_val = KR_val - first * r;
	  KR_val = KR_val / r;
	  KR_val = KR_val + last * powersOfR[ k_kmer - 1 ]; // last * r^k
	  KR_val = KR_val % Prime;
	  return static_cast< u_int64_t >( KR_val);
	}

	/*
	 * This function takes as input a Karp-Rabin value (KR_val)
	 *
	 * target k-mer is IN neighbor of source k-mer
	 */
	u_int64_t update_KRHash_val_IN
	  ( u_int64_t& KR_val_in,       //KR hash of source kmer
	    const unsigned& first,   //character at front of target k-mer
	    const unsigned& last ) { //last character in source k-mer
	  uint128_t KR_val = KR_val_in;
	  KR_val = KR_val - last * powersOfR[ k_kmer - 1 ]; // last * r^k
	  KR_val = KR_val * r;
	  KR_val = KR_val + first * r;
	  KR_val = KR_val % Prime;
	  return static_cast< u_int64_t >( KR_val);
	}

	
	/*
	 * Looks up the minimal perfect hash value, given the Karp-Rabin
	 *
	 */
	u_int64_t perfect_from_KR( const u_int64_t& KR_val ) {
	  return this->bphf->lookup( KR_val );
	}

        /**
         * Build a minimal perfect hash function on the set of integers that our kmers are
         * mapped to via KRH
         */
        void build_minimalPerfectHash(){

            //std::sort(KR_hash_val, KR_hash_val+n_kmer);
            //u_int64_t jj = 0;
            //for (int ii = 1; ii < n_kmer; ii++) {
            //    if (KR_hash_val[ii] != KR_hash_val[jj])
            //        KR_hash_val[++jj] = KR_hash_val[ii];
            //}
            //printf("Found %lli duplicated items from KR_hash_val.  \n", n_kmer-(jj + 1) );

            //auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(KR_hash_val), static_cast<const u_int64_t*>(KR_hash_val+n_kmer));

	  //            bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, data_iterator, nthreads, gammaFactor);

            BOOST_LOG_TRIVIAL(info) << "Building minimal perfect hash function ...";

            std::vector<u_int64_t> KRHash_vec = std::vector<u_int64_t>(this->KRHash.begin(),
               this->KRHash.end());

            // MPHF for our KRHash function values
            this->bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, KRHash_vec, 4, 2.0, true, false);

            BOOST_LOG_TRIVIAL(info) << "Minimal perfect hash function created with "
               << (float) (bphf->totalBitSize())/n_kmer << " bits per element.";
        }

            
};

#endif
