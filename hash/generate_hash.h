#include "BooPHF.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <unordered_set>
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
        u_int64_t k_kmer; //the lengths of the k-mers

        vector<string> kmer_data; // pointer to kmer_data TODO: should read from file

        std::unordered_set<u_int64_t> KRHash; //image of k-mers through Karp-Rabin hash function
        //std::vector<u_int64_t> KRHash_vec; //vector form

        //the image of our k-mers under our Karp-Rabin hash function
        //vector<u_int64_t> KR_hash_val;

        boophf_t* bphf; //MPHF we will generate

        //static double gammaFactor = 2.0; 
        u_int64_t base; // the base for our Karp-Rabin Hash function
        u_int64_t Prime; // the prime for our Karp-Rabin Hash function

        const static unsigned int nthreads = 4; // # of threads for BBHash

        const static short sigma = 4; // alphabet size

        /**
         * Create hash function out of n k-mers of length k
         */
        generate_hash( unordered_set< kmer_t >& kmers , int n , int k ) {
	  std::srand(std::time(NULL));
	  
            this->n_kmer = n; // number of k-mers
            this->k_kmer = k; // length of each k-mer

            //KR_hash_val = (u_int64_t *) calloc(n_kmer, sizeof(u_int64_t)); // space to store KRH image

            //memset(KR_hash_val, 0, n_kmer * sizeof(u_int64_t)); // initialize to 0
        
            printf("Generate hash function ... \n");

            build_KRHash(kmers); // build KR hash function
            build_minimalPerfectHash(); // build minimal perfect hash function
        }


        /**
         * Find the hash value of a k-mer
         */
        u_int64_t get_hash_value(const kmer_t& seq)
        {
            u_int64_t krv = generate_KRHash_val(seq, k_kmer, base, Prime);
            u_int64_t res = this->bphf->lookup(krv);
            return res;
        }

    private:
        // Task4: generate_KRHash_val
        // data is a k-mer
        // k is the length of the k-mer
        // r is the base 
        // P is the prime
        void build_KRHash( unordered_set< kmer_t >& kmers ){
            printf("Build rabin hash function ... \n");

            u_int64_t r; // our base, which we will figure out in loop
            u_int64_t v; // holder for KRH value

            // prime we will mod out by
            u_int64_t P = getPrime(max((u_int64_t)this->sigma, (u_int64_t)k_kmer*n_kmer*n_kmer));

            // Find satisfied base and prime combination
            //memset(KR_hash_val, 0, n_kmer * sizeof(u_int64_t));

            // keep generating new base until we find one that is injective over our k-mers
	    bool f_injective;
	    do
            {
	      f_injective = true; //assume f is injective until evidence otherwise
	      r = randomNumber((u_int64_t) 1, P - 1);     
	      for ( unordered_set< kmer_t >::iterator
		      it1 = kmers.begin(); it1 != kmers.end();
		    ++it1 ) {
		v = generate_KRHash_val( *it1, k_kmer, r, P);

		if (this->KRHash.find(v) == this->KRHash.end())
		  {
		    // this is a new value   
		    this->KRHash.insert(v);
		  }
		else // not injective
		  {
		    printf("Testing base=%d Prime=%d fails...\n", r, P);
		    this->KRHash.clear(); // clear it out and start over
		    f_injective = false;
		    break;
		  }
		
	      }
	    } while (!f_injective);

	    
	    printf("Testing base=%d Prime=%d succeed!\n", r, P);
	    this->base = r;
	    this->Prime = P;

	    return;
	}
        

        /**
         * Given a kmer, find out its KRH using base r and prime P
         */
        u_int64_t generate_KRHash_val(const kmer_t& kmer,
				      const int k, const int r, int P){
            
            u_int64_t val = 0; // what will be the KRH value

            // go through each bp and add value
            for (int i = k - 1;
		 i >= 0; i--)
            {
	      // val += baseNum(kmer.at(i)) * pow(r, i);
	      val += access_kmer( kmer, static_cast<unsigned>(k), static_cast<unsigned>(i)) * pow(r, i);
            }

            val = val % P;

            return u_int64_t(val);
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

            std::vector<u_int64_t> KRHash_vec = std::vector<u_int64_t>(this->KRHash.begin(),
               this->KRHash.end());

            // MPHF for our KRHash function values
            this->bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, KRHash_vec);

            printf("The minimal perfect hashing function is generated. \n");
            printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/n_kmer);
        }

            
};
