#pragma once

#include "BooPHF.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>


using namespace std;

typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
typedef boomphf::mphf<  u_int64_t, hasher_t  > boophf_t;

typedef short kmer_t;

// I leave isPrime and getPrime outside for the usage of others.

// Task 2: Implement the getPrime function.
template <typename T>
int isPrime(T n) // assuming n > 1
{
    T i,root;

    if (n%2 == 0 || n%3 == 0)
        return 0;

    root = (T)sqrt(n);

    for (i=5; i<=root; i+=6)
    {
        if (n%i == 0)
           return 0;
    }

    for (i=7; i<=root; i+=6)
    {
        if (n%i == 0)
           return 0;
    }

    return 1;
}

template <typename T>
T getPrime(T n){

    T c = 1;
    while (1){
        if (isPrime(n+c)) { 
            return n+c;
        }
        c += 1;
    }
    return 0;
}

// Task 3: Implement the randomNumber function.
// Generate a random number between [start, end]. 
// Note that start and end is included.
int randomNumber(const int start, const int end){

    double myRand = rand()/(1.0 + RAND_MAX); 
    int range = end - start + 1;
    int random_num = (myRand * range) + start;
    return random_num;
}

// Task 1-6: The main algorithm
class hash_generator{

    public:
        //data is supposed a 2-d metrix [num_kemr, k]
        hash_generator(kmer_t *data[], int n, int k)
        :n_kmer(n), k_kmer(k)
        {
            nthreads = 4;
            gammaFactor = 2.0;

            kmer_data = data;
            KR_hash_val = (u_int64_t *) calloc(n_kmer, sizeof(u_int64_t));
            memset(KR_hash_val, 0, n_kmer * sizeof(u_int64_t));
        
            printf("Generate hash function ... \n");

            build_minimalPerfectHash();
        }
        u_int64_t get_hash_value(kmer_t * seq)
        {
            u_int64_t krv = generate_KRHash_val(seq, k_kmer, base, Prime);
            u_int64_t res = bphf->lookup(krv);
            return res;
        }
    private:
        // Task4: generate_KRHash_val
        // data is a k-mer
        // k is the length of the k-mer
        // r is the base 
        // P is the prime
        void build_KRHash(){
            printf("Build rabin hash function ... \n");
            int R, P, r;
            u_int64_t v;
            int sigma = 4; // TODO: what is the value of this?
            R = max((u_int64_t)4, (u_int64_t)k_kmer*n_kmer*n_kmer);
            P = getPrime(R);
            // Find satisfied base and prime combination
            memset(KR_hash_val, 0, n_kmer * sizeof(u_int64_t));
            while (1)
            {
                r = randomNumber(0, P-1);     
                // compute rabin hash for each kmer
                for (int i = 0; i < n_kmer; i++)
                { 
                    v = generate_KRHash_val(kmer_data[i], k_kmer, r, P);
                    // printf("generate hash code ... %d\n", v );
                    if (KR_hash_val[i] == 0)
                    {
                        KR_hash_val[i] = v;
                    }
                    else // not injective
                    {
                        printf("Testing base=%d Prime=%d fails", r, P);
                        memset(KR_hash_val, 0, n_kmer * sizeof(u_int64_t));
                        break;
                    }
                }
                printf("Testing base=%d Prime=%d succeed!\n", r, P);
                base = r;
                Prime = P;
                return;
            }
        }

        kmer_t generate_KRHash_val(kmer_t *kmer, const int k, const int r, int P){
            
            u_int64_t val = 0;
            // TODO: The wiki said the exponent of r is in the decresing order, which is different from our paper's.
            // I follow the former for now.
            for (int i = 0; i < k; i++)
            {
                val += kmer[i] * pow(r, k-i-1);
                // printf("generate kmer ... %d\n", kmer[i]);
            }
            // printf("generate mode ... %d\n", val );

            val = val % P;

            return u_int64_t(val);
        }

        void build_minimalPerfectHash(){

            build_KRHash(); // compute KR_hash_val

            std::sort(KR_hash_val, KR_hash_val+n_kmer);
            u_int64_t jj = 0;
            for (int ii = 1; ii < n_kmer; ii++) {
                if (KR_hash_val[ii] != KR_hash_val[jj])
                    KR_hash_val[++jj] = KR_hash_val[ii];
            }
            printf("Found %lli duplicated items from KR_hash_val.  \n", n_kmer-(jj + 1) );

            auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(KR_hash_val), static_cast<const u_int64_t*>(KR_hash_val+n_kmer));

            bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, data_iterator, nthreads, gammaFactor);

            printf("The minimal perfect hashing function is generated. \n");
            printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/n_kmer);
        }

            
    public:
        u_int64_t n_kmer; 
        u_int64_t k_kmer;

        u_int64_t *KR_hash_val = NULL;
        kmer_t **kmer_data = NULL;
        boophf_t *bphf = NULL;

        double gammaFactor; 
        u_int64_t base;
        u_int64_t Prime;
        unsigned int nthreads;
};
