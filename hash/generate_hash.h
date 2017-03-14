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
class HASH{

    public:
        //data is supposed a 2-d metrix [num_kemr, k]
        HASH(kmer_t *data[], int n, int k)
            :n_kmer(n), k_kmer(k)
        {
            nthreads = 4;
            gammaFactor = 2.0;

            kmer_data = data;
            rab_hash_v = (u_int64_t *) calloc(n_kmer, sizeof(u_int64_t));
            memset(rab_hash_v, 0, n_kmer * sizeof(u_int64_t));
        
            printf("Generate hash function ... \n");

            build_minimalPerfectHash();
        }

    private:
        // Task4: rabinHash
        // data is a k-mer
        // k is the length of the k-mer
        // r is the base 
        // P is the prime
        void build_rabinHash(){
            printf("Build rabin hash function ... \n");
            int R, P, r;
            u_int64_t v;
            int sigma = 4; // TODO: what is the value of this?
            R = max((u_int64_t)4, (u_int64_t)k_kmer*n_kmer*n_kmer);
            P = getPrime(R);
            // Find satisfied base and prime combination
            while (1)
            {
                r = randomNumber(0, P-1);     
                // compute rabin hash for each kmer
                for (int i = 0; i < n_kmer; i++)
                { 
                    v = rabinHash(kmer_data[i], k_kmer, r, P);
                    if (rab_hash_v[i] == 0)
                    {
                        rab_hash_v[i] = v;
                    }
                    else // not injective
                    {
                        printf("Testing base=%d Prime=%d fails", r, P);
                        memset(rab_hash_v, 0, n_kmer * sizeof(u_int64_t));
                        break;
                    }
                }
                printf("Testing base=%d Prime=%d succeeds !", r, P);
                base = r;
                Prime = P;
                return;
            }
        }

        kmer_t rabinHash(kmer_t *kmer, const int k, const int r, int P){
            
            u_int64_t val = 0;
            // TODO: The wiki said the exponent of r is in the decresing order, which is different from our paper's.
            // I follow the former for now.
            for (int i = 0; i < k; i++)
            {
                val += kmer[i] * pow(r, k-i-1);
            }
            val = val % P;

            return u_int64_t(val);
        }

        void build_minimalPerfectHash(){

            build_rabinHash();

            std::sort(rab_hash_v, rab_hash_v+n_kmer);
            u_int64_t jj = 0;
            for (int ii = 1; ii < n_kmer; ii++) {
                if (rab_hash_v[ii] != rab_hash_v[jj])
                    rab_hash_v[++jj] = rab_hash_v[ii];
            }
            printf("Found %lli duplicated items ...  \n", n_kmer-(jj + 1) );

            auto data_iterator = boomphf::range(static_cast<const u_int64_t*>(rab_hash_v), static_cast<const u_int64_t*>(rab_hash_v+n_kmer));

            bphf = new boomphf::mphf<u_int64_t, hasher_t>(n_kmer, data_iterator, nthreads, gammaFactor);

            printf("The minimal perfect hashing function is generated ... ");
        }

            
    private:
        u_int64_t n_kmer; 
        u_int64_t k_kmer;
        u_int64_t* rab_hash_v;
        kmer_t **kmer_data;
        double gammaFactor; 
        boophf_t * bphf;
        u_int64_t base;
        u_int64_t Prime;

        unsigned int nthreads;
};
