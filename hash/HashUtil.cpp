#ifndef HASH_UTIL
#define HASH_UTIL

#include <ctime>
#include "../cpp-src/KBFUtil.hpp"
template <typename T>
T getPrime(T);
template <typename T>
bool isPrime(T);
template <typename T>
T randomNumber(const T, const T);
short baseNum(char);
unsigned access_kmer( kmer_t mer,  unsigned i );
void print_kmer( kmer_t mer, unsigned k,ostream& os);
/**Utility functions for generating our hash function*/

/**
 * Return the next prime greater than n 
 */
template <typename T>
T getPrime(T n) {

    T inc = 1;

    while (1){
        if (isPrime(n+inc)) { 
            return n+inc;
        }
        inc += 1;
    }

    return 0;
}

/**
 * Check whether n is a prime. n must be greater than 1.
 */
template <typename T>
bool isPrime(T n) {

    T root, i;

    if (n%2 == 0 || n%3 == 0)
        return false;

    root = sqrt(n);

    for (i=5; i<=root; i+=6)
    {
        if (n%i == 0)
           return false;
    }

    for (i=7; i<=root; i+=6)
    {
        if (n%i == 0)
           return false;
    }

    return true;
}

/**
 * Generate a random number between [start, end]. 
 * Note that start and end is included.
 */
template <typename T>
T randomNumber(const T start, const T end){
    double myRand = std::rand()/(1.0 + RAND_MAX); 
    T range = end - start + 1;
    T random_num = (myRand * range) + start;
    return random_num;
}

/**
 * get number in {0, 1, 2, 3} corresponding to base letter
 */
short baseNum(char base) {

   switch (base) {
      case 'A':
         return 0;
      case 'C':
         return 1;
      case 'G':
         return 2;
      case 'T':
         return 3;
   }

   return -1;
}

unsigned access_kmer( kmer_t mer, unsigned k, unsigned i ) {
  mer = mer >> 2*(k - i - 1);
  kmer_t mask = static_cast<kmer_t>(3);
  mer = mask & mer;
  return static_cast<unsigned>(mer); 
}

void print_kmer( kmer_t mer, unsigned k,ostream& os) {
  for (unsigned i = 0; i < k; ++i) {
    unsigned kk = access_kmer( mer, k, i );
    switch( kk ) {
    case 0:
      os << 'A';
      break;
    case 1:
      os << 'C';
      break;
    case 2:
      os << 'G';
      break;
    case 3:
      os << 'T';
      break;
    }
  }
}

#endif
