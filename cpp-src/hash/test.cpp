#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "generate_hash.h"
typedef short kmer_t;

#include <iostream>
using namespace std;


int main(){

    int n = 1001;

    // int s;

    // s = getPrime(n);

    // printf("the prime is %d\n", s);

    // int k;
    // for (int i; i < 10; i++){
    //     k = randomNumber(0,10);
    //     printf("the random number is %d\n", k);
    // }

    int k = 100;
    kmer_t *data[n];
    for (int i = 0; i < n; i++){
        data[i] = new kmer_t[k];
    }

    for (int i = 0; i < n; i++){
        for (int j = 0; j < k; j++){
            data[i][j] = randomNumber(0, 4);
        }
    }
  
    hash_generator hash_func(data, n, k);
    
    printf("Try query some value\n ");
    for (int i = 0; i < 10; i++)
    {
        u_int64_t h = hash_func.get_hash_value(data[i]);
        
        printf("hash value of ");
        for (int j = 0;  j < k; j++)
        {
            printf("%d", data[i][j]);
        }
        printf(": ");
        printf("hash value %d\n ", h);
    }

    return 1;
}

