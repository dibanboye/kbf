#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>


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
}

// Task 3: Implement the randomNumber function.
// Generate a random number between [start, end]. 
// Note that start and end is included.
template <typename T>
T randomNumber(T start, T end){
    double myRand = rand()/(1.0 + RAND_MAX); 
    printf("my Rand %f", myRand);
    printf("RAND_MAX %d", RAND_MAX);
    T range = end - start + 1;
    T random_num = (myRand * range) + start;
    return random_num;
}

//Task4: rabinHash



