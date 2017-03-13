#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

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

    int c = 1;
    while (1){
        if (isPrime(n+c)) { 
            return n+c;
        }
        c += 1;
    }
}

