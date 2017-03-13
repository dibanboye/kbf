#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "generatehash.h"


int main(){

    int n = 112012;

    int s;

    s = getPrime(n);

    printf("the prime is %d\n", s);

    int k;
    for (int i; i < 10; i++){
        k = randomNumber(0,10);
        printf("the random number is %d\n", k);
    }



    return 1;
}

