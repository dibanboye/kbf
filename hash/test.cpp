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

    printf("the prime is %d", s);

    return 1;


}

