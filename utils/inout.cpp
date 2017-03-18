#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
//#include <boost/math/special_functions/prime.hpp>
#include "../cpp-src/KBFUtil.hpp"
#include "../hash/generate_hash.h"

using namespace std;


//char** ref_mat = get_mat<char>(1000,2)
template<typename T>
T ** get_mat(int row, int col){
   T** ptr = 0;
   ptr = new T* [row];
  for(int id =0; id!=row; ++id){
    ptr[id] = new T[col];
    }
   return ptr;
}

template<typename T>
void clear_mat(T** mat, int row){
  for(int id =0; id!=row; ++id){
     delete [] mat[id];
    }
     delete [] mat;
}

/*
  We may also need to save a ends_ref Mat, in which ends_ref[N,0], ends_ref[N,1] means the starts and ends token of kemer s, 
  where f(s) = N.
  Both ends token should be int value, like: 0,1,2,3,4 repsenting AGCTU.
*/

void gen_ref_mat(unordered_set<kmer_t>& kmers, char ** ref_mat, unsigned k, const generate_hash * hf){
  // ref_mat should be pre_allocated memory.
  //char ** ref_mat = get_mat<char>(kmers.size(), 2);
  for (auto it1 = kmers.begin();	it1 != kmers.end(); ++it1) {
      u_int64_t Ind = hf->get_hash_value( *it1 );
      ref_mat[Ind,0] = get_token_id(  access_kmer( *it1, k, 0 ) );
      ref_mat[Ind,1] = get_token_id(  access_kmer( *it1, k, k-1) );
   }

}

void del_edge(bool** IN, bool** OUT, char** ends_ref, u_int64_t inp, u_int64_t out){

  // delted edge is inp-->out
     u_int64_t out_col = ends_ref[out, 1];
     u_int64_t inp_col = ends_ref[inp, 0];
     IN[out, inp_col] = 0;
     OUT[inp, out_col] =0;
     return; 
}
void add_edge(bool** IN, bool** OUT, char** ends_ref, u_int64_t inp, u_int64_t out){
  // added edge is inp-->out
     u_int64_t out_col = ends_ref[out, 1];
     u_int64_t inp_col = ends_ref[inp, 0];
     IN[out, inp_col] = 1;
     OUT[inp, out_col] =1;
     return; 
}



int main(int argc, char* argv[]) {
    
    return 0;
}
