#ifndef FORMATUTIL_CPP
#define FORMATUTIL_CPP

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <chrono>
#include "KBFUtil.hpp"

using namespace std;

/* Test if a file exists */
bool file_exists(const std::string &name)
{
  if (FILE *file = fopen(name.c_str(), "r"))
  {
    fclose(file);
    return true;
  }
  else
  {
    return false;
  }
}

void write_to_bin(string outfile, unordered_set<kmer_t> &kmers)
{
  BOOST_LOG_TRIVIAL(info) << "Writing binary k-mers to file: " << outfile;
  ofstream binfile(outfile.c_str(), ios::out | ios::binary);

  kmer_t *kdata = new kmer_t[kmers.size()];
  auto it = kmers.begin();
  unsigned i = 0;
  while (it != kmers.end())
  {
    kdata[i] = *it;
    ++it;
    ++i;
  }

  binfile.write((char *)kdata, sizeof(kmer_t) * kmers.size());

  delete[] kdata;

  binfile.close();

  BOOST_LOG_TRIVIAL(info) << "Verifying file " << outfile << "...";
  vector<kmer_t> kmer_vec(kmers.size(), 0);
  ifstream ifile(outfile.c_str(), ios::in | ios::binary);
  ifile.read(reinterpret_cast<char *>(kmer_vec.data()), sizeof(kmer_t) * kmers.size());
  ifile.close();
  auto it2 = kmer_vec.begin();
  it = kmers.begin();
  while (it2 != kmer_vec.end())
  {
    if (*it2 != *it)
    {
      BOOST_LOG_TRIVIAL(fatal) << "File corrupted...";

      exit(1);
    }

    ++it2;
    ++it;
  }

  BOOST_LOG_TRIVIAL(info) << "File has been successfully verified.";
}

void read_from_bin(string filename, unordered_set<kmer_t> &out)
{
  BOOST_LOG_TRIVIAL(info) << "Reading mers from binary file " << filename << " ...";

  ifstream ifile(filename.c_str(), ios::in | ios::binary);
  ifile.seekg(0, ifile.end);
  size_t nbytes = ifile.tellg();
  ifile.seekg(0, ifile.beg);

  BOOST_LOG_TRIVIAL(info) << "File size (Gb): " << nbytes / (1024.0*1024.0*1024.0);
  
  vector<kmer_t> kmer_vec(nbytes / sizeof(kmer_t), 0);

  ifile.read(reinterpret_cast<char *>(kmer_vec.data()), nbytes);
  ifile.close();

  unordered_set<kmer_t> kmers(kmer_vec.begin(), kmer_vec.end());

  out.swap(kmers);

  BOOST_LOG_TRIVIAL(info) << "Input finished, " << kmer_vec.size() << " mers read.";
}

void get_kmers_fasta_or_bin(string fasta_file, unsigned k, unordered_set<kmer_t> &out)
{
  // fasta filename
  string filename = fasta_file;

  // Check if the k-mers from this file have already been stored in binary format
  string kfile = filename.substr(0, filename.find_last_of('.')) + to_string(k) + ".bin";

  if (!file_exists(kfile))
  {
    BOOST_LOG_TRIVIAL(info) << "Binary file for " << k << "-mers does not exist ...";
    // get reads from file
    BOOST_LOG_TRIVIAL(info) << "Parsing fasta ...";
    vector<string> reads = parseFasta(filename);
    out = getKmers(reads, k);

    //save to binary file for future use
    write_to_bin(kfile, out);
  }
  else
  {
    //k-mers are already stored in the binary format!
    read_from_bin(kfile, out);
  }
}

/*
 * This function checks to see if binary information from k, k+1 mers
 * is already present. If so, it simply loads the binary format.
 * If not, it parses the fasta file and writes the k, k+1-mer binary
 * files for future use.
 */

void handle_mers(string fasta_file, unsigned k,
                 unordered_set<kmer_t> &kmers,   //k mers from fasta
                 unordered_set<kmer_t> &edgemers //k + 1 mers from fasta
                 )
{

  get_kmers_fasta_or_bin(fasta_file, k, kmers);
  get_kmers_fasta_or_bin(fasta_file, k + 1, edgemers);
}

#endif
