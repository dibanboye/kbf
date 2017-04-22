#include <random>
#include <math.h>
#include <set>
#include <unordered_set>

/**
 * Functions to help test the FDBG data structure
 */

/**
 * Generate count number of random integers between 0 and max
 * Numbers are not necessarily unique
 */
void randomNumbers(const unsigned& count, const unsigned& max,
   vector<unsigned>& numbers) {

   numbers.clear();

    // get random number generator, rand
   std::random_device rd;
   std::mt19937 rand (rd());

   // we sample from a uniform distribution
   std::uniform_int_distribution<u_int64_t> unif_dist (0, max);

   for (int i = 0; i < count; i++) {

      numbers.push_back(unif_dist(rand));

   }

}

/**
 * Add random edges to the graph
 * count is not the number of edges added, it is the number of kmers inspected in order
 * find edges to add. So there is a maximum of 4*count edges added, and a min of 0.
 * Returns the number of edges that it added
 */
unsigned addRandomEdges(const unsigned& count, FDBG& Graph, unordered_set<kmer_t>& kmers) {

   // Random numbers determining which kmers
   vector<unsigned> randoms;
   randomNumbers(count, Graph.k, randoms);

   // Iterate through kmers skipping by a random number each time   
   unordered_set<kmer_t>::iterator i;
   i = kmers.begin(); 

   // Where we store random kmers
   set<kmer_t> random_kmers;

   for (int m = 0; m < count; m++) {

      //BOOST_LOG_TRIVIAL(debug) << "Skipping " << randoms[m];
      // Skip random number of kmers
      for (int j = 0; (j < randoms[m]) && (i != kmers.end()); j++) {     
         ++i;
         // loop around when we get to the end
         if (i == kmers.end()) {
            i = kmers.begin();
         }
      }

      // Now we have a random kmer
      //BOOST_LOG_TRIVIAL(debug) << "Random kmer " << get_kmer_str(*i, Graph.k);
      random_kmers.insert(*i);
   }

   // Attempt to add edges with these random kmers
   set<kmer_t>::iterator random_kmer_it = random_kmers.begin();
   vector<kmer_t> not_neighbors_in;
   vector<kmer_t> not_neighbors_out;

   unsigned added_edges_count = 0;
   bool added;

   for (random_kmer_it = random_kmers.begin(); random_kmer_it != random_kmers.end();
      ++random_kmer_it) {

      // Look through IN and OUT to see if there exists some possible neighbors
      Graph.getNonExistentEdges(*random_kmer_it, not_neighbors_in, not_neighbors_out);
      //BOOST_LOG_TRIVIAL(debug) << "Looking at the non-existent edges of "
      //   << get_kmer_str(*random_kmer_it, Graph.k);

      // Are any of these actually in our graph
      for (int m=0; m < not_neighbors_in.size(); m++) {
         //BOOST_LOG_TRIVIAL(debug) << "Check if node " << get_kmer_str(not_neighbors_in[m],Graph.k)
         //    << " exists";
         if (kmers.find(not_neighbors_in[m]) != kmers.end()) {
            //BOOST_LOG_TRIVIAL(debug) << "Node " << get_kmer_str(not_neighbors_in[m], Graph.k)
            //   << " exists in our graph";
            added = Graph.dynamicAddEdge(not_neighbors_in[m], *random_kmer_it);
            if (added) {
               added_edges_count++;
            }
         }
      }

      for (int m=0; m < not_neighbors_out.size(); m++) {
         //BOOST_LOG_TRIVIAL(debug) << "Check if node " << get_kmer_str(not_neighbors_out[m],Graph.k)
         //    << " exists";
         if (kmers.find(not_neighbors_out[m]) != kmers.end()) {
            //BOOST_LOG_TRIVIAL(debug) << "Node " << get_kmer_str(not_neighbors_out[m], Graph.k)
             //  << " exists in our graph";

            added = Graph.dynamicAddEdge(*random_kmer_it, not_neighbors_out[m]);
            if (added) {
               added_edges_count++;
            }
         }
      }

   }

   return added_edges_count;
}


/**
 * Remove random edges from the graph
 * count is not the number of edges removed
 * Adds all removed edges to removed
 */
void removeRandomEdges(const unsigned& count, FDBG& Graph,
   unordered_set<kmer_t>& edgemers, unordered_set<kmer_t>& removed) {

   // Random numbers determining which kmers
   vector<unsigned> randoms;
   randomNumbers(count, Graph.k, randoms);

   // Iterate through kmers skipping by a random number each time   
   unordered_set<kmer_t>::iterator i;
   i = edgemers.begin(); 

   // How many we have already removed
   unsigned removed_count = 0;

   kmer_t prefix;
   kmer_t suffix;
   unsigned m = 0;
   while ((removed_count < count) && (removed.size() != edgemers.size())) {

      // Skip random number of edgemers
      for (int j = 0; (j < randoms[m]) && (i != edgemers.end()); j++) {     
         ++i;
         // loop around when we get to the end
         if (i == edgemers.end()) {
            i = edgemers.begin();
         }
      }

      // Now we have a random edgemer
      //BOOST_LOG_TRIVIAL(debug) << "Random edgemer "
      //   << get_kmer_str(*i, Graph.k + 1) << " generated.";

      if (removed.find(*i) == removed.end()) {
         // We have not already removed this edgemer
         Graph.split_edge(*i, prefix, suffix);
         Graph.dynamicRemoveEdge(prefix, suffix);
         removed.insert(*i);
         removed_count++;
      }

      // Use the next random number
      m++;
      if (m >= randoms.size()) {
         // Loop back to the first if we run out
         m = 0;
      }
   }

}

