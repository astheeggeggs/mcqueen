#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#include "util.h"
#include "amino_acids.h"
#include "closest_seqs.h"

// Define a collection of structures.
typedef struct
{
  int seq_idx;
  int hamming_dist;
} index_and_dist;

void hamming_distance_considering_unknown(char **ref_seqs, char *query_seqs,
                                          int num_refs, int num_bases,
                                          int tmp_hamming[num_refs])
{
  int i, r, len_count_ref, len_count_query;
  for(r = 0; r < num_refs; r++)
  {
    tmp_hamming[r] = 0;
    len_count_ref = 0;
    len_count_query = 0;
    for(i = 0; i < num_bases; i++) {
      // Ignore the positions in the query sequence that are unknown.
      if((int)base_to_code(query_seqs[i]) != 4) {
        len_count_query += 1; 
        // Ignore the position in ref sequence that are unknown.
        if((int)base_to_code(ref_seqs[r][i] != 4)) {
          len_count_ref += 1;
          tmp_hamming[r] += (base_to_code(query_seqs[i]) != base_to_code(ref_seqs[r][i]));
        } else {
          printf("unknown nucleotide in the reference sequence\n");
        }
      }
    }
    
    // Gap penalty here is 1/2.
    // Set the gap penalty using the ratios of addition to len_count and
    // multiplication by the factor below.

    tmp_hamming[r]*=2;
    tmp_hamming[r]+=len_count_query-len_count_ref;
  }
}

static int compare_hamming_dsts(const void *p1, const void *p2)
{
  const index_and_dist *w1 = (const index_and_dist *)p1;
  const index_and_dist *w2 = (const index_and_dist *)p2;
  return (w1->hamming_dist - w2->hamming_dist);
}

void knuth_perm_idx_and_dist(index_and_dist *arr, int length)
{
  int i, j;
  for(i = 0; i < length; i++) {
    j = rand_lim(i+1);
    SWAP3(arr[i], arr[j]);
  }
}

void closest_sequences(int num_refs, int tmp_hamming[num_refs], int closest_n, 
                       int closest_refs[closest_n])
{
  assert(num_refs > 0);
  assert(closest_n > 0);

  int i;
  index_and_dist dists[num_refs];

  for(i = 0; i < num_refs; i++) {
    dists[i].seq_idx = i;
    dists[i].hamming_dist = tmp_hamming[i];
  }

  // Initial sort (of the first closest_n elements).
  qsort(dists, num_refs, sizeof(index_and_dist), compare_hamming_dsts);

  // Get start end end
  int start = closest_n-1, end = closest_n-1;
  int largest_hamming = dists[closest_n].hamming_dist;

  while(start > 0 && dists[start-1].hamming_dist == largest_hamming) start--;
  while(end+1 < num_refs && dists[end+1].hamming_dist == largest_hamming) end++;

  int num_last = end - start + 1;

  // Knuth perm
  // Next, randomly permute these indices and choose those to include in
  // closest_refs.
  knuth_perm_idx_and_dist(dists+start, num_last);

  for(i = 0; i < closest_n; i++)
    closest_refs[i] = dists[i].seq_idx;
}
