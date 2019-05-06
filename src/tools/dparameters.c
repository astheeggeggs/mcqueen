// Request decent POSIX version.
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h> // strcasecmp.
#include <strings.h> // strcasecmp.

#include "load_data.h"
#include "closest_seqs.h"
#include "constants.h"
#include "util.h"

FILE *hamming_distance_file = NULL;
const char *ref_path = NULL, *query_path = NULL;

int main(int argc, char **argv)
{
  if(argc - optind != 2) printf("Not passed the correct number of arguments.\n");
  
  ref_path = argv[optind];
  query_path = argv[optind+1];

  char **ref_seqs, **query_seqs;
  int refs_cap, query_cap;
  int q;

  printf("Loading sequences...\n");
  int num_refs = load_seqs(ref_path, &ref_seqs, &refs_cap);
  int num_query = load_seqs(query_path, &query_seqs, &query_cap);
   
  int *hamming = my_calloc(num_refs, sizeof(int), __FILE__, __LINE__);
  hamming_distance_file = fopen("hamming_distance.txt", "w");
  size_t num_bases = strlen(ref_seqs[0]);
  
  fprintf(hamming_distance_file, "sum_hamming:\n");
  Decimal sum_hamming = 0;

  for(q = 0; q < num_query; q++) {
    printf("%i\n", q);
    int closest_refs[1000];
    hamming_distance_considering_unknown(ref_seqs, query_seqs[q],
                                         num_refs, num_bases, hamming);
    closest_sequences(num_refs, hamming, 1000, closest_refs);
    
    int i;
    Decimal current_sum_hamming = 0;

    for(i=0; i<1000; i++){
      current_sum_hamming += hamming[closest_refs[i]]/2;
    }
    current_sum_hamming = current_sum_hamming / 1000;
    sum_hamming += current_sum_hamming;
  }
  
  sum_hamming = sum_hamming / num_query;
  fprintf(hamming_distance_file, DECPRINT"\n", sum_hamming); 
  
  fclose(hamming_distance_file);
  free(hamming);

  return EXIT_SUCCESS;
}