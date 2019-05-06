#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include "util.h"
#include "amino_acids.h"
#include "forward_backward.h"

void initialise_F_numerical_recipes(int num_sites, int closest_n, int total_num_refs,
                                    Decimal codon_to_codon_HLA[NUM_CODONS][num_sites],
                                    char ref_codons[total_num_refs][num_sites],
                                    char ref_triplets[total_num_refs][num_sites],
                                    Decimal R[num_sites],
                                    Decimal F[closest_n][num_sites],
                                    int F_rnorm[num_sites],
                                    int start_site, int end_site,
                                    int closest_refs[closest_n],
                                    int num_bases,
                                    char cons_query_nucls[num_bases])
{
  int r, s = start_site;
  Decimal sum, ref_sum, F_sum; 
  // The code below is for filling the entire F matrix.

  #ifdef USE_TRIPLET
    (void)ref_codons;
    for(r = 0; r < closest_n; r++)
      F[r][0] = codon_to_codon_HLA[(int)triplet_to_codon(ref_triplets[closest_refs[r]][0],&cons_query_nucls[0])][0];
  #else
    (void)ref_triplets;
    (void)cons_query_nucls;
    (void)num_bases;
    for(r = 0; r < closest_n; r++)
    F[r][0] = codon_to_codon_HLA[(int)ref_codons[closest_refs[r]][0]][0];
  #endif

  F_rnorm[0] = 0;

  for(s = start_site+1; s <= end_site; s++) {
    F_sum = 0;
    ref_sum = 0;
    for(r = 0; r < closest_n; r++)
      ref_sum += F[r][s-1];
    ref_sum = ref_sum * R[s] / total_num_refs;
    for(r = 0; r < closest_n; r++)
    {
      sum = ref_sum;
      sum += (1 - R[s]) * F[r][s - 1];
      #ifdef USE_TRIPLET
        sum *= codon_to_codon_HLA[(int) triplet_to_codon(ref_triplets[closest_refs[r]][s],
                                                         &cons_query_nucls[s*3])][s];
        if(isnan(sum))
          printf("site: %i, codon: %i, reference: %i, "DECPRINT"\n", s,
                 (int)triplet_to_codon(ref_triplets[closest_refs[r]][s], &cons_query_nucls[s*3]), r,
                 codon_to_codon_HLA[(int)triplet_to_codon(ref_triplets[closest_refs[r]][s],
                                                          &cons_query_nucls[s*3])][s]);
      #else
        sum *= codon_to_codon_HLA[(int) ref_codons[closest_refs[r]][s]][s];
        if(isnan(sum))
          printf("init hello "DECPRINT"\n", codon_to_codon_HLA[(int)ref_codons[closest_refs[r]][s]][s]);
      #endif  
      F[r][s] = sum;
      F_sum += sum;
    } 
    
    F_rnorm[s] = F_rnorm[s-1];
    
    if(F_sum < BIGI) {
      ++F_rnorm[s];
      for(r = 0; r < closest_n; r++)
        F[r][s] *= BIG;
    }
  }
}

// Update_F_numerical_recipes fills a preinitialised F from start_site
// to end_site.
void update_F_numerical_recipes(int num_sites, int closest_n, int total_num_refs,
                                Decimal codon_to_codon_HLA[NUM_CODONS][num_sites],
                                char ref_codons[total_num_refs][num_sites],
                                char ref_triplets[total_num_refs][num_sites],
                                Decimal R[num_sites],
                                Decimal F_in[closest_n][num_sites],
                                int F_rnorm_in[num_sites],
                                int F_rnorm_out[num_sites],
                                Decimal F_out[closest_n][num_sites],
                                int start_site, int end_site,
                                int closest_refs[closest_n],
                                int num_bases,
                                char cons_query_nucls[num_bases])
{
  int r, s = start_site, prev_F_rnorm;
  Decimal sum, ref_sum, F_sum;

  // The first column of F is the corresponding entries from codon_to_codon_HLA.
  Decimal (*F)[num_sites] = (s == 0 ? F_out : F_in);

  #ifdef USE_TRIPLET
    (void)ref_codons;
    if(s == 0)
    {
      for(r = 0; r < closest_n; r++)
        F_out[r][0] = codon_to_codon_HLA[(int)triplet_to_codon(ref_triplets[closest_refs[r]][0],
                                                               &cons_query_nucls[0])][0];

      F_rnorm_out[0] = 0;
      F_rnorm_in[0] = 0;
      s++;
    }
  #else
    (void)ref_triplets;
    (void)cons_query_nucls;
    (void)num_bases;
    if(s == 0)
    {
      for(r = 0; r < closest_n; r++)
        F_out[r][0] = codon_to_codon_HLA[(int)ref_codons[closest_refs[r]][0]][0];

      F_rnorm_out[0] = 0;
      F_rnorm_in[0] = 0;
      s++;
    }
  #endif

  prev_F_rnorm = F_rnorm_in[s - 1];

  // Now fill the remainder of F.
  // Note: F_norm records the number of times we multiply through by BIG
  // to get the column sum bigger than BIGI.
  for(; s <= end_site; s++)
  {
    F_sum = 0;
    ref_sum = 0;
    for(r = 0; r < closest_n; r++)
      ref_sum += F[r][s-1];
    ref_sum = ref_sum * R[s] / total_num_refs;
    for(r = 0; r < closest_n; r++)
    {
      sum = ref_sum;
      sum += (1 - R[s]) * F[r][s - 1];
      #ifdef USE_TRIPLET
        sum *= codon_to_codon_HLA[(int)triplet_to_codon(ref_triplets[closest_refs[r]][s],
                                                        &cons_query_nucls[s*3])][s];
        if(isnan(sum))
          printf("site: %i, codon: %i, reference: %i, "DECPRINT"\n", s,
                 (int)triplet_to_codon(ref_triplets[closest_refs[r]][s], &cons_query_nucls[s*3]), r,
                 codon_to_codon_HLA[(int)triplet_to_codon(ref_triplets[closest_refs[r]][s],
                                                          &cons_query_nucls[s*3])][s]);
      #else
        sum *= codon_to_codon_HLA[(int)ref_codons[closest_refs[r]][s]][s];
        if(isnan(sum))
          printf("hello "DECPRINT"\n", codon_to_codon_HLA[(int)ref_codons[closest_refs[r]][s]][s]);
      #endif

      F_out[r][s] = sum;
      F_sum += sum;
    } 
    
    F_rnorm_out[s] = prev_F_rnorm;
    
    if(F_sum < BIGI) {
      ++F_rnorm_out[s];
      for(r = 0; r < closest_n; r++)
        F_out[r][s] *= BIG;
    }
    F = F_out;
    prev_F_rnorm = F_rnorm_out[s];
  }
}

void update_B_numerical_recipes(int num_sites, int closest_n, int total_num_refs,
                                Decimal codon_to_codon_HLA[NUM_CODONS][num_sites],
                                char ref_codons[total_num_refs][num_sites],
                                char ref_triplets[total_num_refs][num_sites],
                                Decimal R[num_sites],
                                Decimal B[closest_n][num_sites],
                                int B_rnorm[num_sites],
                                int update,
                                int closest_refs[closest_n],
                                int num_bases,
                                char cons_query_nucls[num_bases])
{
  #ifdef USE_TRIPLET
    (void) ref_codons;
  #else
    (void) ref_triplets;
    (void) cons_query_nucls;
    (void) num_bases;
  #endif

  int r, s, backwards = update;
  Decimal sum, ref_sum, B_sum;
  
  // The final column of B is 1s.
  if(backwards == num_sites - 1)
  {
    for(r = 0; r < closest_n; r++)
      B[r][num_sites - 1] = 1;
    B_rnorm[num_sites - 1] = 0;
    backwards--;
  }
  
  // Now fill the remainder of B.
  // Note: B_norm records the number of times we multiply through by BIG
  // to get the column sum bigger than BIGI.
  for(s = backwards; s >= 0; s--)
  {
    B_sum = 0;
    ref_sum = 0;
    for(r = 0; r < closest_n; r++) {
      #ifdef USE_TRIPLET
        ref_sum += B[r][s + 1] * codon_to_codon_HLA[(int) triplet_to_codon(ref_triplets[closest_refs[r]][s + 1],
                                                                           &cons_query_nucls[(s+1)*3])][s + 1];
      #else
        ref_sum += B[r][s + 1] * codon_to_codon_HLA[(int) ref_codons[closest_refs[r]][s + 1]][s + 1];
      #endif
    }
    ref_sum = ref_sum * R[s + 1] / total_num_refs;
    for(r = 0; r < closest_n; r++)
    {
      sum = ref_sum;
      #ifdef USE_TRIPLET
        sum += (1 - R[s + 1]) * B[r][s + 1] *
               codon_to_codon_HLA[(int) triplet_to_codon(ref_triplets[closest_refs[r]][s + 1],
                                                         &cons_query_nucls[(s+1)*3])][s + 1];
      #else
        sum += (1 - R[s + 1]) * B[r][s + 1] *
               codon_to_codon_HLA[(int)ref_codons[closest_refs[r]][s + 1]][s + 1];
      #endif
      B[r][s] = sum;
      B_sum += sum;
    } 
    
    B_rnorm[s] = B_rnorm[s + 1];
    
    if(B_sum < BIGI) {
       ++B_rnorm[s];
      for(r = 0; r < closest_n; r++)
        B[r][s] *= BIG;
    }
  }
}
