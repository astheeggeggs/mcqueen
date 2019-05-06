#include <stdlib.h>
#include <stdio.h>

#include "amino_acids.h"
#include "constants.h"
#include "util.h"
#include "write_json.h"
#include "assert.h"

FILE *simulated_file = NULL;
FILE *hla_file = NULL;
FILE *mosaic_json = NULL;

void simulate_codon_to_codon_HLA(int num_sites, int num_refs,
                                 int num_sims, int num_hla_types,
                                 Decimal (*codon_to_codon_HLA)[num_sites][NUM_CODONS],
                                 Decimal (*rate_out_of_codons)[num_sites],
                                 int start_site, int end_site,
                                 Decimal mu,
                                 Decimal (*HLA_select_esc)[num_sites],
                                 Decimal HLA_select_rev[num_sites],
                                 char (*hla_types)[num_hla_types],
                                 const char consensus_codons[num_sites],
                                 const Decimal omega[num_sites],
                                 const int NS_TS_sum_from_consensus[num_sites],
                                 const int NS_TV_sum_from_consensus[num_sites],
                                 char (*const type_of_change)[NUM_CODONS][NUM_CODONS],
                                 int (*codon_copied_from)[num_sites],
                                 Decimal (*A)[NUM_CODONS][num_sites])
{
  int sim, s, h, from, to;
  Decimal gamma, HLA_select_esc_prod, HLA_select_rev_prod;

  // To simulate, it is the sequence copying from which is fixed, should consider all
  // possibilities for the codon generated.
  for(sim = 0; sim < num_sims; sim++)
  {
    for(s = start_site; s <= end_site; s++)
    {
      int consensus_codon = consensus_codons[s];
      Decimal rate_out_of_codons_prime[NUM_CODONS];

      HLA_select_esc_prod = 1;
      HLA_select_rev_prod = HLA_select_rev[s];

      for(h = 0; h < num_hla_types; h++) {
        if(hla_types[sim][h] == '1') {
          HLA_select_esc_prod *= HLA_select_esc[h][s];
        }
      }
      // The reversion portion.
      // Set rate_out_of_codons_prime.
      
      // What codon are we copying?
      from = codon_copied_from[sim][s];

      if(from!=consensus_codon) {

        rate_out_of_codons_prime[from]
          = rate_out_of_codons[from][s] +
            (HLA_select_rev_prod - 1) * 
            (omega[s] * (kappa * NS_TS[from][consensus_codon] + 
                         NS_TV[from][consensus_codon]) + 
             kappa * S_TS[from][consensus_codon] + S_TV[from][consensus_codon]);

        // To get the new rate out of codons for reversion,
        // take the standard rate out of codons,
        // add on the extra boost through reversion, but multiply by boolean -
        
        // is it non consensus to consensus?
        // If it is, make sure that it is a non synonymous change using NS_TS and NS_TV.  
        A[sim][from][s] = 1 - (Decimal) num_refs/
                              (num_refs + rate_out_of_codons_prime[from] * 4 * mu);
        
        gamma = (Decimal) phi / rate_out_of_codons_prime[from];

        // Now determine the probability matrix conditional on HLA for the
        // reversion portion.
        Decimal sum_check = 0;
        for(to = 0; to < NUM_CODONS; to++)
        {
          codon_to_codon_HLA[sim][s][to]
            = A[sim][from][s] * ((type_of_change[s][from][to] == NONCON_TO_CON) *
                                (gamma * HLA_select_rev_prod * (omega[s] * 
                                         (kappa * NS_TS[from][to] + NS_TV[from][to]) +
                                          kappa * S_TS[from][to] + S_TV[from][to]) +
                                 (1-phi) * prior_C2[to]) +
                                (type_of_change[s][from][to] == NONCON_TO_NONCON) *
                                (gamma * (omega[s] * (kappa * NS_TS[from][to] + NS_TV[from][to]) +
                                          kappa * S_TS[from][to] + S_TV[from][to]) +
                                          (1-phi) * prior_C2[to])) +
              (from == to) * (type_of_change[s][from][to] == NONCON_TO_NONCON) *
              (1 - A[sim][from][s]);
          sum_check += codon_to_codon_HLA[sim][s][to];
        }
        assert(sum_check > 0.999 && sum_check < 1.0001);

      } else {
        // The codon copied from is consensus at this position.
        // The escape portion.
        // Set rate out of codons_prime.
        rate_out_of_codons_prime[consensus_codon]
          = rate_out_of_codons[consensus_codon][s] +
            (HLA_select_esc_prod - 1) * omega[s] *
            (kappa * NS_TS_sum_from_consensus[s] + NS_TV_sum_from_consensus[s]);

        A[sim][consensus_codon][s]
          = 1 - (Decimal) num_refs /
                (num_refs + rate_out_of_codons_prime[consensus_codon] * 4 * mu);
        
        gamma = (Decimal) phi / rate_out_of_codons_prime[consensus_codon];
        
        Decimal sum_check = 0;
        for(to = 0; to < NUM_CODONS; to++)
        {
          codon_to_codon_HLA[sim][s][to]
            = A[sim][consensus_codon][s] *
              (gamma * (HLA_select_esc_prod * omega[s] *
                        (kappa * NS_TS[consensus_codon][to] +
                        NS_TV[consensus_codon][to]) +
               kappa * S_TS[consensus_codon][to] +
               S_TV[consensus_codon][to]) +
              (1-phi) * prior_C2[to]);

          if(to == consensus_codon)
          {
            codon_to_codon_HLA[sim][s][to]
              = (1 - A[sim][consensus_codon][s]) * (1-(1-phi)*prior_C2[to]) + 
                (1-phi)*prior_C2[to];
          }
          sum_check += codon_to_codon_HLA[sim][s][to];
        }
        assert(sum_check > 0.999 && sum_check < 1.0001);
      }
    }
  }
}

static int bfind_pos(Decimal to_sort[NUM_CODONS])
{
  Decimal where = rand_lim(1);
  int i, j = 0;
  
  for(i = 0; i < NUM_CODONS; i++) {
    if(where > to_sort[i]) {
      j += 1;
    } else {
      break;
    }
  }
  return j;
}

void simulate_sequences(int num_sites, int num_refs, int num_sims,
                        int num_hla_types, Decimal HLA_prevalences[num_hla_types],
                        int ploidy, int n_genes, int n_HLA[n_genes],
                        char ref_codons[num_refs][num_sites], Decimal mu,
                        Decimal (*HLA_select_esc)[num_sites],
                        Decimal HLA_select_rev[num_sites],
                        const char consensus_codons[num_sites],
                        Decimal omega[num_sites],
                        Decimal R[num_sites], 
                        Decimal rate_out_of_codons[num_sites][NUM_CODONS],
                        const int NS_TS_sum_from_consensus[num_sites],
                        const int NS_TV_sum_from_consensus[num_sites],
                        char (*const type_of_change)[NUM_CODONS][NUM_CODONS],
                        Decimal (*sim_codon_to_codon_HLA)[num_sites][NUM_CODONS],
                        Decimal proportion_no_hla,
                        int (*codon_copied_from)[num_sites],
                        int (*simulated_sequences)[num_sites],
                        Decimal (*A)[NUM_CODONS][num_sites],
                        int (*copy_from)[num_sites])
{
  int s, h, p, i, k, sims, codon;

  assert(num_hla_types > 0);
  assert(num_sims > 0);
  assert(ploidy > 0);
  assert(n_genes > 0);

  // First determine the hla_types of each of the sequences to be simulated.
  char (*hla_types)[num_hla_types] = my_malloc(num_sims * sizeof(char[num_hla_types]),
                                               __FILE__,__LINE__);
  
  int cumulative_n_HLA;
  int HLA[ploidy * n_genes];

  // Initialise HLA types.
  for(sims = 0; sims < num_sims; sims++) {          
    for(h = 0; h < num_hla_types; h++) {
      hla_types[sims][h] = '0';
    }
  }

  for(sims = 0; sims < num_sims; sims++)
  {  
    cumulative_n_HLA = 0;
    for(i = 0, k = 0; i < n_genes; i++) {
      for(p = 0; p < ploidy; p++, k++) {
        HLA[k] = cumulative_n_HLA +
                 discrete_sampling_dist(n_HLA[i], &HLA_prevalences[cumulative_n_HLA]);
        hla_types[sims][HLA[k]] = '1';
      }
      cumulative_n_HLA += n_HLA[i];
    }
  }

  // Now randomly allocate a proportion to have no HLA information.
  for(sims = 0; sims < num_sims; sims++) {          
    if(rand_lim(1) < proportion_no_hla) {
      for(h = 0; h < num_hla_types; h++) {
        hla_types[sims][h] = '0';
      }
    }
  }

  // Now save the hla types to a .csv file.
  fprintf(hla_file, "\"\",");
  for(h = 0; h < num_hla_types-1; h++) {
    fprintf(hla_file, "\"%i\",", h+1);
  }
    fprintf(hla_file, "\"%i\"\n", num_hla_types);

  for(sims = 0; sims < num_sims; sims ++) {
    fprintf(hla_file,"\"simulated_seq_%i\"",sims+1);
    for(h = 0; h < num_hla_types; h++) {
      fprintf(hla_file, ",%c", hla_types[sims][h]);
    }
    fprintf(hla_file, "\n");
  }

  Mosaic *true_copy = my_malloc(num_sims * sizeof(Mosaic),__FILE__,__LINE__);
  int *true_seqs = my_malloc(num_sims * num_sites * sizeof(int),__FILE__,__LINE__);
  for(sims = 0; sims < num_sims; sims++) {
    true_copy[sims].which_seqs = &true_seqs[sims*num_sites];
  }

  for(sims = 0; sims < num_sims; sims++) {
    int total_copied_from = 0;

    copy_from[sims][0] = rand_lim(num_refs);
    true_copy[sims].which_seqs[0] = copy_from[sims][0];
    codon_copied_from[sims][0] = (int)ref_codons[copy_from[sims][0]][0];

    // Store copy_from to check how well Viterbi works - it's not needed.
    for(s = 1; s < num_sites; s++) {
      if(rand_lim(1) < R[s]) {
        total_copied_from += 1; 
        copy_from[sims][s] = (int)rand_lim(num_refs);
        true_copy[sims].which_seqs[total_copied_from] = copy_from[sims][s];
        printf("\n%i\n",true_copy[sims].which_seqs[total_copied_from]);
      } else {
        copy_from[sims][s] = copy_from[sims][s-1];
      }
      codon_copied_from[sims][s] = (int)ref_codons[copy_from[sims][s]][s];
    }
    true_copy[sims].num_seqs = total_copied_from+1;
  }

  write_closest_n_json(mosaic_json, true_copy, num_sims);

  free(true_copy);
  free(true_seqs);

  // Create codon_to_codon_HLA for each of the simulated sequences, 
  // dependent on hla type.

  simulate_codon_to_codon_HLA(num_sites, num_refs, num_sims, num_hla_types,
                              sim_codon_to_codon_HLA, rate_out_of_codons,
                              0, num_sites - 1, mu,
                              HLA_select_esc, HLA_select_rev, hla_types,
                              consensus_codons, omega,
                              NS_TS_sum_from_consensus, NS_TV_sum_from_consensus,
                              type_of_change, codon_copied_from, A);

  free(hla_types);

  // codon_to_codon_HLA created for each of the simulated sequences.
  
  // Now need to choose which codon is present at each site to generate the
  // collection of sequences.
  
  // Sum and create a load of uniform random variables.
  for(sims = 0; sims < num_sims; sims++) {
    for(s = 0; s < num_sites; s++) {
      for(codon = 1; codon < NUM_CODONS; codon++) {
        sim_codon_to_codon_HLA[sims][s][codon] += 
          sim_codon_to_codon_HLA[sims][s][codon-1];
      }
    }
  }

  // Now choose codons based on these distributions.
  for(sims = 0; sims < num_sims; sims++) {
    for(s = 0; s < num_sites; s++) {
      simulated_sequences[sims][s] = bfind_pos(sim_codon_to_codon_HLA[sims][s]); 
    }
  }
  
}

void save_simulated_fasta(int num_sims, int num_sites, int (*simulated_sequences)[num_sites])
{
  int s, sim;
  for(sim = 0; sim < num_sims; sim++) {
    fprintf(simulated_file,">simulated_seq_%i\n",sim+1);
    for(s = 0; s < num_sites; s++) {
      fprintf(simulated_file, "%s", code_to_char(simulated_sequences[sim][s]));
    }
    fprintf(simulated_file, "\n");
  }
}
