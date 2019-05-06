#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "constants.h"
#include "util.h"
#include "birth_death_gen.h"
#include "simulate_sequences_gen.h"
#include "amino_acids.h"
#include "assert.h"

gsl_rng *gsl_r;

FILE *simulated_refs_file = NULL;
FILE *simulated_root_file = NULL;
FILE *simulated_queries_file = NULL;
FILE *hla_query_file = NULL;

int create_new_codon(int current_codon, int wildtype_codon, 
                     const int row_sum_type_of_change, 
                     const int type_of_change_matrix[NUM_CODONS][NUM_CODONS])
{
  int c, total_available = 0;
  int new_codon;
  int length_to_codon = row_sum_type_of_change;

  if(type_of_change_matrix[current_codon][wildtype_codon]==1) 
    length_to_codon--;

  int to_codon[length_to_codon];

  for(c = 0; c < NUM_CODONS; c++) {
    if((type_of_change_matrix[current_codon][c] == 1) & (c!=wildtype_codon)) {
      to_codon[total_available] = c;
      total_available++;
    }
  }
  
  assert(total_available==length_to_codon);

  int a = (int) rand_lim(length_to_codon);
  new_codon = to_codon[a];

  return new_codon;
}

void codon_sequence_change(int codon_sequence_length, 
                           Decimal S_portion[NUM_CODONS], Decimal NS_portion[NUM_CODONS],
                           Decimal max_time, Decimal mu,
                           Decimal omega[codon_sequence_length],
                           int wildtype_sequence[codon_sequence_length],
                           int (*codon_sequence_matrix)[codon_sequence_length], 
                           int from, int to,
                           Decimal reversion_selection[codon_sequence_length],
                           Decimal HLA_selection_profile[codon_sequence_length])
{
  int i, j, c;
  Decimal change_time;
  Decimal codon_rates[5];
  int current_codon;
  int new_codon;
  int sub_type;
  int codon_sequence[codon_sequence_length];
  Decimal total_site_codon_rate[NUM_CODONS], rev[NUM_CODONS];

  for(i = 0; i < codon_sequence_length; i++)
    codon_sequence[i] = codon_sequence_matrix[from][i];

  for(i = 0; i < codon_sequence_length; i++) {
    for(c = 0; c < NUM_CODONS; c++) {
      total_site_codon_rate[c] = S_portion[c] + omega[i] * NS_portion[c];
      rev[c] = omega[i] * (kappa * NS_TS[c][wildtype_sequence[i]] + 
                           NS_TV[c][wildtype_sequence[i]]) +
               kappa * S_TS[c][wildtype_sequence[i]] + S_TV[c][wildtype_sequence[i]];
      // Add in the extra reversion contribution to the total rate.
      total_site_codon_rate[c] = (total_site_codon_rate[c] + 
                                 (reversion_selection[i] - 1) * rev[c]) * mu;
    }

    // Rate of movement out of the wildtype site involves the escape parameters.
    total_site_codon_rate[wildtype_sequence[i]] = (S_portion[wildtype_sequence[i]] +
                                                  omega[i] * HLA_selection_profile[i] * 
                                                  NS_portion[wildtype_sequence[i]]) * mu;
    
    Decimal exp_parameter = (Decimal) 1 / total_site_codon_rate[codon_sequence[i]];
    change_time = gsl_ran_exponential(gsl_r, exp_parameter);

    while(change_time < max_time) {
      current_codon = codon_sequence[i];

      // Choose what kind of substitution it is...
      if(current_codon == wildtype_sequence[i]) {
        // Given that the sequence is wildtype at this position, there are 4 choices.
        codon_rates[0] = kappa * beta_S[current_codon]; // S_TS.
        codon_rates[1] = beta_V[current_codon]; // S_TV.
        codon_rates[2] = HLA_selection_profile[i] * omega[i] * kappa * alpha_S[current_codon]; // NS_TS - escape.
        codon_rates[3] = HLA_selection_profile[i] * omega[i] * alpha_V[current_codon]; // NS_TV - escape;

        Decimal dummy_total_rate = (codon_rates[0]+codon_rates[1]+codon_rates[2]+codon_rates[3]) * mu;
        
        assert(((dummy_total_rate - total_site_codon_rate[current_codon]) < 0.001) &&
               ((dummy_total_rate - total_site_codon_rate[current_codon]) > -0.001));
        
        for(j = 0; j < 4; j++)
          codon_rates[j] = mu * codon_rates[j] / total_site_codon_rate[current_codon];

        sub_type = discrete_sampling_dist(4, codon_rates);

        // Now randomly assign what the codon should be.
        if(sub_type == 0) {
          // S_TS
          new_codon = create_new_codon(current_codon, wildtype_sequence[i], 
                                       beta_S[current_codon], S_TS);
        } else if(sub_type == 1) {
          // S_TV
          new_codon = create_new_codon(current_codon, wildtype_sequence[i], 
                                       beta_V[current_codon], S_TV);
        } else if(sub_type == 2) {
          // NS_TS;
          new_codon = create_new_codon(current_codon, wildtype_sequence[i],
                                       alpha_S[current_codon], NS_TS);
        } else {
          assert(sub_type == 3);
          // NS_TV;
          new_codon = create_new_codon(current_codon, wildtype_sequence[i], 
                                       alpha_V[current_codon], NS_TV);
        }

      } else {
        // The codon is not wildtype at this position, so reversion can occur.
        codon_rates[0] = kappa * (beta_S[current_codon] - S_TS[current_codon][wildtype_sequence[i]]);
        codon_rates[1] = beta_V[current_codon] - S_TV[current_codon][wildtype_sequence[i]];
        codon_rates[2] = omega[i] * kappa * (alpha_S[current_codon] - 
                                             NS_TS[current_codon][wildtype_sequence[i]]);
        codon_rates[3] = omega[i] * (alpha_V[current_codon] - NS_TV[current_codon][wildtype_sequence[i]]);
        codon_rates[4] = reversion_selection[i] * rev[current_codon];
        
        Decimal dummy_total_rate = (codon_rates[0]+codon_rates[1]+codon_rates[2]+codon_rates[3]+codon_rates[4]) * mu;

        assert(((dummy_total_rate - total_site_codon_rate[current_codon]) < 0.001) &&
               ((dummy_total_rate - total_site_codon_rate[current_codon]) > -0.001));
        
        for(j = 0; j < 5; j++)
          codon_rates[j] = mu * codon_rates[j] / total_site_codon_rate[current_codon];

        sub_type = discrete_sampling_dist(5, codon_rates);

        // Now randomly assign what the codon should be.
        if(sub_type == 0) {
          // S_TS.
          new_codon = create_new_codon(current_codon, wildtype_sequence[i], 
                                       beta_S[current_codon], S_TS);
        } else if(sub_type == 1) {
          // S_TV.
          new_codon = create_new_codon(current_codon, wildtype_sequence[i],
                                       beta_V[current_codon], S_TV);
        } else if(sub_type == 2) {
          // NS_TS - but not a reversion to wildtype.
          new_codon = create_new_codon(current_codon, wildtype_sequence[i],
                                       alpha_S[current_codon], NS_TS);
        } else if(sub_type == 3) {
          // NS_TV - but not a reversion to wildtype.
          new_codon = create_new_codon(current_codon, wildtype_sequence[i],
                                       alpha_V[current_codon], NS_TV);
        } else {
          assert(sub_type == 4);
          // Reversion.
          new_codon = wildtype_sequence[i];
        }
      }

      codon_sequence[i] = new_codon;
      exp_parameter = (Decimal) 1 / total_site_codon_rate[new_codon];
      change_time += gsl_ran_exponential(gsl_r, exp_parameter);
    
    }
    codon_sequence_matrix[to][i] = codon_sequence[i];
  }

}

void pass_HLA(int ploidy, int n_genes, int starting_node,
              Tree *tree, int number_of_leaves, int total_n_HLA,
              int n_HLA[n_genes], Decimal HLA_prevalences[total_n_HLA])
{
  int parent_nodes[number_of_leaves];
  int daughter_nodes[number_of_leaves];
  
  int parent_nodes_length = 1;
  parent_nodes[0] = starting_node;
  
  int n, p, k, i, d;
  int next_generation_exists = 1;
  int cumulative_n_HLA;
  int HLA[ploidy * n_genes];

  assert(n_genes > 0);
  assert(ploidy > 0);

  while(next_generation_exists == 1) {
    d = 0;
    for(n = 0; n < parent_nodes_length; n++) {
      // Determine the HLA types of one of the daughters.
      cumulative_n_HLA = 0;
      for(i = 0, k = 0; i < n_genes; i++) {
        for(p = 0; p < ploidy; p++, k++) {
          HLA[k] = cumulative_n_HLA + 
                   discrete_sampling_dist(n_HLA[i], &HLA_prevalences[cumulative_n_HLA]);
        }
        cumulative_n_HLA += n_HLA[i];
      }

      if(tree[parent_nodes[n]].daughter_nodes[0] != -1)
      {
        daughter_nodes[d] = tree[parent_nodes[n]].daughter_nodes[0];
        // Either there's a seen coalescence.
        if(tree[parent_nodes[n]].daughter_nodes[1] != -1)
        { 
          d++;
          daughter_nodes[d] = tree[parent_nodes[n]].daughter_nodes[1];
          // Choose the daughter to get the parental HLA.
          int which_daughter = (int) rand_lim(2);
          for(p = 0; p < (ploidy * n_genes); p++) 
          {
            // Fill in the details.
            tree[daughter_nodes[d - which_daughter]].HLAs[p] = tree[parent_nodes[n]].HLAs[p];
            tree[daughter_nodes[d - (1 - which_daughter)]].HLAs[p] = HLA[p];
          }
          d++;
        // Or there's an unseen coalescence.
        } else {
          // In which case, there's always a transmission.
          for(p = 0; p < (ploidy * n_genes); p++)
            tree[daughter_nodes[d]].HLAs[p] = HLA[p];
          d++;
        }
      }
    }
    // Set the new parent nodes as the old daughter nodes.
    for(n = 0; n < d; n++) parent_nodes[n] = daughter_nodes[n];
    parent_nodes_length = d;

    if(d == 0) next_generation_exists = 0; 
  }
}

void pass_codon_sequence_change(int codon_sequence_length, int ploidy,
                                int n_genes, int total_n_HLA,
                                int starting_node, Decimal mu,
                                int (*codon_sequence_matrix)[codon_sequence_length],
                                Tree *tree,
                                int number_of_leaves,
                                Decimal S_portion[NUM_CODONS], Decimal NS_portion[NUM_CODONS],
                                Decimal (*HLA_selection_profiles)[codon_sequence_length],
                                int wildtype_sequence[codon_sequence_length],
                                Decimal omega[codon_sequence_length],
                                Decimal reversion_selection[codon_sequence_length])
{
  // First, create a list of parent nodes - the maximum that this can be is the number of leaves.
  int parent_nodes[number_of_leaves];
  int daughter_nodes[number_of_leaves];

  int parent_nodes_length = 1;
  parent_nodes[0] = starting_node;
  
  int i, p, s, h, d, daughter;
  int next_generation_exists = 1;
  Decimal daughter_node_time;
  Decimal HLA_selection_profile[codon_sequence_length];
  char HLAs_to_check[total_n_HLA];

  for(s = 0; s < codon_sequence_length; s++) HLA_selection_profile[s] = 1;

  while(next_generation_exists == 1) {
    d = 0;
    for(p = 0; p < parent_nodes_length; p++) {
      for(daughter = 0; daughter < 2; daughter++) {
        if(tree[parent_nodes[p]].daughter_nodes[daughter] != -1) {
          daughter_nodes[d] = tree[parent_nodes[p]].daughter_nodes[daughter];
          daughter_node_time = tree[daughter_nodes[d]].node_time - 
                               tree[parent_nodes[p]].node_time;
          // DEV: testing to see what happens when the leaf nodes are
          // set at a fixed length.
          if(daughter_nodes[d] < number_of_leaves) daughter_node_time = 3;

          for(h = 0; h < total_n_HLA; h++) HLAs_to_check[h] = '0';
          for(i = 0; i < (ploidy * n_genes); i++) HLAs_to_check[tree[daughter_nodes[d]].HLAs[i]] = '1';

          for(s = 0; s < codon_sequence_length; s++) {
            HLA_selection_profile[s] = 1;
            for(h = 0; h < total_n_HLA; h++) {
              if(HLAs_to_check[h] == '1') HLA_selection_profile[s] *= HLA_selection_profiles[h][s];
            }
          }

          codon_sequence_change(codon_sequence_length,
                                S_portion, NS_portion,
                                daughter_node_time, mu,
                                omega, wildtype_sequence,
                                codon_sequence_matrix,
                                parent_nodes[p], daughter_nodes[d],
                                reversion_selection,
                                HLA_selection_profile);
          d++;
        }
      }
    }
    
    // Set the new parent nodes as the old daughter nodes.
    for(p = 0; p < d; p++) parent_nodes[p] = daughter_nodes[p];
    parent_nodes_length = d;

    if(d == 0) next_generation_exists = 0; 
  }
}

void save_simulated_ref_and_query_fasta(int num_queries, int num_refs, int num_sims, 
                                        int all_sequences[num_sims], int codon_sequence_length, 
                                        int (*simulated_sequences)[codon_sequence_length],
                                        Tree *tree, int ploidy, int n_genes)
{
  int s, sim, h;
  // Knuth sampling to get queries.
  knuth_sample(num_queries, num_sims, all_sequences);

  // Write the references to a .fasta file.
  for(sim = 0; sim < num_refs; sim++) {
    // Make sure that that the sampled sequence is a leaf.
    assert(tree[all_sequences[sim]].daughter_nodes[0] == -1);
    assert(tree[all_sequences[sim]].daughter_nodes[1] == -1);
    fprintf(simulated_refs_file, ">simulated_seq_%i\n", all_sequences[sim]+1);
    for(s = 0; s < codon_sequence_length; s++) {
      fprintf(simulated_refs_file, "%s", code_to_char(simulated_sequences[all_sequences[sim]][s]));
    }
    fprintf(simulated_refs_file, "\n");
  }

  // Write the queries to a .fasta file.
  for(sim = num_refs; sim < num_sims; sim++) {
    fprintf(simulated_queries_file, ">simulated_seq_%i_HLA", all_sequences[sim]+1);
    for(h = 0; h < (ploidy * n_genes); h++) fprintf(simulated_queries_file, "_%i", tree[all_sequences[sim]].HLAs[h]);
    fprintf(simulated_queries_file, "\n");
    for(s = 0; s < codon_sequence_length; s++) {
      fprintf(simulated_queries_file, "%s", code_to_char(simulated_sequences[all_sequences[sim]][s]));
    }
    fprintf(simulated_queries_file, "\n");
  }
}

