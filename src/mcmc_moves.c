#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <signal.h> // abort()

#include "util.h"
#include "amino_acids.h"
#include "forward_backward.h"
#include "mcmc_moves.h"

int rec_accept = 0, rec_reject = 0;
int sel_accept = 0, sel_reject = 0;
int merge_accept = 0, merge_reject = 0;
int split_accept = 0, split_reject = 0;
int grow_accept = 0, grow_reject = 0;
int wcoeff_accept = 0, wcoeff_reject = 0;
int mut_accept = 0, mut_reject = 0;

// Generate codon_to_codon_HLA - the probability of switching from codon to
// codon in the presence of an HLA.
//
// Only updates HLA which_HLA, unless which_HLA == -1 when we update all
// codon_to_codon_HLA[num_query][from_codon][sites]
// start_site to end_site (inclusive).
//
// do_esc_or_rev is DO_HLA_ESC, DO_HLA_REV or DO_HLA_BOTH.
//
// DEV: swap codon_to_codon_HLA[num_queries][NUM_CODONS][num_sites]
//      to codon_to_codon_HLA[num_queries][num_sites][NUM_CODONS]
//
void create_codon_to_codon_HLA(int num_sites, int total_num_refs,
                               int num_query, int num_hla_types,
                               int num_ref_site_codons[num_sites], 
                               char ref_site_codons[num_sites][NUM_CODONS],
                               Decimal (*codon_to_codon_HLA)[NUM_CODONS][num_sites],
                               Decimal (*rate_out_of_codons)[num_sites],
                               int start_site, int end_site, int which_HLA,
                               Decimal mu,
                               Decimal (*HLA_select_esc)[num_sites],
                               Decimal (*HLA_select_rev)[num_sites],
                               char *const* hla_types,
                               const char consensus_codons[num_sites],
                               const Decimal omega[num_sites],
                               const int NS_TS_sum_from_consensus[num_sites],
                               const int NS_TV_sum_from_consensus[num_sites],
                               char (*const type_of_change)[NUM_CODONS][NUM_CODONS],
                               char (*const query_codons)[num_sites],
                               // tmp variables
                               Decimal (*A)[NUM_CODONS][num_sites],
                               int do_esc_or_rev)
{
  int q, s, h, c1, from;
  Decimal gamma, HLA_select_esc_prod, HLA_select_rev_prod;

  for(q = 0; q < num_query; q++)
  {
    // If which_HLA is -1 run through everything,
    // if which_HLA is present; evaluate the change.
    if(do_esc_or_rev == DO_HLA_ESC && which_HLA != DO_ALL_HLAS &&
       hla_types[q][which_HLA] == '0') continue;

    for(s = start_site; s <= end_site; s++)
    {
      int consensus_codon = consensus_codons[s];
      Decimal rate_out_of_codons_prime[NUM_CODONS];

      HLA_select_esc_prod = 1;
      HLA_select_rev_prod = HLA_select_rev[0][s];

      for(h = 0; h < num_hla_types; h++) {
        if(hla_types[q][h] == '1') {
          HLA_select_esc_prod *= HLA_select_esc[h][s];
        }
      }

      if(do_esc_or_rev == DO_HLA_BOTH || do_esc_or_rev == DO_HLA_REV)
      {
        // The reversion portion:
        // Set rate_out_of_codons_prime.
        for(from = 0; from < num_ref_site_codons[s]; from++)
        {
          c1 = ref_site_codons[s][from];
          if(c1 != consensus_codon) {

            rate_out_of_codons_prime[c1]
              = rate_out_of_codons[c1][s] +
                (HLA_select_rev_prod - 1) *
                (omega[s] * (kappa * NS_TS[c1][consensus_codon] + 
                             NS_TV[c1][consensus_codon]) +
                 kappa * S_TS[c1][consensus_codon] + S_TV[c1][consensus_codon]);

            // To get the new rate out of codons for reversion,
            // take the standard rate out of codons,
            // add on the extra boost through reversion, but multiply by boolean -
            // is it non consensus to consensus?
            // If it is consensus, make sure that it is a non synonymous change 
            // using NS_TS and NS_TV.
          
            A[q][c1][s] = 1 - total_num_refs/
                              (total_num_refs + rate_out_of_codons_prime[c1] * 4 * mu);
           
            gamma = (Decimal) phi / rate_out_of_codons_prime[c1];

            // Now determine the probability matrix conditional on HLA for the
            // reversion portion.
            int to = query_codons[q][s];

            codon_to_codon_HLA[q][c1][s]
              = A[q][c1][s] * ((type_of_change[s][c1][to] == NONCON_TO_CON) *
                              (gamma * HLA_select_rev_prod * (omega[s] *
                                        (kappa * NS_TS[c1][to] + NS_TV[c1][to]) +
                                        kappa * S_TS[c1][to] + S_TV[c1][to]) +
                              (1-phi) * prior_C2[to]) +
                              (type_of_change[s][c1][to] == NONCON_TO_NONCON) *
                              (gamma * (omega[s] * (kappa * NS_TS[c1][to] + NS_TV[c1][to]) +
                                        kappa * S_TS[c1][to] + S_TV[c1][to]) +
                                        (1-phi) * prior_C2[to])) +
                (c1==to) * (type_of_change[s][c1][to] == NONCON_TO_NONCON) *
                (1 - A[q][c1][s]);       
          }
        }
      }
      if(do_esc_or_rev == DO_HLA_BOTH || do_esc_or_rev == DO_HLA_ESC)
      {
        // The escape portion:
        // Set rate out of codons_prime.
        int to = query_codons[q][s];

        rate_out_of_codons_prime[consensus_codon]
          = rate_out_of_codons[consensus_codon][s] +
            (HLA_select_esc_prod - 1) * omega[s] *
            (kappa * NS_TS_sum_from_consensus[s] + NS_TV_sum_from_consensus[s]);

        A[q][consensus_codon][s]
          = 1 - (total_num_refs) /
                (total_num_refs + rate_out_of_codons_prime[consensus_codon] * 4 * mu);

        if((int)consensus_codon == 15 || (int)consensus_codon == 35)
        {
          gamma = phi / (kappa * NS_TS_sum_from_consensus[s] + NS_TV_sum_from_consensus[s]);
          codon_to_codon_HLA[q][consensus_codon][s]
            = A[q][consensus_codon][s] *
              (gamma * (kappa * NS_TS[consensus_codon][to] + NS_TV[consensus_codon][to]) +
               (1-phi) * prior_C2[to]);
        }
        else
        {
          gamma = (Decimal) phi / rate_out_of_codons_prime[consensus_codon];

          codon_to_codon_HLA[q][consensus_codon][s]
            = A[q][consensus_codon][s] *
              (gamma * (HLA_select_esc_prod * omega[s] *
                        (kappa * NS_TS[consensus_codon][to] + NS_TV[consensus_codon][to]) +
               kappa * S_TS[consensus_codon][to] + S_TV[consensus_codon][to]) +
               (1-phi) * prior_C2[to]);
        }

        if(to == (int)consensus_codon)
        {
            codon_to_codon_HLA[q][consensus_codon][s]
              = (1 - A[q][consensus_codon][s]) * (1-(1-phi)*prior_C2[to]) + 
                (1-phi)*prior_C2[to];
        }

      }
    }
  }
}

void recombination_move(Decimal temperature,
                        int num_sites, int closest_n, int total_num_refs, 
                        int num_query, Decimal *R,
                        Decimal (*F)[closest_n][num_sites],
                        Decimal (*B)[closest_n][num_sites],
                        int F_rnorm[num_query][num_sites],
                        int B_rnorm[num_query][num_sites],
                        Decimal (*codon_to_codon_HLA)[NUM_CODONS][num_sites],
                        char ref_codons[total_num_refs][num_sites],
                        char ref_triplets[total_num_refs][num_sites],
                        Decimal (*tmp_F)[closest_n][num_sites],
                        Decimal *curr_llk_across_queries,
                        int closest_refs[num_query][closest_n],
                        int tmp_F_rnorm[num_query][num_sites],
                        int num_bases,
                        char cons_query_nucls[num_query][num_bases])
{
  status("recombination move\n");
  
  // R[i] is recombination probability just before codon i.
  // Determine between which sites to alter the recombination rate
  int where_rec = rand_lim(num_sites - 1) + 1;
  // where_rec = 1;
  int q, r;
  Decimal ref_sum, sum_llk_across_queries = 0;

  Decimal old_rec_rate = - total_num_refs * log(1 - R[where_rec]);
  // Decimal old_rec_rate = - log(1 - R[where_rec]);
  Decimal rec_rate_prime = old_rec_rate * exp(rand_lim(4) - 2);
  Decimal rec_prime = 1 - exp(-rec_rate_prime / total_num_refs);
  // Decimal rec_prime = 1 - exp(-rec_rate_prime);

  // Accept or reject rec_prime.
  // First store the old recombination probability.
  Decimal rec_old = R[where_rec];

  // Then enter the new recombination probability in its place.
  R[where_rec] = rec_prime;
  Decimal llk_ratio;
  
  if(sample_prior == false) {
    for(q = 0; q < num_query; q++)
    {
      update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                 codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                 F[q], F_rnorm[q], tmp_F_rnorm[q], tmp_F[q],
                                 where_rec, where_rec, closest_refs[q],
                                 num_bases, cons_query_nucls[q]);
      ref_sum = 0;
      for(r = 0; r < closest_n; r++) ref_sum += tmp_F[q][r][where_rec] * B[q][r][where_rec];
      sum_llk_across_queries += log(ref_sum) - (tmp_F_rnorm[q][where_rec] + 
                                                B_rnorm[q][where_rec])* log(BIG);
    }

    // Determine the acceptance ratio.
    status("current log-likelihood: "DECPRINT", proposed log_likelihood: "
           DECPRINT"\n", *curr_llk_across_queries, sum_llk_across_queries);
  
    llk_ratio = (1 / temperature) * (sum_llk_across_queries - *curr_llk_across_queries);

  } else {
    llk_ratio = 0;
  }

  if(isnan(llk_ratio)){
    die("recombination_move llk_ratio: NaN");
  } 
  
  status("lambda_rec: "DECPRINT"\n",lambda_rec);
  Decimal acceptance_ratio = llk_ratio -
                             lambda_rec * (rec_rate_prime - old_rec_rate) +
                             log(rec_rate_prime / old_rec_rate);
  status("likelihood ratio contribution: "DECPRINT"\nprior_contribution: "
         DECPRINT"\nnon symmetric contribution: "DECPRINT"\n",
         llk_ratio,lambda_rec * (rec_rate_prime - old_rec_rate),
         (Decimal) log(rec_rate_prime / old_rec_rate));

  Decimal alpha = exp(MIN2(0, acceptance_ratio));

  if(alpha * RAND_MAX > rand())
  {
    // Accept the move.
    rec_accept++;
    R[where_rec] = rec_prime;

    if(sample_prior == false) {
      *curr_llk_across_queries = sum_llk_across_queries;

      // Update F and B.
      for(q = 0; q < num_query; q++)
      {
        for(r = 0; r < closest_n; r++)
          F[q][r][where_rec] = tmp_F[q][r][where_rec];
        F_rnorm[q][where_rec] = tmp_F_rnorm[q][where_rec];

        if(where_rec < num_sites-1) {
          update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                     codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                     F[q], F_rnorm[q], F_rnorm[q],
                                     F[q], where_rec+1, num_sites-1, closest_refs[q],
                                     num_bases, cons_query_nucls[q]);
        }

        update_B_numerical_recipes(num_sites, closest_n, total_num_refs,
                                   codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                   B[q], B_rnorm[q], where_rec, closest_refs[q],
                                   num_bases, cons_query_nucls[q]);
      }
    }
  }
  else
  {
    // Reject the move.
    rec_reject++;
    // Replace the old recombination probability.
    R[where_rec] = rec_old;
  }
}

void selection_move(Decimal temperature,
                    int num_sites, int closest_n, int total_num_refs,
                    int num_query, int num_hla_types,
                    int num_ref_site_codons[num_sites], 
                    char ref_site_codons[num_sites][NUM_CODONS],
                    Decimal mu, Decimal *R,
                    Decimal (*F)[closest_n][num_sites],
                    Decimal (*B)[closest_n][num_sites],
                    int F_rnorm[num_query][num_sites],
                    int B_rnorm[num_query][num_sites],
                    Decimal (*codon_to_codon_HLA)[NUM_CODONS][num_sites],
                    Decimal (*rate_out_of_codons)[num_sites],
                    Decimal (*HLA_select_esc)[num_sites],
                    Decimal (*HLA_select_rev)[num_sites],
                    char *const* hla_types,
                    const char consensus_codons[num_sites],
                    char (*ref_codons)[num_sites],
                    char (*ref_triplets)[num_sites],
                    char (*query_codons)[num_sites],
                    Decimal omega[num_sites],
                    const int NS_TS_sum_from_consensus[num_sites],
                    const int NS_TV_sum_from_consensus[num_sites],
                    char (*const type_of_change)[NUM_CODONS][NUM_CODONS],
                    // tmp variables.
                    Decimal (*A)[NUM_CODONS][num_sites],
                    Decimal (*tmp_codon_to_codon_HLA)[NUM_CODONS][num_sites],
                    Decimal (*tmp_F)[closest_n][num_sites],
                    Decimal *curr_llk_across_queries,
                    int closest_refs[num_query][closest_n],
                    int tmp_F_rnorm[num_query][num_sites],
                    Decimal (*A_tmp)[NUM_CODONS][num_sites],
                    int num_bases, char cons_query_nucls[num_query][num_bases])
                    
{
  status("selection move\n");
  
  // Omega[i] is the selection coefficient at site i. 
  int rand_site = rand_lim(num_sites);
   
  Decimal omega_prime = omega[rand_site] * exp(rand_lim(2)-1);

  // Accept or reject omega_prime.
  // First store the old selection parameter.
  Decimal prev_omega = omega[rand_site];

  // Then enter the new rate in its place.
  omega[rand_site] = omega_prime;

  Decimal tmp_codon_rate[NUM_CODONS];
  Decimal llk_ratio;
  Decimal ref_sum, sum_llk_across_queries = 0;
  int c1, q, r, from;
  
  // Store the current rate_out_of_codons.
  // DEV: Merge these loops?
  if(sample_prior == false) {
    for(from = 0; from < num_ref_site_codons[rand_site]; from++)
    {
      c1 = ref_site_codons[rand_site][from];
      tmp_codon_rate[c1] = rate_out_of_codons[c1][rand_site];
      rate_out_of_codons[c1][rand_site]
        = kappa * beta_S[c1] + beta_V[c1] +
          omega_prime * (kappa * alpha_S[c1] + alpha_V[c1]);
    }

    for(q = 0; q < num_query; q++) {
      for(from = 0; from < num_ref_site_codons[rand_site]; from++) {
        c1 = ref_site_codons[rand_site][from];
        tmp_codon_to_codon_HLA[q][c1][rand_site]
          = codon_to_codon_HLA[q][c1][rand_site];
        A_tmp[q][c1][rand_site] = A[q][c1][rand_site];
      }
    }

    // Determine new probability matrix.
    create_codon_to_codon_HLA(num_sites, total_num_refs,
                              num_query, num_hla_types,
                              num_ref_site_codons, ref_site_codons,
                              codon_to_codon_HLA, rate_out_of_codons,
                              rand_site, rand_site, DO_ALL_HLAS, mu,
                              HLA_select_esc, HLA_select_rev, hla_types, 
                              consensus_codons, omega, 
                              NS_TS_sum_from_consensus, NS_TV_sum_from_consensus,
                              type_of_change, query_codons, A, DO_HLA_BOTH);

    // Evaluate the likelihood based on this change and decide whether to 
    // accept or reject the move.

    for(q = 0; q < num_query; q++)
    {
      update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                 codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                 F[q], F_rnorm[q], tmp_F_rnorm[q], tmp_F[q],
                                 rand_site, rand_site, closest_refs[q],
                                 num_bases, cons_query_nucls[q]);

      ref_sum = 0;
      for(r = 0; r < closest_n; r++)
        ref_sum += tmp_F[q][r][rand_site] * B[q][r][rand_site];
      
      if(isnan(log(ref_sum)))
        die("log(ref_sum) is NaN");
  
      sum_llk_across_queries += log(ref_sum) - (tmp_F_rnorm[q][rand_site] + 
                                                B_rnorm[q][rand_site])* log(BIG);
    }

    // Determine the acceptance ratio.
    status("current log-likelihood: "DECPRINT", proposed log_likelihood: "
           DECPRINT"\n", *curr_llk_across_queries, sum_llk_across_queries); 
  
    llk_ratio = (1 / temperature) * (sum_llk_across_queries - *curr_llk_across_queries);

    if(isnan(llk_ratio)){
      die("selection_move llk_ratio: NaN");
    }
  } else {
    llk_ratio = 0;
  }

  Decimal prior_diff;
  switch(omega_prior)
  {
    case EXPONENTIAL_PRIOR:
      prior_diff = - lambda_sel * (omega_prime - prev_omega);
      break;
    case GAMMA_PRIOR:
      prior_diff = (omega_gamma_k_shape - 1) * log(omega_prime / prev_omega) +
                   (prev_omega - omega_prime) / omega_gamma_theta_scale;
      break;
    case LOG_NORMAL_PRIOR:
      prior_diff = log(prev_omega / omega_prime) -
                   (pow(log(omega_prime / omega_mu_sel), 2) -
                    pow(log(prev_omega / omega_mu_sel), 2)) /
                   (2 * pow(omega_sigma_sel, 2));
      break;
    case FLAT_PRIOR:
      prior_diff = 0;
      break;
    default:
      die("No prior detected!");
  }
  
  Decimal acceptance_ratio = llk_ratio + prior_diff + log(omega_prime / prev_omega);
  Decimal alpha = exp(MIN2(0, acceptance_ratio));

  if(alpha * RAND_MAX > rand())
  {
    // Accept the move.
    sel_accept++;
    if(sample_prior == false) {
      *curr_llk_across_queries = sum_llk_across_queries;
    
      // Update F and B.
      for(q = 0; q < num_query; q++)
      {
        for(r = 0; r < closest_n; r++)
          F[q][r][rand_site] = tmp_F[q][r][rand_site];
        F_rnorm[q][rand_site] = tmp_F_rnorm[q][rand_site];

        if(rand_site < num_sites-1)
        {
          update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                     codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                     F[q], F_rnorm[q], F_rnorm[q], F[q],
                                     rand_site+1, num_sites-1, closest_refs[q],
                                     num_bases, cons_query_nucls[q]);
        }
          update_B_numerical_recipes(num_sites, closest_n, total_num_refs,
                                     codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                     B[q], B_rnorm[q], rand_site, closest_refs[q],
                                     num_bases, cons_query_nucls[q]);
      }
    }
  }
  else
  {
    // Reject the move.
    sel_reject++;

    if(sample_prior == false) {
      // Reset the values.
      for(from = 0; from < num_ref_site_codons[rand_site]; from++) {
        c1 = ref_site_codons[rand_site][from];
        rate_out_of_codons[c1][rand_site] = tmp_codon_rate[c1];
      }

      for(q = 0; q < num_query; q++) {
        for(from = 0; from < num_ref_site_codons[rand_site]; from++) {
          c1 = ref_site_codons[rand_site][from];
          codon_to_codon_HLA[q][c1][rand_site]
             = tmp_codon_to_codon_HLA[q][c1][rand_site];
          A[q][c1][rand_site] = A_tmp[q][c1][rand_site];
        }
      }
    }
    
    // Replace the old selection parameter.
    omega[rand_site] = prev_omega;

  }
}

Decimal evaluate_prior_diff_merge(int HLA_coeff_prior,
                                  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal mu_sel_norm_prior, Decimal sigma_sel_norm_prior,
                                  Decimal gamma_k_shape, Decimal gamma_theta_scale,
                                  Decimal old_coeff0, Decimal old_coeff1,
                                  Decimal new_coeff, Decimal GAMMA_CONST)
{
  Decimal prior_diff; 

  switch(HLA_coeff_prior)
  {
    case LOG_NORMAL_PRIOR:
      // Log-normal prior ratio (and binomial for choosing splitting sites).
      prior_diff = log(((1 - prior_add_window) * old_coeff0 * old_coeff1 * 
                       sqrt(2 * M_PI) * sigma_sel_prior) / (new_coeff * prior_add_window)) 
                   + (- pow(log(new_coeff / mu_sel_prior), 2) + 
                      pow(log(old_coeff0 / mu_sel_prior), 2) +
                      pow(log(old_coeff1 / mu_sel_prior), 2)) / (2 * pow(sigma_sel_prior, 2));
      break;
    case NORMAL_PRIOR:
      // Normal prior ratio (and binomial for choosing splitting sites).
      prior_diff = log(sigma_sel_norm_prior * sqrt(2 * M_PI)) + 
                   (pow(old_coeff0 - mu_sel_norm_prior, 2) + 
                    pow(old_coeff1 - mu_sel_norm_prior, 2) - 
                    pow(new_coeff - mu_sel_norm_prior, 2)) / (2 * pow(sigma_sel_norm_prior,2)) + 
                   log((1 - prior_add_window)/prior_add_window);
      break;
    case FLAT_PRIOR:
      // Flat prior (and binomial for choosing splitting sites).
      prior_diff = log((1 - prior_add_window) / prior_add_window);
      break;
    case GAMMA_PRIOR:
      // Gamma prior ratio (and binomial for choosing splitting sites).
      prior_diff = (gamma_k_shape - 1) * (log (new_coeff) - log(old_coeff0 * old_coeff1)) + 
                   log((1 - prior_add_window)/prior_add_window) + 
                   (old_coeff0 + old_coeff1 - new_coeff)/gamma_theta_scale + GAMMA_CONST;
      break;
    default: die("Unknown prior!");
  }
  return(prior_diff);
}

Decimal evaluate_prior_diff_split(int HLA_coeff_prior,
                                  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal mu_sel_norm_prior, Decimal sigma_sel_norm_prior,
                                  Decimal gamma_k_shape, Decimal gamma_theta_scale,
                                  Decimal old_coeff0,
                                  Decimal new_coeff0, Decimal new_coeff1,
                                  Decimal GAMMA_CONST)
{
  Decimal prior_diff;

  switch(HLA_coeff_prior)
  {
    case LOG_NORMAL_PRIOR:
      // Log-normal prior ratio (and binomial for choosing splitting sites).
      prior_diff = log((prior_add_window * old_coeff0) / 
                     ((1 - prior_add_window) * new_coeff0 * new_coeff1 * sqrt(2 * M_PI) * 
                     sigma_sel_prior)) 
                   + (- pow(log(new_coeff0 / mu_sel_prior), 2) - 
                      pow(log(new_coeff1 / mu_sel_prior), 2) + 
                      pow(log(old_coeff0 / mu_sel_prior), 2)) / (2 * pow(sigma_sel_prior, 2));
      break;
    case NORMAL_PRIOR:
      // Normal prior ratio (and binomial for choosing splitting sites).
      prior_diff = - log(sigma_sel_norm_prior * sqrt(2 * M_PI)) + 
                   (pow(old_coeff0 - mu_sel_norm_prior, 2) - 
                    pow(new_coeff0 - mu_sel_norm_prior, 2) + 
                    pow(new_coeff1 - mu_sel_norm_prior, 2)) / (2 * pow(sigma_sel_norm_prior,2)) +
                   log(prior_add_window/(1 - prior_add_window));
      break;
    case FLAT_PRIOR:
      // Flat prior (and binomial for choosing splitting sites).
      prior_diff = log((prior_add_window) / (1 - prior_add_window));
      break;
    case GAMMA_PRIOR:          
      // Gamma prior ratio (and binomial for choosing splitting sites).
      prior_diff = (gamma_k_shape - 1) * (log (new_coeff0 * new_coeff1) - log(old_coeff0)) + 
                   log(prior_add_window/(1 - prior_add_window)) - 
                   (new_coeff0 + new_coeff1 - old_coeff0)/gamma_theta_scale - GAMMA_CONST;
      break;
    default: die("Unknown prior!");
  } 
  return(prior_diff);
}

Decimal evaluate_prior_diff_coeff(int HLA_coeff_prior,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal mu_sel_norm_prior, Decimal sigma_sel_norm_prior,
                                  Decimal gamma_k_shape, Decimal gamma_theta_scale,
                                  Decimal old_coeff0, Decimal new_coeff)
{
  Decimal prior_diff;

  switch(HLA_coeff_prior)
  {
    case LOG_NORMAL_PRIOR:
      // Log-normal prior ratio.
      prior_diff = log(old_coeff0 / new_coeff) -
                   (pow(log(new_coeff / mu_sel_prior), 2) -
                    pow(log(old_coeff0 / mu_sel_prior), 2)) /
                    (2 * pow(sigma_sel_prior, 2));
      break;
    case NORMAL_PRIOR:
      // Normal prior ratio.
      prior_diff = (pow((old_coeff0 - mu_sel_norm_prior), 2) - pow((new_coeff - mu_sel_norm_prior),2)) / 
                   (2 * pow(sigma_sel_norm_prior, 2));
      break;
    case FLAT_PRIOR:
      // Flat prior ratio.
      prior_diff = 0;
      break;
    case GAMMA_PRIOR:
      // Gamma prior ratio.
      prior_diff = (gamma_k_shape - 1) * log(new_coeff / old_coeff0) + 
                   (old_coeff0 - new_coeff) / gamma_theta_scale;
      break;
    default: die("Unknown prior!");
  }
  return(prior_diff);
}

void window_move(Decimal temperature,
                 int num_sites, int closest_n, int total_num_refs,
                 int num_query, int num_hla_types,
                 int num_ref_site_codons[num_sites], 
                 char ref_site_codons[num_sites][NUM_CODONS],
                 int do_esc_or_rev, int window_action,
                 Decimal mu, Decimal *R,
                 Decimal (*F)[closest_n][num_sites],
                 Decimal (*B)[closest_n][num_sites],
                 int F_rnorm[num_query][num_sites],
                 int B_rnorm[num_query][num_sites],
                 Decimal (*codon_to_codon_HLA)[NUM_CODONS][num_sites],
                 Decimal (*rate_out_of_codons)[num_sites],
                 Decimal (*HLA_select_esc)[num_sites],
                 Decimal (*HLA_select_rev)[num_sites],
                 char *const* hla_types,
                 const char consensus_codons[num_sites],
                 char (*ref_codons)[num_sites],
                 char (*ref_triplets)[num_sites],
                 char (*query_codons)[num_sites],
                 const Decimal omega[num_sites],
                 const int NS_TS_sum_from_consensus[num_sites],
                 const int NS_TV_sum_from_consensus[num_sites],
                 char (*const type_of_change)[NUM_CODONS][NUM_CODONS],
                 Decimal (*A)[NUM_CODONS][num_sites],
                 HLASelection *hla_windows,
                 // tmp variables.
                 Decimal (*tmp_codon_to_codon_HLA)[NUM_CODONS][num_sites],
                 Decimal (*tmp_F)[closest_n][num_sites],
                 Decimal *curr_llk_across_queries,
                 int closest_refs[num_query][closest_n],
                 int tmp_F_rnorm[num_query][num_sites],
                 Decimal (*A_tmp)[NUM_CODONS][num_sites],
                 int num_bases, char cons_query_nucls[num_query][num_bases],
                 Decimal GAMMA_CONST)
{
  int h, w, s, q, r, c1, from;

  // Pick the HLA type at random.
  h = rand_lim(num_hla_types);

  status("DO ESC OR REV: %i\n", do_esc_or_rev);

  // Reversion windows not done by HLA anymore.
  // Keep the structure but only use first reversion window set.
  HLAWindowSet *wset = (do_esc_or_rev == DO_HLA_ESC ? hla_windows[h].esc
                                                    : hla_windows[0].rev);

  SelectionWindow *windows = wset->windows;
  Decimal *hla_selection = (do_esc_or_rev == DO_HLA_ESC ? HLA_select_esc[h]
                                                        : HLA_select_rev[0]);
  Decimal prior_add_window = (do_esc_or_rev == DO_HLA_ESC ? prior_add_window_esc
                                                          : prior_add_window_rev);
  Decimal HLA_coeff_prior = (do_esc_or_rev == DO_HLA_ESC ? HLA_coeff_prior_esc
                                                         : HLA_coeff_prior_rev);
  Decimal mu_sel_prior = (do_esc_or_rev == DO_HLA_ESC ? mu_sel_prior_esc
                                                      : mu_sel_prior_rev);
  Decimal sigma_sel_prior = (do_esc_or_rev == DO_HLA_ESC ? sigma_sel_prior_esc
                                                         : sigma_sel_prior_rev);
  Decimal mu_sel_norm_prior = (do_esc_or_rev == DO_HLA_ESC ? mu_sel_norm_prior_esc
                                                           : mu_sel_norm_prior_rev);
  Decimal sigma_sel_norm_prior = (do_esc_or_rev == DO_HLA_ESC ? sigma_sel_norm_prior_esc
                                                              : sigma_sel_norm_prior_rev);
  Decimal gamma_k_shape = (do_esc_or_rev == DO_HLA_ESC ? gamma_k_shape_esc
                                                       : gamma_k_shape_rev);
  Decimal gamma_theta_scale = (do_esc_or_rev == DO_HLA_ESC ? gamma_theta_scale_esc
                                                           : gamma_theta_scale_rev);

  Decimal prob_split = 1, prob_merge = 1;
  if(window_action == DO_MERGE_OR_SPLIT)
  {
    // Choose between merge and split.
    int transition_pts = wset->num_windows - 1;
    prob_split = MIN2(1, (Decimal)(num_sites - wset->num_windows) /
                                 (Decimal)wset->num_windows *
                                  prior_add_window / (1-prior_add_window));
    status("PROB SPLIT: "DECPRINT"\n", prob_split);
    prob_merge = MIN2(1, (Decimal)transition_pts / (Decimal)(num_sites - transition_pts) *
                                 (1-prior_add_window) / prior_add_window);

    double urand = rand_lim(1);

    if(urand < prob_split / (prob_split + prob_merge))
      window_action = DO_SPLIT;
    else
      window_action = DO_MERGE;
  }

  int old_start_site;
  Decimal old_coeff0, old_coeff1;
  
  // Perform the move.
  switch(window_action)
  {
    case DO_MERGE:
      status("merge move\n");
      status("NUMBER OF WINDOWS: %i\n", wset->num_windows);
      w = selec_window_merge_rnd(wset, &old_start_site, &old_coeff0, &old_coeff1);
      break;
    case DO_SPLIT:
      status("split move\n");
      status("NUMBER OF WINDOWS: %i\n", wset->num_windows);
      w = selec_window_split_rnd(wset, num_sites, &old_coeff0);
      break;
    case DO_GROW:
      status("grow move\n");
      w = selec_window_grow_rnd(wset, &old_start_site);
      if(w == -1) {
        grow_reject++;
        return;
      }
      break;
    case DO_COEFF:
      status("HLA coeff move\n");
      w = selec_window_change_coeff_rnd(wset, &old_coeff0);
      break;
    default:
      die("No action!");
  }

  int start_site, end_site;

  if(window_action != DO_GROW) {
    start_site = windows[w].start;
    end_site = start_site + windows[w].length - 1;

    for(s = start_site; s <= end_site; s++) {
      hla_selection[s] = windows[w].coeff;
    }

    if(window_action == DO_SPLIT) {
      end_site += windows[w+1].length;
      for(; s <= end_site; s++) {
        hla_selection[s] = windows[w+1].coeff;
      }
    }
  } else {
    // The site which is moved.
    int new_start_site = windows[w+1].start;
    Decimal new_selection;
    if(old_start_site < new_start_site) {
      start_site = old_start_site;
      end_site = new_start_site -1;
      new_selection = windows[w].coeff;
    } else {
      start_site = new_start_site;
      end_site = old_start_site - 1;
      new_selection = windows[w+1].coeff;
    }
    for(s = start_site; s <= end_site; s++) {
      hla_selection[s] = new_selection;
    }
  } 

  // Set up tmp variables.
  for(q = 0; q < num_query; q++) {
    for(s = start_site; s <= end_site; s++) {
      for(from = 0; from < num_ref_site_codons[s]; from++)
      { 
        c1 = ref_site_codons[s][from];
        tmp_codon_to_codon_HLA[q][c1][s] = codon_to_codon_HLA[q][c1][s];
        A_tmp[q][c1][s] = A[q][c1][s];
      }
    }
  }

  Decimal ref_sum, sum_llk_across_queries = 0;
  Decimal llk_ratio;

  if(sample_prior == false) {
    // Determine new probability matrix.
    create_codon_to_codon_HLA(num_sites, total_num_refs,
                              num_query, num_hla_types, 
                              num_ref_site_codons, ref_site_codons,
                              codon_to_codon_HLA, rate_out_of_codons,
                              start_site, end_site, h, mu,
                              HLA_select_esc, HLA_select_rev, hla_types,
                              consensus_codons, omega,
                              NS_TS_sum_from_consensus, NS_TV_sum_from_consensus,
                              type_of_change, query_codons, A,
                              do_esc_or_rev);

    // Evaluate the likelihood based on this change and decide whether to 
    // accept or reject the move.

    for(q = 0; q < num_query; q++)
    {
      update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                 codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                 F[q], F_rnorm[q], tmp_F_rnorm[q], tmp_F[q],
                                 start_site, end_site, closest_refs[q],
                                 num_bases, cons_query_nucls[q]);

      ref_sum = 0;
      for(r = 0; r < closest_n; r++)
        ref_sum += tmp_F[q][r][end_site] * B[q][r][end_site];
      
      if(isnan(log(ref_sum)))
        die("log(ref_sum) is NaN");

      sum_llk_across_queries += log(ref_sum) - (tmp_F_rnorm[q][end_site] + 
                                B_rnorm[q][end_site])* log(BIG);
    }

    // Determine the acceptance ratio.
    status("current log-likelihood: "DECPRINT", proposed log-likelihood: "DECPRINT"\n",
           *curr_llk_across_queries, sum_llk_across_queries);
  
    llk_ratio = (1 / temperature) * (sum_llk_across_queries - *curr_llk_across_queries);

    if(isnan(llk_ratio)){
      switch(window_action)
      {
        case DO_MERGE:
          die("window_move merge llk_ratio: NaN");
          break;
        case DO_SPLIT:
          die("window_move split llk_ratio: NaN");
          break;
        case DO_GROW:
          die("window_move grow llk_ratio: NaN");
          break;
        case DO_COEFF:
          die("window_move coeff llk_ratio: NaN");
          break;
        default:
          die("window_move llk_ratio: NaN");
      }
    }
  } else {
    llk_ratio = 0;
  }

  Decimal acceptance_ratio;
  Decimal coeff_sum, prior;
  Decimal prob_birth_in_dim_minus_one, prob_death, prob_birth_over_death;
  Decimal prob_death_in_dim_plus_one, prob_birth, prob_death_over_birth;
  Decimal prior_diff, new_coeff;
  Decimal new_coeff0, new_coeff1;

  switch(window_action)
  {
    case DO_MERGE:
      // Note: the move has been made, so the number of windows is one less than the 
      // 'current' number of windows. 
      new_coeff = wset->windows[w].coeff;
      prior = prior_add_window / (1 - prior_add_window);

      status("PROB MERGE: "DECPRINT"\n", prob_merge);

      prob_birth_in_dim_minus_one = ((num_sites - (Decimal)(wset->num_windows)) / 
                                    (Decimal)(wset->num_windows)) * prior;
      prob_death = MIN2(1, 1 / prob_birth_in_dim_minus_one);

      status("THE SAME AS PROB DEATH?: "DECPRINT"\n", prob_death);

      prob_birth_in_dim_minus_one = MIN2(1, prob_birth_in_dim_minus_one);
      prob_birth_over_death = prob_birth_in_dim_minus_one / prob_death;

      prior_diff = evaluate_prior_diff_merge(HLA_coeff_prior, prior_add_window,
                                             mu_sel_prior, sigma_sel_prior,
                                             mu_sel_norm_prior, sigma_sel_norm_prior,
                                             gamma_k_shape, gamma_theta_scale,                                    
                                             old_coeff0, old_coeff1, new_coeff,
                                             GAMMA_CONST);
      
      acceptance_ratio = llk_ratio + prior_diff +
                         log(prob_birth_over_death *  
                             ((wset->num_windows) * new_coeff) / 
                             ((num_sites - (wset->num_windows)) * pow(old_coeff0 + old_coeff1, 2)));
      break;
    case DO_SPLIT:
      // Note: the move has been made, so the number of windows is one more than the 
      // 'current' number of windows. 
      new_coeff0 = wset->windows[w].coeff;
      new_coeff1 = wset->windows[w+1].coeff;
      status("NUMBER OF WINDOWS: %i\n", wset->num_windows);

      coeff_sum = new_coeff0 + new_coeff1;
      prior = (1 - prior_add_window) / prior_add_window;

      status("PROB SPLIT: "DECPRINT"\n", prob_split);

      prob_death_in_dim_plus_one = (Decimal)(wset->num_windows - 1) / 
                                   (Decimal)(num_sites - (wset->num_windows - 1)) * prior;
      prob_birth = MIN2(1, 1 / prob_death_in_dim_plus_one);
      
      status("THE SAME AS PROB BIRTH?: "DECPRINT"\n", prob_birth);

      prob_death_in_dim_plus_one  = MIN2(1, prob_death_in_dim_plus_one);
      prob_death_over_birth = prob_death_in_dim_plus_one / prob_birth;

      prior_diff = evaluate_prior_diff_split(HLA_coeff_prior, prior_add_window,
                                             mu_sel_prior, sigma_sel_prior,
                                             mu_sel_norm_prior, sigma_sel_norm_prior,
                                             gamma_k_shape, gamma_theta_scale, 
                                             old_coeff0, new_coeff0, new_coeff1,
                                             GAMMA_CONST);
      
      acceptance_ratio = llk_ratio + prior_diff +
                         log(prob_death_over_birth *
                             ((num_sites - (wset->num_windows - 1)) * pow(coeff_sum, 2)) / 
                             ((wset->num_windows - 1) * old_coeff0));
     break;
    case DO_GROW:
      acceptance_ratio = llk_ratio;
      break;
    case DO_COEFF:
      new_coeff = wset->windows[w].coeff;

      prior_diff = evaluate_prior_diff_coeff(HLA_coeff_prior,
                                             mu_sel_prior, sigma_sel_prior,
                                             mu_sel_norm_prior, sigma_sel_norm_prior,
                                             gamma_k_shape, gamma_theta_scale, 
                                             old_coeff0, new_coeff);

      acceptance_ratio = llk_ratio + prior_diff +
                         log(new_coeff / old_coeff0);
      break;
  }

  Decimal alpha = exp(MIN2(0, acceptance_ratio));

  if(alpha * RAND_MAX > rand())
  {
    // Accept the move.
    switch(window_action) {
      case DO_MERGE: merge_accept++; break;
      case DO_SPLIT: split_accept++; break;
      case DO_GROW: 
      grow_accept++; 
      selec_window_grow_accept(wset, wset->windows[w+1].start, old_start_site, num_sites);
      break;
      case DO_COEFF: wcoeff_accept++; break;
    }
    
    if(sample_prior == false) {
      *curr_llk_across_queries = sum_llk_across_queries;
      // Update F and B.
      for(q = 0; q < num_query; q++)
      {
        for(s = start_site; s <= end_site; s++)
        {
          for(r = 0; r < closest_n; r++)
            F[q][r][s] = tmp_F[q][r][s];

          F_rnorm[q][s] = tmp_F_rnorm[q][s];
        }

        if(end_site + 1 < num_sites)
        {
          update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                     codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                     F[q], F_rnorm[q], F_rnorm[q], F[q],
                                     end_site+1, num_sites-1, closest_refs[q],
                                     num_bases, cons_query_nucls[q]);
        }
      
        update_B_numerical_recipes(num_sites, closest_n, total_num_refs,
                                   codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                   B[q], B_rnorm[q], end_site, closest_refs[q],
                                   num_bases, cons_query_nucls[q]);
      }
    }
  }
  else
  {
    // Reject the move.
    switch(window_action) {
      case DO_MERGE: merge_reject++; break;
      case DO_SPLIT: split_reject++; break;
      case DO_GROW: grow_reject++; break;
      case DO_COEFF: wcoeff_reject++; break;
    }

    if(sample_prior == false) {
      // Reset the values.
      for(q = 0; q < num_query; q++) {
        for(s = start_site; s <= end_site; s++) {
          for(from = 0; from < num_ref_site_codons[s]; from++) {
            c1 = ref_site_codons[s][from];
            codon_to_codon_HLA[q][c1][s] = tmp_codon_to_codon_HLA[q][c1][s];
            A[q][c1][s] = A_tmp[q][c1][s];
          }
        }
      }
    }

    // Undo the window move which was applied.  
    switch(window_action) {
      case DO_MERGE:
        selec_window_merge_undo(wset, w, old_start_site, old_coeff0, old_coeff1);
        break;
      case DO_SPLIT:
        selec_window_split_undo(wset, w, old_coeff0);
        break;
      case DO_GROW:
        selec_window_grow_undo(wset, w, old_start_site);
        break;
      case DO_COEFF:
        selec_window_change_coeff_undo(wset, w, old_coeff0);
        break;
    }

    start_site = windows[w].start;
    end_site = start_site + windows[w].length - 1;

    for(s = start_site; s <= end_site; s++)
      hla_selection[s] = windows[w].coeff;

    if(window_action == DO_MERGE || window_action == DO_GROW)
    {
      start_site = windows[w+1].start;
      end_site = start_site + windows[w+1].length - 1;
      for(s = start_site; s <= end_site; s++)
        hla_selection[s] = windows[w+1].coeff;
    }
  }
}

void mutation_move(Decimal temperature,
                   int num_sites, int closest_n, int total_num_refs, 
                   int num_query, int num_ref_site_codons[num_sites], 
                   char ref_site_codons[num_sites][NUM_CODONS],
                   Decimal *mu, Decimal *R,
                   Decimal (**F_old)[closest_n][num_sites],
                   Decimal (*B)[closest_n][num_sites],
                   int (**F_rnorm_old)[num_sites],
                   int (*B_rnorm)[num_sites],
                   Decimal (**codon_to_codon_HLA_old)[NUM_CODONS][num_sites],
                   char (*ref_codons)[num_sites],
                   char (*ref_triplets)[num_sites],
                   char (*query_codons)[num_sites],
                   Decimal (**A_old)[NUM_CODONS][num_sites],
                   // tmp variables.
                   Decimal (**A_tmp)[NUM_CODONS][num_sites],
                   Decimal (**codon_to_codon_HLA_tmp)[NUM_CODONS][num_sites],
                   Decimal (**F_tmp)[closest_n][num_sites],
                   Decimal *curr_llk_across_queries,
                   int closest_refs[num_query][closest_n],
                   int (**F_rnorm_tmp)[num_sites],
                   int num_bases, char cons_query_nucls[num_query][num_bases])
{
  status("mutation move\n");
  int q, r, s, c1, from, to;

  // Pointers to tmp variables, these switch if we accept the move.
  Decimal (*A)[NUM_CODONS][num_sites] = *A_old;
  Decimal (*tmp_A)[NUM_CODONS][num_sites] = *A_tmp;
  Decimal (*tmp_codon_to_codon_HLA)[NUM_CODONS][num_sites] = *codon_to_codon_HLA_tmp;
  Decimal (*codon_to_codon_HLA)[NUM_CODONS][num_sites] = *codon_to_codon_HLA_old;

  // Note: rate_out_of_codons_prime has a different meaning within this move.
  // Here it is rate_out_of_codons after it has been modified by the proposed
  // mu_prime - nothing to do with HLA. 
  // In previous moves it referred to the rate out of codons after it had been 
  // modified by the individuals HLA profile.

  Decimal rate_out_of_codons_prime;
  Decimal mu_prime = *mu * exp((rand_lim(1) - 0.5)/4);
  Decimal ref_sum, sum_llk_across_queries = 0;
  Decimal llk_ratio;

  if(sample_prior == false) {
    for(q = 0; q < num_query; q++)
    {
      for(s = 0; s < num_sites; s++)
      {
        to = query_codons[q][s];
        for(from = 0; from < num_ref_site_codons[s]; from++)
        {
          c1 = ref_site_codons[s][from];
          if(A[q][c1][s] == 0) {
            rate_out_of_codons_prime = 0;
          } else {
            rate_out_of_codons_prime = (total_num_refs / (1 - A[q][c1][s]) - total_num_refs) / (4 * *mu);
          }
          tmp_A[q][c1][s] = 1 - total_num_refs/(total_num_refs + rate_out_of_codons_prime * 4 * mu_prime);

          if(tmp_A[q][c1][s] > 1) status("A_tmp: "DECPRINT"\n", tmp_A[q][c1][s]);
          if(to != c1) {
            if(A[q][c1][s] == 0) {
              tmp_codon_to_codon_HLA[q][c1][s] = 0;
            } else {
              tmp_codon_to_codon_HLA[q][c1][s] = codon_to_codon_HLA[q][c1][s] *
                                                 tmp_A[q][c1][s] / A[q][c1][s];
            }
          } else {
            tmp_codon_to_codon_HLA[q][c1][s] = 1 - tmp_A[q][c1][s] + 
                                               (1 - phi) * prior_C2[c1] * tmp_A[q][c1][s];
          }
        }
      }
    }

    Decimal (*F)[closest_n][num_sites] = *F_old;
    int (*F_rnorm)[num_sites] = *F_rnorm_old;

    Decimal (*tmp_F)[closest_n][num_sites] = *F_tmp;
    int (*tmp_F_rnorm)[num_sites] = *F_rnorm_tmp;

    // Evaluate the likelihood based on this change and decide whether to 
    // accept or reject the move.

    for(q = 0; q < num_query; q++)
    {
      update_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                 tmp_codon_to_codon_HLA[q], ref_codons, ref_triplets, R,
                                 F[q], F_rnorm[q], tmp_F_rnorm[q], tmp_F[q],
                                 0, num_sites-1, closest_refs[q],
                                 num_bases, cons_query_nucls[q]);
      ref_sum = 0;

      for(r = 0; r < closest_n; r++)
        ref_sum += tmp_F[q][r][num_sites - 1];

      sum_llk_across_queries += log(ref_sum) - (tmp_F_rnorm[q][num_sites-1])* log(BIG);
    }

    // Determine the acceptance ratio.
    status("current log-likelihood: "DECPRINT", proposed log-likelihood: "
           DECPRINT"\n", *curr_llk_across_queries, sum_llk_across_queries);

    llk_ratio = (1 / temperature) * (sum_llk_across_queries - *curr_llk_across_queries);

    if(isnan(llk_ratio))
      die("mutation_move llk_ratio: NaN");

  } else {
    llk_ratio = 0;
  }

  Decimal acceptance_ratio = llk_ratio - mu_prior * (mu_prime - *mu) + log(mu_prime / *mu);
  Decimal alpha = exp(MIN2((Decimal)0, acceptance_ratio));

  if(alpha * RAND_MAX > rand())
  {
    // Accept the move.
    mut_accept++;
    *mu = mu_prime;
    if(sample_prior == false) {
      // Take new values.
      Decimal (*swap_A)[NUM_CODONS][num_sites];
      Decimal (*swap_F)[closest_n][num_sites];
      Decimal (*swap_codon_to_codon_HLA)[NUM_CODONS][num_sites];
      int (*swap_F_rnorm)[num_sites];

      // Swap the pointers.
      SWAP(*A_old, *A_tmp, swap_A);
      SWAP(*F_old, *F_tmp, swap_F);
      SWAP(*codon_to_codon_HLA_old, *codon_to_codon_HLA_tmp, swap_codon_to_codon_HLA);
      SWAP(*F_rnorm_old, *F_rnorm_tmp, swap_F_rnorm);

      for(q = 0; q < num_query; q++)
      {
        update_B_numerical_recipes(num_sites, closest_n, total_num_refs,
                                   (*codon_to_codon_HLA_old)[q], ref_codons,
                                   ref_triplets, R, B[q], B_rnorm[q], num_sites-1, 
                                   closest_refs[q],
                                   num_bases, cons_query_nucls[q]);
      }
      // Update the log-likelihood.
      *curr_llk_across_queries = sum_llk_across_queries;
    }
  }
  else
  {
    // Reject the move.
    mut_reject++;
  }
}
