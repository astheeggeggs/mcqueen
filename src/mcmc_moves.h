#ifndef MCMC_MOVES_H_
#define MCMC_MOVES_H_

#include "selection_window.h"

#define DO_HLA_ESC 0
#define DO_HLA_REV 1
#define DO_HLA_BOTH 2
#define DO_ALL_HLAS -1

// Window move definitions.
#define DO_MERGE_OR_SPLIT 0
#define DO_MERGE 1
#define DO_SPLIT 2
#define DO_GROW 3
#define DO_COEFF 4

int rec_accept, rec_reject;
int sel_accept, sel_reject;
int merge_accept, merge_reject;
int split_accept, split_reject;
int grow_accept, grow_reject;
int wcoeff_accept, wcoeff_reject;
int mut_accept, mut_reject;

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
                               int do_esc_or_rev);

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
                        char cons_query_nucls[num_query][num_bases]);

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
                    int num_bases, char cons_query_nucls[num_query][num_bases]);

Decimal evaluate_prior_diff_merge(int HLA_coeff_prior,
                                  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal mu_sel_norm_prior, Decimal sigma_sel_norm_prior,
                                  Decimal gamma_k_shape, Decimal gamma_theta_scale,
                                  Decimal old_coeff0, Decimal old_coeff1,
                                  Decimal new_coeff, Decimal GAMMA_CONST);

Decimal evaluate_prior_diff_split(int HLA_coeff_prior,
                                  Decimal prior_add_window,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal mu_sel_norm_prior, Decimal sigma_sel_norm_prior,
                                  Decimal gamma_k_shape, Decimal gamma_theta_scale,
                                  Decimal old_coeff0,
                                  Decimal new_coeff1, Decimal new_coeff2,
                                  Decimal GAMMA_CONST);

Decimal evaluate_prior_diff_coeff(int HLA_coeff_prior,
                                  Decimal mu_sel_prior, Decimal sigma_sel_prior,
                                  Decimal mu_sel_norm_prior, Decimal sigma_sel_norm_prior,
                                  Decimal gamma_k_shape, Decimal gamma_theta_scale,
                                  Decimal old_coeff0, Decimal new_coeff);

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
                 Decimal GAMMA_CONST);

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
                   int num_bases, char cons_query_nucls[num_query][num_bases]);

#endif /* MCMC_MOVES_H_ */
