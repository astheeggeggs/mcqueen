#ifndef CONSTS_H_
#define CONSTS_H_

#include <stdio.h>
#include <stdbool.h>

#include "amino_acids.h"

// Define the accuracy of the type "Decimal";

// typedef float Decimal;
// #define DECPRINT "%f"
// #define DECPRINTJSON "%.10f"
// #define DECPRINTTIME "%.2f"

typedef double Decimal;
#define DECPRINT "%f"
#define DECPRINTJSON "%.10f"
#define DECPRINTTIME "%.2f"

// typedef long double Decimal;
// #define DECPRINT "%Lf"
// #define DECPRINTJSON "%.10Lf"
// #define DECPRINTTIME "%.2Lf"

#define CON_TO_CON 1
#define CON_TO_NONCON 2
#define NONCON_TO_CON 3
#define NONCON_TO_NONCON 4

#define LOG_NORMAL_PRIOR 0
#define NORMAL_PRIOR 1
#define GAMMA_PRIOR 2
#define FLAT_PRIOR 3
#define EXPONENTIAL_PRIOR 4

// Extern constants.

extern int HLA_coeff_prior_esc;
extern int HLA_coeff_prior_rev;
extern int omega_prior;
extern bool sample_prior;

// Large numbers for normalisation when calculating F.
extern const Decimal BIG;
extern const Decimal BIGI;

// Initial mutation rate, mu.
extern const Decimal init_mutation_rate;

// Initial recombination probability between all sites, R[i].
extern const Decimal init_R;

// Parameter of exponential prior on rate of recombination.
extern const Decimal p_rec;
extern Decimal lambda_rec;

// Parameter of exponential prior for selection coeff.
extern const Decimal lambda_sel;

// Parameter of log normal prior for selection coeff.
extern const Decimal omega_mu_sel;
// -> Std dev.
extern const Decimal omega_sigma_sel;

// Parameter of gamma prior for selection coeff.
// -> Shape parameter
extern const Decimal omega_gamma_k_shape;
// -> Scale parameter
extern const Decimal omega_gamma_theta_scale;

// Initial HLA escape selection parameter across all sites.
extern const Decimal init_HLA_sel;

// Initial HLA reversion selection parameter across all sites.
extern const Decimal init_HLA_rev;

// Prior for all selection parameters (log-normal).
// -> Mean.
extern const Decimal mu_sel_prior_esc;
extern const Decimal mu_sel_prior_rev;
// -> Std dev.
extern const Decimal sigma_sel_prior_esc;
extern const Decimal sigma_sel_prior_rev;

// Prior for all HLA selection parameters - if choosing truncated normal
// -> Mean.
extern const Decimal mu_sel_norm_prior_esc;
extern const Decimal mu_sel_norm_prior_rev;
// -> Std dev.
extern const Decimal sigma_sel_norm_prior_esc;
extern const Decimal sigma_sel_norm_prior_rev;

// Prior for all HLA selection parameters - if choosing a gamma.
// -> Shape parameter
extern const Decimal gamma_k_shape_esc;
extern const Decimal gamma_k_shape_rev;
// -> Scale parameter
extern const Decimal gamma_theta_scale_esc;
extern const Decimal gamma_theta_scale_rev;

// Selection windows per HLA type.
extern const int init_num_selec_windows_esc;
extern const int init_num_selec_windows_rev;

// Parameter of prior on selection window expansion.
extern const Decimal selection_window_geom_param;

// Success probability in a binomial trial for adding a window.
extern const Decimal prior_add_window_esc;
extern const Decimal prior_add_window_rev;

// Parameter of exponential prior for mu.
extern const Decimal mu_prior;

extern const int num_codons;

// Synonymous transitions.
extern const int beta_S[NUM_CODONS];

// Synonymous transverions.
extern const int beta_V[NUM_CODONS];

// Non-synonymous transitions.
extern const int alpha_S[NUM_CODONS];

// Non-synonymous transversions.
extern const int alpha_V[NUM_CODONS];

// Transition tranversion ratio.
extern Decimal kappa;

// Probabilility of 2 or more base pair changes in a codon.
extern Decimal phi;

// Prior for C1 - the prior probability of observing codon C1 - 
// just the prevalence from a large African dataset.
extern Decimal prior_C1[NUM_CODONS];
extern Decimal *prior_C2;

// Matrix of non-synonymous transitions (1 if NS_TS, 0 if not).
extern const int NS_TS[NUM_CODONS][NUM_CODONS];

// Matrix of non-synonymous transversions (1 if NS_TV, 0 if not).
extern const int NS_TV[NUM_CODONS][NUM_CODONS];

// Matrix of synonymous transitions (1 if S_TS, 0 if not).
extern const int S_TS[NUM_CODONS][NUM_CODONS];

// Matrix of non-synonymous transversions (1 if S_TV, 0 if not).
extern const int S_TV[NUM_CODONS][NUM_CODONS];

// Matrix of two (or more) step transitions (1 if two_step, 0 if not).
extern const int two_step[NUM_CODONS][NUM_CODONS];

typedef struct
{
  int index_of_sequence;
  Decimal likelihood_of_sequence;
} index_and_likelihood;

#endif /* CONSTS_H_ */
