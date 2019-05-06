#ifndef LOAD_FROM_JSON_H_
#define LOAD_FROM_JSON_H_

#include "constants.h"
#include "selection_window.h"
#include "write_json.h"

void load_from_json(const char *path, Decimal *llk, Decimal *mu,
                    Decimal *R, Decimal *omega,
                    HLASelection *hla_selection,
                    int num_sites, int num_hla_types,
                    int num_query, int closest_n, int closest_refs[num_query][closest_n]);

void load_closest_n_from_json(const char *path, int num_sims, int num_sites,
                              Mosaic *true_copy);

void load_fixed_parameters_from_json(const char *path, Decimal *kappa_new, Decimal *phi_new,
                                     Decimal *codon_prior_new);

void load_lengths_for_simulation_from_json(const char *path, Decimal *kappa, Decimal *mu,
                                           int *codon_sequence_length, int *total_n_HLA,
                                           int *ploidy, int *n_genes);

void load_parameters_for_simulation_from_json(const char *path, int codon_sequence_length,
                                              Decimal omega[codon_sequence_length], Decimal R[codon_sequence_length],
                                              Decimal reversion_selection[codon_sequence_length],
                                              int total_n_HLA, int n_genes, int n_HLA[n_genes],
                                              Decimal HLA_prevalences[total_n_HLA],
                                              Decimal HLA_selection_profiles[total_n_HLA][codon_sequence_length]);

#endif /* LOAD_FROM_JSON_H_ */
