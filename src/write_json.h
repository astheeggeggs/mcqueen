#ifndef WRITE_JSON_H_
#define WRITE_JSON_H_

#include "constants.h"
#include "selection_window.h"

typedef struct
{
  int *which_seqs;
  int num_seqs;
} Mosaic;

void write_json(FILE *final_state_json,
	            int chain_i, Decimal curr_llk_across_queries, Decimal mu, int num_sites,
	            Decimal omega[num_sites], Decimal R[num_sites], HLASelection *hla_windows,
	            int num_hla_types,
	            int num_query, int closest_n, int closest_refs[num_query][closest_n]);

void write_summary_json(FILE *summary_sim_json,
                        Decimal mu, int num_sites, int ploidy,
                        int n_genes, int num_HLA[n_genes],
                        int total_n_HLA,
                        Decimal HLA_prevalences[total_n_HLA],
                        Decimal omega[num_sites], Decimal R[num_sites],
                        Decimal reversion_selection[num_sites],
                        Decimal HLA_selection_profiles[total_n_HLA][num_sites]);

void write_closest_n_json(FILE *closest_n_json, Mosaic *true_copy, int num_sims);

#endif /* WRITE_JSON_H_ */