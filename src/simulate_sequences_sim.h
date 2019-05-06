#ifndef SIMULATE_SEQUENCES_H_
#define SIMULATE_SEQUENCES_H_

extern FILE *simulated_file;
extern FILE *hla_file;
extern FILE *mosaic_json;

void simulate_codon_to_codon_HLA(int num_sites, int num_refs,
                                 int num_sims, int num_hla_types,
                                 Decimal (*codon_to_codon_HLA)[num_sites][NUM_CODONS],
                                 Decimal (*rate_out_of_codons)[num_sites],
                                 int start_site, int end_site,
                                 Decimal mu,
                                 Decimal (*HLA_select_esc)[num_sites],
                                 Decimal (*HLA_select_rev)[num_sites],
                                 char (*hla_types)[num_hla_types],
                                 const char consensus_codons[num_sites],
                                 const Decimal omega[num_sites],
                                 const int NS_TS_sum_from_consensus[num_sites],
                                 const int NS_TV_sum_from_consensus[num_sites],
                                 char (*const type_of_change)[NUM_CODONS][NUM_CODONS],
                                 int (*codon_copied_from)[num_sites],
                                 Decimal (*A)[NUM_CODONS][num_sites]);

int bfind_pos(Decimal to_sort[NUM_CODONS]);

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
                        int (*copy_from)[num_sites]);

void save_simulated_fasta(int num_sims, int num_sites, int (*simulated_sequences)[num_sites]);

#endif /* SIMULATE_SEQUENCES_H_ */