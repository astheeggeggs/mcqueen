#ifndef SIMULATE_SEQS_GEN_
#define SIMULATE_SEQS_GEN_

#include "amino_acids.h"

extern FILE *simulated_refs_file;
extern FILE *simulated_root_file;
extern FILE *simulated_queries_file;
extern FILE *hla_query_file;

int create_new_codon(int current_codon, int wildtype_codon, 
                     const int row_sum_type_of_change,
                     const int type_of_change_matrix[NUM_CODONS][NUM_CODONS]);

void codon_sequence_change(int codon_sequence_length, 
                           Decimal S_portion[NUM_CODONS], Decimal NS_portion[NUM_CODONS],
                           Decimal max_time, Decimal mu,
                           Decimal omega[codon_sequence_length],
                           int wildtype_sequence[codon_sequence_length],
                           int (*codon_sequence_matrix)[codon_sequence_length], 
                           int from, int to,
                           Decimal reversion_selection[codon_sequence_length],
                           Decimal HLA_selection_profile[codon_sequence_length]);

void pass_HLA(int ploidy, int n_genes, int starting_node,
              Tree *tree, int number_of_leaves, int total_n_HLA,
              int n_HLA[n_genes], Decimal HLA_prevalences[total_n_HLA]);

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
                                Decimal reversion_selection[codon_sequence_length]);

void save_simulated_ref_and_query_fasta(int num_queries, int num_refs, int num_sims, 
                                        int all_sequences[num_sims], int codon_sequence_length, 
                                        int (*simulated_sequences)[codon_sequence_length],
                                        Tree *tree, int ploidy, int n_genes);

#endif /* SIMULATE_SEQS_GEN_ */
