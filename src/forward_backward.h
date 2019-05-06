#ifndef FORWARD_BACKWARD_H_
#define FORWARD_BACKWARD_H_

#define USE_TRIPLET 1

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
                                    char cons_query_nucls[num_bases]);

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
                                char cons_query_nucls[num_bases]);

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
                                char cons_query_nucls[num_bases]);

#endif /* FORWARD_BACKWARDS_H_ */
