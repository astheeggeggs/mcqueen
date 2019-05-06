#ifndef AMINO_ACIDS_H_
#define AMINO_ACIDS_H_

#define ALPHABET_SIZE 4
#define TRIPLET_ALPHABET_SIZE 5

extern const char alphabet[ALPHABET_SIZE];

#define NUM_CODONS (ALPHABET_SIZE*ALPHABET_SIZE*ALPHABET_SIZE)
#define NUM_TRIPLETS (TRIPLET_ALPHABET_SIZE*TRIPLET_ALPHABET_SIZE*TRIPLET_ALPHABET_SIZE)

#define is_base(x) (dnachar_to_code[(int)(x)] != 4)

extern const char dnachar_to_code[128];

extern const char rna_to_codon[NUM_CODONS];

extern const char codon_to_dna[NUM_CODONS][4];

char triplet_to_codon(char triplet, const char consensus[3]);

#define base_to_code(x) dnachar_to_code[(int)x]

#define amino_to_code(codon) \
        (base_to_code((codon)[0])*ALPHABET_SIZE*ALPHABET_SIZE + \
         base_to_code((codon)[1])*ALPHABET_SIZE + base_to_code((codon)[2]))

#define triplet_to_code(triplet) \
        (base_to_code((triplet)[0])*TRIPLET_ALPHABET_SIZE*TRIPLET_ALPHABET_SIZE + \
         base_to_code((triplet)[1])*TRIPLET_ALPHABET_SIZE + base_to_code((triplet)[2]))

#define code_to_char(x) codon_to_dna[(int)x]
#define code_to_amino(x) rna_to_codon[(int)x]

void generate_codons(char codons[NUM_CODONS][4]);

// char triplet_to_codon(char triplet, const char consensus[3]);

void generate_consensus(char *const* seqs, int num_seqs, int len, char *consensus);

void generate_consensus_subset(char *const* seqs, int closest_n, int len, char *consensus, 
                               int closest_refs[closest_n]);

void map_gaps_to_consensus(char **seqs, int num_seqs, int len,
                           const char *consensus);

void list_ref_codons_at_all_sites(int num_sites, int num_refs,
                                  int num_ref_codons[num_sites],
                                  char ref_site_codons[num_sites][NUM_CODONS],
                                  char ref_codons[num_refs][num_sites]);

#endif /* AMINO_ACIDS_H_ */
