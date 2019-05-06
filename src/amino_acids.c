#include "amino_acids.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>

// Exported.
const char alphabet[] = "TCAG";

// T->0, C->1, A->2, G->3, N->4.
const char dnachar_to_code[128] = { 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 2, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 2, 4, 1, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
                                    4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

// ttt ... ggg.
const char rna_to_codon[NUM_CODONS] =
{'F','F','L','L','S','S','S','S',
 'Y','Y','*','*','C','C','*','W',
 'L','L','L','L','P','P','P','P',
 'H','H','Q','Q','R','R','R','R',
 'I','I','I','M','T','T','T','T',
 'N','N','K','K','S','S','R','R',
 'V','V','V','V','A','A','A','A',
 'D','D','E','E','G','G','G','G'};

const char codon_to_dna[NUM_CODONS][4] =
{"TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", 
 "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", 
 "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", 
 "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", 
 "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", 
 "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", 
 "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", 
 "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"};

static const char triplet_codons[125] =
{   0,  1,  2,  3, -1,  4,  5,  6,  7, -1,  8,  9, 10, 11, -1, 12,
   13, 14, 15, -1, -1, -1, -1, -1, -1, 16, 17, 18, 19, -1, 20, 21,
   22, 23, -1, 24, 25, 26, 27, -1, 28, 29, 30, 31, -1, -1, -1, -1,
   -1, -1, 32, 33, 34, 35, -1, 36, 37, 38, 39, -1, 40, 41, 42, 43,
   -1, 44, 45, 46, 47, -1, -1, -1, -1, -1, -1, 48, 49, 50, 51, -1,
   52, 53, 54, 55, -1, 56, 57, 58, 59, -1, 60, 61, 62, 63, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
   -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};

char triplet_to_codon(char triplet, const char consensus[3])
{
  int codon = triplet_codons[(int)triplet];

  if(codon == -1)
  {
    int a = triplet%5, b = (triplet/5)%5, c = (triplet/25)%5;

    if(a == 4) a = consensus[2];
    if(b == 4) b = consensus[1];
    if(c == 4) c = consensus[0];

    codon = c*16+b*4+a;
  }

  return codon;
}

void generate_codons(char codons[NUM_CODONS][4])
{
  int x = 0, i, j, k;
  for(i = 0; i < ALPHABET_SIZE; i++) {
    for(j = 0; j < ALPHABET_SIZE; j++) {
      for(k = 0; k < ALPHABET_SIZE; k++) {
        codons[x][0] = alphabet[i];
        codons[x][1] = alphabet[j];
        codons[x][2] = alphabet[k];
        codons[x][3] = '\0';
        x++;
      }
    }
  }
}

static inline void increment_base_count(int counts[4], char c)
{
  switch(c) {
    case 'T': case 't': counts[0]++; break;
    case 'C': case 'c': counts[1]++; break;
    case 'A': case 'a': counts[2]++; break;
    case 'G': case 'g': counts[3]++; break;
  }
}

void generate_consensus(char *const* seqs, int num_seqs, int len, char *consensus)
{
  // Counts of T, C, A, G at a given position.
  int base_counts[4];
  int i, s;
  // Loop over positions.
  for(i = 0; i < len; i++)
  {
    memset(base_counts, 0, sizeof(int) * 4);

    for(s = 0; s < num_seqs; s++)
      increment_base_count(base_counts, seqs[s][i]);
    
    char con_base = 'T';
    int count = base_counts[0];
    if(base_counts[1] > count) { con_base = 'C'; count = base_counts[1];}
    if(base_counts[2] > count) { con_base = 'A'; count = base_counts[2];}
    if(base_counts[3] > count) { con_base = 'G'; count = base_counts[3];}
    if(count == 0) { printf("Error: consensus is a sequence ambiguity at site %i\n", i);}
    consensus[i] = con_base;
  }
  consensus[len] = '\0';
}

void generate_consensus_subset(char *const* seqs, int closest_n, int len, char *consensus, 
                               int closest_refs[closest_n])
{
  // Counts of T, C, A, G at a given position.
  int base_counts[4];
  int i, s;
  // Loop over positions.
  for(i = 0; i < len; i++)
  {
    memset(base_counts, 0, sizeof(int) * 4);

    for(s = 0; s < closest_n; s++)
      increment_base_count(base_counts, seqs[closest_refs[s]][i]);

    char con_base = 'T';
    int count = base_counts[0];
    if(base_counts[1] > count) { con_base = 'C'; count = base_counts[1];}
    if(base_counts[2] > count) { con_base = 'A'; count = base_counts[2];}
    if(base_counts[3] > count) { con_base = 'G'; count = base_counts[3];}
    if(count == 0) { printf("Error: consensus for query is a sequence ambiguity at site %i\n", i);}
    consensus[i] = con_base;
  }
  consensus[len] = '\0';
}

void map_gaps_to_consensus(char **seqs, int num_seqs, int len,
                           const char *consensus)
{
  int i, s;
  for(i = 0; i < len; i++) {
    for(s = 0; s < num_seqs; s++) {
      if(!is_base(seqs[s][i]))
        seqs[s][i] = consensus[i];
    }
  }
}

void list_ref_codons_at_all_sites(int num_sites, int num_refs,
                                  int num_ref_codons[num_sites],
                                  char ref_site_codons[num_sites][NUM_CODONS],
                                  char ref_codons[num_refs][num_sites])
{
  int s, r, c, i;
  char codons_seen[NUM_CODONS];

  for(s = 0; s < num_sites; s++)
  {
    memset(codons_seen, 0, sizeof(char) * NUM_CODONS);

    for(r = 0; r < num_refs; r++)
      codons_seen[(int)ref_codons[r][s]] = 1;

    for(c = 0, i = 0; c < NUM_CODONS; c++)
      if(codons_seen[c] == 1)
        ref_site_codons[s][i++] = c;

    num_ref_codons[s] = i;
  }
}
