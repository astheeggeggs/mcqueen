// request decent POSIX version.
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h> // seed random.
#include <string.h>
#include <strings.h> // strcasecmp.
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <limits.h>
#include <assert.h>
#include <getopt.h>

#include "constants.h"
#include "util.h"
#include "load_data.h"
#include "amino_acids.h"
#include "selection_window.h"
#include "load_from_json.h"
#include "write_json.h"
#include "simulate_sequences_sim.h"

#define DEFAULT_NUM_SIMS 500
#define DEFAULT_PROPORTION_NO_HLA 0.2

FILE *summary_json_file = NULL;

int main(int argc, char **argv)
{
  init_rand();
  setup_gsl();

  int ploidy, n_genes, num_hla_types, num_sims = DEFAULT_NUM_SIMS;
  Decimal mu, proportion_no_hla = DEFAULT_PROPORTION_NO_HLA;
  char *cons_option_path = NULL;

  static struct option longopts[] =
  {
    {"n_sims",            required_argument, NULL, 's'},
    {"consensus_fasta",   required_argument, NULL, 'c'},
    {"prop_no_HLA",       required_argument, NULL, 'n'},
    {NULL, 0, NULL, 0}
  };

  char shortopts[] = "s:c:n:p";

  int optional;
  while((optional = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch (optional) {
      case 's': num_sims = atoi(optarg); break;
      case 'c': cons_option_path = optarg; break;
      case 'n': proportion_no_hla = atof(optarg); break;
      default: die("Unknown option: %c", optional);
    }
  }

  if(argc - optind != 3) {
    fprintf(stderr,
            "usage: dsim [options] <parameters.json> <ref_seqs.fasta> <outdir>\n"
            "  options:\n"
            "    --n_sims            -s <S> Set number of simulated sequences [default: %i].\n"
            "    --consensus_fasta   -c <C> Set a consensus sequence [default: consensus of ref_seqs.fasta].\n"
            "    --prop_no_HLA       -n <N> Set a proportion of the query sequences to not have HLA information "
                                "[default: "DECPRINT"].\n",
            (int) DEFAULT_NUM_SIMS, (Decimal) DEFAULT_PROPORTION_NO_HLA);
    exit(EXIT_FAILURE);
  }
  
  const char *json_passed = argv[optind];
  const char *ref_path = argv[optind+1];
  const char *out_dir = argv[optind+2];
  // Set the sequence file from which to define consensus as the reference set (as default).
  char *cons_path = argv[optind+1];
  // If an optional argument is passed, then set the consensus filepath as this instead.
  if(cons_option_path != NULL)
    cons_path = cons_option_path;

  if(!mkpath(out_dir, 0755)) die("Cannot create output dir: %s", out_dir);

  time_t rawtime;
  struct tm * timeinfo;

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  char date[20];
  sprintf(date, "%d_%d_%d_%d_%d_%d", timeinfo->tm_year+1900, 
          timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
          timeinfo->tm_min, timeinfo->tm_sec);

  char simulated_sequences_filename[PATH_MAX+1];
  char summary_json_filename[PATH_MAX+1];
  char hla_filename[PATH_MAX+1];
  char mosaic_json_filename[PATH_MAX+1];

  sprintf(simulated_sequences_filename, "%s/%s_simulated_sequences.fasta", out_dir, date);  
  sprintf(summary_json_filename, "%s/%s_summary_dsim.json", out_dir, date);
  sprintf(hla_filename, "%s/%s_simulated_hla.csv", out_dir, date);
  sprintf(mosaic_json_filename, "%s/%s_mosaic.json", out_dir, date);
  
  printf("Writing simulated sequences to: %s.\n", simulated_sequences_filename);
  printf("Writing summary .json to: %s.\n", summary_json_filename);
  printf("Writing simulated HLAs to: %s.\n", hla_filename);
  printf("Writing mosaic of sequences copied from to: %s.\n", mosaic_json_filename);

  simulated_file = fopen(simulated_sequences_filename, "w");  
  hla_file = fopen(hla_filename, "w");  
  summary_json_file = fopen(summary_json_filename, "w"); 
  mosaic_json = fopen(mosaic_json_filename, "w");  

  if(simulated_file == NULL)
    die("Cannot open output file: simulated_sequences.fasta.");
  if(hla_file == NULL)
    die("Cannot open output file: simulated_hla.csv.");
  if(summary_json_file == NULL)
    die("Cannot open output file: summary_dsim.json.");
  if(mosaic_json == NULL)
    die("Cannot open output file: mosaic.json.");

  char **ref_seqs, **cons_seqs;
  int refs_cap, cons_cap;

  printf("Loading sequences...\n");
  int num_refs = load_seqs(ref_path, &ref_seqs, &refs_cap);
  int num_cons = load_seqs(cons_path, &cons_seqs, &cons_cap);

  printf("Loaded %i reference sequences.\n", num_refs);
  printf("Loaded %i sequences to determine consensus.\n", num_cons);

  // Check all sequences have the same length.
  size_t len = strlen(ref_seqs[0]);

  if(len % 3 != 0) die("Sequences contain partial codons [length mod 3 != 0].");

  int r, s, i;
  int c1, c2, c;

  for(r = 1; r < num_refs; r++)
    if(strlen(ref_seqs[r]) != len)
      die("Ref sequences aren't all the same length.");
  
  for(c = 0; c < num_cons; c++) {
    if(strlen(cons_seqs[c]) != len) {
      die("Sequences from which to derive the consensus sequence aren't all "
          "the same length.");
    }
  }

  printf("Mapping gaps to consensus...\n");
  char consensus[len+1];

  generate_consensus(cons_seqs, num_cons, len, consensus);
  
  printf("Consensus:\n%s\n", consensus);
  map_gaps_to_consensus(ref_seqs, num_refs, len, consensus);

  int num_sites = len / 3;
  assert(num_sites > 0);

  char (*ref_codons)[num_sites] = my_malloc(num_refs * sizeof(char[num_sites]),
                                            __FILE__,__LINE__);

  printf("Number of references: %i.\n", num_refs);

  for(r = 0; r < num_refs; r++)
    for(s = 0; s < num_sites; s++)
      ref_codons[r][s] = amino_to_code(ref_seqs[r]+s*3);

  // Free sequence data - only work with codons from now on.
  for(r = 0; r < num_refs; r++) free(ref_seqs[r]);
  free(ref_seqs);

  for(c = 0; c < num_cons; c++) free(cons_seqs[c]);
  free(cons_seqs);

  int *num_ref_site_codons = my_malloc(num_sites * sizeof(int),__FILE__,__LINE__);
  char (*ref_site_codons)[NUM_CODONS] = my_malloc(num_sites * sizeof(char[NUM_CODONS]),
                                                  __FILE__,__LINE__);
  list_ref_codons_at_all_sites(num_sites, num_refs,
                               num_ref_site_codons, ref_site_codons,
                               ref_codons);

  printf("Reference codons: ");
  for(s = 0; s < num_sites; s++) {
    if(s > 0) printf(" | ");
    printf("%i", (int) ref_site_codons[s][0]);
    for(i = 1; i < num_ref_site_codons[s]; i++)
      printf(",%i", (int) ref_site_codons[s][i]);
  }
  printf("\n");

  if(json_passed == NULL) die("No .json file passed.");

  int num_sites_copy = num_sites;

  load_lengths_for_simulation_from_json(json_passed, &kappa, &mu,
                                        &num_sites, &num_hla_types, &ploidy, &n_genes);

  if(num_sites != num_sites_copy) die("Mismatch between the length of passed parameters and "
                                       "the length of the reference sequences: [%i != %i]." ,
                                       num_sites_copy, num_sites);
  
  printf("mu "DECPRINT".\n", mu);
  printf("kappa "DECPRINT".\n", kappa);
  // HLA prevalences, and the number of HLA types in each gene under consideration.  
  Decimal HLA_prevalences[num_hla_types];
  int num_HLA[n_genes];

  // Selection coeff in absence of HLA, per site, and recombination probabilities.
  Decimal *omega = my_malloc(num_sites * sizeof(Decimal),__FILE__,__LINE__);
  Decimal *R = my_malloc(num_sites * sizeof(Decimal),__FILE__,__LINE__);
  
  // HLA_select_esc - HLA selection coefficients associated to escape.
  // HLA_select_rev - coefficients associated to reversion.
  Decimal (*HLA_select_esc)[num_sites] = my_malloc(num_hla_types * sizeof(Decimal[num_sites]),
                                                   __FILE__,__LINE__);
  Decimal *HLA_select_rev = my_malloc(num_sites * sizeof(Decimal),__FILE__,__LINE__);

  load_parameters_for_simulation_from_json(json_passed, num_sites, omega, R, HLA_select_rev,
                                           num_hla_types, n_genes, num_HLA,
                                           HLA_prevalences, HLA_select_esc);

  printf("omega:\n");
  for (i = 0; i < num_sites; i++) printf(DECPRINT" ", omega[i]);
  printf("\n");

  printf("R:\n");
  for (i=0; i<num_sites; i++) printf(DECPRINT" ", R[i]);
  printf("\n");

  printf("HLA_select_rev:\n");
  for (i=0; i<num_sites; i++) printf(DECPRINT" ", HLA_select_rev[i]);
  printf("\n");

  printf("n_genes: %i.\n", n_genes);
  printf("num_HLA: %i.\n", num_hla_types);

  // Set 'lambda'.
  Decimal (*rate_out_of_codons)[num_sites] = my_malloc(NUM_CODONS * sizeof(Decimal[num_sites]),
                                                       __FILE__,__LINE__);

  for(c1 = 0; c1 < NUM_CODONS; c1++)
  {
    Decimal v1 = kappa * beta_S[c1] + beta_V[c1];
    Decimal v2 = kappa * alpha_S[c1] + alpha_V[c1];

    for(s = 0; s < num_sites; s++) {
      rate_out_of_codons[c1][s] = v1 + omega[s] * v2;
    }
  }

  // Determine the consensus codon sequence.
  char *consensus_codons = my_malloc(num_sites * sizeof(char),__FILE__,__LINE__);
  printf("Consensus codons:\n");
  for(s = 0; s < num_sites; s++) {
    consensus_codons[s] = amino_to_code(consensus+s*3);
    printf("%4i ", consensus_codons[s]);
  }
  printf("\n");

  // Print the first reference sequence for comparison.
  printf("Reference codons 1:\n");
  for(s = 0; s < num_sites; s++) {
    printf("%4i ", ref_codons[1][s]);
  }
  printf("\n");

  // Number of possible non-synonymous transitions and transversions
  // from the consensus codon at each site. 
  int NS_TS_sum_from_consensus[num_sites];
  int NS_TV_sum_from_consensus[num_sites];

  for(s = 0; s < num_sites; s++) {
    NS_TS_sum_from_consensus[s] = 0;
    NS_TV_sum_from_consensus[s] = 0;
    for(c1 = 0; c1 < NUM_CODONS; c1++) {
      NS_TS_sum_from_consensus[s] += NS_TS[(int)consensus_codons[s]][c1];
      NS_TV_sum_from_consensus[s] += NS_TV[(int)consensus_codons[s]][c1];
    }
  }

  char (*type_of_change)[NUM_CODONS][NUM_CODONS]
    = my_malloc(num_sites * sizeof(char[NUM_CODONS][NUM_CODONS]),__FILE__,__LINE__);

  // Determine which sites in the queries are consensus to consensus,
  // consensus to non-consensus, non-consensus to consensus and non-consensus
  // to non-consensus.
  for(s = 0; s < num_sites; s++) {
    int consensus_codon = consensus_codons[s];
    for(c1 = 0; c1 < NUM_CODONS; c1++) {
      for(c2 = 0; c2 < NUM_CODONS; c2++)
        type_of_change[s][c1][c2] = NONCON_TO_NONCON;

      type_of_change[s][c1][consensus_codon] = NONCON_TO_CON;
    }
  }

  for(s = 0; s < num_sites; s++) {
    int consensus_codon = consensus_codons[s];
    for(c1 = 0; c1 < NUM_CODONS; c1++)
      type_of_change[s][consensus_codon][c1] = CON_TO_NONCON;

    type_of_change[s][consensus_codon][consensus_codon] = CON_TO_CON;
  }

  // Start of simulation:
  Decimal (*sim_codon_to_codon_HLA)[num_sites][NUM_CODONS]
    = my_malloc(num_sims * sizeof(Decimal[num_sites][NUM_CODONS]),__FILE__,__LINE__);
  Decimal (*A)[NUM_CODONS][num_sites]
    = my_malloc(num_sims * sizeof(Decimal[NUM_CODONS][num_sites]),__FILE__,__LINE__);

  int (*codon_copied_from)[num_sites] = my_malloc(num_sims * sizeof(int[num_sites]),
                                                  __FILE__,__LINE__);
  int (*copy_from)[num_sites] = my_malloc(num_sims * sizeof(int[num_sites]),
                                          __FILE__,__LINE__);
  int (*simulated_sequences)[num_sites] = my_malloc(num_sims * sizeof(int[num_sites]),
                                                    __FILE__,__LINE__ );

  // Now write this information to a .json and then free the variables.
  printf("Writing the parameters to a summary .json file.\n");

  write_summary_json(summary_json_file, mu, num_sites, ploidy,
                     n_genes, num_HLA, num_hla_types, HLA_prevalences,
                     omega, R, HLA_select_rev, HLA_select_esc);

  simulate_sequences(num_sites, num_refs, num_sims,
                     num_hla_types, HLA_prevalences,
                     ploidy, n_genes, num_HLA,
                     ref_codons,  mu,
                     HLA_select_esc, HLA_select_rev,
                     consensus_codons, omega, R, 
                     rate_out_of_codons,
                     NS_TS_sum_from_consensus,
                     NS_TV_sum_from_consensus,
                     type_of_change,
                     sim_codon_to_codon_HLA,
                     proportion_no_hla,
                     codon_copied_from,
                     simulated_sequences,
                     A, copy_from);
  
  save_simulated_fasta(num_sims, num_sites, simulated_sequences);

  printf("Simulated %i sequences from reference panel of size %i\n.", num_sims, num_refs);

  free(rate_out_of_codons);
  free(ref_codons);
  free(num_ref_site_codons);
  free(ref_site_codons);
  free(consensus_codons);
  free(type_of_change);
  free(sim_codon_to_codon_HLA);
  free(HLA_select_esc);
  free(HLA_select_rev);
  free(R);
  free(omega);
  free(codon_copied_from);
  free(copy_from);
  free(simulated_sequences);
  free(A);

  // Close the simulated sequence file
  fclose(simulated_file);
  fclose(hla_file);
  fclose(summary_json_file);
  fclose(mosaic_json);

  clearup_gsl();
  return EXIT_SUCCESS;
}
