// Request decent POSIX version.
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h> // Seed random.
#include <string.h>
#include <strings.h> // strcasecmp.
#include <math.h>
#include <unistd.h>
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
#include "closest_seqs.h"
#include "forward_backward.h"
#include "mcmc_moves.h"

#define DEFAULT_CHAIN_LEN 50000
#define DEFAULT_SAMPLE_LEN 100
#define DEFAULT_CLOSEST_N 100

// Files opened by dmodel_parse_cmdline()
FILE *log_file = NULL;
FILE *window_file = NULL;
FILE *priors_file = NULL;
FILE *summary_file = NULL;
FILE *final_state_json = NULL;
FILE *interim_final_state_json = NULL;

FILE *hamming_file = NULL;
FILE *hamming_distance_file = NULL;

// Variables set by the command line by dmodel_parse_cmdline()
char interim_json_path[PATH_MAX+1];

int chain_length = DEFAULT_CHAIN_LEN, sample_rate = DEFAULT_SAMPLE_LEN;
int runtime_limit_hrs = 0, closest_n = DEFAULT_CLOSEST_N;
const char *ref_path = NULL, *query_path = NULL, *hla_path = NULL;
const char *cons_path, *out_dir = NULL;

char *cons_option_path = NULL, *resume_state_path = NULL;
char *json_parameters_path = NULL, *true_mosaic_sequences = NULL;
bool runtime_limit_hrs_set = false;
bool number_of_states_set = false;
bool separate_reference_set = true;
bool no_HLA_inference = false;
bool no_rev_inference = false;
bool no_recombination = false;
bool no_simulated_annealing = false;
// Possible optimisations:
// 1) accept/reject window move before moving window.
// 2) re-order matrix dimensions.
// 3) tweaking values of BIG and BIGI.

// Search 'DEV' for other issues.
// when loading a json file, check the closest_n too. 
// think about making sure that the furthest away of the closest_n are the same sequences.

static void print_usage()
{
  fprintf(stderr,
          "usage: dmodel [options] <ref_seqs.fasta> <query.fasta> <outdir>\n"
          "  options:\n"
          "    --HLA_inference              -H <HLA.csv> Pass a .csv file with HLA information.\n"
          "    --chain_length               -n <N> Set chain length [default: %i].\n"
          "    --sample_every               -s <S> Sample every <S> steps [default: %i].\n"
          "    --runtime_hours              -t <T> Run for T hours.\n"
          "    --consensus_fasta            -c <C> Consensus set as consensus of this .fasta file [default: <ref_seqs.fasta>].\n"
          "    --resume_json                -r <state.json> Resume state.\n"
          "    --closest_n                  -l <L> Closest L sequences [default: %i].\n"
          "    --parameters_json            -p <parameters.json> a .json file containing the codon usage, transition transversion ratio and\n"
          "                                    proportion of codon changes which involve at least two nucleotide substitutions.\n"
          "    --true_mosaic_fasta          -m <mosaic.json> True mosaic of sequences from simulation.\n"
          "    --separate_reference_fasta   -q Separate reference set [default: %s].\n"
          "    --no_HLA_inference           -i No HLA inference [default: %s].\n"
          "    --no_rev_inference           -j No reversion inference [default: %s].\n"
          "    --no_recombination_inference -o No recombination [default: %s].\n"
          "    --omega_gamma                -v Gamma prior on omega coefficients [default: exponential prior].\n"
          "    --omega_log_normal           -e Log Normal prior on omega coefficients [default: exponential prior].\n"
          "    --omega_uniform              -w Improper uniform prior on omega coefficients [default: exponential prior].\n"
          "    --HLA_coeff_esc_gamma        -g Gamma prior on HLA coefficients [default: log normal prior].\n"
          "    --HLA_coeff_esc_normal       -k Normal prior on HLA coefficients [default: log normal prior].\n"
          "    --HLA_coeff_esc_uniform      -f Flat prior on HLA coefficients [default: log normal prior].\n"
          "    --coeff_rev_gamma            -a Gamma prior on HLA coefficients [default: log normal prior].\n"
          "    --coeff_rev_normal           -d Normal prior on HLA coefficients [default: log normal prior].\n"
          "    --coeff_rev_uniform          -u Flat prior on HLA coefficients [default: log normal prior].\n"
          "    --sample_prior               -b Sample all parameters from the specified priors.\n"
          "    --no_simulated_annealing     -z Don't perform simulated annealing during the initial portion of the MCMC.\n",
          DEFAULT_CHAIN_LEN, DEFAULT_SAMPLE_LEN, DEFAULT_CLOSEST_N, separate_reference_set ? "true" : "false",
          no_HLA_inference ? "true" : "false", no_rev_inference ? "true" : "false",
          no_recombination ? "true" : "false");
  exit(EXIT_FAILURE);
}

// Also opens output files
static void dmodel_parse_cmdline(int argc, char **argv)
{

  static struct option longopts[] =
  {
    {"HLA_inference",              required_argument, NULL, 'H'},
    {"chain_length",               required_argument, NULL, 'n'},
    {"sample_every",               required_argument, NULL, 's'},
    {"runtime_hours",              required_argument, NULL, 't'},
    {"consensus_fasta",            required_argument, NULL, 'c'},
    {"resume_json",                required_argument, NULL, 'r'},
    {"closest_n",                  required_argument, NULL, 'l'},
    {"parameters_json",            required_argument, NULL, 'p'},
    {"true_mosaic_fasta",          required_argument, NULL, 'm'},
    {"separate_reference_fasta",   no_argument,       NULL, 'q'},
    {"no_HLA_inference",           no_argument,       NULL, 'i'},
    {"no_rev_inference",           no_argument,       NULL, 'j'}, 
    {"no_recombination_inference", no_argument,       NULL, 'o'},
    {"omega_gamma",                no_argument,       NULL, 'v'},
    {"omega_log_normal",           no_argument,       NULL, 'e'},
    {"omega_uniform",              no_argument,       NULL, 'w'},
    {"HLA_coeff_esc_gamma",        no_argument,       NULL, 'g'},
    {"HLA_coeff_esc_normal",       no_argument,       NULL, 'k'},
    {"HLA_coeff_esc_uniform",      no_argument,       NULL, 'f'},
    {"HLA_coeff_rev_gamma",        no_argument,       NULL, 'a'},
    {"HLA_coeff_rev_normal",       no_argument,       NULL, 'd'},
    {"HLA_coeff_rev_uniform",      no_argument,       NULL, 'u'},
    {"sample_prior",               no_argument,       NULL, 'b'},
    {"no_simulated_annealing",     no_argument,       NULL, 'z'},
    {NULL, 0, NULL, 0}
  };

  char shortopts[] = "H:n:s:t:c:r:l:p:m:qijovewgkfadubz";

  int optional;
  while((optional = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(optional) {
      case 'H': hla_path = optarg; break;
      case 'n': 
        chain_length = atoi(optarg);
        number_of_states_set = true;
        break;
      case 's': sample_rate = atoi(optarg); break;
      case 't':
        runtime_limit_hrs = atoi(optarg);
        runtime_limit_hrs_set = true;
      break;
      case 'c': cons_option_path = optarg; break;
      case 'r': resume_state_path = optarg; break;
      case 'l': closest_n = atoi(optarg); break;
      case 'p': json_parameters_path = optarg; break;
      case 'm': true_mosaic_sequences = optarg; break;
      case 'q': separate_reference_set = false; break;
      case 'i': no_HLA_inference = true; break;
      case 'j': no_rev_inference = true; break;
      case 'o': no_recombination = true; break;
      case 'v': omega_prior = GAMMA_PRIOR; break;
      case 'e': omega_prior = LOG_NORMAL_PRIOR; break;
      case 'w': omega_prior = FLAT_PRIOR; break;
      case 'g': HLA_coeff_prior_esc = GAMMA_PRIOR; break;
      case 'k': HLA_coeff_prior_esc = NORMAL_PRIOR; break;
      case 'f': HLA_coeff_prior_esc = FLAT_PRIOR; break;
      case 'a': HLA_coeff_prior_rev = GAMMA_PRIOR; break;
      case 'd': HLA_coeff_prior_rev = NORMAL_PRIOR; break;
      case 'u': HLA_coeff_prior_rev = FLAT_PRIOR; break;
      case 'b': sample_prior = true; break;
      case 'z': no_simulated_annealing = true; break;
      default: die("Unknown option: %c", optional);
    }
  }

  if(argc - optind != 3) print_usage();
  if(runtime_limit_hrs_set == true && number_of_states_set == true)
    die("Set both a time limit and number of states");

  ref_path = argv[optind];
  query_path = argv[optind+1];
  out_dir = argv[optind+2];
  // Set the sequence file from which to define consensus as the reference set (as default).
  cons_path = argv[optind];

  // If an optional argument is passed, then set the consensus filepath as this instead.
  if(cons_option_path != NULL)
    cons_path = cons_option_path;

  if(!mkpath(out_dir, 0755)) die("Cannot create output dir: %s", out_dir);

  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char date[100];
  sprintf(date, "%d_%d_%d_%d_%d_%d", timeinfo->tm_year+1900, 
          timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
          timeinfo->tm_min, timeinfo->tm_sec);

  char log_path[PATH_MAX+1], window_path[PATH_MAX+1];
  char priors_path[PATH_MAX+1], summary_path[PATH_MAX+1];
  char hamming_filepath[PATH_MAX+1], hamming_distance_filepath[PATH_MAX+1];
  char json_path[PATH_MAX+1];

  sprintf(log_path, "%s/%s_log.txt", out_dir, date);  
  sprintf(window_path, "%s/%s_window.txt", out_dir, date);
  sprintf(priors_path, "%s/%s_priors.txt", out_dir, date);
  sprintf(summary_path, "%s/%s_summary.txt", out_dir, date);
  sprintf(json_path, "%s/%s_final_state.json", out_dir, date);
  sprintf(interim_json_path, "%s/%s_interim_final_state.json", out_dir, date);
  sprintf(hamming_filepath, "%s/%s_hamming.txt", out_dir, date);
  sprintf(hamming_distance_filepath, "%s/%s_hamming_distance.txt", out_dir, date);
  
  while(futil_file_exists(log_path))
  {
    printf("This filename already exists, wait for a second to change the folder name.\n");
    sleep(1);
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(date, "%d_%d_%d_%d_%d_%d\n", timeinfo->tm_year+1900, 
            timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);

    sprintf(log_path, "%s/%s_log.txt", out_dir, date);  
    sprintf(window_path, "%s/%s_window.txt", out_dir, date);
    sprintf(priors_path, "%s/%s_priors.txt", out_dir, date);
    sprintf(summary_path, "%s/%s_summary.txt", out_dir, date);
    sprintf(json_path, "%s/%s_final_state.json", out_dir, date);
    sprintf(interim_json_path, "%s/%s_interim_final_state.json", out_dir, date);
    sprintf(hamming_filepath, "%s/%s_hamming.txt", out_dir, date);
    sprintf(hamming_distance_filepath, "%s/%s_hamming_distance.txt", out_dir, date);
  }

  printf("Writing logs to: %s\n", log_path);
  if(hla_path != NULL || no_rev_inference == false)
    printf("Writing windows to: %s\n", window_path);
  printf("Writing priors to: %s\n", priors_path);
  printf("Writing summary to: %s\n", summary_path);
  printf("Writing final state to: %s\n", json_path);
  printf("Writing interim final state to: %s\n", interim_json_path);

  log_file = fopen(log_path, "w");
  if(hla_path != NULL || no_rev_inference == false) {
    if(hla_path != NULL) printf("Performing HLA selection inference.\n");
    if(no_rev_inference) printf("Performing reversion inference.\n");
    window_file = fopen(window_path, "w");
  } else {
    printf("No HLA.csv file passed, and not performing reversion inference.\n");
  }

  priors_file = fopen(priors_path, "w");  
  summary_file = fopen(summary_path, "w");
  final_state_json = fopen(json_path, "w");
  hamming_file = fopen(hamming_filepath, "w");
  hamming_distance_file = fopen(hamming_distance_filepath, "w");

  if(log_file == NULL)
    die("Cannot open output file: %s", log_path);
  if(window_file == NULL && (hla_path != NULL || no_rev_inference == false))
    die("Cannot open output file: %s", window_path);
  if(priors_file == NULL)
    die("Cannot open output file: %s", priors_path);
  if(summary_file == NULL)
    die("Cannot open output file: %s", summary_path);

  // Add this information to the summary file.
  fprintf(summary_file, "chain length: %i\n""sample every: %i\n"
                        "reference sequence file: %s\n"
                        "query sequence file: %s"
                        "%s%s\n"
                        "consensus sequence file: %s\n", 
          chain_length, sample_rate, ref_path, query_path,
          hla_path != NULL ? "\nHLA csv file: " : "",
          hla_path != NULL ? hla_path : "", cons_path);
}

static void memory_allocation(int num_refs, int num_query, int num_cons, int closest_n,
                                int num_sites, int num_bases, int num_hla_types)
{
  // Approximate the total amount of memory allocation.
  size_t total_memory_alloc = 0;
  // ref_seqs.
  total_memory_alloc += num_refs * num_bases * sizeof(char);
  // query_seqs.
  total_memory_alloc += num_query * num_bases * sizeof(char);
  // cons_seqs.
  total_memory_alloc += num_cons * num_bases * sizeof(char);
  // hla_types.
  total_memory_alloc += sizeof(char) * num_query * (num_hla_types +1);
  // closest_refs.
  total_memory_alloc += num_query * sizeof(int[closest_n]);
  // hamming.
  total_memory_alloc += num_refs;
  // omega.
  total_memory_alloc += num_sites * sizeof(Decimal);
  // R.
  total_memory_alloc += num_sites * sizeof(Decimal);
  // HLASelection.
  total_memory_alloc += num_hla_types * sizeof(HLASelection);
  // ref_triplets.
  total_memory_alloc += num_refs * sizeof(char[num_sites]);
  // ref_codons.
  total_memory_alloc += num_refs * sizeof(char[num_sites]);
  // query_codons.
  total_memory_alloc += num_query * sizeof(char[num_sites]);
  // cons_query_nucls.
  total_memory_alloc += num_query * sizeof(char[num_bases]);
  
  printf("Memory allocation before freeing sequence allocation: %zu bytes = %zu GB\n", 
         total_memory_alloc,total_memory_alloc/(1<<30));
  // (1<<30) is 2^30 in binary.
  total_memory_alloc -= (num_refs + num_query + num_cons) * num_bases * sizeof(char);

  // num_ref_site_codons.
  total_memory_alloc += num_sites * sizeof(int);
  // ref_site_codons.
  total_memory_alloc += num_sites * sizeof(char[NUM_CODONS]);
  // F.
  total_memory_alloc += num_query * sizeof(Decimal[closest_n][num_sites]);
  // B.
  total_memory_alloc += num_query * sizeof(Decimal[closest_n][num_sites]);
  // F_rnorm.
  total_memory_alloc += num_query * sizeof(int[num_sites]);
  // B_rnorm.
  total_memory_alloc += num_query * sizeof(int[num_sites]);
  // tmp_F.
  total_memory_alloc += num_query * sizeof(Decimal[closest_n][num_sites]);
  // tmp_codon_to_codon_HLA.
  total_memory_alloc += num_query * sizeof(Decimal[NUM_CODONS][num_sites]);
  // tmp_F_rnorm.
  total_memory_alloc += num_query * sizeof(int[num_sites]);
  // consensus_codons.
  total_memory_alloc += num_sites * sizeof(char);
  // codon_to_codon_HLA.
  total_memory_alloc += num_query * sizeof(Decimal[NUM_CODONS][num_sites]);
  // type_of_change.
  total_memory_alloc += num_sites * sizeof(char[NUM_CODONS][NUM_CODONS]);
  // A.
  total_memory_alloc += num_query * sizeof(Decimal[NUM_CODONS][num_sites]);
  // tmp_A.
  total_memory_alloc += num_query * sizeof(Decimal[NUM_CODONS][num_sites]);
  // rate_out_of_codons.
  total_memory_alloc += NUM_CODONS * sizeof(Decimal[num_sites]);
  // HLA_select_esc.
  total_memory_alloc += num_hla_types * sizeof(Decimal[num_sites]);
  // HLA_select_rev.
  total_memory_alloc += num_hla_types * sizeof(Decimal[num_sites]);

  printf("Total memory allocation: %zu bytes\n = %zu GB\n",total_memory_alloc,
         total_memory_alloc/(1<<30));
  // (1<<30) is 2^30 in binary.
}

void fill_closest_n(int closest_n, int (*closest)[closest_n], bool skip_first,
                    char **ref_seqs, char **query_seqs, int num_refs, int num_query,
                    int num_bases, int tmp_hamming[num_refs])
{
  // If want 100 closest and skip first, find 101 closest
  // int q, r, n = skip_first ? closest_n + 1 : closest_n;
  int q, n = skip_first ? closest_n + 1 : closest_n;

  int (*tmp_closest_refs)[n] = closest;
  if(skip_first) tmp_closest_refs = my_malloc(num_query * n * sizeof(int), __FILE__, __LINE__);
  assert(n <= num_refs);

  for(q = 0; q < num_query; q++) {
    hamming_distance_considering_unknown(ref_seqs, query_seqs[q],
                                         num_refs, num_bases, tmp_hamming);
    
    // Write hamming distances to hamming distance file
    // for(r = 0; r < num_refs; r++)
    //   fprintf(hamming_distance_file, DECPRINT" ",(Decimal)tmp_hamming[r]/2);
    // fprintf(hamming_distance_file, "\n");

    closest_sequences(num_refs, tmp_hamming, n, tmp_closest_refs[q]);
  }

  if(skip_first) {
    for(q = 0; q < num_query; q++) {
      memcpy(closest[q], tmp_closest_refs[q] + 1,
             closest_n * sizeof(closest[0][0]));
    }
    free(tmp_closest_refs);
  }
}

int main(int argc, char **argv)
{
  init_rand();
  setup_gsl();

  Decimal lgam = lgamma(gamma_k_shape_esc);
  Decimal gamma_const_esc = gamma_k_shape_esc * log(gamma_theta_scale_esc) + lgam;

  lgam = lgamma(gamma_k_shape_rev);
  Decimal gamma_const_rev = gamma_k_shape_rev * log(gamma_theta_scale_rev) + lgam;
  
  dmodel_parse_cmdline(argc, argv);

  char **ref_seqs, **query_seqs, **cons_seqs;
  int refs_cap, query_cap, cons_cap;
  int q, r, s, h, c;
  int c1, c2;

  printf("Loading sequences...\n");
  int num_refs = load_seqs(ref_path, &ref_seqs, &refs_cap);
  int num_query = load_seqs(query_path, &query_seqs, &query_cap);
  int num_cons = load_seqs(cons_path, &cons_seqs, &cons_cap);

  assert(num_refs > 0);
  assert(num_query > 0);
  assert(num_cons > 0);

  printf("Loaded %i reference sequences.\n", num_refs);
  printf("Loaded %i query sequences.\n", num_query);
  printf("Loaded %i sequences to determine consensus.\n", num_cons);

  // Check all sequences have the same length.
  size_t num_bases = strlen(ref_seqs[0]);

  if(num_bases % 3 != 0)
    die("Sequences contain partial codons [%zu mod 3 != 0].", num_bases);

  for(r = 1; r < num_refs; r++)
    if(strlen(ref_seqs[r]) != num_bases)
      die("Reference sequences aren't all the same length.");

  for(q = 0; q < num_query; q++)
    if(strlen(query_seqs[q]) != num_bases)
      die("Query sequences aren't all the same length.");

  for(c = 0; c < num_cons; c++) {
    if(strlen(cons_seqs[c]) != num_bases) {
      die("Sequences from which to derive the consensus sequence aren't all "
          "the same length.");
    }
  }
  
  char **hla_types;
  char **bools = NULL, *data = NULL;
  int num_hla_types = 0;

  if(hla_path  != NULL){
    printf("Loading HLA types...\n");
    num_hla_types = load_hla_csv(hla_path, &hla_types, num_query);
    assert(num_hla_types > 0);
    printf("Loaded %i HLA types.\n", num_hla_types);
  } else {
    bools = my_malloc(sizeof(char*) * num_query, __FILE__, __LINE__);
    data = my_malloc(sizeof(char) * num_query * (1+1), __FILE__, __LINE__);
    for(q = 0; q < num_query; q++){
      bools[q] = data + q + 2;
      bools[q][0] = '0';
      bools[q][1] = '\0';
    }
    hla_types = bools;
  }

  // Load in optional .json viral paramters.
  if(json_parameters_path != NULL)
    load_fixed_parameters_from_json(json_parameters_path, &kappa, &phi, prior_C1);
  
  printf("kappa:"DECPRINT"\n", kappa);
  printf("phi:"DECPRINT"\n", phi);
  
  lambda_rec = - log(1 - p_rec) * num_refs;
  lambda_rec = 1 / lambda_rec;
  printf("lambda_rec:"DECPRINT"\n", lambda_rec);
  
  for(c = 0; c < NUM_CODONS; c++) printf(DECPRINT" ", prior_C1[c]);
  printf("\n");

  // Each codon is a site.
  int num_sites = num_bases / 3;
  assert(num_sites > 0);

  memory_allocation(num_refs, num_query, num_cons, closest_n, num_sites, num_bases, 
                    num_hla_types);

  printf("Loaded %i sites.\n", num_sites);

  printf("Mapping gaps to consensus...\n");
  char consensus[num_bases+1];

  // Set the consensus sequence - default is the consensus of the references.
  // Otherwise it's the consensus of another sequence file that is passed.
  generate_consensus(cons_seqs, num_cons, num_bases, consensus);
  printf("Consensus:\n%s\n", consensus);

  int (*closest_refs)[closest_n] = my_malloc(num_query * sizeof(int[closest_n]),
                                           __FILE__, __LINE__);

  int *tmp_hamming = my_calloc(num_refs, sizeof(int), __FILE__, __LINE__);

  fill_closest_n(closest_n, closest_refs, !separate_reference_set,
                 ref_seqs, query_seqs, num_refs, num_query, num_bases, tmp_hamming);

  // Print closest refs for query sample 0
  printf("Example - closest references to query sequence 1:\n");
  for(r = 0; r < closest_n; r++) printf("%i ", closest_refs[0][r]);
  printf("\n");
  
  if(true_mosaic_sequences != NULL)
  {
    Mosaic *true_copy = my_malloc(num_query * sizeof(Mosaic), __FILE__, __LINE__);
    int *true_seqs = my_malloc(num_query * num_sites * sizeof(int), __FILE__, __LINE__);

    for(q = 0; q < num_query; q++) {
      true_copy[q].which_seqs = &true_seqs[q*num_sites];
      true_copy[q].num_seqs = 0;
    }

    load_closest_n_from_json(true_mosaic_sequences, num_query, num_sites, true_copy);

    for(q = 0; q < num_query; q++)
    { 
      // Permute the recombinant sequences.
      knuth_perm(true_copy[q].num_seqs, true_copy[q].which_seqs);
    
      // Some will have n closest relatives included, some will have n+1.
      // Determine floor(closest_n/true_copy[q].num_seqs).
      int n_closest_relatives = closest_n / true_copy[q].num_seqs;
      int remain = 0;
    
      if(closest_n % true_copy[q].num_seqs != 0)
        remain = closest_n % true_copy[q].num_seqs;
    
      int m;
      int m_total = 0;

      for(r = 0; r < (true_copy[q].num_seqs - remain); r++) {
        hamming_distance_considering_unknown(ref_seqs, ref_seqs[true_copy[q].which_seqs[r]],
                                             num_refs, num_bases, tmp_hamming);
        int mosaic_closest_refs[n_closest_relatives];
        closest_sequences(num_refs, tmp_hamming, n_closest_relatives, mosaic_closest_refs);
      
        // have to put these in the right place.
        for(m = 0; m < n_closest_relatives; m++, m_total++)
          closest_refs[q][m_total] = mosaic_closest_refs[m];
      }
      
      for(r = (true_copy[q].num_seqs - remain); r < true_copy[q].num_seqs; r++) {
        hamming_distance_considering_unknown(ref_seqs, ref_seqs[true_copy[q].which_seqs[r]],
                                             num_refs, num_bases, tmp_hamming);
        int mosaic_closest_refs[n_closest_relatives+1];
        closest_sequences(num_refs, tmp_hamming, (n_closest_relatives+1), mosaic_closest_refs);
        
        // have to put these in the right place.
        for(m = 0; m < (n_closest_relatives+1); m++, m_total++)
          closest_refs[q][m_total] = mosaic_closest_refs[m];
      }
    }
    // Free the Mosaic object, it's no longer needed.
    free(true_copy);
    free(true_seqs);
  }

  // Done with Hamming.
  fclose(hamming_distance_file);
  free(tmp_hamming);

   // Selection coeff in absence of HLA, per site.
  Decimal *omega = my_malloc(num_sites * sizeof(Decimal), __FILE__, __LINE__);
  Decimal init_omega = 1.0;
  for(s = 0; s < num_sites; s++) omega[s] = init_omega;

  // Generate R (population recombination probability)
  // R[i] is probability just before codon i
  Decimal *R = my_malloc(num_sites * sizeof(Decimal), __FILE__, __LINE__);
  for(s = 0; s < num_sites; s++) R[s] = init_R;

  if(no_recombination == true)
    for(s = 0; s < num_sites; s++) R[s] = 0;

  Decimal mu = init_mutation_rate;

  // Set up selection windows.
  if (hla_path == NULL){
    // We need to set up the reversion coefficients to 1,
    // as these get used to generate codon_to_codon_HLA.
    num_hla_types = 1; 
  }

  HLASelection *hla_windows = my_malloc(num_hla_types * sizeof(HLASelection),
                                        __FILE__, __LINE__);
  selec_window_hla_alloc(num_sites, num_hla_types, hla_windows);

  Decimal resume_llk = 0;
  int i;

  if(resume_state_path != NULL)
  {
    // Load last state in resume_state_path.
    for(i = 0; i < closest_n; i++)
      printf("%i ", closest_refs[(num_query-1)][i]);
    printf("\n\n");
    load_from_json(resume_state_path, &resume_llk, &mu, R, omega, hla_windows,
                   num_sites, num_hla_types, num_query, closest_n, closest_refs);
    for(i = 0; i < closest_n; i++)
      printf("%i ", closest_refs[(num_query-1)][i]);
    printf("\n\n");
  }
  else
  { 
    Decimal log_mu_sel_prior_esc = log(mu_sel_prior_esc);
    Decimal log_mu_sel_prior_rev = log(mu_sel_prior_rev);

    for(h = 0; h < num_hla_types; h++) {
      selec_windows_randomise(hla_windows[h].esc, num_sites, init_num_selec_windows_esc,
                              log_mu_sel_prior_esc, sigma_sel_prior_esc);
      selec_windows_randomise(hla_windows[h].rev, num_sites, init_num_selec_windows_rev,
                              log_mu_sel_prior_rev, sigma_sel_prior_rev);
    }
    
    // Set the escape windows to 1 if there's no HLA inference.
    if(hla_path == NULL || no_HLA_inference == true) {
      for(h = 0; h < num_hla_types; h++)
        for(i = 0; i < init_num_selec_windows_esc; i++)
          hla_windows[h].esc->windows[i].coeff = 1;
    }
    
    // Set the reversion windows to 1 if there's no reversion inference.
    if(no_rev_inference == true) {
      for(i = 0; i < init_num_selec_windows_rev; i++)
        hla_windows[0].rev->windows[i].coeff = 1;
    }
    
    // Uncomment to initialise HLA coefficients at 1 (used in simulations).
    // for(h = 0; h < num_hla_types; h++) {
    //   for(i = 0; i < init_num_selec_windows_esc; i++) {
    //     hla_windows[h].esc->windows[i].coeff = 1;
    //   }
    // }
    
    // Set the unused reversion coefficients to 1.
    for(h = 1; h < num_hla_types; h++) {
      for(i = 0; i < init_num_selec_windows_rev; i++) {
        hla_windows[h].rev->windows[i].coeff = 1;
      }
    }
    
    // for(h=0; h < num_hla_types; h++){
    //   selec_window_print(hla_windows[h].esc);
    //   printf("windows esc - %i\n",hla_windows[h].esc->num_windows);
    //   selec_window_print(hla_windows[h].rev);
    //   printf("windows rev - %i\n",hla_windows[h].rev->num_windows);
    // }
  }

  char (*ref_triplets)[num_sites] = my_malloc(num_refs * sizeof(char[num_sites]),
                                              __FILE__, __LINE__);
  char (*ref_codons)[num_sites] = my_malloc(num_refs * sizeof(char[num_sites]),
                                            __FILE__, __LINE__);
  char (*query_codons)[num_sites] = my_malloc(num_query * sizeof(char[num_sites]),
                                              __FILE__, __LINE__);
  // Determine a nucleotide consensus for each query sequence.
  char (*cons_query_nucls)[num_bases] = my_malloc(num_query * sizeof(char[num_bases]),
                                                  __FILE__, __LINE__);

  // At this point we fill in the missing data.
  char consensus_dummy[num_bases+1];
  for(q = 0; q < num_query; q++) {
    // Create a consensus based on the sequences that are closest to the query sequence.
    generate_consensus_subset(ref_seqs, closest_n, num_bases, consensus_dummy, closest_refs[q]);
    // Then need to set missing nucleotides in the queries equal to this consensus.
    for(s = 0; s < (int)num_bases; s++) {
      if(!is_base(query_seqs[q][s])) {
        query_seqs[q][s] = consensus_dummy[s];
        if(!is_base(consensus_dummy[s])) {
          die("The consensus is missing data for the seqs closest to this query.");
        }
      }
    }
    // Fill cons_query_nucls with the consensus nucleotides of the closest_n to the queries.
    for(s = 0; s < (int)num_bases; s++) {
      cons_query_nucls[q][s] = base_to_code(consensus_dummy[s]);
      status("%4i ", cons_query_nucls[q][s]);
    }
  }
  
  // Have to think a bit more carefully - ambiguous nucleotides still exist in the ref_seqs. 
  for(r = 0; r < num_refs; r++)
    for(s = 0; s < num_sites; s++)
      ref_triplets[r][s] = triplet_to_code(ref_seqs[r]+s*3);

  // Here, also generate the ref_codons in the usual way, and use this to determine the collection of codons
  // at each site (used in codon to codon HLA). Everywhere else, use the ref_triplets and triplet to codon.
  // triplet to codon, either just returns the equivalent encoding in the 4 letter alphabet TCAG, or, if theres
  // a gap, uses the query specific consensus to fill it in.

  map_gaps_to_consensus(ref_seqs, num_refs, num_bases, consensus);

  for(r = 0; r < num_refs; r++)
    for(s = 0; s < num_sites; s++)
      ref_codons[r][s] = amino_to_code(ref_seqs[r]+s*3);

  for(q = 0; q < num_query; q++)
    for(s = 0; s < num_sites; s++)
      query_codons[q][s] = amino_to_code(query_seqs[q]+s*3);

  // Free sequence data - only work with codons from now on.
  for(r = 0; r < num_refs; r++) free(ref_seqs[r]);
  for(q = 0; q < num_query; q++) free(query_seqs[q]);
  for(c = 0; c < num_cons; c++) free(cons_seqs[c]);
  free(ref_seqs);
  free(query_seqs);
  free(cons_seqs);

  int *num_ref_site_codons = my_malloc(num_sites * sizeof(int), __FILE__, __LINE__);
  char (*ref_site_codons)[NUM_CODONS] = my_malloc(num_sites * sizeof(char[NUM_CODONS]),
                                                  __FILE__, __LINE__);
  list_ref_codons_at_all_sites(num_sites, num_refs,
                               num_ref_site_codons, ref_site_codons,
                               ref_codons);

  printf("Reference codons: ");

  for(s = 0; s < num_sites; s++) {
    if(s > 0) printf(" | ");
    printf("%i", (int)ref_site_codons[s][0]);
    for(i = 1; i < num_ref_site_codons[s]; i++)
      printf(",%i", (int)ref_site_codons[s][i]);
  }
  printf("\n");

  int sample_i = 1;

  // Important change here - change the number of references.
  // This is set to the closest_n sequences to each query.
  // Note: if we are to consider variable numbers of sequences 
  // for each query, this will have to be coded differently.

  // DEV: changed num_refs to total num_refs and reallocated num_refs as closest_n
  // From here, closest_n and num_refs are the same.

  int total_num_refs = num_refs;

  // Note: Indexed as: x[query]...
  Decimal (*F)[closest_n][num_sites] = my_malloc(num_query * sizeof(Decimal[closest_n][num_sites]),
                                                 __FILE__, __LINE__);
  Decimal (*B)[closest_n][num_sites] = my_malloc(num_query * sizeof(Decimal[closest_n][num_sites]),
                                                 __FILE__, __LINE__);
  int (*F_rnorm)[num_sites] = my_malloc(num_query * sizeof(int[num_sites]), __FILE__, __LINE__);
  int (*B_rnorm)[num_sites] = my_malloc(num_query * sizeof(int[num_sites]), __FILE__, __LINE__);

  // Temp output
  Decimal (*tmp_F)[closest_n][num_sites] = my_malloc(num_query * sizeof(Decimal[closest_n][num_sites]),
                                                     __FILE__, __LINE__);
  Decimal (*tmp_codon_to_codon_HLA)[NUM_CODONS][num_sites]
    = my_malloc(num_query * sizeof(Decimal[NUM_CODONS][num_sites]), __FILE__, __LINE__);
  int (*tmp_F_rnorm)[num_sites] = my_malloc(num_query * sizeof(int[num_sites]), __FILE__, __LINE__);

  // Determine the consensus codon sequence.
  char *consensus_codons = my_malloc(num_sites * sizeof(char), __FILE__, __LINE__);
  printf("Consensus codons:\n");
  for(s = 0; s < num_sites; s++) {
    consensus_codons[s] = amino_to_code(consensus+s*3);
    printf("%4i ", consensus_codons[s]);
  }
  printf("\n");

  for(s = 0; s < num_sites; s++) {
    if(s != num_sites -1) {
      // fprintf(hamming_file, "%s", codon_to_dna[(int)consensus_codons[s]]);
      fprintf(summary_file, "%s", codon_to_dna[(int)consensus_codons[s]]);
    }
  }
  // fprintf(hamming_file, "%s\n", codon_to_dna[(int)consensus_codons[num_sites-1]]);
  fprintf(summary_file, "%s\n", codon_to_dna[(int)consensus_codons[num_sites-1]]);

  // Print the first reference sequence for comparison.
  printf("Reference codons 1:\n");
  for(s = 0; s < num_sites; s++) {
    printf("%4i ", ref_codons[1][s]);
  }
  printf("\n");

  // For the hamming_file, print the number of each reference sequence closest to each query sequence 
  // This can then be read into R and used to create plots.
  // for(q = 0; q < num_query; q++) {
  //   fprintf(hamming_file, "%i\n", closest_refs[q][0]);
  // }

  fclose(hamming_file);
  fclose(summary_file);

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

  Decimal (*codon_to_codon_HLA)[NUM_CODONS][num_sites]
    = my_calloc(num_query, sizeof(Decimal[NUM_CODONS][num_sites]), __FILE__, __LINE__);

  // type_of_change[num_sites][NUM_CODONS][NUM_CODONS]
  char (*type_of_change)[NUM_CODONS][NUM_CODONS]
    = my_malloc(num_sites * sizeof(char[NUM_CODONS][NUM_CODONS]), __FILE__, __LINE__);

  // A[num_query][NUM_CODONS][num_sites]
  Decimal (*A)[NUM_CODONS][num_sites]
    = my_malloc(num_query * sizeof(Decimal[NUM_CODONS][num_sites]), __FILE__, __LINE__);
  Decimal (*tmp_A)[NUM_CODONS][num_sites]
    = my_malloc(num_query * sizeof(Decimal[NUM_CODONS][num_sites]), __FILE__, __LINE__);

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

  // Set 'Lambda'.
  Decimal (*rate_out_of_codons)[num_sites] = my_malloc(NUM_CODONS * sizeof(Decimal[num_sites]),
                                                       __FILE__, __LINE__);

  for(c1 = 0; c1 < NUM_CODONS; c1++)
  {
    Decimal v1 = kappa * beta_S[c1] + beta_V[c1];
    Decimal v2 = kappa * alpha_S[c1] + alpha_V[c1];

    for(s = 0; s < num_sites; s++) {
      rate_out_of_codons[c1][s] = v1 + omega[s] * v2;
    }
  }
 
  // Generate HLA_selection_esc - HLA selection coefficients associated to escape.
  Decimal (*HLA_select_esc)[num_sites] = my_malloc(num_hla_types * sizeof(Decimal[num_sites]),
                                                   __FILE__, __LINE__);
  // Generate HLA_selection_rev - HLA selection coefficients associated to reversion.
  Decimal (*HLA_select_rev)[num_sites] = my_malloc(num_hla_types * sizeof(Decimal[num_sites]),
                                                   __FILE__, __LINE__);

  // Copy window coeff into array per site.
  for(h = 0; h < num_hla_types; h++) {
    update_HLA_selection(hla_windows[h].esc, num_sites, HLA_select_esc[h]);
    update_HLA_selection(hla_windows[h].rev, num_sites, HLA_select_rev[h]);
  }

  Decimal ref_sum, curr_llk_across_queries = 0;

  // Initial setup codon_to_codon_HLA.
  Decimal sum_prior = 0;
  for(i = 0; i < (NUM_CODONS-1); i++) sum_prior += prior_C1[i];
  prior_C1[NUM_CODONS-1] = 1 - sum_prior;

  if(sample_prior == false)
  {
    create_codon_to_codon_HLA(num_sites, total_num_refs, 
                              num_query, num_hla_types,
                              num_ref_site_codons, ref_site_codons,  
                              codon_to_codon_HLA, rate_out_of_codons,
                              0, num_sites-1, DO_ALL_HLAS, mu,
                              HLA_select_esc, HLA_select_rev, hla_types,
                              consensus_codons, omega,
                              NS_TS_sum_from_consensus, NS_TV_sum_from_consensus,
                              type_of_change, query_codons, A,
                              DO_HLA_BOTH);

    // Initialise F and B.
    for(q = 0; q < num_query; q++)
    {
      initialise_F_numerical_recipes(num_sites, closest_n, total_num_refs,
                                     codon_to_codon_HLA[q], ref_codons,
                                     ref_triplets, R, F[q], F_rnorm[q],
                                     0, num_sites-1, closest_refs[q],
                                     num_bases, cons_query_nucls[q]); // this line is for the ambiguous codons.
    
      update_B_numerical_recipes(num_sites, closest_n, total_num_refs,
                                 codon_to_codon_HLA[q], ref_codons, ref_triplets,
                                 R, B[q], B_rnorm[q], num_sites-1,
                                 closest_refs[q],
                                 num_bases, cons_query_nucls[q]);

      ref_sum = 0;
      for(r = 0; r < closest_n; r++) ref_sum += F[q][r][num_sites-1];
      curr_llk_across_queries += log(ref_sum) - F_rnorm[q][num_sites-1] * log(BIG);
    }
    printf("Total log-likelihood using initialise F: "DECPRINT".\n", curr_llk_across_queries);
  }
  
  // Run a check if we resumed a previous run.
  if(resume_state_path != NULL)
    printf("Likelihood from loaded .json file: "DECPRINT".\n", resume_llk);

  // Begin MCMC:
  time_t start, end;
  Decimal diff;
  time(&start);
  int hours_passed = 1;
  
  // Define relative probabilities of performing each move.
  int esc_mergesplit, rev_mergesplit, esc_grow, rev_grow;
  int esc_wind_coeff, rev_wind_coeff;
  int mut, rec, sel;

  // Turn on and off moves for easy debugging
  bool rec_on = true, sel_on = true;
  bool esc_mergesplit_on = true, rev_mergesplit_on = true;
  bool esc_grow_on = true, rev_grow_on = true;
  bool esc_wind_coeff_on = true, rev_wind_coeff_on = true;
  bool mut_on = true;

  if(no_HLA_inference == true) {
    esc_mergesplit_on = false;
    esc_grow_on = false;
    esc_wind_coeff_on = false;
  }

  if(no_rev_inference == true) {
    rev_mergesplit_on = false;
    rev_grow_on = false;
    rev_wind_coeff_on = false;   
  }
  
  if(no_recombination == true) rec_on = false;

  if(hla_path == NULL) {
    esc_mergesplit_on = false;
    esc_grow_on = false;
    esc_wind_coeff_on = false;
  }

  // Danny used this multiplier.
  // for escapes we multiply by num_hla_types.
  int danny_prob = 5;

  rec = (rec_on ? 20 : 0);
  sel = (sel_on ? 20 : 0) + rec;
  esc_mergesplit = (esc_mergesplit_on ? num_hla_types*danny_prob : 0) + sel;
  rev_mergesplit = (rev_mergesplit_on ? danny_prob : 0) + esc_mergesplit;
  esc_grow = (esc_grow_on ? num_hla_types * danny_prob : 0) + rev_mergesplit;
  rev_grow = (rev_grow_on ? danny_prob : 0) + esc_grow;
  esc_wind_coeff = (esc_wind_coeff_on ? num_hla_types * danny_prob : 0) + rev_grow;
  rev_wind_coeff = (rev_wind_coeff_on ? danny_prob : 0) + esc_wind_coeff;
  mut = (mut_on ? 1 : 0 ) + rev_wind_coeff;

  int which_move = mut;
  
  // Headers for the log_file.
  fprintf(log_file,"state posterior mu ");
  for(s = 0; s < num_sites; s++) {
    fprintf(log_file, "omega[%i] ", s);
  }
  for(s = 1; s < num_sites; s++) {
    fprintf(log_file, "R[%i] ", s);
  }
  fprintf(log_file, "\n");
  if (hla_path != NULL) 
    fprintf(window_file, "number of HLAs: %i\n", num_hla_types);
  if (hla_path == NULL && no_rev_inference == false) 
    fprintf(window_file, "reversion inference only\n");
  
  fprintf(priors_file, "Closest sequences: %i\n"
                       "Mutation\n"
                       "mu_prior: "DECPRINT"\n"
                       "Recombination\n"
                       "lambda_rec: "DECPRINT"\n"
                       "Omega\n",
                       closest_n, mu_prior, lambda_rec);

  switch(omega_prior){
    case EXPONENTIAL_PRIOR:
      fprintf(priors_file, "omega prior: exponential\n"
                           "lambda_sel: "DECPRINT"\n", lambda_sel);
      break;
    case GAMMA_PRIOR:
      fprintf(priors_file, "omega prior: gamma\n"
                           "omega_gamma_k_shape: "DECPRINT" "
                           "omega_gamma_theta_scale: "DECPRINT"\n",
                           omega_gamma_k_shape, omega_gamma_theta_scale);
      break;
    case LOG_NORMAL_PRIOR:
      fprintf(priors_file, "omega prior: log normal\n"
                           "omega_mu_sel: "DECPRINT" "
                           "omega_sigma_sel: "DECPRINT"\n",
                           omega_mu_sel, omega_sigma_sel);
      break;
    case FLAT_PRIOR:
      fprintf(priors_file, "omega prior: uniform\n");
      break;
    default:
      die("No prior detected!");
  }

  fprintf(priors_file, "Escape windows\n");

  switch(HLA_coeff_prior_esc)
  {
    case LOG_NORMAL_PRIOR:
      // Print priors to the file priors.txt.
      printf("Log Normal prior for escape coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: log normal\n"
                           "log_normal_mu_esc: "DECPRINT" log_normal_sigma_esc: "DECPRINT"\n",
                           mu_sel_prior_esc, sigma_sel_prior_esc);
      break;
    case NORMAL_PRIOR:
      printf("Truncated Normal prior for escape coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: truncated normal\n"
                           "normal_mu_esc: "DECPRINT" normal_sigma_esc: "DECPRINT"\n",
                           mu_sel_norm_prior_esc, sigma_sel_norm_prior_esc);
      break;
    case GAMMA_PRIOR:
      printf("Gamma prior for escape coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: gamma\n"
                           "gamma_k_shape_esc: "DECPRINT" gamma_theta_scale_esc: "DECPRINT"\n",
                           gamma_k_shape_esc, gamma_theta_scale_esc);
      break;
    case FLAT_PRIOR:
      printf("Flat prior for escape coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: uniform\n");
      break;
    default: 
      die("No prior detected!");
  }

  fprintf(priors_file, "prior_add_window_esc: "DECPRINT"\n"
                       "Reversion windows\n", prior_add_window_esc);
  
    switch(HLA_coeff_prior_rev)
  {
    case LOG_NORMAL_PRIOR:
      // Print priors to the file priors.txt.
      printf("Log Normal prior for reversion coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: log normal\n"
                           "log_normal_mu_rev: "DECPRINT" log_normal_sigma_rev: "DECPRINT"\n",
                           mu_sel_prior_rev, sigma_sel_prior_rev);
      break;
    case NORMAL_PRIOR:
      printf("Truncated Normal prior for reversion coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: truncated normal\n"
                           "normal_mu_rev: "DECPRINT" normal_sigma_rev: "DECPRINT"\n",
                           mu_sel_norm_prior_rev, sigma_sel_norm_prior_rev);
      break;
    case GAMMA_PRIOR:
      printf("Gamma prior for reversion coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: gamma\n"
                           "gamma_k_shape_rev: "DECPRINT" gamma_theta_scale_rev: "DECPRINT"\n",
                           gamma_k_shape_rev, gamma_theta_scale_rev);
      break;
    case FLAT_PRIOR:
      printf("Flat prior for reversion coefficients selected.\n");
      fprintf(priors_file, "coefficient prior: uniform\n");
      break;
    default: 
      die("No prior detected!");
  }

  fprintf(priors_file, "prior_add_window_rev: "DECPRINT"\n"
                       "selection_window_geom_param: "DECPRINT"\n",
                       prior_add_window_rev, selection_window_geom_param);

  // Close the prior file
  fclose(priors_file);

  int chain_i;
  Decimal initial_temperature = 500;
  if(no_simulated_annealing == true) initial_temperature = 1;
  Decimal temperature;

  for(chain_i = 0; chain_i < chain_length; chain_i++, sample_i++)
  {
    status("chain: %i\n", chain_i);
    int move;
    int choose_move = rand_lim(which_move);
    move = choose_move < rec ? 0 :
           (choose_move < sel ? 1 : 
            (choose_move < esc_mergesplit ? 2 :
             (choose_move < rev_mergesplit ? 3 :
              (choose_move < esc_grow ? 4 :
               (choose_move < rev_grow ? 5 :
                (choose_move < esc_wind_coeff ? 6 :
                 (choose_move < rev_wind_coeff ? 7 : 8)))))));

    int waction;

    // Determine temperature
    temperature = 1 + (initial_temperature-1) * pow(0.999, chain_i);

    // Fill in log and window file if required.
    if(sample_i == sample_rate)
    {
      sample_i = 0;
      // Print state to log.file.
      fprintf(log_file,"%i "DECPRINT" "DECPRINT" ", chain_i, curr_llk_across_queries,mu);
      for(s = 0; s < num_sites; s++) {
        fprintf(log_file, DECPRINT" ", omega[s]);
      }
      for(s = 1; s < num_sites; s++) {
        fprintf(log_file, DECPRINT" ", R[s]);
      }
      fprintf(log_file, "\n");
      
      if (hla_path != NULL){
        fprintf(window_file, "%i\n", chain_i);
      
        selec_window_print_to_file(hla_windows[0].rev, window_file);

        for(h = 0; h < num_hla_types; h++) 
          selec_window_print_to_file(hla_windows[h].esc, window_file);
      } else if (no_rev_inference == false){
        fprintf(window_file, "%i\n", chain_i);
        selec_window_print_to_file(hla_windows[0].rev, window_file);
      }
    }

    switch(move)
    {
      case 0:
        recombination_move(temperature, num_sites, closest_n, total_num_refs,
                           num_query, R,
                           F, B,
                           F_rnorm, B_rnorm,
                           codon_to_codon_HLA, ref_codons, ref_triplets,
                           tmp_F,
                           &curr_llk_across_queries,
                           closest_refs,
                           tmp_F_rnorm,
                           num_bases, cons_query_nucls);
        break;
      case 1:
        selection_move(temperature, num_sites, closest_n, total_num_refs,
                       num_query, num_hla_types, 
                       num_ref_site_codons, ref_site_codons,
                       mu, R, F, B,
                       F_rnorm, B_rnorm,
                       codon_to_codon_HLA, rate_out_of_codons,
                       HLA_select_esc, HLA_select_rev,
                       hla_types, consensus_codons, ref_codons, ref_triplets, query_codons,
                       omega, NS_TS_sum_from_consensus, NS_TV_sum_from_consensus,
                       type_of_change, A, tmp_codon_to_codon_HLA,
                       tmp_F,
                       &curr_llk_across_queries,
                       closest_refs, tmp_F_rnorm, tmp_A,
                       num_bases, cons_query_nucls);
        break;
      case 2: // esc DO_MERGE_OR_SPLIT
      case 3: // rev DO_MERGE_OR_SPLIT
      case 4: // esc DO_GROW
      case 5: // rev DO_GROW
      case 6: // esc DO_COEFF
      case 7: // rev DO_COEFF
        waction = move <= 3 ? DO_MERGE_OR_SPLIT
                            : (move <= 5 ? DO_GROW : DO_COEFF);

        window_move(temperature, num_sites, closest_n, total_num_refs,
                    num_query, num_hla_types,
                    num_ref_site_codons, ref_site_codons,
                    move & 1 ? DO_HLA_REV : DO_HLA_ESC, waction, mu, R,
                    F, B,
                    F_rnorm, B_rnorm,
                    codon_to_codon_HLA, rate_out_of_codons,
                    HLA_select_esc, HLA_select_rev, hla_types,
                    consensus_codons, ref_codons, ref_triplets, query_codons, omega,
                    NS_TS_sum_from_consensus, NS_TV_sum_from_consensus,
                    type_of_change, A, hla_windows,
                    tmp_codon_to_codon_HLA, tmp_F,
                    &curr_llk_across_queries, 
                    closest_refs, tmp_F_rnorm, tmp_A,
                    num_bases, cons_query_nucls,
                    move & 1 ? gamma_const_rev : gamma_const_esc);
        break;
      case 8:
        mutation_move(temperature, num_sites, closest_n, total_num_refs,
                      num_query, 
                      num_ref_site_codons, ref_site_codons,
                      &mu, R, &F, B,
                      &F_rnorm, B_rnorm,
                      &codon_to_codon_HLA, ref_codons, ref_triplets, query_codons,
                      &A, &tmp_A, &tmp_codon_to_codon_HLA,
                      &tmp_F,
                      &curr_llk_across_queries,
                      closest_refs,
                      &tmp_F_rnorm,
                      num_bases, cons_query_nucls);
        break;
    }
    
    if(runtime_limit_hrs_set) {
      time(&end);
      diff = difftime(end,start);
      int diff_hours = (int) (diff / (60 * 60));
      if(diff_hours >= runtime_limit_hrs) {
        write_json(final_state_json, chain_i, curr_llk_across_queries, mu, num_sites,
                   omega, R, hla_windows, num_hla_types, num_query, closest_n, closest_refs);
        break;
      } else if(chain_i == (chain_length - 1)) {
        chain_length *= 2;
      } else if(diff_hours >= hours_passed) {
        interim_final_state_json = fopen(interim_json_path, "w");
        write_json(interim_final_state_json, chain_i, curr_llk_across_queries, mu, num_sites,
                   omega, R, hla_windows, num_hla_types, num_query, closest_n, closest_refs);
        fclose(interim_final_state_json);
        hours_passed+=1;
      }
    } else if(chain_i == (chain_length-1)) {
      write_json(final_state_json, chain_i, curr_llk_across_queries, mu, num_sites,
                 omega, R, hla_windows, num_hla_types, num_query, closest_n, closest_refs);
      break;
    } else {
      time(&end);
      diff = difftime(end,start);
      int diff_hours = (int) (diff / (60 * 60));
      if(diff_hours >= hours_passed) {
        interim_final_state_json = fopen(interim_json_path, "w");
        write_json(interim_final_state_json, chain_i, curr_llk_across_queries, mu, num_sites,
                   omega, R, hla_windows, num_hla_types, num_query, closest_n, closest_refs);
        fclose(interim_final_state_json);
        hours_passed += 1;
      }
    }
  }
  printf("Proposed %i moves.\n",chain_i);

  // Print out the various moves accepted and rejected during the MCMC.
  printf("rec_accept: %i rec_reject: %i\n", rec_accept, rec_reject);
  printf("sel_accept: %i sel_reject: %i\n", sel_accept, sel_reject);
  printf("merge_accept: %i merge_reject: %i\n", merge_accept, merge_reject);
  printf("split_accept: %i split_reject: %i\n", split_accept, split_reject);
  printf("grow_accept: %i grow_reject: %i\n", grow_accept, grow_reject);
  printf("window_coeff_accept: %i window_coeff_reject: %i\n", wcoeff_accept, wcoeff_reject);
  printf("mut_accept: %i mut_reject: %i\n", mut_accept, mut_reject);

  time(&end);
  diff = difftime(end,start);
  // Print the computation time in seconds.
  printf("time taken: "DECPRINTTIME"\n", diff);
  
  #ifdef USE_TRIPLET
    printf("USED REF TRIPLET\n");
  #else
    printf("USED REF CODON\n");
  #endif
  
  // Free momory.
  free(closest_refs);
  free(rate_out_of_codons);
  free(ref_codons);
  free(ref_triplets);
  free(query_codons);
  free(num_ref_site_codons);
  free(ref_site_codons);
  free(consensus_codons);
  free(cons_query_nucls);
  free(type_of_change);
  free(codon_to_codon_HLA);
  free(tmp_codon_to_codon_HLA);
  free(HLA_select_esc);
  free(HLA_select_rev);
  free(A);
  free(tmp_A);
  free(R);
  free(omega);
  free(F);
  free(tmp_F);
  free(B);
  free(F_rnorm);
  free(B_rnorm);
  free(tmp_F_rnorm);

  selec_window_hla_dealloc(hla_windows);
  free(hla_windows);

  if(hla_path != NULL){
    free(hla_types[0]);
    free(hla_types);
  } else {
    free(data);
    free(bools);
  }
  
  // Close the log and window files.
  fclose(log_file);
  if (hla_path != NULL || no_rev_inference == false)
    fclose(window_file);
  fclose(final_state_json);

  clearup_gsl();

  return EXIT_SUCCESS;
}
