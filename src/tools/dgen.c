// Request decent POSIX version.
#define _XOPEN_SOURCE 700

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h> // Seed random.
#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <string.h>
#include <getopt.h>

#include "util.h"
#include "load_data.h"
#include "constants.h"
#include "amino_acids.h"
#include "load_from_json.h"
#include "birth_death_gen.h"
#include "simulate_sequences_sim.h"
#include "simulate_sequences_gen.h"
#include "write_tree.h"

Decimal past_sampling = -1;
Decimal query_fraction = 0.01;
int N, M, num_queries = -1;
Decimal lambda, mu_tree;

bool write_tree_to_file = false;
bool write_newick_tree_to_file = false;

FILE *tree_data_file = NULL;
FILE *newick_tree_data_file = NULL;
FILE *summary_file = NULL;
FILE *json_summary_file = NULL;

const char *cons_path = NULL, *root_path = NULL, *out_dir = NULL;
char *json_parameters_path = NULL;

static void print_usage()
{
  fprintf(stderr,
          "usage: dgen [options] <.json file> <N> <M> <lambda> <mu_tree> <outdir>\n"
          "  options:\n"
          "    --prop_past_sampling -s <S> Set sampling proportion for historic sampling [default: M/N].\n"
          "    --consensus_fasta    -c <C> Pass .fasta file from which to generate the consensus sequence.\n"
          "    --root_fasta         -r <R> Pass .fasta file giving the sequence a the root node.\n"
          "    --write_newick_tree  -n Write a Newick tree file.\n"
          "    --write_tree_matrix  -g Write a tree file.\n"
          "    --num_queries        -q provide the number of queries.\n"
          "    --prop_query         -f provide the query fraction [default: "DECPRINT"].\n", query_fraction);
  exit(EXIT_FAILURE);
}

static void dgen_parse_cmdline(int argc, char **argv)
{
  // DEV: could make it optional to write the HLA types at the internal nodes.
  // DEV: talk to Gil about generating his approximation to the coalescent 
  // with recombination - different trees as you track along the viral sequence.
  static struct option longopts[] =
  {
    {"prop_past_sampling", required_argument, NULL, 's'},
    {"consensus_fasta",    required_argument, NULL, 'c'},
    {"root_fasta",         required_argument, NULL, 'r'},
    {"write_newick_tree",  no_argument,       NULL, 'n'},
    {"write_tree_matrix",  no_argument,       NULL, 'g'},
    {"num_queries",        required_argument, NULL, 'q'},
    {"prop_query",         required_argument, NULL, 'f'},
    {NULL, 0, NULL, 0}
  };

  char shortopts[] = "s:j:c:r:ngq:f:";
  
  int optional;
  while((optional = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch (optional) {
      case 's': past_sampling = atof(optarg); break;
      case 'c': cons_path = optarg; break;
      case 'r': root_path = optarg; break;
      case 'n': write_newick_tree_to_file = true; break;
      case 'g': write_tree_to_file = true; break;
      case 'q': num_queries = atoi(optarg); break;
      case 'f': query_fraction = atof(optarg); break;
      default: die("Unknown option: %c", optional);
    }
  }

  if(argc - optind != 6) print_usage();
  
  json_parameters_path = argv[optind];
  printf("%s\n", json_parameters_path);
  N = atoi(argv[optind+1]);
  M = atoi(argv[optind+2]);
  
  if(M > N)
    die("Sample size larger than the infected population size: [M=%i, N=%i]", M, N);

  if(past_sampling < 0) past_sampling = (Decimal) M/N;
  if(past_sampling > 1) 
    die("Past sampling proportion is greater than 1: [M=%i, N=%i]", M, N);

  lambda = atof(argv[optind+3]);
  mu_tree = atof(argv[optind+4]);

  out_dir = argv[optind+5];

  if(!mkpath(out_dir, 0755)) die("Cannot create output dir: %s", out_dir);
  
  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char date[100];
  sprintf(date, "%d_%d_%d_%d_%d_%d", timeinfo->tm_year+1900, 
          timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
          timeinfo->tm_min, timeinfo->tm_sec);

  char output_data_path[PATH_MAX+1];
  char newick_data_path[PATH_MAX+1];
  char summary_path[PATH_MAX+1];
  char json_summary_path[PATH_MAX+1];
  char simulated_refs_path[PATH_MAX+1];
  char simulated_root_path[PATH_MAX+1];
  char simulated_queries_path[PATH_MAX+1];
  char hla_query_path[PATH_MAX+1];

  sprintf(output_data_path, "%s/%s_tree.txt", out_dir, date);
  sprintf(newick_data_path, "%s/%s_newick.tre", out_dir, date);
  sprintf(summary_path, "%s/%s_run_summary.txt", out_dir, date);
  sprintf(json_summary_path, "%s/%s_summary.json", out_dir, date);
  sprintf(simulated_refs_path, "%s/%s_simulated_birth_death_refs.fasta", out_dir, date);
  sprintf(simulated_root_path, "%s/%s_simulated_birth_death_root.fasta", out_dir, date);
  sprintf(simulated_queries_path, "%s/%s_simulated_birth_death_queries.fasta", out_dir, date);
  sprintf(hla_query_path, "%s/%s_simulated_hla_queries_birth_death.csv", out_dir, date);
  
  while(futil_file_exists(output_data_path))
  {
    printf("This filename already exists, wait for a second to change the folder name.\n");
    sleep(1);
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(date, "%d_%d_%d_%d_%d_%d\n", timeinfo->tm_year+1900, 
            timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);

    sprintf(output_data_path, "%s/%s_tree.txt", out_dir, date);
    sprintf(newick_data_path, "%s/%s_newick.tre", out_dir, date);
    sprintf(summary_path, "%s/%s_summary.txt", out_dir, date);
    sprintf(json_summary_path, "%s/%s_summary_dgen.json", out_dir, date);
    sprintf(simulated_refs_path, "%s/%s_simulated_birth_death_refs.fasta", out_dir, date);
    sprintf(simulated_root_path, "%s/%s_simulated_birth_death_root.fasta", out_dir, date);
    sprintf(simulated_queries_path, "%s/%s_simulated_birth_death_queries.fasta", out_dir, date);
    sprintf(hla_query_path, "%s/%s_simulated_hla_queries_birth_death.csv", out_dir, date);  
  }
  
  if(write_tree_to_file == true) printf("Writing tree information to: %s\n", output_data_path);
  if(write_newick_tree_to_file == true) printf("Writing Newick tree information to: %s\n", newick_data_path);
  printf("Writing summary information to: %s\n", summary_path);
  printf("Writing passed .json information to: %s\n", summary_path);
  printf("Writing simulated reference sequences to: %s\n", simulated_refs_path);
  printf("Writing simulated root sequence to: %s\n", simulated_root_path);
  printf("Writing simulated query sequences to: %s\n", simulated_queries_path);
  printf("Writing simulated HLA query data to: %s\n", hla_query_path);
  
  if(write_tree_to_file == true) tree_data_file = fopen(output_data_path, "w");
  if(write_newick_tree_to_file == true) newick_tree_data_file = fopen(newick_data_path, "w");
  simulated_refs_file = fopen(simulated_refs_path, "w");
  simulated_root_file = fopen(simulated_root_path, "w");
  simulated_queries_file = fopen(simulated_queries_path, "w");
  summary_file = fopen(summary_path, "w");
  json_summary_file = fopen(json_summary_path, "w");
  hla_query_file = fopen(hla_query_path, "w");

  if(tree_data_file == NULL && write_tree_to_file == true)
    die("Cannot open output file: %s", output_data_path);
  if(newick_tree_data_file == NULL && write_newick_tree_to_file == true)
    die("Cannot open output file: %s", newick_data_path);
  if(summary_file == NULL) die("Cannot open output file: %s", summary_path);
  if(json_summary_file == NULL) die("Cannot open output file: %s", json_summary_path);
  if(hla_query_file == NULL) die("Cannot open output file: %s", hla_query_path);  
  if(simulated_refs_file == NULL) die("Cannot open output file: %s", simulated_refs_path);
  if(simulated_root_file == NULL) die("Cannot open output file: %s", simulated_root_path);
  if(simulated_queries_file == NULL) die("Cannot open output file: %s", simulated_queries_path);

  // Add this information to the summary file.
  fprintf(summary_file, "Total lineages at present: %i\n"
                        "Sampled lineages at present: %i\n"
                        "Historical sampling proportion: "DECPRINT"\n"
                        "Birth rate in the genealogy (lambda): "DECPRINT"\n"
                        "Death rate in the genealogy (mu_tree): "DECPRINT"\n", 
          N, M, past_sampling, lambda, mu_tree);
}

static void memory_allocation(int num_cons, int num_root, int codon_sequence_length,
                              int max_internal_node_index,
                              int max_leaf_node_index, int total_nodes,
                              int ploidy, int n_genes, int total_n_HLA,
                              int leaf_node_index)
{
  // Approximate the total amount of memory allocation.
  size_t total_memory_alloc = 0;
  // cons_seqs.
  if(num_cons > 0) total_memory_alloc += num_cons * codon_sequence_length * 3 * sizeof(char);
  // root_seqs.
  if(num_root > 0) total_memory_alloc += num_root * codon_sequence_length * 3 * sizeof(char);
  // internal_node_times.
  total_memory_alloc += max_internal_node_index * sizeof(Decimal);
  // leaf_node_times.
  total_memory_alloc += max_leaf_node_index * sizeof(Decimal);
  // tree.
  total_memory_alloc += (total_nodes+1) * sizeof(Tree);
  // HLAs_in_tree.
  total_memory_alloc += (total_nodes+1)*ploidy*n_genes * sizeof(int);
  // codon_sequence_matrix.
  total_memory_alloc += total_nodes * sizeof(int[codon_sequence_length]);
  // HLA_selection_profiles.
  total_memory_alloc += total_n_HLA * sizeof(Decimal[codon_sequence_length]);
  // hla_types.
  total_memory_alloc += leaf_node_index * sizeof(int[total_n_HLA]);
  
  printf("Memory allocation before freeing sequence allocation: %zu bytes = %zu GB.\n", 
         total_memory_alloc,total_memory_alloc/(1<<30));
  // (1<<30) is 2^30 in binary.
  if(num_cons > 0) total_memory_alloc -= num_cons * codon_sequence_length * 3 * sizeof(char);
  if(num_root > 0) total_memory_alloc -= num_root * codon_sequence_length * 3 * sizeof(char);

  printf("Total memory allocation: %zu bytes = %zu GB.\n",total_memory_alloc,
         total_memory_alloc/(1<<30));
}

int main(int argc, char **argv)
{
  init_rand();
  setup_gsl_dgen();

  dgen_parse_cmdline(argc, argv);

  int i, j, h, p, k, c;
  
  int cons_cap, num_cons = -1, root_cap, num_root = -1;
  char **cons_seqs, **root_seqs;

  int internal_node_index = 0, leaf_node_index = M;
  int max_internal_node_index = 2000000;
  int max_leaf_node_index = 2000000;
  int sum_seen_or_unseen = 0;

  int ploidy, codon_sequence_length = -1, n_genes, total_n_HLA;
  Decimal mu;

  if(json_parameters_path == NULL) die("No .json file passed.");

  load_lengths_for_simulation_from_json(json_parameters_path, &kappa, &mu,
                                        &codon_sequence_length, &total_n_HLA, &ploidy, &n_genes);

  if(ploidy < 1) die("Ploidy is less than 1.");
  if(n_genes < 1) die("Number of genes is less than 1.");
  
  if(cons_path != NULL) {
    printf("Loading sequences to obtain consensus...\n");
    num_cons = load_seqs(cons_path, &cons_seqs, &cons_cap);
    assert(num_cons > 0);
    printf("Loaded %i sequences to determine consensus.\n", num_cons);
    
    codon_sequence_length = strlen(cons_seqs[0]);
    printf("Codon_sequence_length: %i\n",codon_sequence_length/3);
    
    if(codon_sequence_length % 3 != 0)
      die("Sequences contain partial codons [%i mod 3 != 0].", codon_sequence_length);

    for(c = 0; c < num_cons; c++) {
      if((int) strlen(cons_seqs[c]) != codon_sequence_length) {
        die("Sequences from which to derive the consensus sequence aren't all "
            "the same length.");
      }
    }
    codon_sequence_length = codon_sequence_length/3;
  }

  if(root_path != NULL) {
    printf("Loading sequences to obtain root...\n");
    num_root = load_seqs(root_path, &root_seqs, &root_cap);
    printf("Loaded %i sequences to determine root.\n", num_root);
    
    if(cons_path == NULL)
      die("Did not pass a file to find the consensus sequence.");
    
    if((int) (strlen(root_seqs[0])/3) != codon_sequence_length)
      die("Sequences used to determine the root are different lengths to those used for the consensus.");

    for(c = 0; c < num_root; c++) {
      if((int) strlen(root_seqs[c]) != 3*codon_sequence_length) {
        die("Sequences from which to derive the root sequence aren't all "
            "the same length.");
      }
    }
  }

  Decimal *internal_node_times = my_malloc(max_internal_node_index * sizeof(Decimal) , __FILE__, __LINE__);
  Decimal *leaf_node_times = my_malloc(max_leaf_node_index * sizeof(Decimal), __FILE__, __LINE__);
  int *seen_or_unseen = my_malloc(max_internal_node_index * sizeof(int), __FILE__, __LINE__);

  birth_death_simulation_backwards(max_internal_node_index, max_leaf_node_index,
                                   internal_node_times, 
                                   leaf_node_times,
                                   &internal_node_index, &leaf_node_index,
                                   seen_or_unseen,
                                   N, M, lambda, mu_tree, past_sampling);

  for(i = 0; i < internal_node_index; i++) sum_seen_or_unseen += seen_or_unseen[i];
  
  int total_nodes = (2 * leaf_node_index) - 1 + internal_node_index - sum_seen_or_unseen;
  Tree *tree = my_malloc((total_nodes+1) * sizeof(Tree), __FILE__, __LINE__);
  // Now malloc the memory that this points to.
  int *HLAs_in_tree = my_malloc((total_nodes+1) * ploidy * n_genes * sizeof(int), __FILE__, __LINE__);

  for(i = 0; i < total_nodes; i++) 
    tree[i].HLAs = &HLAs_in_tree[i * ploidy * n_genes];

  construct_birth_death_tree(leaf_node_index, internal_node_index,
                             leaf_node_times, internal_node_times,
                             M, seen_or_unseen,
                             tree);

  // Reverse the direction that time is measured in the tree.
  // DEV: Don't need to do this, waste of computation - sort.
  // DEV: The parent times are wrong when there are unseen nodes.
  for(i = 0; i < total_nodes; i++)
    tree[i].node_time = tree[total_nodes-1].node_time - tree[i].node_time;
   
  int root_node = tree[total_nodes-1].node;
  
  if(write_newick_tree_to_file == true) {
    write_newick_tree(newick_tree_data_file, tree, root_node, 1);
    fclose(newick_tree_data_file);
  }

  Decimal S_portion[NUM_CODONS];
  Decimal NS_portion[NUM_CODONS];

  for(c = 0; c < NUM_CODONS; c++) {
    S_portion[c] = kappa * beta_S[c] + beta_V[c];
    NS_portion[c] = kappa * alpha_S[c] + alpha_V[c];
  }
  
  int n_HLA[n_genes];

  printf("Total number of HLA types: %i.\n", total_n_HLA);

  Decimal HLA_prevalences[total_n_HLA];
  int wildtype_sequence[codon_sequence_length];
  Decimal *R = my_malloc(codon_sequence_length * sizeof(Decimal), __FILE__, __LINE__);
  Decimal *omega = my_malloc(codon_sequence_length * sizeof(Decimal), __FILE__, __LINE__);
  Decimal *reversion_selection = my_malloc(codon_sequence_length * sizeof(Decimal), __FILE__, __LINE__);

  memory_allocation(num_cons, num_root, codon_sequence_length,
                    max_internal_node_index, max_leaf_node_index, 
                    total_nodes, ploidy, n_genes, total_n_HLA,
                    leaf_node_index);

  int (*codon_sequence_matrix)[codon_sequence_length] = my_malloc(total_nodes *
                                                                  sizeof(int[codon_sequence_length]),
                                                                  __FILE__, __LINE__);
  Decimal (*HLA_selection_profiles)[codon_sequence_length] = my_malloc(total_n_HLA * sizeof(Decimal[codon_sequence_length]),
                                                                       __FILE__, __LINE__);

  load_parameters_for_simulation_from_json(json_parameters_path, codon_sequence_length,
                                           omega, R, reversion_selection, total_n_HLA,
                                           n_genes, n_HLA, HLA_prevalences,
                                           HLA_selection_profiles);
  
  Decimal sum_check;
  for(i = 0, k = 0; i < n_genes; i++) {
    sum_check = 0;
    for(h = 0; h < n_HLA[i]; h++, k++) {
      sum_check += HLA_prevalences[k];
    }
    if(sum_check > 1.00001 || sum_check < 0.9999) die("HLA prevalences for gene %i do not sum to 1\n", i+1);
  }
  
  if(cons_path != NULL) {
    printf("Mapping gaps to consensus...\n");
    // Set the consensus sequence - the consensus of the optional sequence file 
    // that is passed.
    char wildtype_sequence_dummy[3*codon_sequence_length+1];
    generate_consensus(cons_seqs, num_cons, 3*codon_sequence_length, wildtype_sequence_dummy);
    printf("Wildtype sequence:\n%s\n", wildtype_sequence_dummy);
    // By default, set the root as the wildtype sequence.
    for(i = 0; i < codon_sequence_length; i++)
      wildtype_sequence[i] = (int) amino_to_code(wildtype_sequence_dummy+i*3);

    if(root_path == NULL) {
      for(i = 0; i < codon_sequence_length; i++)
        codon_sequence_matrix[root_node][i] = wildtype_sequence[i];
    } else {
      printf("Mapping gaps to root...\n");
      char root_sequence_dummy[3*codon_sequence_length+1];
      generate_consensus(root_seqs, num_root, 3*codon_sequence_length, root_sequence_dummy);
      printf("Root sequence:\n%s\n", root_sequence_dummy);
      
      for(i = 0; i < codon_sequence_length; i++)
        codon_sequence_matrix[root_node][i] = (int) amino_to_code(root_sequence_dummy+i*3);
      printf("Number of root sequences: %i.\n", num_root);
      for(c = 0; c < num_root; c++) free(root_seqs[c]);
      free(root_seqs);
    }
    printf("Number of consensus sequences: %i.\n", num_cons);
    for(c = 0; c < num_cons; c++) free(cons_seqs[c]);
    free(cons_seqs);
  
  } else {
    for(i = 0; i < codon_sequence_length; i++) {
      // Sample the root sequence according to the HIV codon usage information.
      codon_sequence_matrix[root_node][i] = discrete_sampling_dist(NUM_CODONS, prior_C1);
      // As default, set the root node to the consensus sequence.  
      wildtype_sequence[i] = codon_sequence_matrix[root_node][i];
    }
  }
  
  // No matter what is read in, there is no recombination simulated - so make sure it's set to 0.
  for(i = 0; i < codon_sequence_length; i++) R[i] = 0;

  write_summary_json(json_summary_file,
                     mu, codon_sequence_length, ploidy,
                     n_genes, n_HLA, total_n_HLA,
                     HLA_prevalences,
                     omega, R, reversion_selection,
                     HLA_selection_profiles);

  free(R);
  
  fprintf(simulated_root_file, ">root_sequence\n");
  for(i = 0; i < codon_sequence_length; i++)
    fprintf(simulated_root_file, "%s", code_to_char(codon_sequence_matrix[root_node][i]));
  fprintf(simulated_root_file, "\n");

  int root_HLA[ploidy * n_genes];
  int cumulative_n_HLA = 0;

  for(i = 0, k = 0; i < n_genes; i++) {
    for(p = 0; p < ploidy; p++, k++) {
      root_HLA[k] = cumulative_n_HLA + 
                    discrete_sampling_dist(n_HLA[i], &HLA_prevalences[cumulative_n_HLA]);
      tree[root_node].HLAs[k] = root_HLA[k];
    }
    cumulative_n_HLA = cumulative_n_HLA + n_HLA[i];
  }

  printf("Passing HLA information...\n");
  pass_HLA(ploidy, n_genes, root_node,
                           tree, leaf_node_index, total_n_HLA,
                           n_HLA, HLA_prevalences);
  printf("Passed HLA information\n");
  
  // printf("Printing the tree\n");
  // for(i = 0; i < total_nodes; i++) {
  //   printf("%i %i %i "DECPRINT" %i", tree[i].node, tree[i].daughter_nodes[0],
  //          tree[i].daughter_nodes[1], tree[i].node_time,
  //          tree[i].seen_or_unseen);
  //   for(j = 0; j < (ploidy * n_genes); j++) {
  //     printf(" %i", tree[i].HLAs[j]);
  //   }
  //   printf("\n");
  // }

  if(write_tree_to_file == true) {
    write_tree(tree_data_file, tree, root_node, ploidy, n_genes);
    fclose(tree_data_file);
  }

  printf("Passing sequence information...\n");
  
  pass_codon_sequence_change(codon_sequence_length, ploidy,
                             n_genes, total_n_HLA, 
                             root_node, mu,
                             codon_sequence_matrix,
                             tree,
                             leaf_node_index,
                             S_portion, NS_portion,
                             HLA_selection_profiles,
                             wildtype_sequence,
                             omega, reversion_selection);

  printf("Passed sequence information\n"
         "Now generating .fasta files of reference and query sequences, and\n"
         "a .csv file of the HLA information associated to the query sequences.\n");

  if(num_queries < 0) {
    // Set the number of query sequences.
    num_queries = (int) (query_fraction * leaf_node_index);
    printf("Number of queries: %i.\n", num_queries);
  } else {
    printf("Number of queries: %i.\n", num_queries);
  }

  if(num_queries > leaf_node_index) die("Number of query sequences larger than the number of leaves");
  int *all_sequences = my_malloc(leaf_node_index * sizeof(int), __FILE__, __LINE__);
  int num_refs = leaf_node_index - num_queries;

  for(i = 0; i < leaf_node_index; i++) all_sequences[i] = i;

  save_simulated_ref_and_query_fasta(num_queries, num_refs, leaf_node_index,
                                     all_sequences, codon_sequence_length, codon_sequence_matrix,
                                     tree, ploidy, n_genes);

  // Now save the hla types to a .csv file.
  fprintf(hla_query_file, "\"\",");
  for(h = 0; h < total_n_HLA-1; h++) fprintf(hla_query_file, "\"%i\",", h+1);
  fprintf(hla_query_file, "\"%i\"\n", total_n_HLA);
  
  // Write the HLA types of the leaves to a file.
  int (*hla_types)[total_n_HLA] = my_malloc(leaf_node_index * sizeof(int[total_n_HLA]), __FILE__, __LINE__);

  for(i = 0; i < leaf_node_index; i++)
  {
    for(h = 0; h < total_n_HLA; h++)
      hla_types[i][h] = 0;
    for(j = 0; j < (n_genes * ploidy); j++)
      hla_types[i][tree[i].HLAs[j]] = 1;
  }

  // Write the query HLA types to a .csv file.
  for(i = num_refs; i < leaf_node_index; i++)
  {
    fprintf(hla_query_file,"\"simulated_seq_%i_HLA", all_sequences[i]+1);
    for(h = 0; h < (ploidy * n_genes); h++) fprintf(hla_query_file, "_%i", tree[all_sequences[i]].HLAs[h]);
    fprintf(hla_query_file, "\"");
    for(h = 0; h < total_n_HLA; h++) {
      fprintf(hla_query_file, ",%i", hla_types[all_sequences[i]][h]);
    }
    fprintf(hla_query_file, "\n");
  }

  free(hla_types); 
  free(internal_node_times);
  free(leaf_node_times);
  free(seen_or_unseen);
  free(codon_sequence_matrix);
  free(HLA_selection_profiles);
  free(all_sequences);
  free(omega);
  free(reversion_selection);

  free(tree[0].HLAs);
  free(tree);

  fclose(summary_file);
  fclose(json_summary_file);
  fclose(simulated_refs_file);
  fclose(simulated_root_file);
  fclose(simulated_queries_file);
  fclose(hla_query_file);

  clearup_gsl_dgen();
  return EXIT_SUCCESS;
}
