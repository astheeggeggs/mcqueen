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
#include "birth_death_gen.h"
#include "write_tree.h"
#include "load_data.h"

Decimal past_sampling = -1;
Decimal query_fraction = 0.01;
int N, M, n_trees, num_queries = -1;
Decimal lambda, mu_tree;

FILE *summary_file = NULL;
const char *out_dir = NULL;

static void print_usage()
{
  fprintf(stderr,
          "usage: dscaling [options] <N> <M> <lambda> <mu_tree> <n_trees> <outdir>\n"
          "  options:\n"
          "    prop_past_sampling  -s <S> Set sampling proportion for historic sampling [default: M/N].\n"
          "    num_queries         -q <Q> provide the number of queries.\n"
          "    prop_query          -f <F> provide the query fraction.\n");
  exit(EXIT_FAILURE);
}

static void dgen_parse_cmdline(int argc, char **argv)
{
  // DEV: include the codon sequence length as an optional argument.
  // Also, all the other scalar parameters?
  // DEV: want to be able to save the tree in newick format too?
  //  - to make sure that everthing is correct.
  // DEV: could make it optional to write the HLA types at the internal nodes.
  // DEV: talk to Gil about generating the his approximation to the coalescent. 
  // with recombination - different trees as you track along the viral sequence.
  static struct option longopts[] =
  {
    {"prop_past_sampling", required_argument, NULL, 's'},
    {"num_queries",        required_argument, NULL, 'q'},
    {"prop_query",         required_argument, NULL, 'f'},
    {NULL, 0, NULL, 0}
  };

  char shortopts[] = "s:q:f:";
  
  int optional;
  while((optional = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1) {
    switch(optional) {
      case 's': past_sampling = atof(optarg); break;
      case 'q': num_queries = atoi(optarg); break;
      case 'f': query_fraction = atof(optarg); break;
      default: die("Unknown option: %c", optional);
    }
  }
 
  if(argc - optind != 6) print_usage();

  N = atoi(argv[optind]);
  M = atoi(argv[optind+1]);

  if(past_sampling < 0) past_sampling = (Decimal) M/N;
  if(past_sampling > 1) 
    die("Past sampling proportion is greater than 1: [M = %i, N = %i].", M, N);

  lambda = atof(argv[optind+2]);
  mu_tree = atof(argv[optind+3]);
  n_trees = atof(argv[optind+4]);

  out_dir = argv[optind+5];

  if(!mkpath(out_dir, 0755)) die("Cannot create output dir: %s.", out_dir);
  
  time_t rawtime;
  struct tm * timeinfo;

  time(&rawtime);
  timeinfo = localtime(&rawtime);
  char date[100];
  sprintf(date, "%d_%d_%d_%d_%d_%d", timeinfo->tm_year+1900, 
          timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
          timeinfo->tm_min, timeinfo->tm_sec);

  char summary_path[PATH_MAX+1];

  sprintf(summary_path, "%s/%s_summary.txt", out_dir, date);
  
  while(futil_file_exists(summary_path))
  {
    printf("This filename already exists, wait for a second to change the folder name.\n");
    sleep(1);
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    sprintf(date, "%d_%d_%d_%d_%d_%d\n", timeinfo->tm_year+1900, 
            timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);

    sprintf(summary_path, "%s/%s_summary.txt", out_dir, date);
  }
  
  printf("Writing summary information to: %s.\n", summary_path);
  
  summary_file = fopen(summary_path, "w");

  if(summary_file == NULL) die("Cannot open output file: %s", summary_path);

  // Add this information to the summary file.
  fprintf(summary_file, "Total lineages at present: %i\n"
                        "Sampled lineages at present: %i\n"
                        "Historical sampling proportion: "DECPRINT"\n"
                        "Birth rate in the genealogy (lambda): "DECPRINT"\n"
                        "Death rate in the genealogy (mu_tree): "DECPRINT"\n", 
          N, M, past_sampling, lambda, mu_tree);

}

int main(int argc, char **argv)
{
  init_rand();
  setup_gsl_dgen();

  dgen_parse_cmdline(argc, argv);
  
  int i, j, k;

  int internal_node_index, leaf_node_index;
  int max_internal_node_index = 2000000;
  int max_leaf_node_index = 2000000;
  int sum_seen_or_unseen;

  // Note: codon_sequence_length is overwritten by the length of the sequences in the sequence file
  // passed to determine the consensus.

  Decimal *internal_node_times = my_malloc(max_internal_node_index * sizeof(Decimal), __FILE__, __LINE__);
  Decimal *leaf_node_times = my_malloc(max_leaf_node_index * sizeof(Decimal), __FILE__, __LINE__);
  int *seen_or_unseen = my_malloc(max_internal_node_index * sizeof(int), __FILE__, __LINE__);
  Decimal *branch_lengths = my_malloc(max_leaf_node_index * sizeof(Decimal), __FILE__, __LINE__);
  Decimal *parent_time = my_malloc(max_leaf_node_index * sizeof(Decimal), __FILE__, __LINE__);
  int *nodes_passed_through = my_malloc(max_leaf_node_index* sizeof(int), __FILE__, __LINE__);

  k = 0;

  for(j = 0; j < n_trees; j++)
  {
    internal_node_index = 0;
    leaf_node_index = M;
    sum_seen_or_unseen = 0;
    birth_death_simulation_backwards(max_internal_node_index, max_leaf_node_index,
  	                                 internal_node_times, 
  	                                 leaf_node_times,
                                     &internal_node_index, &leaf_node_index,
                                     seen_or_unseen,
                                     N, M, lambda, mu_tree, past_sampling);

    for(i = 0; i < internal_node_index; i++) sum_seen_or_unseen += seen_or_unseen[i];
    
    int total_nodes = (2 * leaf_node_index) - 1 + internal_node_index - sum_seen_or_unseen;
    Tree *tree = my_malloc((total_nodes+1) * sizeof(Tree), __FILE__, __LINE__);

    construct_birth_death_tree(leaf_node_index, internal_node_index,
  	                           leaf_node_times, internal_node_times,
  	                           M, seen_or_unseen,
  	                           tree);

    int root_node = tree[total_nodes-1].node;

    find_branch_length(tree, tree[root_node].node, parent_time, &branch_lengths[k],
                       &nodes_passed_through[k]);
    k += leaf_node_index;

    free(tree);
    printf("%i trees checked.\n", j+1);
  }
  
  Decimal mean_nodes_passed_through = 0;
  for(j = 0; j < k; j++) mean_nodes_passed_through += nodes_passed_through[j];

  mean_nodes_passed_through = 2 * (Decimal) mean_nodes_passed_through / k;
  printf("mean number of nodes passed through: "DECPRINT"\n", mean_nodes_passed_through);
  fprintf(summary_file, "mean number of nodes passed through: "DECPRINT"\n", mean_nodes_passed_through);
  
  Decimal total_branch = 0;
  for(i = 0; i < k; i++) total_branch += branch_lengths[i];
  total_branch = total_branch / k;

  printf("mean branch length: "DECPRINT".\n", total_branch);
  fprintf(summary_file, "mean branch length: "DECPRINT"\n", total_branch);
  
  free(internal_node_times);
  free(leaf_node_times);
  free(seen_or_unseen);
  free(branch_lengths);
  free(parent_time);
  free(nodes_passed_through);

  fclose(summary_file);
 
  clearup_gsl_dgen();

  return EXIT_SUCCESS;
}