#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <unistd.h>

#include "constants.h"
#include "util.h"
#include "birth_death_gen.h"
#include "assert.h"

gsl_rng *gsl_r;

void setup_gsl_dgen()
{
  gsl_rng_env_setup();
  const gsl_rng_type *T = gsl_rng_default;
  gsl_r = gsl_rng_alloc(T);
  Decimal seed_rng = time (NULL) * getpid();
  gsl_rng_set (gsl_r, seed_rng);
}

void clearup_gsl_dgen()
{
  gsl_rng_free(gsl_r);
}

void birth_death_simulation_backwards(int max_internal_node_index,
                                      int max_leaf_node_index,
                                      Decimal *internal_node_times, 
                                      Decimal *leaf_node_times,
                                      int *internal_node_index,
                                      int *leaf_node_index,
                                      int *seen_or_unseen,
                                      int N, int M, Decimal lambda, Decimal mu, 
                                      Decimal past_sampling)
{
  int l;
  Decimal p_death = mu / (lambda + mu);
  Decimal current_ratio;
  Decimal current_time = 0;
  Decimal choose_birth_death, choose_sampled_coal, choose_sampled_death;
  *leaf_node_index = M;
  Decimal exp_parameter;

  // Pre-allocate present leaf times at 0.
  for(l = 0; l < M; l++) leaf_node_times[l] = 0;
 
  while(M > 0) {
    exp_parameter = (Decimal) 1 / (N * (lambda + mu));
    current_time += gsl_ran_exponential(gsl_r, exp_parameter);
    choose_birth_death = rand_lim(1);
    if(choose_birth_death > p_death) { // Birth.;
      choose_sampled_coal = rand_lim(1);
      current_ratio = (Decimal) M / N;
      if(choose_sampled_coal < (current_ratio * current_ratio)) {
        // Coalescence event in the sampled tree.
        M = M - 1;
        internal_node_times[*internal_node_index] = current_time;
        seen_or_unseen[*internal_node_index] = 1;
        *internal_node_index += 1;
      } else if(choose_sampled_coal < current_ratio) {
        // Record the time;
        // Coalescence between an unseen lineage and a sampled lineage.
        internal_node_times[*internal_node_index] = current_time;
        seen_or_unseen[*internal_node_index] = 0;
        *internal_node_index += 1;
      }
      N = N - 1;
    } else { // Death.
      choose_sampled_death = rand_lim(1);
      if(choose_sampled_death < past_sampling) {
        // Record the time:
        // This leaf is sampled in the past.
        leaf_node_times[*leaf_node_index] = current_time;
        *leaf_node_index += 1;
        M = M + 1;
      }
      N = N + 1;
    }
    if(*leaf_node_index > max_leaf_node_index) die("Not enough space for leaves\n");
    if(*internal_node_index > max_internal_node_index) die("Not enough space for internal nodes\n");
  }
}

void construct_birth_death_tree(int leaf_node_index, int internal_node_index,
                                Decimal leaf_node_times[leaf_node_index],
                                Decimal internal_node_times[internal_node_index],
                                int M, int *seen_or_unseen,
                                Tree *tree)
                                                      
{
  Decimal current_time = 0;
  int leaf_nodes[leaf_node_index];
  int l, i, j;
  int parent_node = leaf_node_index-1;
  int number_of_lineages = M;
  int max_available_nodes = leaf_node_index;
  int *available_nodes = my_malloc(max_available_nodes * sizeof(int), __FILE__, __LINE__);
  int min_leaf_index = M;
  int unseen_node_pos;
  int unseen_node;
  int coalescing_nodes[2];

  assert(leaf_node_index > 0);
  assert(M > 0);
  assert(internal_node_index > 0);

  for(l = 0; l < leaf_node_index; l++) leaf_nodes[l] = l;
  for(i = 0; i < M; i++) available_nodes[i] = i;

  for(l = 0; l < leaf_node_index; l++) {
    tree[l].node = leaf_nodes[l];
    tree[l].daughter_nodes[0] = -1;
    tree[l].daughter_nodes[1] = -1;
    tree[l].node_time = leaf_node_times[l];
  }
  // Set the final root node as unseen for the algorithm below to finish.
  seen_or_unseen[internal_node_index -1] = 0;
  
  for(i = 0; i < internal_node_index; i++) {
    // Do I add some leaf nodes to the collection of available nodes?
    for(j = min_leaf_index; j < leaf_node_index; j++) {
      if(leaf_node_times[j] <= internal_node_times[i] && 
          leaf_node_times[j] > current_time) {
        // Node at index j is available.
        available_nodes[number_of_lineages] = leaf_nodes[j];
        number_of_lineages++;
      }
    }

    if(seen_or_unseen[i] == 0) {
      // Unseen transmission
      unseen_node_pos = rand_lim(number_of_lineages);
      unseen_node = available_nodes[unseen_node_pos];
      parent_node++;
      available_nodes[unseen_node_pos] = parent_node;

      tree[i+leaf_node_index].node = parent_node;
      tree[i+leaf_node_index].daughter_nodes[0] = unseen_node;
      tree[i+leaf_node_index].daughter_nodes[1] = -1;
      tree[i+leaf_node_index].node_time = internal_node_times[i];

    } else {
      // Seen transmission
      // Sample without replacement.
      knuth_sample(2, number_of_lineages, available_nodes);

      coalescing_nodes[0] = available_nodes[number_of_lineages-1];
      coalescing_nodes[1] = available_nodes[number_of_lineages-2];

      parent_node++;
      available_nodes[number_of_lineages-2] = parent_node;

      // Need to rearrange the available nodes - shuffle stuff.
      tree[i+leaf_node_index].node = parent_node;
      tree[i+leaf_node_index].daughter_nodes[0] = coalescing_nodes[0];
      tree[i+leaf_node_index].daughter_nodes[1] = coalescing_nodes[1];
      tree[i+leaf_node_index].node_time = internal_node_times[i];

      number_of_lineages--;
    }
    current_time = internal_node_times[i];
  }

  // Add in the seen or unseen vector as a column in the matrix.
  for(i=0; i<leaf_node_index; i++) tree[i].seen_or_unseen = -1;

  for(i=leaf_node_index; i<(leaf_node_index + internal_node_index); i++)
    tree[i].seen_or_unseen = seen_or_unseen[i-leaf_node_index];

  seen_or_unseen[internal_node_index -1] = 1;

  free(available_nodes);
}
