#ifndef BIRTH_DEATH_GEN_H_
#define BIRTH_DEATH_GEN_H_

#include <stdio.h>
#include "constants.h"

typedef struct
{
  int node, daughter_nodes[2];
  Decimal node_time;
  int seen_or_unseen;
  int *HLAs;
} Tree;

void setup_gsl_dgen();

void clearup_gsl_dgen();

void birth_death_simulation_backwards(int max_internal_node_index,
	                                  int max_leaf_node_index,
	                                  Decimal *internal_node_times, 
	                                  Decimal *leaf_node_times,
	                                  int *internal_node_index,
	                                  int *leaf_node_index,
	                                  int *seen_or_unseen,
	                                  int N, int M, Decimal lambda, Decimal mu, 
	                                  Decimal past_sampling);


void construct_birth_death_tree(int leaf_node_index, int internal_node_index,
  	                            Decimal leaf_node_times[leaf_node_index],
  	                            Decimal internal_node_times[internal_node_index],
  	                            int M, int *seen_or_unseen,
                                Tree *tree);

#endif /* BIRTH_DEATH_GEN_H_ */
