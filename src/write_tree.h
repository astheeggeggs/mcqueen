#ifndef WRITE_TREE_
#define WRITE_TREE_ 

#include "birth_death_gen.h"

void write_newick_node(FILE *out, Tree *tree, int node, int depth,
                       bool oneline, Decimal *parent_time);

void find_branch_length(Tree *tree, int node,
                        Decimal *parent_time, Decimal *branch_lengths,
                        int *nodes_passed_through);

void write_newick_tree(FILE *out, Tree *tree, int root_node, bool oneline);

void write_tree(FILE *out, Tree *tree, int root_node, int ploidy, int n_genes);

#endif /* WRITE_TREE_ */
