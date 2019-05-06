#include <stdio.h>
#include <stdbool.h>

#include "write_tree.h"
#include "constants.h"
#include "birth_death_gen.h"

// write out the newick notation of a tree
void write_newick_node(FILE *out, Tree *tree, int node, int depth,
                       bool oneline, Decimal *parent_time)
{
  int i;
  if (tree[node].daughter_nodes[0] < 0) { // Leaf node.
    
    if (!oneline){
      for (i=0; i<depth; i++){
        fprintf(out, "  ");
      }
    }
    fprintf(out, "%i:"DECPRINT, tree[node].node, tree[node].node_time - parent_time[node]);
  
  } else { // Internal node.
    
    if (tree[node].daughter_nodes[1] != -1){
      if (oneline){
        fprintf(out, "(");
      } else {
        
        for (i=0; i<depth; i++){
          fprintf(out, "  ");
        }
        fprintf(out, "(\n");

      }
        
      // Daughter 1.
      parent_time[tree[node].daughter_nodes[0]] = tree[node].node_time;
      write_newick_node(out, tree, tree[node].daughter_nodes[0], depth+1, oneline, parent_time);
      
      if (oneline){
        fprintf(out, ",");
      } else {
        fprintf(out, ",\n");
      }
        
      // Daughter 2.
      parent_time[tree[node].daughter_nodes[1]] = tree[node].node_time;
      write_newick_node(out, tree, tree[node].daughter_nodes[1], depth+1, oneline, parent_time);
        
      if (!oneline) {
        fprintf(out, "\n");            
        for (i=0; i<depth; i++){
          fprintf(out, "  ");
        }
      }
      fprintf(out, ")");
        
      if (depth > 0){
        fprintf(out, "%i:"DECPRINT, tree[node].node, tree[node].node_time - parent_time[node]);
      } else {
        fprintf(out, "%i", tree[node].node);
      }
    } else {
      parent_time[tree[node].daughter_nodes[0]] = parent_time[tree[node].node];
      // Then skip to the daughter of the first node.
      write_newick_node(out, tree, tree[node].daughter_nodes[0], depth, oneline, parent_time);
    }
  }
}

void find_branch_length(Tree *tree, int node,
                        Decimal *parent_time, Decimal *branch_lengths,
                        int *nodes_passed_through)
{

  if (tree[node].daughter_nodes[0] < 0) { // Leaf node.
    branch_lengths[tree[node].node] =  -(tree[node].node_time - parent_time[node]);
  } else { // Internal node.
    
    if (tree[node].daughter_nodes[1] != -1){
        
      // Daughter 1.
      parent_time[tree[node].daughter_nodes[0]] = tree[node].node_time;
      nodes_passed_through[tree[node].daughter_nodes[0]] = 1;
      find_branch_length(tree, tree[node].daughter_nodes[0], parent_time, branch_lengths,
                         nodes_passed_through);
        
      // Daughter 2.
      parent_time[tree[node].daughter_nodes[1]] = tree[node].node_time;
      nodes_passed_through[tree[node].daughter_nodes[1]] = 1;
      find_branch_length(tree, tree[node].daughter_nodes[1], parent_time, branch_lengths,
                         nodes_passed_through);
        
      branch_lengths[tree[node].node] = -(tree[node].node_time - parent_time[node]);

    } else {
      parent_time[tree[node].daughter_nodes[0]] = parent_time[tree[node].node];
      nodes_passed_through[tree[node].daughter_nodes[0]] = nodes_passed_through[tree[node].node] + 1;
      
      // Then skip to the daughter of the first node.
      find_branch_length(tree, tree[node].daughter_nodes[0], parent_time, branch_lengths,
                         nodes_passed_through);
    }
  }
}

// write out the newick notation of a tree
void write_newick_tree(FILE *out, Tree *tree, int root_node, bool oneline)
{
  Decimal parent_time[root_node];
  write_newick_node(out, tree, tree[root_node].node, 0, oneline, parent_time);
  if (oneline)
    fprintf(out, ";");
  else
    fprintf(out, ";\n");
}

void write_tree(FILE *out, Tree *tree, int root_node, int ploidy, int n_genes){
  int i,h;
  fprintf(out, "node, daughter_node_1, daughter_node_2, node_time,"
               " seen_or_unseen, ");
  for (h=0; h<(ploidy*n_genes)-1; h++)
    fprintf(out, "HLA[%i], ", h);
  fprintf(out, "HLA[%i]\n", h);

  for (i=0; i<root_node; i++){
    fprintf(out, "%i, %i, %i, "DECPRINT", %i,", tree[i].node,
            tree[i].daughter_nodes[0], tree[i].daughter_nodes[1],
            tree[i].node_time, tree[i].seen_or_unseen);
    for (h=0; h<(ploidy*n_genes)-1; h++)
      fprintf(out, "%i, ", tree[i].HLAs[h]);
    fprintf(out, "%i\n", tree[i].HLAs[h]);
  }
}
