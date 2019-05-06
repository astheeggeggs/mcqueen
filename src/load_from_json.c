#include "load_from_json.h"
#include "util.h"
#include "constants.h"
#include "selection_window.h"

#include "cJSON.h"
#include "string_buffer.h"
#include "write_json.h"
#include "amino_acids.h"

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <zlib.h>

// Definition to check existence of classes within the .json file.
#define json_assert(x,msg,...) ((x) ? (void)0 : call_die(__FILE__,__LINE__,msg,##__VA_ARGS__))

static void load_windows(cJSON *w, HLAWindowSet *hlaset, int num_sites)
{
  cJSON *start, *coeff;
  SelectionWindow tmp_window, *windows = hlaset->windows;
  int i, starti;

  for(i = 0; w != NULL; w = w->next, i++) {
    start = cJSON_GetObjectItem(w,"start");
    json_assert(start != NULL && start->type == cJSON_Number, "No valid 'start'.");
    coeff = cJSON_GetObjectItem(w,"coeff");
    json_assert(coeff != NULL && coeff->type == cJSON_Number, "No valid 'coeff'.");
    starti = start->valueint;

    printf(" %i start: %i coeff: %f\n", i, starti, coeff->valuedouble);

    json_assert(starti < num_sites, "Windows greater than num sites.");
    json_assert(windows[starti].start == starti, "Prog error: [%i vs %i].",
                windows[starti].start, starti);

    SWAP(windows[i], windows[starti], tmp_window);
    windows[i].coeff = coeff->valuedouble;
  }

  hlaset->num_windows = i;
  selec_windows_sort(hlaset);
  selec_windows_set_lengths(hlaset, num_sites);
}

void load_from_json(const char *path, Decimal *llk, Decimal *mu,
                    Decimal *R, Decimal *omega,
                    HLASelection *hla_selection,
                    int num_sites, int num_hla_types,
                    int num_query, int closest_n,
                    int closest_refs[num_query][closest_n])
{
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open file %s.", path);

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 8192); // 2^13 bytes = 8KB

  // Read in the whole file.
  while(strbuf_gzreadline(&sbuf, gz)) {}
  gzclose(gz);

  cJSON *root, *state, *statenum, *mu_obj, *llk_obj, *sites, *site, *obj;
  cJSON *hlas, *hla, *esc, *rev;
  int s, h, num_states;

  json_assert(sbuf.end > 0, "Empty json file.");

  root = cJSON_Parse(sbuf.b);
  state = cJSON_GetObjectItem(root,"chain");
  json_assert(state != NULL && state->type == cJSON_Array, "No 'chain' object.");
  json_assert(state->child != NULL, "No states in chain.");
  state = state->child;

  // Get the last state.
  for(num_states = 1; state->next != NULL; state = state->next, num_states++) {}
  printf("We got %i states - loading last one.\n", num_states);

  statenum = cJSON_GetObjectItem(state,"state");
  json_assert(statenum != NULL && statenum->type == cJSON_Number, "No 'state' object.");
  printf("Got state number: %i\n", statenum->valueint);

  mu_obj = cJSON_GetObjectItem(state,"mu");
  json_assert(mu_obj != NULL && mu_obj->type == cJSON_Number, "No valid 'mu' value.");
  *mu = mu_obj->valuedouble;
  
  llk_obj = cJSON_GetObjectItem(state,"llkhood");
  json_assert(llk_obj != NULL && llk_obj->type == cJSON_Number, "No valid 'llkhood' value.");
  *llk = llk_obj->valuedouble;

  sites = cJSON_GetObjectItem(state,"sites");
  json_assert(sites != NULL && sites->type == cJSON_Array && sites->child != NULL,
              "No 'sites' object.");
  site = sites->child;

  // Loop over sites.
  for(s = 0; site != NULL; site = site->next, s++)
  {
    obj = cJSON_GetObjectItem(site,"omega");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'omega'.");
    omega[s] = obj->valuedouble;
    obj = cJSON_GetObjectItem(site,"R");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'R'.");
    R[s] = obj->valuedouble;

    printf("omega[%i]:"DECPRINT"; R[%i]:"DECPRINT"\n", s, omega[s], s, R[s]);
  }

  json_assert(num_sites == s, "Mismatch in sites: [%i vs %i].", num_sites, s);

  hlas = cJSON_GetObjectItem(state,"hlas");
  json_assert(hlas != NULL && hlas->type == cJSON_Array && hlas->child != NULL,
              "No 'hlas' object");
  hla = hlas->child;

  for(h = 0; hla != NULL; hla = hla->next, h++)
  {
    obj = cJSON_GetObjectItem(hla,"id");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'id'.");

    printf("HLA: %i\n", obj->valueint);

    printf("esc:\n");
    esc = cJSON_GetObjectItem(hla,"esc");
    json_assert(esc != NULL && esc->type == cJSON_Array && esc->child != NULL,
                "No valid 'esc' escape windows.");

    load_windows(esc->child, hla_selection[h].esc, num_sites);

    printf("rev:\n");
    rev = cJSON_GetObjectItem(hla,"rev");
    json_assert(rev != NULL && rev->type == cJSON_Array && rev->child != NULL,
                "No valid 'esc' escape windows.");

    load_windows(rev->child, hla_selection[h].rev, num_sites);
  }

  json_assert(num_hla_types == h, "Mismatch in HLAs: [%i vs %i].", num_hla_types, h);

  // Now pull out the sequences chosen for each query sequence.
  cJSON *query, *closest_seqs, *seq;
  query = cJSON_GetObjectItem(root,"query");
  json_assert(query != NULL && query->type == cJSON_Array, "No 'query' object.");
  json_assert(query->child != NULL, "No states in chain.");
  query = query->child;
  int q, i;

  for(q = 0; query != NULL; query = query->next, q++)
  {
    printf("loading references closest to query %i.\n", q+1);

    closest_seqs = cJSON_GetObjectItem(query, "closest_n");
    json_assert(closest_seqs != NULL && closest_seqs->type == cJSON_Array,
                "No 'closest_seqs' object");
    json_assert(closest_seqs->child != NULL, "No entries in closest seq.");
    seq = closest_seqs->child;

    for(i = 0; seq != NULL && i < closest_n; seq = seq->next, i++)
    {
      json_assert(seq->type == cJSON_Number, "Closest seq is not number.");
      printf("%i ", seq->valueint);
      closest_refs[q][i] = seq->valueint;
    }

    json_assert(i == closest_n && seq == NULL, "Mismatch in closest N.");

    printf("\n");
  }

  json_assert(q == num_query, "Mismatch in number of queries: [%i vs %i].",
              q, num_query);

  cJSON_Delete(root);
  strbuf_dealloc(&sbuf);
}

void load_closest_n_from_json(const char *path, int num_sims, int num_sites,
                              Mosaic *true_copy)
{
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open file %s.", path);

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 8192); // 2^13 bytes = 8KB

  // Read in the whole file.
  while(strbuf_gzreadline(&sbuf, gz)) {}
  gzclose(gz);

  cJSON *root, *true_sequences, *sim, *seq_num;
  int n_sims, num_closest;

  json_assert(sbuf.end > 0, "Empty json file.");

  root = cJSON_Parse(sbuf.b);
  sim = cJSON_GetObjectItem(root,"true_sequences");
  json_assert(sim != NULL, "Shit is NULL son.");
  json_assert(sim != NULL && sim->type == cJSON_Array, "No 'true_sequences' object.");
  json_assert(sim->child != NULL, "No sequences in file.");
  sim = sim->child;
  
  printf("num sims: %i\n",num_sims);
  int n_seqs;

  // Get all the mosaics
  for(n_sims = 0; sim->next != NULL; sim = sim->next, n_sims++)
  {
    seq_num = cJSON_GetObjectItem(sim,"number");
    json_assert(seq_num != NULL && seq_num->type == cJSON_Number, "No 'number' object.");

    true_sequences = cJSON_GetObjectItem(sim,"closest_n");
    json_assert(true_sequences != NULL && true_sequences->type == cJSON_Array, "No 'closest_n' object.");

    num_closest = cJSON_GetArraySize(true_sequences);
    json_assert(num_closest <= num_sites, "num_seqs too big.");
    true_copy[n_sims].num_seqs = num_closest;
    
    for(n_seqs = 0 ; n_seqs < num_closest ; n_seqs++)
    {
      cJSON *subitem = cJSON_GetArrayItem(true_sequences, n_seqs);
      json_assert(subitem != NULL && subitem->type == cJSON_Number, "No 'num_seqs' object.");
      true_copy[n_sims].which_seqs[n_seqs] = subitem->valueint;
    }

  }

  seq_num = cJSON_GetObjectItem(sim,"number");
  json_assert(seq_num != NULL && seq_num->type == cJSON_Number, "No 'number' object.");
  
  true_sequences = cJSON_GetObjectItem(sim,"closest_n");
  json_assert(true_sequences != NULL && true_sequences->type == cJSON_Array, "No 'closest_n' object.");

  num_closest = cJSON_GetArraySize(true_sequences);
  true_copy[n_sims].num_seqs = num_closest;
    
  for(n_seqs = 0 ; n_seqs < num_closest ; n_seqs++)
  {
    cJSON * subitem = cJSON_GetArrayItem(true_sequences, n_seqs);
    true_copy[n_sims].which_seqs[n_seqs] = subitem->valueint;
  }
  printf("We got %i simulated sequence mosaics\n", n_sims);

  cJSON_Delete(root);
  strbuf_dealloc(&sbuf);
}

void load_fixed_parameters_from_json(const char *path, Decimal *kappa, Decimal *phi,
                                     Decimal *codon_prior)
{
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open file %s.", path);

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 8192); // 2^13 bytes = 8KB

  // Read in the whole file.
  while(strbuf_gzreadline(&sbuf, gz)) {}
  gzclose(gz);

  cJSON *parameters, *kappa_obj, *phi_obj, *codons, *codon, *obj;
  int c;

  json_assert(sbuf.end > 0, "Empty json file.");

  parameters = cJSON_Parse(sbuf.b);
  json_assert(parameters != NULL, "No 'parameters' object.");

  kappa_obj = cJSON_GetObjectItem(parameters,"kappa");
  json_assert(kappa_obj != NULL && kappa_obj->type == cJSON_Number, "No valid 'kappa' value.");
  *kappa = kappa_obj->valuedouble;

  phi_obj = cJSON_GetObjectItem(parameters,"phi");
  json_assert(phi_obj != NULL && phi_obj->type == cJSON_Number, "No valid 'phi' value.");
  *phi = phi_obj->valuedouble;
  
  codons = cJSON_GetObjectItem(parameters,"codons");
  json_assert(codons != NULL && codons->type == cJSON_Array && codons->child != NULL,
              "No 'codons' object.");
  codon = codons->child;

  // Loop over codons.
  for(c = 0; codon != NULL; codon = codon->next, c++)
  {
    obj = cJSON_GetObjectItem(codon,"codon_prior");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'codon_prior'.");
    codon_prior[c] = obj->valuedouble;
    printf("codon_prior[%i]:"DECPRINT"\n", c, codon_prior[c]);
  }

  json_assert(NUM_CODONS == c, "Mismatch in number of codon types: [%i vs %i].", NUM_CODONS, c);

  cJSON_Delete(parameters);
  strbuf_dealloc(&sbuf);
}

void load_lengths_for_simulation_from_json(const char *path, Decimal *kappa, Decimal *mu,
                                           int *codon_sequence_length, int *total_n_HLA,
                                           int *ploidy, int *n_genes)
{
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open file %s.", path);

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 8192); // 2^13 bytes = 8KB

  // Read in the whole file.
  while(strbuf_gzreadline(&sbuf, gz)) {}
  gzclose(gz);

  cJSON *parameters, *kappa_obj, *mu_obj, *seq_length_obj, *total_n_HLA_obj;
  cJSON *ploidy_obj, *n_genes_obj;

  json_assert(sbuf.end > 0, "Empty json file.");

  parameters = cJSON_Parse(sbuf.b);
  json_assert(parameters != NULL, "No 'parameters' object.");

  kappa_obj = cJSON_GetObjectItem(parameters, "kappa");
  json_assert(kappa_obj != NULL && kappa_obj->type == cJSON_Number, "No valid 'kappa' value.");
  *kappa = kappa_obj->valuedouble;

  mu_obj = cJSON_GetObjectItem(parameters, "mu");
  json_assert(mu_obj != NULL && mu_obj->type == cJSON_Number, "No valid 'mu' value.");
  *mu = mu_obj->valuedouble;

  seq_length_obj = cJSON_GetObjectItem(parameters, "codon_sequence_length");
  json_assert(seq_length_obj != NULL && seq_length_obj->type == cJSON_Number,
              "No valid 'codon_sequence_length' value.");
  *codon_sequence_length = seq_length_obj->valueint;

  total_n_HLA_obj = cJSON_GetObjectItem(parameters, "total_n_HLA");
  json_assert(total_n_HLA_obj != NULL && total_n_HLA_obj->type == cJSON_Number,
              "No valid 'total_n_HLA' value.");
  *total_n_HLA = total_n_HLA_obj->valueint;

  ploidy_obj = cJSON_GetObjectItem(parameters, "ploidy");
  json_assert(ploidy_obj != NULL && ploidy_obj->type == cJSON_Number, "No valid 'ploidy' value.");
  *ploidy = ploidy_obj->valueint;

  n_genes_obj = cJSON_GetObjectItem(parameters, "n_genes");
  json_assert(n_genes_obj != NULL && n_genes_obj->type == cJSON_Number, "No valid 'n_genes' value.");
  *n_genes = n_genes_obj->valueint;


  printf("kappa: "DECPRINT".\n", *kappa);
  printf("mu: "DECPRINT".\n", *mu);
  printf("ploidy: %i.\n", *ploidy);
  printf("n_genes: %i.\n", *n_genes);
  printf("total_n_HLA: %i.\n", *total_n_HLA);
  printf("codon_sequence_length: %i.\n", *codon_sequence_length);

  cJSON_Delete(parameters);
  strbuf_dealloc(&sbuf);
}

void load_parameters_for_simulation_from_json(const char *path, int codon_sequence_length,
                                              Decimal omega[codon_sequence_length], Decimal R[codon_sequence_length],
                                              Decimal reversion_selection[codon_sequence_length],
                                              int total_n_HLA, int n_genes, int n_HLA[n_genes],
                                              Decimal HLA_prevalences[total_n_HLA],
                                              Decimal HLA_selection_profiles[total_n_HLA][codon_sequence_length])
{
  gzFile gz = gzopen(path, "r");
  if(gz == NULL) die("Cannot open file %s.", path);

  StrBuf sbuf;
  strbuf_alloc(&sbuf, 8192); // 2^13 bytes = 8KB

  // Read in the whole file.
  while(strbuf_gzreadline(&sbuf, gz)) {}
  gzclose(gz);

  cJSON *parameters, *obj, *sites, *site;
  cJSON *HLA_prevalences_obj, *HLA_prevalence_obj, *HLA, *n_HLA_obj;
  int s, h;

  json_assert(sbuf.end > 0, "Empty json file.");

  parameters = cJSON_Parse(sbuf.b);
  json_assert(parameters != NULL, "No 'parameters' object.");

  sites = cJSON_GetObjectItem(parameters, "sites");
  json_assert(sites != NULL && sites->type == cJSON_Array && sites->child != NULL,
              "No 'sites' object");
  site = sites->child;
  
  HLA_prevalences_obj = cJSON_GetObjectItem(parameters, "HLA_prevalences");
  json_assert(HLA_prevalences_obj != NULL && HLA_prevalences_obj->type == cJSON_Array && HLA_prevalences_obj->child != NULL,
              "No 'HLA_prevalences' object.");
  HLA_prevalence_obj = HLA_prevalences_obj->child;

  for(h = 0; HLA_prevalence_obj != NULL; HLA_prevalence_obj = HLA_prevalence_obj->next, h++)
  {
    json_assert(HLA_prevalence_obj != NULL && HLA_prevalence_obj->type == cJSON_Number,
                "No valid 'HLA_prevalences'.");
    HLA_prevalences[h] = HLA_prevalence_obj->valuedouble;
    printf("HLA_prevalences[%i]: "DECPRINT".\n", h, HLA_prevalences[h]);
  }
  
  n_HLA_obj = cJSON_GetObjectItem(parameters, "n_HLA");
  json_assert(n_HLA_obj != NULL && n_HLA_obj->type == cJSON_Array && n_HLA_obj->child != NULL,
              "No 'n_HLA' object.");
  n_HLA_obj = n_HLA_obj->child;
  int n_HLA_check = 0;
  for(h = 0; n_HLA_obj != NULL; n_HLA_obj = n_HLA_obj->next, h++)
  {
    json_assert(n_HLA_obj != NULL && n_HLA_obj->type == cJSON_Number, "No valid 'n_HLA'.");
    n_HLA[h] = n_HLA_obj->valueint;
    n_HLA_check += n_HLA[h];
    printf("n_HLA[%i]: %i.\n", h, n_HLA[h]); 
  }

  json_assert(n_HLA_check == total_n_HLA, "Mismatch in the number of HLA types: [%i vs %i].", n_HLA_check, total_n_HLA);
  
  // Loop over sites.
  for(s = 0; site != NULL; site = site->next, s++)
  {
    obj = cJSON_GetObjectItem(site,"omega");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'omega'.");
    omega[s] = obj->valuedouble;

    obj = cJSON_GetObjectItem(site,"R");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'R'.");
    R[s] = obj->valuedouble;

    obj = cJSON_GetObjectItem(site,"reversion");
    json_assert(obj != NULL && obj->type == cJSON_Number, "No valid 'reversion_selection'.");
    reversion_selection[s] = obj->valuedouble;
    
    obj = cJSON_GetObjectItem(site,"HLA");
    HLA = obj->child;
    json_assert(obj != NULL && obj->type == cJSON_Array, "No valid 'HLA'.");
    for(h = 0; HLA != NULL; HLA = HLA->next, h++) {
      HLA_selection_profiles[h][s] = HLA->valuedouble;
    }
  }

  json_assert(s == codon_sequence_length, "Mismatch in sequence length: [%i vs %i].", s, codon_sequence_length);

  cJSON_Delete(parameters);
  strbuf_dealloc(&sbuf);
}
