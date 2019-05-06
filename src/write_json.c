#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include "constants.h"
#include "selection_window.h"
#include "write_json.h"

void write_json(FILE *final_state_json,
	              int chain_i, Decimal curr_llk_across_queries, Decimal mu, int num_sites,
	              Decimal omega[num_sites], Decimal R[num_sites], HLASelection *hla_windows,
	              int num_hla_types,
                int num_query, int closest_n, int closest_refs[num_query][closest_n])
{
  int i, h, w, q, n;
  fprintf(final_state_json, "{\n \"chain\":[{\n    \"state\":%i,\n"
                            "    \"llkhood\": "DECPRINTJSON",\n"
                            "    \"mu\": "DECPRINTJSON",\n"
                            "    \"sites\": [",
                            chain_i, curr_llk_across_queries, mu);
  for(i = 0; i < num_sites; i++) {
  	if (i != 0) {
  		fprintf(final_state_json, "              ");
  	}
    fprintf(final_state_json, "{\"omega\":"DECPRINTJSON", \"R\":"DECPRINTJSON"}", omega[i], R[i]);
  	fprintf(final_state_json, "%s\n", i != (num_sites-1) ? "," : "],");
  }
  fprintf(final_state_json, "    \"hlas\": [{\n");
  for(h = 0; h < num_hla_types; h++) {
  	if(h != 0) {
  		fprintf(final_state_json, "    {\n");
  	}
  	fprintf(final_state_json, "      \"id\": %i,\n"
  		                      "      \"esc\": [", h+1);
  	for(w = 0; w < hla_windows[h].esc->num_windows; w++) {
  	  if(w != 0) {
  	  	fprintf(final_state_json, "                      ");
  	  }
  	  fprintf(final_state_json,"{\"start\":%i, \"len\":%i, \"coeff\":"DECPRINTJSON"}",
  			  hla_windows[h].esc->windows[w].start, hla_windows[h].esc->windows[w].length,
  		      hla_windows[h].esc->windows[w].coeff);
  	  fprintf(final_state_json, "%s\n", w != (hla_windows[h].esc->num_windows - 1) ? "," : "],");
  	}
    fprintf(final_state_json, "      \"rev\": [");
  	for(w = 0; w < hla_windows[h].rev->num_windows; w++) {
  	  if(w != 0) {
  	  	fprintf(final_state_json, "                      ");
  	  }
  	  fprintf(final_state_json,"{\"start\":%i, \"len\":%i, \"coeff\":"DECPRINTJSON"}",
  			  hla_windows[h].rev->windows[w].start, hla_windows[h].rev->windows[w].length,
  		      hla_windows[h].rev->windows[w].coeff);
  	  fprintf(final_state_json, "%s\n", w != (hla_windows[h].rev->num_windows - 1) ? "," : "]");
  	}
  	fprintf(final_state_json, "%s\n", h != (num_hla_types-1) ? "    }," : "    }]");
  }
  fprintf(final_state_json, "  }],\n \"query\":[{\n");
  for(q = 0; q < num_query; q++) {
    fprintf(final_state_json, "    \"q\":%i,\n"
                              "    \"closest_n\": [", q);
    for(n = 0; n < closest_n; n++) {
      fprintf(final_state_json, "%i", closest_refs[q][n]);
      fprintf(final_state_json, "%s", n != (closest_n-1) ? "," : "]");
    }
    fprintf(final_state_json, "%s\n", q != (num_query-1) ? "\n  },\n    {" : "\n  }]");
  }
  fprintf(final_state_json, "}");
}

void write_summary_json(FILE *summary_sim_json,
                        Decimal mu, int num_sites, int ploidy,
                        int n_genes, int num_HLA[n_genes],
                        int total_n_HLA,
                        Decimal HLA_prevalences[total_n_HLA],
                        Decimal omega[num_sites], Decimal R[num_sites],
                        Decimal reversion_selection[num_sites],
                        Decimal HLA_selection_profiles[total_n_HLA][num_sites])
{
  int i, h, s;
  fprintf(summary_sim_json, "{\n    \"kappa\": "DECPRINTJSON",\n"
                               "    \"HLA_prevalences\": [",
                          kappa);
  for (h=0; h<total_n_HLA-1; h++){
    fprintf(summary_sim_json, DECPRINTJSON",", HLA_prevalences[h]);
  }
  fprintf(summary_sim_json, DECPRINTJSON"],\n"
                            "    \"n_HLA\": [", HLA_prevalences[h]);
  for (i=0; i<n_genes-1; i++) fprintf(summary_sim_json, "%i,", num_HLA[i]);
  fprintf(summary_sim_json, "%i],\n"
                            "    \"mu\": "DECPRINTJSON",\n"
                            "    \"codon_sequence_length\": %i,\n"
                            "    \"total_n_HLA\": %i,\n"
                            "    \"ploidy\": %i,\n"
                            "    \"n_genes\": %i,\n"
                            "    \"sites\": [{\"omega\": "DECPRINTJSON", \"R\": "DECPRINTJSON", \"reversion\": "DECPRINTJSON", \"HLA\":[",
                          num_HLA[i], mu, num_sites, total_n_HLA, ploidy, n_genes, omega[0], R[0], reversion_selection[0]);
  for (h=0; h<total_n_HLA-1;h++) fprintf(summary_sim_json, DECPRINTJSON",", HLA_selection_profiles[h][0]);
  fprintf(summary_sim_json, DECPRINTJSON"]}", HLA_selection_profiles[h][0]);

  for (s=1; s<num_sites; s++){
    fprintf(summary_sim_json, ",\n              {\"omega\": "DECPRINTJSON", \"R\": "DECPRINTJSON", \"reversion\": "DECPRINTJSON", \"HLA\":[",
                              omega[s], R[s], reversion_selection[s]);
    for (h=0; h<total_n_HLA-1; h++) fprintf(summary_sim_json, DECPRINTJSON",", HLA_selection_profiles[h][s]);
    fprintf(summary_sim_json, DECPRINTJSON"]}", HLA_selection_profiles[h][s]);
  }
  fprintf(summary_sim_json, "]\n}");
}

void write_closest_n_json(FILE *closest_n_json, Mosaic *true_copy, int num_sims)
{
  int sims, i;
  fprintf(closest_n_json, "{\n  \"true_sequences\":[");

  for (sims = 0; sims < num_sims; sims++)
  {
    fprintf(closest_n_json, "{\n    \"number\": %i,\n", sims+1);
    fprintf(closest_n_json, "    \"closest_n\": [");
    
    for (i = 0; i < true_copy[sims].num_seqs; i++){
      fprintf(closest_n_json, "%i",true_copy[sims].which_seqs[i]);
      fprintf(closest_n_json, "%s", i != (true_copy[sims].num_seqs-1) ? "," : "]\n");
    }

    fprintf(closest_n_json, "%s", sims != (num_sims-1) ? "    },\n    " : "  }]\n}\n");
  }
}
