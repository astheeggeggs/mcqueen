#ifndef SELECTION_WINDOW_H_
#define SELECTION_WINDOW_H_

#include <stdio.h>
#include "constants.h"

typedef struct
{
  int start, length;
  Decimal coeff;
} SelectionWindow;

typedef struct
{
  SelectionWindow *windows;
  int num_windows;
} HLAWindowSet;

typedef struct
{
  HLAWindowSet *esc, *rev;
} HLASelection;

void setup_gsl();
void clearup_gsl();

void selec_window_hla_alloc(int num_sites, int num_hla_types,
                            HLASelection *HLAs);

void selec_window_hla_dealloc(HLASelection *HLAs);

void selec_windows_set_lengths(HLAWindowSet *hla_set, int num_sites);
void selec_windows_sort(HLAWindowSet *hla_set);

// Returns index of first (left-hand) window.
int selec_window_split_rnd(HLAWindowSet *hla_set, int num_sites,
                           Decimal *old_coeff);

void selec_window_split_undo(HLAWindowSet *hla_set, int windex, Decimal old_coeff);

// Returns index of new bigger window.
int selec_window_merge_rnd(HLAWindowSet *hla_set, int *old_start_site,
                           Decimal *old_coeff0, Decimal *old_coeff1);

void selec_window_merge_undo(HLAWindowSet *hla_set, int windex, int old_start_site,
                             Decimal old_coeff0, Decimal old_coeff1);

// Returns index of left hand window
// or -1 if move failed.
int selec_window_grow_rnd(HLAWindowSet *hla_set, int *old_start_site);
void selec_window_grow_undo(HLAWindowSet *hla_set, int windex, int old_start_site);

void selec_window_grow_accept(HLAWindowSet *hla_set, int new_start_site, 
                              int old_start_site, int num_sites);

int selec_window_change_coeff_rnd(HLAWindowSet *hla_set, Decimal *old_coeff0);
void selec_window_change_coeff_undo(HLAWindowSet *hla_set,
                                    int windex, Decimal old_coeff);

void selec_windows_randomise(HLAWindowSet *hla_set, int num_hla_types, int num_windows,
                             Decimal log_mu_sel_prior, Decimal sigma_sel_prior);

void selec_windows_set(HLAWindowSet *hla_set, int num_sites, int num_windows, int start[num_windows], 
                       Decimal coeff[num_windows]);

void selec_window_print(HLAWindowSet *hla_set);

void selec_window_print_to_file(HLAWindowSet *hla_set, FILE *print_to_file);

void update_HLA_selection(const HLAWindowSet *winset, int num_sites, Decimal hla_select[num_sites]);

#endif /* SELECTION_WINDOW_H_ */
