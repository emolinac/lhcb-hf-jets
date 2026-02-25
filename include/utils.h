#ifndef UTILS_ALGORITHMS_H
#define UTILS_ALGORITHMS_H

#include "TRandom3.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TRandom3.h"

double get_hwhm(TH1F* h);

double get_median_from_cumulative(TH1F* h);

void determine_log10binning(int Nbins, double x_i, double x_f, double* binning);

void determine_eqsizebinning(int Nbins, double x_i, double x_f, double* binning);

void apply_unfolded_weights(TH3F* h_rl_jetpt_weight, TH2F* h_rl_jetpt);

void apply_unfolded_weights(TH3D* h_rl_jetpt_weight, TH2D* h_rl_jetpt);

void apply_jet_weight_to_npairs(TH3D* h_npair, TH1F* h_purity_jet, TH1F* h_efficiency_jet);

void project_nominal_phase_space(TH2D* h_2d, TH1F* h_1d);

void substract_stat_error(TH1F* hnominal, TH1F* hsyst);

void get_tau_from_uoflow_eec(TH1F* h_eec, TH1F* h_tau);

void get_tau_binning_from_eec_binning(double* tau_binning, double* eec_binning, double average_pt2_jet);

void normalize_by_njets(TH1F* h, double h_njet_content, double h_njet_error);

void normalize_by_njets(TH1D* h, double h_njet_content, double h_njet_error);

void square_root_bins(TH1F* h);

void regularize_correction_factors(TH2F* h);

void regularize_correction_factors(TH2D* h);

void regularize_correction_factors(TH3F* h);

void regularize_correction_factors(TH3D* h);

void set_histo_with_systematics(TH1F* hrelerror, TH1F* hnominal, TH1F* hsystematic, int syst_index, bool print_table = false);

void set_histoa_errors_as_histob_content(TH1F* hnominal, TH1F* hnominalerror);

void set_histo_sqrt_content(TH1F* h);

void set_histo_null_errors(TH1F* h);

void set_shift_histo(TH2F* href, TH2F* hshift, TRandom3* rndm);

void set_shift_histo(TH2D* href, TH2D* hshift, TRandom3* rndm);

void smear_pseudodata(TH1D* hpseudodata, TH1D* hrealdata, TRandom3* rndm);

void smear_pseudodata(TH1F* hpseudodata, TH1F* hrealdata, TRandom3* rndm);

void smear_pseudodata(TH3D* hpseudodata, TH3D* hrealdata, TRandom3* rndm);

void smear_pseudodata(TH3F* hpseudodata, TH3F* hrealdata, TRandom3* rndm);

void smear_pseudodata(TH2D* hpseudodata, TH2D* hrealdata, TRandom3* rndm);

void smear_pseudodata(TH2F* hpseudodata, TH2F* hrealdata, TRandom3* rndm);

void set_data_ntuple_branches(TNtuple* ntuple, float* event_weight, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror);

void set_data_ntuple_branches(TNtuple* ntuple, float* event_weight, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror, float* eq_charge);

void set_unfolding_ntuple_branches(TNtuple* ntuple, float* R_L_reco, float* R_L_truth, float* jet_pt_reco, float* jet_pt_truth, float* weight_pt_reco, float* weight_pt_truth);

void set_unfolding_jet_ntuple_branches(TNtuple* ntuple, float* jet_pt_reco, float* jet_pt_truth);

void set_unity_content(TH1F* h);

double get_jes_jer_factor(const double jet_pt, TRandom3 *myRNG);

double weight(double h1_E, double h2_E, double jet_E);

void smooth_nominal_phase_space(TH1F* h_to_smooth, TH1F* h_nominal_phase_space);

#endif