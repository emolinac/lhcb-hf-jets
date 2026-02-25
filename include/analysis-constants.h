#ifndef ANALYSIS_CONSTANTS_H
#define ANALYSIS_CONSTANTS_H

// Masses (GeV)
const double rho_mass      = 0.77526;  // PDG 2023
const double omega_mass    = 0.78266; // PDG 2023
const double eta_mass      = 0.547862; // PDG 2023
const double etaprime_mass = 0.95778;  // PDG 2023
const double kaonp_mass    = 0.493677; // PDG 2023
const double kaonm_mass    = 0.493677; // PDG 2023
const double kaon_mass     = 0.497611; // PDG 2023
const double pi_mass       = 0.134977;
const double phi_mass      = 1.019455;
const double mass_res      = 0.008; // Mass resolution parameter (see src-resolution)

// Analysis
const int nominal_niter = 4;
const int niter_ct = 1; // if we are separating by polarities then 1 is correct!
const int reg_par_window = 3;

// Visual
const double std_marker_size  = 1.0;
const int    std_marker_style = 8;
const int    std_line_width   = 3;

const int std_marker_style_jet_pt[]  = {24,25,27,42,46,40,45};
const int corr_marker_style_jet_pt[] = {20,21,33,43,47,41,45};

const int std_marker_color_jet_pt[]  = {868,797,618,41,46,38,418,861,617,1};
const int corr_marker_color_jet_pt[] = {868,797,618,41,46,38,418,861,617,1};

const double text_size_correction_plots = 0.015;

#endif