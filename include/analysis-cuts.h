#ifndef ANALYSIS_CUTS_H
#define ANALYSIS_CUTS_H

#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "analysis-constants.h"
#include "analysis-binning.h"

bool apply_jet_cuts(double jet_eta, double jet_pt);

bool apply_muon_cuts(double mu_pt);

bool apply_kaon_cuts(double kaon_pt);

bool apply_chargedparticle_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta);

bool apply_chargedparticle_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet);

bool apply_chargedparticle_momentum_cuts(double charge, double p, double pt, double eta);

bool apply_chargedparticle_momentum_cuts(double charge, double p, double pt, double eta, double deltaR_h_jet);

bool apply_particle_cuts(double p, double pt, double chi2ndf, double probnnghost, double eta);

bool apply_particle_cuts(double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet);

bool apply_particle_momentum_cuts(double p, double pt, double eta);

bool apply_particle_momentum_cuts(double p, double pt, double eta, double deltaR_h_jet);

#endif
