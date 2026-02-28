#include "analysis-cuts.h"
#include "TCut.h"
#include "TMath.h"
#include "TString.h"
#include "analysis-constants.h"
#include "analysis-binning.h"

// Jet cuts
const double jet_eta_min = 2.5;
const double jet_eta_max = 4.0;
const double jet_pt_min  = 15;

const double jet_radius = 0.5;

// Z boson cuts
const double muon_pt_min  = 20.;
const double lhcb_eta_min = 2;
const double lhcb_eta_max = 4.5;
const double dimuon_mass_min = 60;
const double dimuon_mass_max = 120;
const double muon_trackprob_min = 0.001;

// Track cuts
const double track_chi2ndf_max     = 3;
const double track_p_min           = 4;
const double track_p_max           = 1000;
const double track_pt_min          = 0.25;
const double track_probnnghost_max = 0.5;

bool apply_jet_cuts(double jet_eta, double jet_pt)
{
        if (jet_eta < jet_eta_min || jet_eta > jet_eta_max) 
                return false;
        
        if (jet_pt < unfolding_jet_pt_binning[0])
                return false;

        return true;    
}

bool apply_muon_cuts(double deltaR_mu_jet, double mu_pt, double mu_eta)
{
        if (deltaR_mu_jet < jet_radius) 
                return false; 
        
        if (mu_pt < muon_pt_min)
                return false;

        if (mu_eta < lhcb_eta_min || mu_eta > lhcb_eta_max) 
                return false;

        return true;
}

bool apply_chargedparticle_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (chi2ndf > track_chi2ndf_max)
                return false;

        if (probnnghost > track_probnnghost_max)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        return true;
}

bool apply_chargedparticle_cuts(double charge, double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (chi2ndf > track_chi2ndf_max)
                return false;

        if (probnnghost > track_probnnghost_max)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        if (deltaR_h_jet > jet_radius)
                return false;

        return true;
}

bool apply_chargedparticle_momentum_cuts(double charge, double p, double pt, double eta)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        return true;
}

bool apply_chargedparticle_momentum_cuts(double charge, double p, double pt, double eta, double deltaR_h_jet)
{
        if (charge == 0)
                return false;

        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        if (deltaR_h_jet > jet_radius)
                return false;

        return true;
}

bool apply_particle_cuts(double p, double pt, double chi2ndf, double probnnghost, double eta)
{
        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (chi2ndf > track_chi2ndf_max)
                return false;

        if (probnnghost > track_probnnghost_max)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        return true;
}

bool apply_particle_cuts(double p, double pt, double chi2ndf, double probnnghost, double eta, double deltaR_h_jet)
{
        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (chi2ndf > track_chi2ndf_max)
                return false;

        if (probnnghost > track_probnnghost_max)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        if (deltaR_h_jet > jet_radius)
                return false;

        return true;
}

bool apply_particle_momentum_cuts(double p, double pt, double eta)
{
        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        return true;
}

bool apply_particle_momentum_cuts(double p, double pt, double eta, double deltaR_h_jet)
{
        if (p < track_p_min || p > track_p_max)
                return false;

        if (pt < track_pt_min)
                return false;

        if (eta < lhcb_eta_min || eta > lhcb_eta_max)
                return false;

        if (deltaR_h_jet > jet_radius)
                return false;

        return true;
}
