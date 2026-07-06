#include "../Settings.h"
#include "../../Helpers_IC.h"

#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.cpp"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/utils.cpp"
#include "../../include/utils.h"
#include "../../include/utils-visual.cpp"
#include "../../include/utils-visual.h"

void macro_determine_binning()
{
        double binning[10 + 1];
        double binning_log[8+1];
        double binning_corr_log[nbin_rl_altlogbin + 1];
        double binning_weight[nbin_weight + 1];
        double binning_ptprod[nbin_ptprod + 1];
        double binning_h_pt[nbin_h_pt + 1];

        double binning_tau_log[nbin_tau_logbin + 1];

        determine_eqsizebinning(10, rl_logmin, rl_logmax, binning);
        determine_log10binning(nbin_h_pt, track_pt_min, 50., binning_h_pt);
        determine_log10binning(nbin_weight, weight_absmin, weight_absmax, binning_weight);
        determine_log10binning(nbin_ptprod, ptprod_absmin, ptprod_absmax, binning_ptprod);
        determine_log10binning(8, rl_logmin, rl_logmax, binning_log);
        determine_log10binning(nbin_tau_logbin, tau_min, tau_max, binning_tau_log);
        determine_log10binning(nbin_rl_altlogbin, rl_logmin, rl_logmax, binning_corr_log);
        
        // Pt prod
        std::cout<<"const double h_pt_binning[] = {track_pt_min";
        
        for (int i = 1 ; i < nbin_h_pt ; i++)
                std::cout<<", "<<binning_h_pt[i];

        std::cout<<", 50};"<<std::endl;

        // Pt prod
        std::cout<<"const double ptprod_binning[] = {ptprod_absmin";
        
        for (int i = 1 ; i < nbin_ptprod ; i++)
                std::cout<<", "<<binning_ptprod[i];

        std::cout<<", ptprod_absmax};"<<std::endl;

        // Weight
        std::cout<<"const double weight_binning[] = {weight_absmin";
        
        for (int i = 1 ; i < nbin_weight ; i++)
                std::cout<<", "<<binning_weight[i];

        std::cout<<", weight_absmax};"<<std::endl;

        // rl binning
        std::cout<<"const double linear_rl_binning[]              = {rl_logmin";
        
        for (int i = 1 ; i < 10 ; i++)
                std::cout<<", "<<binning[i];

        std::cout<<", rl_logmax};"<<std::endl;

        // rl log bin
        std::cout<<"const double rl_nominal_binning[]           = {rl_logmin";
        
        for (int i = 1 ; i < 8 ; i++)
                std::cout<<", "<<binning_log[i];

        std::cout<<", rl_logmax};"<<std::endl;

        // rl alternative logbinning
        std::cout<<"const double rl_altlogbinning[]           = {rl_logmin";
        
        for (int i = 1 ; i < nbin_rl_altlogbin ; i++)
                std::cout<<", "<<binning_corr_log[i];

        std::cout<<", rl_logmax};"<<std::endl;

        // tau logbinning
        std::cout<<"const double tau_nominal_binning[]          = {tau_min";
        
        for (int i = 1 ; i < nbin_tau_logbin ; i++)
                std::cout<<", "<<binning_tau_log[i];
        
        std::cout<<", tau_max};"<<std::endl;

        // TFile* f = new TFile((output_folder + namef_ntuple_mc_eec).c_str());
        // TNtuple* ntuple = (TNtuple*) f->Get(name_ntuple_mc.c_str());

        // TH1F* h = new TH1F("h","", nbin_ptprod, ptprod_binning);
        // // TH1F* h = new TH1F("h","", nbin_weight, weight_binning);

        // ntuple->Project("h","h1_pt*h2_pt");

        // h->Draw();

        // TH1F* h_cumulative = (TH1F*) h->GetCumulative();

        // h_cumulative->Draw();
}