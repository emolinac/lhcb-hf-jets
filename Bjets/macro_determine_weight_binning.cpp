#include "Settings.h"
#include "../Helpers_IC.h"

#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

void macro_determine_weight_binning()
{
        const int nbin_weight_test = 30;

        double binning_weight[nbin_weight_test + 1];
        
        determine_log10binning(nbin_weight_test, weight_absmin, weight_absmax, binning_weight);
        
        // Weight
        std::cout<<"const double weight_binning[] = {weight_absmin";
        
        for (int i = 1 ; i < nbin_weight_test ; i++)
                std::cout<<", "<<binning_weight[i];

        std::cout<<", weight_absmax};"<<std::endl;

        TFile* f = new TFile((output_folder + "ntuple_bjets_mcreco.root").c_str());
        TTree* ntuple = (TTree*) f->Get("BTree");

        TH1F* h = new TH1F("h","", nbin_weight_test, binning_weight);

        ntuple->Project("h","pair_weight");

        h->Draw();

        // TH1F* h_cumulative = (TH1F*) h->GetCumulative();

        // h_cumulative->Draw();
}