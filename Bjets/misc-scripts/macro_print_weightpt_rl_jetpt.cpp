#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "../Settings.h"

#include "../../Helpers_IC.h"

#include "../../include/analysis-constants.h"
#include "../../include/analysis-binning.h"
#include "../../include/analysis-cuts.cpp"
#include "../../include/analysis-cuts.h"
#include "../../include/directories.h"
#include "../../include/names.h"
#include "../../include/TBJetsMC.h"
#include "../../include/TBJetsMC.C"
#include "../../include/utils.cpp"
#include "../../include/utils.h"
#include "../../include/utils-visual.cpp"
#include "../../include/utils-visual.h"

using namespace std;

TCut wta_cut = "WTA_reco_dist<0.005";

void macro_print_weightpt_rl_jetpt()
{
        TFile* fin = new TFile((output_folder + "ntuple_bjets_data.root").c_str());

        TNtuple* ntuple = (TNtuple*) fin->Get("BTree");

        TH2D* h[nbin_jet_pt];

        TCanvas* c = new TCanvas("c","",1500,500);
        c->Draw();
        c->Divide(3,1);
        
        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                h[bin] = new TH2D(Form("h%i",bin),"",nbin_rl_nominal, rl_nominal_binning,nbin_weight,weight_binning);
                h[bin]->SetTitle(Form("%.0f<p_{T,jet}<%.0f;R_{L};w",jet_pt_binning[bin],jet_pt_binning[bin + 1]));
                
                ntuple->Project(Form("h%i",bin),"pair_weight:pair_rl",wta_cut + Form("jet_pt>%f&&jet_pt<%f", jet_pt_binning[bin],jet_pt_binning[bin + 1]));

                c->cd(bin + 1);
                
                h[bin]->Draw("col");
                gPad->SetLogx(1);
                gPad->SetLogy(1);
                gPad->SetLogz(1);
                h[bin]->Smooth();
        }

        c->Print("../plots/weightpt_rl_jet_pt_distributions.pdf");
}