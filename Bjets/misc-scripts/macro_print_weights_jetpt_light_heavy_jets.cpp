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
#include "../../include/TBJetsMC.h"
#include "../../include/TBJetsMC.C"
#include "../../include/utils.cpp"
#include "../../include/utils.h"
#include "../../include/utils-visual.cpp"
#include "../../include/utils-visual.h"

using namespace std;

TCut wta_cut = "WTA_reco_dist<0.005";

void macro_print_weights_jetpt_light_heavy_jets(double jet_pt_min = 20, double jet_pt_max = 100)
{
        TFile *f_light = new TFile((output_folder + "ntuple_zjets.root").c_str());
        TFile *f_heavy = new TFile((output_folder + "ntuple_bjets_data.root").c_str());

        TTree* t_light = (TNtuple*) f_light->Get("ntuple_data_pair");
        TTree* t_heavy = (TTree*) f_heavy->Get("BTree");

        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        TH1F* h_light = new TH1F("h_light", "", nbin_weight, weight_binning);
        TH1F* h_heavy = new TH1F("h_heavy", "", nbin_weight, weight_binning);

        set_histogram_style(h_light, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size);
        set_histogram_style(h_heavy, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[2], std_marker_size);

        t_light->Project("h_light","weight_pt", Form("jet_pt>%f&&jet_pt<%f", jet_pt_min, jet_pt_max));
        t_heavy->Project("h_heavy","pair_weight", wta_cut + Form("jet_pt>%f&&jet_pt<%f", jet_pt_min, jet_pt_max));

        h_light->Scale(1./h_light->Integral());
        h_heavy->Scale(1./h_heavy->Integral());

        hs->Add(h_light, "P E1");
        hs->Add(h_heavy, "P E1");

        l->SetHeader(Form("%.0f< #it{p}_{T, jet}<%.0f GeV", jet_pt_min, jet_pt_max));
        l->AddEntry(h_light, "Z-tagged jet","p");
        l->AddEntry(h_heavy, "B-tagged jet","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";weights;norm. data");

        l->Draw("same");

        gPad->SetLogx(1);
}
