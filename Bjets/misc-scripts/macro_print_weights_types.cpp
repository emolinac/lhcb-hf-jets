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

void macro_print_weights_types()
{
        TFile *f_mc     = new TFile((output_folder + "ntuple_bjets_mc.root").c_str());
        TFile *f_mcreco = new TFile((output_folder + "ntuple_bjets_mcreco.root").c_str());
        TFile *f_data   = new TFile((output_folder + "ntuple_bjets_data.root").c_str());

        TTree* t_mc     = (TTree*) f_mc->Get("BTree");
        TTree* t_mcreco = (TTree*) f_mcreco->Get("BTree");
        TTree* t_data   = (TTree*) f_data->Get("BTree");

        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        TH1F* h_mc = new TH1F("h_mc", "", nbin_weight, weight_binning);
        TH1F* h_mcreco = new TH1F("h_mcreco", "", nbin_weight, weight_binning);
        TH1F* h_data = new TH1F("h_data", "", nbin_weight, weight_binning);

        set_histogram_style(h_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size);
        set_histogram_style(h_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size);
        set_histogram_style(h_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[2], std_marker_size);

        t_mc->Project("h_mc","pair_weight",Form("jet_pt>%f&&jet_pt<%f", 20., 100.));
        t_mcreco->Project("h_mcreco","pair_weight",Form("jet_pt>%f&&jet_pt<%f", 20., 100.));
        t_data->Project("h_data","pair_weight",Form("jet_pt>%f&&jet_pt<%f", 20., 100.));

        h_mc->Scale(1./h_mc->Integral());
        h_mcreco->Scale(1./h_mcreco->Integral());
        h_data->Scale(1./h_data->Integral());

        hs->Add(h_mc, "P E1");
        hs->Add(h_mcreco, "P E1");
        hs->Add(h_data, "P E1");

        l->AddEntry(h_mc,"h_mc","p");
        l->AddEntry(h_mcreco,"h_mcreco","p");
        l->AddEntry(h_data,"h_data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";weights;norm.");

        // hs->SetMaximum(1.4);
        // hs->SetMinimum(0.6);

        l->Draw("same");

        gPad->SetLogx(1);
}
