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

void macro_print_jet_kinematics()
{
        TFile *f_mc     = new TFile((output_folder + "ntuple_bjets_mc.root").c_str());
        TFile *f_mcreco = new TFile((output_folder + "ntuple_bjets_mcreco.root").c_str());
        TFile *f_data   = new TFile((output_folder + "ntuple_bjets_data.root").c_str());

        TTree* t_mc     = (TTree*) f_mc->Get("BTree");
        TTree* t_mcreco = (TTree*) f_mcreco->Get("BTree");
        TTree* t_data   = (TTree*) f_data->Get("BTree");

        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        TH1F* h_jetpt_mc     = new TH1F("h_jetpt_mc"    , "", 50,12.5, 100);
        TH1F* h_jetpt_mcreco = new TH1F("h_jetpt_mcreco", "", 50,12.5, 100);
        TH1F* h_jetpt_data   = new TH1F("h_jetpt_data"  , "", 50,12.5, 100);

        TH1F* h_y_mc     = new TH1F("h_y_mc"    , "", 50,2.5, 4);
        TH1F* h_y_mcreco = new TH1F("h_y_mcreco", "", 50,2.5, 4);
        TH1F* h_y_data   = new TH1F("h_y_data"  , "", 50,2.5, 4);

        set_histogram_style(h_jetpt_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_jetpt_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_jetpt_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_y_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_y_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_y_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        t_mc->Project("h_jetpt_mc","jet_pt",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_jetpt_mcreco","jet_pt",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_jetpt_data","jet_pt",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        t_mc->Project("h_y_mc","jet_rap",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_y_mcreco","jet_rap",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_y_data","jet_rap",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        h_jetpt_mc->Scale(1./h_jetpt_mc->Integral());
        h_jetpt_mcreco->Scale(1./h_jetpt_mcreco->Integral());
        h_jetpt_data->Scale(1./h_jetpt_data->Integral());

        h_y_mc->Scale(1./h_y_mc->Integral());
        h_y_mcreco->Scale(1./h_y_mcreco->Integral());
        h_y_data->Scale(1./h_y_data->Integral());

        TCanvas* c = new TCanvas("c", "", 800, 600);
        c->Draw();

        // Plotting jet pt
        hs->Add(h_jetpt_mcreco, "P E1");
        hs->Add(h_jetpt_data, "P E1");

        // l->AddEntry(h_jetpt_mc,"h_jetpt_mc","p");
        l->AddEntry(h_jetpt_mcreco,"mc(reco)","p");
        l->AddEntry(h_jetpt_data,"data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{p}_{T, jet}(GeV);Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(1);
        
        c->Print("./plots/kinematics_jet_pt.pdf");

        // Plotting jet rapidity
        hs = new THStack();
        hs->Add(h_y_mcreco, "P E1");
        hs->Add(h_y_data, "P E1");

        // l->AddEntry(h_jetpt_mc,"h_jetpt_mc","p");
        l = new TLegend();
        l->AddEntry(h_y_mcreco, "mc(reco)","p");
        l->AddEntry(h_y_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{y}_{jet};Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(0);

        c->Print("./plots/kinematics_jet_rapidity.pdf");
}
