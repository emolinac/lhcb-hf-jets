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

void macro_print_hf_kinematics()
{
        TFile *f_mc     = new TFile((output_folder + "ntuple_bjets_mc.root").c_str());
        TFile *f_mcreco = new TFile((output_folder + "ntuple_bjets_mcreco.root").c_str());
        TFile *f_data   = new TFile((output_folder + "ntuple_bjets_data.root").c_str());

        TTree* t_mc     = (TTree*) f_mc->Get("BTree");
        TTree* t_mcreco = (TTree*) f_mcreco->Get("BTree");
        TTree* t_data   = (TTree*) f_data->Get("BTree");

        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        TH1F* h_hfpt_mc     = new TH1F("h_hfpt_mc"    , "", 50,12.5, 100);
        TH1F* h_hfpt_mcreco = new TH1F("h_hfpt_mcreco", "", 50,12.5, 100);
        TH1F* h_hfpt_data   = new TH1F("h_hfpt_data"  , "", 50,12.5, 100);

        TH1F* h_y_mc     = new TH1F("h_y_mc"    , "", 50,2.5, 4);
        TH1F* h_y_mcreco = new TH1F("h_y_mcreco", "", 50,2.5, 4);
        TH1F* h_y_data   = new TH1F("h_y_data"  , "", 50,2.5, 4);

        TH1F* h_mass_mc     = new TH1F("h_mass_mc"    , "", 50, 5.15, 5.55);
        TH1F* h_mass_mcreco = new TH1F("h_mass_mcreco", "", 50, 5.15, 5.55);
        TH1F* h_mass_data   = new TH1F("h_mass_data"  , "", 50, 5.15, 5.55);

        TH1F* h_jpsi_mass_mc     = new TH1F("h_jpsi_mass_mc"    , "", 50, 2.99, 3.2);
        TH1F* h_jpsi_mass_mcreco = new TH1F("h_jpsi_mass_mcreco", "", 50, 2.99, 3.2);
        TH1F* h_jpsi_mass_data   = new TH1F("h_jpsi_mass_data"  , "", 50, 2.99, 3.2);

        TH1F* h_kaon_mass_mc     = new TH1F("h_kaon_mass_mc"    , "", 50, 0.47, 0.52);
        TH1F* h_kaon_mass_mcreco = new TH1F("h_kaon_mass_mcreco", "", 50, 0.47, 0.52);
        TH1F* h_kaon_mass_data   = new TH1F("h_kaon_mass_data"  , "", 50, 0.47, 0.52);

        set_histogram_style(h_hfpt_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_hfpt_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_hfpt_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_y_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_y_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_y_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_mass_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_mass_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_mass_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_jpsi_mass_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_jpsi_mass_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_jpsi_mass_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_kaon_mass_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_kaon_mass_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_kaon_mass_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        t_mc->Project("h_hfpt_mc","HF_pt",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_hfpt_mcreco","HF_pt",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_hfpt_data","HF_pt",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        TString hf_y = "0.5 * log((HF_e + HF_pz)/(HF_e - HF_pz))";
        t_mc->Project("h_y_mc",hf_y,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_y_mcreco",hf_y,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_y_data",hf_y,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        TString hf_mass = "TMath::Sqrt((HF_e*HF_e) - (HF_px*HF_px) - (HF_py*HF_py) - (HF_pz*HF_pz))";
        t_mc->Project("h_mass_mc",hf_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_mass_mcreco",hf_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_mass_data",hf_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        TString kaon_mass = "TMath::Sqrt((K_e*K_e) - (K_px*K_px) - (K_py*K_py) - (K_pz*K_pz))";
        // t_mc->Project("h_kaon_mass_mc",kaon_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_kaon_mass_mcreco",kaon_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_kaon_mass_data",kaon_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        TString jpsi_mass = "TMath::Sqrt((Jpsi_e*Jpsi_e) - (Jpsi_px*Jpsi_px) - (Jpsi_py*Jpsi_py) - (Jpsi_pz*Jpsi_pz))";
        // t_mc->Project("h_jpsi_mass_mc",jpsi_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_jpsi_mass_mcreco",jpsi_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_jpsi_mass_data",jpsi_mass,Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        h_hfpt_mc->Scale(1./h_hfpt_mc->Integral());
        h_hfpt_mcreco->Scale(1./h_hfpt_mcreco->Integral());
        h_hfpt_data->Scale(1./h_hfpt_data->Integral());

        h_y_mc->Scale(1./h_y_mc->Integral());
        h_y_mcreco->Scale(1./h_y_mcreco->Integral());
        h_y_data->Scale(1./h_y_data->Integral());

        h_mass_mc->Scale(1./h_mass_mc->Integral());
        h_mass_mcreco->Scale(1./h_mass_mcreco->Integral());
        h_mass_data->Scale(1./h_mass_data->Integral());

        h_jpsi_mass_mc->Scale(1./h_jpsi_mass_mc->Integral());
        h_jpsi_mass_mcreco->Scale(1./h_jpsi_mass_mcreco->Integral());
        h_jpsi_mass_data->Scale(1./h_jpsi_mass_data->Integral());

        h_kaon_mass_mc->Scale(1./h_kaon_mass_mc->Integral());
        h_kaon_mass_mcreco->Scale(1./h_kaon_mass_mcreco->Integral());
        h_kaon_mass_data->Scale(1./h_kaon_mass_data->Integral());

        TCanvas* c = new TCanvas("c", "", 800, 600);
        c->Draw();

        // Plotting jet pt
        hs->Add(h_hfpt_mcreco, "P E1");
        hs->Add(h_hfpt_data, "P E1");

        // l->AddEntry(h_hfpt_mc,"h_hfpt_mc","p");
        l->AddEntry(h_hfpt_mcreco,"mc(reco)","p");
        l->AddEntry(h_hfpt_data,"data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{p}_{T, HF}(GeV);Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(1);
        
        c->Print("./plots/kinematics_hf_pt.pdf");

        // Plotting jet rapidity
        hs = new THStack();
        hs->Add(h_y_mcreco, "P E1");
        hs->Add(h_y_data, "P E1");

        // l->AddEntry(h_hfpt_mc,"h_hfpt_mc","p");
        l = new TLegend();
        l->AddEntry(h_y_mcreco, "mc(reco)","p");
        l->AddEntry(h_y_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{y}_{HF};Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(0);

        c->Print("./plots/kinematics_hf_rapidity.pdf");

        // Plotting hf mass
        hs = new THStack();
        hs->Add(h_mass_mcreco, "P E1");
        hs->Add(h_mass_data, "P E1");

        // l->AddEntry(h_hfpt_mc,"h_hfpt_mc","p");
        l = new TLegend(0.7,0.7,0.9,0.9);
        l->AddEntry(h_mass_mcreco, "mc(reco)","p");
        l->AddEntry(h_mass_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{M}_{HF}(GeV);Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(0);

        c->Print("./plots/kinematics_hf_mass.pdf");

        // Plotting jpsi mass
        hs = new THStack();
        hs->Add(h_jpsi_mass_mcreco, "P E1");
        hs->Add(h_jpsi_mass_data, "P E1");

        // l->AddEntry(h_hfpt_mc,"h_hfpt_mc","p");
        l = new TLegend(0.7,0.7,0.9,0.9);
        l->AddEntry(h_jpsi_mass_mcreco, "mc(reco)","p");
        l->AddEntry(h_jpsi_mass_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{M}_{J/#{Psi}}(GeV);Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(0);

        c->Print("./plots/kinematics_jpsi_mass.pdf");

        // Plotting hf mass
        hs = new THStack();
        hs->Add(h_kaon_mass_mcreco, "P E1");
        hs->Add(h_kaon_mass_data, "P E1");

        // l->AddEntry(h_hfpt_mc,"h_hfpt_mc","p");
        l = new TLegend(0.7,0.7,0.9,0.9);
        l->AddEntry(h_kaon_mass_mcreco, "mc(reco)","p");
        l->AddEntry(h_kaon_mass_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{M}_{Kaon}(GeV);Normalized distributions");

        l->Draw("same");

        gPad->SetLogy(1);

        c->Print("./plots/kinematics_kaon_mass.pdf");
}
