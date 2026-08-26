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

void macro_print_pairs_kinematics()
{
        TFile *f_mc     = new TFile((output_folder + "ntuple_bjets_mc.root").c_str());
        TFile *f_mcreco = new TFile((output_folder + "ntuple_bjets_mcreco.root").c_str());
        TFile *f_data   = new TFile((output_folder + "ntuple_bjets_data.root").c_str());

        TTree* t_mc     = (TTree*) f_mc->Get("BTree");
        TTree* t_mcreco = (TTree*) f_mcreco->Get("BTree");
        TTree* t_data   = (TTree*) f_data->Get("BTree");

        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        TH1F* h_rl_mc     = new TH1F("h_rl_mc"    , "", nbin_rl_nominal, rl_nominal_binning);
        TH1F* h_rl_mcreco = new TH1F("h_rl_mcreco", "", nbin_rl_nominal, rl_nominal_binning);
        TH1F* h_rl_data   = new TH1F("h_rl_data"  , "", nbin_rl_nominal, rl_nominal_binning);

        TH1F* h_weights_mc     = new TH1F("h_weights_mc"    , "", nbin_weight, weight_binning);
        TH1F* h_weights_mcreco = new TH1F("h_weights_mcreco", "", nbin_weight, weight_binning);
        TH1F* h_weights_data   = new TH1F("h_weights_data"  , "", nbin_weight, weight_binning);

        TH1F* h_eec_mc     = new TH1F("h_eec_mc"    , "", nbin_rl_nominal, rl_nominal_binning);
        TH1F* h_eec_mcreco = new TH1F("h_eec_mcreco", "", nbin_rl_nominal, rl_nominal_binning);
        TH1F* h_eec_data   = new TH1F("h_eec_data"  , "", nbin_rl_nominal, rl_nominal_binning);

        set_histogram_style(h_rl_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_rl_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_rl_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_weights_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_weights_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_weights_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        set_histogram_style(h_eec_mc, corr_marker_color_jet_pt[0], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);
        set_histogram_style(h_eec_mcreco, corr_marker_color_jet_pt[1], std_line_width-1, corr_marker_style_jet_pt[1], std_marker_size + 0.5);
        set_histogram_style(h_eec_data, corr_marker_color_jet_pt[2], std_line_width-1, corr_marker_style_jet_pt[0], std_marker_size + 0.5);

        t_mc->Project("h_rl_mc","pair_rl",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_rl_mcreco","pair_rl",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_rl_data","pair_rl",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        TString hf_y = "0.5 * log((HF_e + HF_pz)/(HF_e - HF_pz))";
        t_mc->Project("h_weights_mc","pair_weight",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_mcreco->Project("h_weights_mcreco","pair_weight",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));
        t_data->Project("h_weights_data","pair_weight",Form("jet_pt>%f&&jet_pt<%f", 12.5, 100.));

        // EEC needs special filling
        TChain *BTree = new TChain("BTree", "B-jets Tree Variables");
        BTree->Add((output_folder + "ntuple_bjets_mcreco.root").c_str());

        vector<float> *pair_rl = 0, *pair_weight = 0;
        float jet_pt;
        BTree->SetBranchAddress("pair_rl"        , &pair_rl);
        BTree->SetBranchAddress("pair_weight"    , &pair_weight);
        BTree->SetBranchAddress("jet_pt"    , &jet_pt);
        
        for (int ev = 0; ev < BTree->GetEntries(); ev++) {
                BTree->GetEntry(ev);

                if(jet_pt < 12.5)
                        continue;

                if (!pair_rl->empty()) {
                        ULong_t vector_size = pair_rl->size();

                        float *rl_info         = pair_rl->data();
                        float *weight_info     = pair_weight->data();
                        
                        for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                h_eec_mcreco->Fill(rl_info[vector_index], weight_info[vector_index]);
                        }
                }
        }

        TChain *BTree_data = new TChain("BTree", "B-jets Tree Variables");
        BTree_data->Add((output_folder + "ntuple_bjets_data.root").c_str());

        vector<float> *pair_rl_data = 0, *pair_weight_data = 0;
        float jet_pt_data;
        BTree_data->SetBranchAddress("pair_rl"    , &pair_rl_data);
        BTree_data->SetBranchAddress("pair_weight", &pair_weight_data);
        BTree_data->SetBranchAddress("jet_pt"     , &jet_pt_data);
        
        for (int ev = 0; ev < BTree_data->GetEntries(); ev++) {
                BTree_data->GetEntry(ev);

                if (jet_pt_data < 12.5)
                        continue;

                if (!pair_rl_data->empty()) {
                        ULong_t vector_size = pair_rl_data->size();

                        float *rl_info         = pair_rl_data->data();
                        float *weight_info     = pair_weight_data->data();
                        
                        for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                h_eec_data->Fill(rl_info[vector_index], weight_info[vector_index]);
                        }
                }
        }
        
        // -------------------------

        // TString hf_mass = "TMath::Sqrt((HF_e*HF_e) - (HF_px*HF_px) - (HF_py*HF_py) - (HF_pz*HF_pz))";
        // t_mc->Project("h_eec_mc","pair_rl",Form("pair_weight*(jet_pt>%f&&jet_pt<%f)", 12.5, 100.));
        // t_mcreco->Project("h_eec_mcreco","pair_rl",Form("pair_weight*(jet_pt>%f&&jet_pt<%f)", 12.5, 100.));
        // t_data->Project("h_eec_data","pair_rl",Form("pair_weight*(jet_pt>%f&&jet_pt<%f)", 12.5, 100.));

        h_rl_mc->Scale(1./h_rl_mc->Integral());
        h_rl_mcreco->Scale(1./h_rl_mcreco->Integral());
        h_rl_data->Scale(1./h_rl_data->Integral());

        h_weights_mc->Scale(1./h_weights_mc->Integral());
        h_weights_mcreco->Scale(1./h_weights_mcreco->Integral());
        h_weights_data->Scale(1./h_weights_data->Integral());

        h_eec_mc->Scale(1./h_eec_mc->Integral());
        h_eec_mcreco->Scale(1./h_eec_mcreco->Integral());
        h_eec_data->Scale(1./h_eec_data->Integral());

        TCanvas* c = new TCanvas("c", "", 800, 600);
        c->Draw();

        // Plotting jet pt
        hs->Add(h_rl_mcreco, "P E1");
        hs->Add(h_rl_data, "P E1");

        // l->AddEntry(h_rl_mc,"h_rl_mc","p");
        l->AddEntry(h_rl_mcreco,"mc(reco)","p");
        l->AddEntry(h_rl_data,"data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{R_{L}};Normalized distributions");

        l->Draw("same");

        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/kinematics_pairs_rl.pdf");

        // Plotting jet rapidity
        hs = new THStack();
        hs->Add(h_weights_mcreco, "P E1");
        hs->Add(h_weights_data, "P E1");

        // l->AddEntry(h_rl_mc,"h_rl_mc","p");
        l = new TLegend();
        l->AddEntry(h_weights_mcreco, "mc(reco)","p");
        l->AddEntry(h_weights_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";weights;Normalized distributions");

        l->Draw("same");

        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/kinematics_pairs_weights.pdf");

        // Plotting hf mass
        hs = new THStack();
        hs->Add(h_eec_mcreco, "P E1");
        hs->Add(h_eec_data, "P E1");

        // l->AddEntry(h_rl_mc,"h_rl_mc","p");
        l = new TLegend();
        l->AddEntry(h_eec_mcreco, "mc(reco)","p");
        l->AddEntry(h_eec_data, "data","p");
        
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{R_{L}};Normalized distributions");

        l->Draw("same");

        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/kinematics_eec_rl.pdf");
}
