#include "./include/analysis-constants.h"
#include "./include/analysis-binning.h"
#include "./include/analysis-cuts.cpp"
#include "./include/analysis-cuts.h"
#include "./include/utils.cpp"
#include "./include/utils.h"
#include "./include/utils-visual.cpp"
#include "./include/utils-visual.h"
#include "./include/directories.h"

void macro_print_mc_eec()
{
        TFile* fin = new TFile((output_folder + "ntuple_mc_bjets.root").c_str());

        TTree* tree = (TTree*) fin->Get("BTree");

        TNtuple* ntuple_jet  = (TNtuple*) fin->Get("ntuple_true_jet");
        TNtuple* ntuple_pair = (TNtuple*) fin->Get("ntuple_true_pair");

        TH1F* h1 = new TH1F("h1","",nbin_rl_nominal, rl_nominal_binning);
        TH1F* h2 = new TH1F("h2","",nbin_rl_nominal, rl_nominal_binning);
        TH1F* h3 = new TH1F("h3","",nbin_rl_nominal, rl_nominal_binning);

        TH1F* h1_wta = new TH1F("h1_wta","",nbin_rl_nominal, rl_nominal_binning);
        TH1F* h2_wta = new TH1F("h2_wta","",nbin_rl_nominal, rl_nominal_binning);
        TH1F* h3_wta = new TH1F("h3_wta","",nbin_rl_nominal, rl_nominal_binning);

        TH1F* h_jet     = new TH1F("h_jet", "", nbin_jet_pt_corrections, jet_pt_binning);
        TH1F* h_jet_wta = new TH1F("h_jet_wta", "", nbin_jet_pt_corrections, jet_pt_binning);
        
        ntuple_jet->Project("h_jet", "jet_pt");
        ntuple_pair->Project("h1", "R_L", "weight_pt*(jet_pt>20&&jet_pt<30)");
        ntuple_pair->Project("h2", "R_L", "weight_pt*(jet_pt>30&&jet_pt<50)");
        ntuple_pair->Project("h3", "R_L", "weight_pt*(jet_pt>50&&jet_pt<100)");
        
        ntuple_jet->Project("h_jet_wta", "jet_pt", "wta_distance<0.005");
        ntuple_pair->Project("h1_wta", "R_L", "weight_pt*(jet_pt>20&&jet_pt<30&&wta_distance<0.005)");
        ntuple_pair->Project("h2_wta", "R_L", "weight_pt*(jet_pt>30&&jet_pt<50&&wta_distance<0.005)");
        ntuple_pair->Project("h3_wta", "R_L", "weight_pt*(jet_pt>50&&jet_pt<100&&wta_distance<0.005)");
        
        set_histogram_style(h1, corr_marker_color_jet_pt[0], std_line_width, corr_marker_style_jet_pt[0], std_marker_size);
        set_histogram_style(h2, corr_marker_color_jet_pt[1], std_line_width, corr_marker_style_jet_pt[1], std_marker_size);
        set_histogram_style(h3, corr_marker_color_jet_pt[2], std_line_width, corr_marker_style_jet_pt[2], std_marker_size);

        set_histogram_style(h1_wta, std_marker_color_jet_pt[0], std_line_width, std_marker_style_jet_pt[0], std_marker_size);
        set_histogram_style(h2_wta, std_marker_color_jet_pt[1], std_line_width, std_marker_style_jet_pt[1], std_marker_size);
        set_histogram_style(h3_wta, std_marker_color_jet_pt[2], std_line_width, std_marker_style_jet_pt[2], std_marker_size);

        h1->Scale(1./h_jet->GetBinContent(1), "width");
        h2->Scale(1./h_jet->GetBinContent(2), "width");
        h3->Scale(1./h_jet->GetBinContent(3), "width");

        h1_wta->Scale(1./h_jet_wta->GetBinContent(1), "width");
        h2_wta->Scale(1./h_jet_wta->GetBinContent(2), "width");
        h3_wta->Scale(1./h_jet_wta->GetBinContent(3), "width");

        THStack* hs = new THStack();
        hs->Add(h1,"PE1");
        hs->Add(h2,"PE1");
        hs->Add(h3,"PE1");
        hs->Add(h1_wta,"PE1");
        hs->Add(h2_wta,"PE1");
        hs->Add(h3_wta,"PE1");

        hs->SetMaximum(0.35);

        TCanvas* c = new TCanvas();
        c->Draw();

        hs->Draw("NOSTACK");
        hs->SetTitle(";R_{L};#Sigma_{EEC}");

        gPad->SetLogx(1);

        c->Print("./mc_eec.pdf");
}