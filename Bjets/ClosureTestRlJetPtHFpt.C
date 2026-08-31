#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"

#include "../Helpers_IC.h"

#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/TBJetsMC.h"
#include "../include/TBJetsMC.C"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

using namespace std;

int niter_npair = 2;

void ClosureTestRlJetPtHFpt(int NumEvts = -1,
                 int NumIters = 4,  
                 bool DoShapeClosure = false)
{
        gStyle->SetPaintTextFormat("3.3f");
        gROOT->ForceStyle();
        gStyle->SetOptStat(0);

        TRandom3 *myRNG = new TRandom3(0);

        const int nRuns =  DoShapeClosure ? 1 : 1;

        bool smear_by_data = (!DoShapeClosure);  
        
        const int fixSmear = 42;

        // int NumIters = 1;
        // int RegIDS = 5;

        /////////////////////   Get Files /////////////////////////////////

        TFile *file_reco = new TFile((output_folder + "bjets_simpleobservable_mcreco.root").c_str(), "READ"); 
        TFile *file_data = new TFile((output_folder + "bjets_simpleobservable_data.root").c_str(), "READ");
        TFile *file_truth = new TFile((output_folder + "bjets_simpleobservable_mc.root").c_str(), "READ"); 
        TFile *file_unfold = new TFile((output_folder + "bjets_corrections.root").c_str(), "READ"); 

        TFile *file_write = new TFile((output_folder + "bjets_closuretest_npairs_rljetpthfpt.root").c_str(), "RECREATE");

        std::cout << "############################## Unfolding 1D Jet distribution ##############################" << std::endl;

        // 1D JET CORRECTION
        /////////////////////   Get histograms /////////////////////////////////

        TH1D *h1_jetpt_reco = (TH1D *)file_reco->Get("Jet_pT");
        TH1D *h1_jetpt_data = (TH1D *)file_data->Get("Jet_pT");  

        std::cout<<"MCReco jet pt entries = "<<h1_jetpt_reco->Integral()<<std::endl;
        std::cout<<"Data   jet pt entries = "<<h1_jetpt_data->Integral()<<std::endl;

        // h1_jetpt_reco->Draw();
        TH1D *h1_jetpt_pseudodata = (TH1D *)h1_jetpt_reco->Clone("jetpt_pseudodata");

        /////////////////////   Get Truth histograms /////////////////////////////////
        TH1D *h1_jetpt_truth = (TH1D *)file_truth->Get("Jet_pT");

        /////////////////////   Get Purity & Efficiency Hists /////////////////////////////////
        TH1F *h1_purity_jetpt     = (TH1F *)file_unfold->Get("purity_jetpt");
        TH1F *h1_efficiency_jetpt = (TH1F *)file_unfold->Get("efficiency_jetpt");

        TH1F *h1_purity_HFpt     = (TH1F *)file_unfold->Get("purity_HFpt");
        TH1F *h1_efficiency_HFpt = (TH1F *)file_unfold->Get("efficiency_HFpt");

        /////////////////////   Get Response Matrices /////////////////////////////////
        RooUnfoldResponse *response_jetpt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_jetpt");

        TH2 *h2_response_jetpt = (TH2 *)response_jetpt->Hresponse();

        // response_jetpt->UseOverflow();

        h1_jetpt_pseudodata->Write("jetpt_uncorr");

        // Multiply by purity
        h1_jetpt_pseudodata->Multiply(h1_purity_jetpt);

        h1_jetpt_pseudodata->Write("jetpt_pur");
        
        // Unfold
        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_pseudodata, NumIters);
        
        h1_jetpt_pseudodata = (TH1D *)unfold_jetpt.Hreco();

        h1_jetpt_pseudodata->Write("jetpt_pur_unf");

        // Divide by efficiency
        h1_jetpt_pseudodata->Divide(h1_efficiency_jetpt);

        h1_jetpt_pseudodata->Write("jetpt_pur_unf_eff");

        TH1D *h1_jetpt_pseudodata_ratio = (TH1D *)h1_jetpt_pseudodata->Clone("h1_jetpt_pseudodata_ratio");
        
        h1_jetpt_pseudodata_ratio->Divide(h1_jetpt_truth);

        h1_jetpt_reco->Write("h1_jetpt_reco");
        h1_jetpt_truth->Write("h1_jetpt_truth");
        h1_jetpt_pseudodata->Write("h1_jetpt_pseudodata");
        h1_jetpt_pseudodata_ratio->Write("pseudodata_to_truth_jetpt");  
        h1_purity_jetpt->Write("purity_jetpt");
        h1_efficiency_jetpt->Write("efficiency_jetpt");
        h2_response_jetpt->Write("response_jetpt");

        std::cout << "############################## Unfolding 3D Jet distribution ##############################" << std::endl;

        /////////////////////   Get histograms /////////////////////////////////
        TH3D *h3_jet_reco = (TH3D *)file_reco->Get("h3_HFptetajetpt");
        TH3D *h3_jet_data = (TH3D *)file_data->Get("h3_HFptetajetpt");  

        std::cout<<"MCReco jet pt entries = "<<h3_jet_reco->Integral()<<std::endl;
        std::cout<<"Data   jet pt entries = "<<h3_jet_data->Integral()<<std::endl;

        // h1_jetpt_reco->Draw();
        TH3D *h3_jet_pseudodata = (TH3D *)h3_jet_reco->Clone("h3_jet_pseudodata");

        /////////////////////   Get Truth histograms /////////////////////////////////
        TH3D *h3_jet_truth = (TH3D *)file_truth->Get("h3_HFptetajetpt");

        /////////////////////   Get Purity & Efficiency Hists /////////////////////////////////
        TH3F *h3_purity_HFptetajetpt     = (TH3F *)file_unfold->Get("purity_HFptetajetpt");
        TH3F *h3_efficiency_HFptetajetpt = (TH3F *)file_unfold->Get("efficiency_HFptetajetpt");

        /////////////////////   Get Response Matrices /////////////////////////////////
        RooUnfoldResponse *response_jet = (RooUnfoldResponse *)file_unfold->Get("Roo_response_HFptetajetpt");

        // response_jetpt->UseOverflow();

        h3_jet_pseudodata->Multiply(h3_purity_HFptetajetpt);

        RooUnfoldBayes unfold_jet(response_jet, h3_jet_pseudodata, NumIters);
        
        h3_jet_pseudodata = (TH3D *)unfold_jet.Hreco();

        h3_jet_pseudodata->Divide(h3_efficiency_HFptetajetpt);

        TH1F* h1_jet_3dcorr = new TH1F("h1_jet_3dcorr", "", ptbinsize, pt_binedges);

        h1_jet_3dcorr = (TH1F*) h3_jet_pseudodata->ProjectionZ();

        h1_jet_3dcorr->Write("h1_jet_3dcorr");

        TH1D *h1_jet_pseudodata_ratio = (TH1D *)h1_jet_3dcorr->Clone("h1_jet_pseudodata_ratio");
        
        h1_jet_pseudodata_ratio->Divide(h1_jetpt_truth);

        h1_jet_pseudodata_ratio->Write("h1_jet_3dcorr_ratio");

        std::cout << "############################## Unfolding 3D Npair distribution ##############################" << std::endl;
        
        // Get all the necessary histograms
        TH3D *h3_rl_jetpt_weight_truth      = (TH3D*) file_truth->Get("h_npair_mc");
        TH3D *h3_rl_jetpt_weight            = (TH3D*) file_reco->Get("h3_rl_jetpt_HFpt");
        TH3D *h3_rl_jetpt_weight_pseudodata = (TH3D*) h3_rl_jetpt_weight->Clone("h3_rl_jetpt_weight_pseudodata");
        
        TH3D *h3_eff_rl_jetpt_weight    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_HFpt");
        TH3D *h3_purity_rl_jetpt_weight = (TH3D *)file_unfold->Get("purity_rl_jetpt_HFpt");
        
        RooUnfoldResponse *response_rl_jetpt_weight = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair_HFpt");
        
        // Correct the distributions
        h3_rl_jetpt_weight_pseudodata->Multiply(h3_purity_rl_jetpt_weight);
        
        RooUnfoldBayes unfold_rl_jetpt_weight(response_rl_jetpt_weight, h3_rl_jetpt_weight_pseudodata, niter_npair);

        h3_rl_jetpt_weight_pseudodata = (TH3D *)unfold_rl_jetpt_weight.Hreco();
        
        h3_rl_jetpt_weight_pseudodata->Divide(h3_eff_rl_jetpt_weight);

        // Calculate EECs from the npair 3D distributions
        TH2D *h2_eec_pseudodata = new TH2D("h2_eec_pseudodata", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_truth      = new TH2D("h2_eec_truth"     , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec            = new TH2D("h2_eec"           , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        h2_eec_pseudodata = (TH2D*) h3_rl_jetpt_weight_pseudodata->Project3D("yx");
        h2_eec_truth      = (TH2D*) h3_rl_jetpt_weight_truth->Project3D("yx");
        h2_eec            = (TH2D*) h3_rl_jetpt_weight->Project3D("yx");

        // Project EECs into 1D histograms
        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hmcreco_eec[ptbinsize]; 
        TH1F* h_eec_mcreco_truth_ratio[ptbinsize]; 

        TH1F* hmcreco_eec_jet3dcorr[ptbinsize]; 
        TH1F* h_eec_mcreco_jet3dcorr_truth_ratio[ptbinsize]; 

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                if (h1_jetpt_truth->GetBinContent(bin + 1) == 0)
                        continue;

                hmc_eec[bin]                  = new TH1F(Form("hmc_eec%i",bin)    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hmcreco_eec[bin]              = new TH1F(Form("hmcreco_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                h_eec_mcreco_truth_ratio[bin] = new TH1F(Form("pseudodata_to_truth_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                hmcreco_eec_jet3dcorr[bin]              = new TH1F(Form("hmcreco_jet3dcorr_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                h_eec_mcreco_jet3dcorr_truth_ratio[bin] = new TH1F(Form("pseudodata_jet3dcorr_to_truth_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_eec[bin]              , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hmcreco_eec[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hmcreco_eec_jet3dcorr[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_eec_truth     , hmc_eec[bin]    , bin + 1);
                project_nominal_phase_space(h2_eec_pseudodata, hmcreco_eec[bin], bin + 1);
                project_nominal_phase_space(h2_eec_pseudodata, hmcreco_eec_jet3dcorr[bin], bin + 1);

                hmc_eec[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hmcreco_eec[bin]->Scale(1./h1_jetpt_pseudodata->GetBinContent(bin + 1),"width");
                hmcreco_eec_jet3dcorr[bin]->Scale(1./h1_jet_3dcorr->GetBinContent(bin + 1),"width");

                h_eec_mcreco_truth_ratio[bin]->Divide(hmcreco_eec[bin], hmc_eec[bin]);
                h_eec_mcreco_jet3dcorr_truth_ratio[bin]->Divide(hmcreco_eec_jet3dcorr[bin], hmc_eec[bin]);

                hmc_eec[bin]->Write();
                hmcreco_eec[bin]->Write();
                h_eec_mcreco_truth_ratio[bin]->Write();

                hmcreco_eec_jet3dcorr[bin]->Write();
                h_eec_mcreco_jet3dcorr_truth_ratio[bin]->Write();
        }
        
        file_write->Close();
}
