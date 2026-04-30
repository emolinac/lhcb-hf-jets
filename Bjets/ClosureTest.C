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

void ClosureTest(int NumEvts = -1,
                 int dataset1 = 91599,
                 int dataset2 = 91599,
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

        TFile *file_write = new TFile((output_folder + "bjets_closuretest.root").c_str(), "RECREATE");

        /////////////////////   Get histograms /////////////////////////////////

        TH1D *h1_jetpt_reco = (TH1D *)file_reco->Get("Jet_pT");
        TH1D *h1_jetpt_data = (TH1D *)file_data->Get("Jet_pT");  
        h1_jetpt_reco->Draw();
        TH1D *h1_jetpt_final = (TH1D *)h1_jetpt_reco->Clone("jetpt_final");

        /////////////////////   Get Truth histograms /////////////////////////////////
        TH1D *h1_jetpt_truth = (TH1D *)file_truth->Get("Jet_pT");

        /////////////////////   Get Purity & Efficiency Hists /////////////////////////////////
        TH1D *h1_purity_jetpt = (TH1D *)file_unfold->Get("purity_jetpt");
        TH1D *h1_efficiency_jetpt = (TH1D *)file_unfold->Get("efficiency_jetpt");

        /////////////////////   Get Response Matrices /////////////////////////////////
        RooUnfoldResponse *response_jetpt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_jetpt");

        TH2 *h2_response_jetpt = (TH2 *)response_jetpt->Hresponse();

        response_jetpt->UseOverflow();

        // Multiply by purity
        h1_jetpt_final->Multiply(h1_purity_jetpt);

        // Unfold
        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_final, NumIters);
        
        h1_jetpt_final = (TH1D *)unfold_jetpt.Hreco();
        
        // Divide by efficiency
        h1_jetpt_final->Divide(h1_efficiency_jetpt);

        TH1D *h1_jetpt_final_ratio = (TH1D *)h1_jetpt_final->Clone("h1_jetpt_final_ratio");
        
        h1_jetpt_final_ratio->Divide(h1_jetpt_truth);

        h1_jetpt_reco->Write("h1_jetpt_reco");
        h1_jetpt_truth->Write("h1_jetpt_truth");
        h1_jetpt_final->Write("h1_jetpt_final");
        h1_jetpt_final_ratio->Write("pseudodata_to_truth_jetpt");  

        int binlow_jet  = h1_jetpt_truth->FindBin(ptMin);
        int binhigh_jet = h1_jetpt_truth->FindBin(ptMax);

        double Njets_reco  = h1_jetpt_reco->Integral(binlow_jet, binhigh_jet);
        double Njets_truth = h1_jetpt_truth->Integral(binlow_jet, binhigh_jet);
        double Njets_final = h1_jetpt_final->Integral(binlow_jet, binhigh_jet);
                
        ////////////////////////////////////
        // Smearing the jet pt distribution
        ///////////////////////////////////
        TH1D *h1_jetpt_closure_error;
        for (int i = 0; i < nRuns; i++) {
                TH1D *h1_jetpt_smear = (TH1D *)h1_jetpt_reco->Clone(Form("jetpt_smear%d", i));
                
                if (smear_by_data)
                        SmearObservables(h1_jetpt_smear, h1_jetpt_data, myRNG);
                
                // Multiply by purity
                h1_jetpt_smear->Multiply(h1_purity_jetpt);
                
                RooUnfoldBayes unfold_jetpt_smear(response_jetpt, h1_jetpt_smear, NumIters);      
                
                h1_jetpt_smear = (TH1D *)unfold_jetpt_smear.Hreco();

                // Divide by efficiency
                h1_jetpt_smear->Divide(h1_efficiency_jetpt);
                
                TH1D *h1_jetpt_ratio_smear = (TH1D *)h1_jetpt_reco->Clone(Form("h1_jetpt_ratio_smear%d", i));
                //TH1D *h1_jetpt_pull_smear = (TH1D *)h1_jetpt_reco->Clone(Form("h1_jetpt_pull_smear%d", i));
                h1_jetpt_ratio_smear->Divide(h1_jetpt_smear, h1_jetpt_truth);
                h1_jetpt_ratio_smear->Write();
                //h1_jetpt_pull_smear->Write(); 
        }

        std::cout << "############################## Unfolding 3D Npair distribution ##############################" << std::endl;
        
        // Get all the necessary histograms
        TH3D *h3_rl_jetpt_weight_truth = (TH3D*) file_truth->Get("h_npair_mc");
        TH3D *h3_rl_jetpt_weight       = (TH3D*) file_reco->Get("h3_rl_jetpt_weight");
        TH3D *h3_rl_jetpt_weight_data  = (TH3D*) file_data->Get("h3_rl_jetpt_weight");
        TH3D *h3_rl_jetpt_weight_final = (TH3D*) h3_rl_jetpt_weight->Clone("h3_rl_jetpt_weight_final");
        
        TH3D *h3_eff_rl_jetpt_weight    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_weight");
        TH3D *h3_purity_rl_jetpt_weight = (TH3D *)file_unfold->Get("purity_rl_jetpt_weight");
        
        RooUnfoldResponse *response_rl_jetpt_weight = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair");
        
        // Correct the distributions
        h3_rl_jetpt_weight_final->Multiply(h3_rl_jetpt_weight_final, h3_purity_rl_jetpt_weight);
        
        RooUnfoldBayes unfold_rl_jetpt_weight(response_rl_jetpt_weight, h3_rl_jetpt_weight_final, NumIters);
        
        h3_rl_jetpt_weight_final = (TH3D *)unfold_rl_jetpt_weight.Hreco();
        
        h3_rl_jetpt_weight_final->Divide(h3_rl_jetpt_weight_final, h3_eff_rl_jetpt_weight);

        // Calculate EECs from the npair 3D distributions
        TH2D *h2_eec_final = new TH2D("h2_eec_final", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_truth = new TH2D("h2_eec_truth", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec       = new TH2D("h2_eec", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        apply_unfolded_weights(h3_rl_jetpt_weight_final, h2_eec_final);
        apply_unfolded_weights(h3_rl_jetpt_weight_truth, h2_eec_truth);
        apply_unfolded_weights(h3_rl_jetpt_weight      , h2_eec);

        // Project EECs into 1D histograms
        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hmcreco_eec[ptbinsize]; 
        TH1F* h_eec_mcreco_truth_ratio[ptbinsize]; 

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                if (h1_jetpt_truth->GetBinContent(bin + 1) == 0)
                        continue;

                hmc_eec[bin]                  = new TH1F(Form("hmc_eec%i",bin)    , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hmcreco_eec[bin]              = new TH1F(Form("hmcreco_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                h_eec_mcreco_truth_ratio[bin] = new TH1F(Form("pseudodata_to_truth_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_eec[bin]    , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hmcreco_eec[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_eec_truth, hmc_eec[bin]    , bin + 1);
                project_nominal_phase_space(h2_eec_final, hmcreco_eec[bin], bin + 1);

                hmc_eec[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hmcreco_eec[bin]->Scale(1./h1_jetpt_final->GetBinContent(bin + 1),"width");

                h_eec_mcreco_truth_ratio[bin]->Divide(hmcreco_eec[bin], hmc_eec[bin]);

                hmc_eec[bin]->Write();
                hmcreco_eec[bin]->Write();
                h_eec_mcreco_truth_ratio[bin]->Write();
        }
        
        // // Smearing the Observables
        // TH1F* hmcreco_eec_smeared[ptbinsize][nRuns]; 
        // TH1F* h_eec_mcreco_truth_ratio_smeared[ptbinsize][nRuns]; 
        // TH2D* h2_eec_smeared[nRuns];

        // for (int i = 0; i < nRuns; i++) {
        //         TH3D *h3_rl_jetpt_weight_smear = (TH3D *)h3_rl_jetpt_weight->Clone(Form("rl_jetpt_weight_smeared%d", i));
                
        //         if (smear_by_data)
        //                 SmearObservables(h3_rl_jetpt_weight_smear, h3_rl_jetpt_weight_data, myRNG);
                
        //         // Correct the smeared pseudodata
        //         h3_rl_jetpt_weight_smear->Multiply(h3_rl_jetpt_weight_smear, h3_purity_rl_jetpt_weight);
        //         RooUnfoldBayes unfold_rl_jetpt_weight_smear(response_rl_jetpt_weight, h3_rl_jetpt_weight_smear, NumIters);

        //         h3_rl_jetpt_weight_smear = (TH3D *)unfold_rl_jetpt_weight_smear.Hreco();
                
        //         h3_rl_jetpt_weight_smear->Divide(h3_rl_jetpt_weight_smear, h3_eff_rl_jetpt_weight);
                
        //         // Estimate the EEC for this smeared iteration
        //         h2_eec_smeared[i] = new TH2D(Form("h2_eec_smeared%i", i), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        //         apply_unfolded_weights(h3_rl_jetpt_weight_smear, h2_eec_smeared[i]);

        //         // Estimate the mcreco/truth ratio in 1D for the smeared thing
        //         for (int bin = 0 ; bin < ptbinsize ; bin++) {
        //                 if (h1_jetpt_truth->GetBinContent(bin + 1) == 0)
        //                         continue;

        //                 hmcreco_eec_smeared[bin][i]              = new TH1F(Form("hmcreco_eec%i_smeared%i",bin,i)            , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
        //                 h_eec_mcreco_truth_ratio_smeared[bin][i] = new TH1F(Form("pseudodata_to_truth_eec%i_smeared%i",bin,i), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        
        //                 set_histogram_style(hmcreco_eec_smeared[bin][i], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

        //                 project_nominal_phase_space(h2_eec_smeared[i], hmcreco_eec_smeared[bin][i], bin + 1);

        //                 hmcreco_eec_smeared[bin][i]->Scale(1./h1_jetpt_final->GetBinContent(bin + 1),"width");

        //                 h_eec_mcreco_truth_ratio_smeared[bin][i]->Divide(hmcreco_eec_smeared[bin][i], hmc_eec[bin]);

        //                 hmcreco_eec_smeared[bin][i]->Write();
        //                 h_eec_mcreco_truth_ratio_smeared[bin][i]->Write();
        //         }
        // }
        
        //file_write->Write();
        file_write->Close();
}
