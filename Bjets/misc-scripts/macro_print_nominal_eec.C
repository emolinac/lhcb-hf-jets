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

int niter_npair = 4;

void macro_print_nominal_eec()
{
        /////////////////////   Get Files /////////////////////////////////

        TFile *file_data = new TFile((output_folder + "bjets_3duncorrecteddistributions_data.root").c_str(), "READ");
        TFile *file_truth = new TFile((output_folder + "bjets_3duncorrecteddistributions_mc.root").c_str(), "READ"); 
        TFile *file_unfold = new TFile((output_folder + "bjets_corrections.root").c_str(), "READ"); 

        TFile *file_write = new TFile((output_folder + "bjets_nominal_eec.root").c_str(), "RECREATE");

        /////////////////////   Get histograms /////////////////////////////////

        TH1D *h1_jetpt_data = (TH1D *)file_data->Get("Jet_pT");  

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

        h1_jetpt_data->Write("jetpt_uncorr");

        // Multiply by purity
        h1_jetpt_data->Multiply(h1_purity_jetpt);

        h1_jetpt_data->Write("jetpt_pur");
        
        // Unfold
        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_data, 4);
        
        h1_jetpt_data = (TH1D *)unfold_jetpt.Hreco();

        h1_jetpt_data->Write("jetpt_pur_unf");

        // Divide by efficiency
        h1_jetpt_data->Divide(h1_efficiency_jetpt);

        h1_jetpt_data->Write("jetpt_pur_unf_eff");

        h1_jetpt_truth->Write("h1_jetpt_truth");
        h1_jetpt_data->Write("h1_jetpt_data");
        h1_purity_jetpt->Write("purity_jetpt");
        h1_efficiency_jetpt->Write("efficiency_jetpt");
        h2_response_jetpt->Write("response_jetpt");

        std::cout << "############################## Unfolding 3D Npair distribution ##############################" << std::endl;
        
        // Get all the necessary histograms
        TH3D *h3_rl_jetpt_weight_truth = (TH3D*) file_truth->Get("h_npair_mc");
        TH3D *h3_rl_jetpt_weight_data  = (TH3D*) file_data->Get("h3_rl_jetpt_weight");
        
        TH3D *h3_eff_rl_jetpt_weight    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_weight");
        TH3D *h3_purity_rl_jetpt_weight = (TH3D *)file_unfold->Get("purity_rl_jetpt_weight");
        
        RooUnfoldResponse *response_rl_jetpt_weight = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair");
        
        // Correct the distributions
        h3_rl_jetpt_weight_data->Multiply(h3_purity_rl_jetpt_weight);
        
        RooUnfoldBayes unfold_rl_jetpt_weight(response_rl_jetpt_weight, h3_rl_jetpt_weight_data, 7);

        h3_rl_jetpt_weight_data = (TH3D *)unfold_rl_jetpt_weight.Hreco();
        
        h3_rl_jetpt_weight_data->Divide(h3_eff_rl_jetpt_weight);

        std::cout << "############################## Unfolding 3D Npair with HF distribution ##############################" << std::endl;
        
        // Get all the necessary histograms
        TH3D *h3_rl_jetpt_weight_truth_whf      = (TH3D*) file_truth->Get("h_npair_whf_mc");
        TH3D *h3_rl_jetpt_weight_data_whf       = (TH3D*) file_data->Get("h3_rl_jetpt_weight_whf");
        
        TH3D *h3_eff_rl_jetpt_weight_whf    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_weight_whf");
        TH3D *h3_purity_rl_jetpt_weight_whf = (TH3D *)file_unfold->Get("purity_rl_jetpt_weight_whf");
        
        RooUnfoldResponse *response_rl_jetpt_weight_whf = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair_whf");
        
        // Correct the distributions
        h3_rl_jetpt_weight_data_whf->Multiply(h3_purity_rl_jetpt_weight_whf);
        
        RooUnfoldBayes unfold_rl_jetpt_weight_whf(response_rl_jetpt_weight_whf, h3_rl_jetpt_weight_data_whf, niter_npair);

        h3_rl_jetpt_weight_data_whf = (TH3D *)unfold_rl_jetpt_weight_whf.Hreco();
        
        h3_rl_jetpt_weight_data_whf->Divide(h3_eff_rl_jetpt_weight_whf);

        std::cout << "############################## Unfolding 3D Npair without HF distribution ##############################" << std::endl;
        
        // Get all the necessary histograms
        TH3D *h3_rl_jetpt_weight_truth_wohf      = (TH3D*) file_truth->Get("h_npair_wohf_mc");
        TH3D *h3_rl_jetpt_weight_data_wohf       = (TH3D*) file_data->Get("h3_rl_jetpt_weight_wohf");
        
        TH3D *h3_eff_rl_jetpt_weight_wohf    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_weight_wohf");
        TH3D *h3_purity_rl_jetpt_weight_wohf = (TH3D *)file_unfold->Get("purity_rl_jetpt_weight_wohf");
        
        RooUnfoldResponse *response_rl_jetpt_weight_wohf = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair_wohf");
        
        // Correct the distributions
        h3_rl_jetpt_weight_data_wohf->Multiply(h3_purity_rl_jetpt_weight_wohf);
        
        RooUnfoldBayes unfold_rl_jetpt_weight_wohf(response_rl_jetpt_weight_wohf, h3_rl_jetpt_weight_data_wohf, niter_npair);

        h3_rl_jetpt_weight_data_wohf = (TH3D *)unfold_rl_jetpt_weight_wohf.Hreco();
        
        h3_rl_jetpt_weight_data_wohf->Divide(h3_eff_rl_jetpt_weight_wohf);

        // Calculate EECs from the npair 3D distributions
        TH2D *h2_eec_data  = new TH2D("h2_eec_data" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_truth = new TH2D("h2_eec_truth", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        apply_unfolded_weights(h3_rl_jetpt_weight_data  , h2_eec_data);
        apply_unfolded_weights(h3_rl_jetpt_weight_truth , h2_eec_truth);

        TH2D *h2_eec_data_whf  = new TH2D("h2_eec_data_whf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_truth_whf = new TH2D("h2_eec_truth_whf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        apply_unfolded_weights(h3_rl_jetpt_weight_data_whf  , h2_eec_data_whf);
        apply_unfolded_weights(h3_rl_jetpt_weight_truth_whf , h2_eec_truth_whf);

        TH2D *h2_eec_data_wohf  = new TH2D("h2_eec_data_wohf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_truth_wohf = new TH2D("h2_eec_truth_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        apply_unfolded_weights(h3_rl_jetpt_weight_data_wohf  , h2_eec_data_wohf);
        apply_unfolded_weights(h3_rl_jetpt_weight_truth_wohf , h2_eec_truth_wohf);

        // Calculate Npair 2D from the npair 3D distributions
        TH2D *h2_npair_data  = new TH2D("h2_npair_data" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_npair_truth = new TH2D("h2_npair_truth", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        h2_npair_data  = (TH2D*) h3_rl_jetpt_weight_data->Project3D("yx");
        h2_npair_truth = (TH2D*) h3_rl_jetpt_weight_truth->Project3D("yx");
        
        TH2D *h2_npair_data_whf  = new TH2D("h2_npair_data_whf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_npair_truth_whf = new TH2D("h2_npair_truth_whf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        h2_npair_data_whf  = (TH2D*) h3_rl_jetpt_weight_data_whf->Project3D("yx");
        h2_npair_truth_whf = (TH2D*) h3_rl_jetpt_weight_truth_whf->Project3D("yx");
        
        TH2D *h2_npair_data_wohf  = new TH2D("h2_npair_data_wohf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_npair_truth_wohf = new TH2D("h2_npair_truth_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        h2_npair_data_wohf  = (TH2D*) h3_rl_jetpt_weight_data_wohf->Project3D("yx");
        h2_npair_truth_wohf = (TH2D*) h3_rl_jetpt_weight_truth_wohf->Project3D("yx");
        
        // Project EECs into 1D histograms
        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hdata_eec[ptbinsize]; 
        TH1F* h_eec_data_truth_ratio[ptbinsize]; 

        TH1F* hdata_eec_corr[ptbinsize]; 
        TH1F* h_eec_data_truth_ratio_corr[ptbinsize]; 

        TH1F* hmc_eec_whf[ptbinsize]; 
        TH1F* hdata_eec_whf[ptbinsize]; 
        TH1F* h_eec_data_truth_ratio_whf[ptbinsize]; 

        TH1F* hmc_eec_wohf[ptbinsize]; 
        TH1F* hdata_eec_wohf[ptbinsize]; 
        TH1F* h_eec_data_truth_ratio_wohf[ptbinsize]; 

        TH1F* hmc_npair[ptbinsize]; 
        TH1F* hdata_npair[ptbinsize]; 
        TH1F* hdata_npair_corr[ptbinsize]; 
        TH1F* h_npair_data_truth_ratio[ptbinsize]; 

        TH1F* hmc_npair_whf[ptbinsize]; 
        TH1F* hdata_npair_whf[ptbinsize]; 
        TH1F* h_npair_data_truth_ratio_whf[ptbinsize]; 

        TH1F* hmc_npair_wohf[ptbinsize]; 
        TH1F* hdata_npair_wohf[ptbinsize]; 
        TH1F* h_npair_data_truth_ratio_wohf[ptbinsize]; 

        std::cout<<"Calculating per bin"<<std::endl;
        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                if (h1_jetpt_truth->GetBinContent(bin + 1) == 0)
                        continue;

                // EECs
                hmc_eec[bin]   = new TH1F(Form("hmc_eec%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec[bin] = new TH1F(Form("hdata_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec_corr[bin] = new TH1F(Form("hdata_eec_corr%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_eec[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec_corr[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_eec_truth , hmc_eec[bin]  , bin + 1);
                project_nominal_phase_space(h2_eec_data  , hdata_eec[bin], bin + 1);

                hmc_eec[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_eec[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");

                hmc_eec_whf[bin]   = new TH1F(Form("hmc_eec_whf%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec_whf[bin] = new TH1F(Form("hdata_eec_whf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_eec_whf[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec_whf[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_eec_truth_whf , hmc_eec_whf[bin]  , bin + 1);
                project_nominal_phase_space(h2_eec_data_whf  , hdata_eec_whf[bin], bin + 1);

                hmc_eec_whf[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_eec_whf[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");

                hmc_eec_wohf[bin]   = new TH1F(Form("hmc_eec_wohf%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec_wohf[bin] = new TH1F(Form("hdata_eec_wohf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_eec_wohf[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec_wohf[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_eec_truth_wohf , hmc_eec_wohf[bin]  , bin + 1);
                project_nominal_phase_space(h2_eec_data_wohf  , hdata_eec_wohf[bin], bin + 1);

                hmc_eec_wohf[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_eec_wohf[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");

                hdata_eec_corr[bin]->Add(hdata_eec_whf[bin], hdata_eec_wohf[bin], 1, 1);

                // Npairs
                hmc_npair[bin]   = new TH1F(Form("hmc_npair%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_npair[bin] = new TH1F(Form("hdata_npair%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_npair_corr[bin] = new TH1F(Form("hdata_npair_corr%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_npair[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_npair[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_npair_corr[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_npair_truth , hmc_npair[bin]  , bin + 1);
                project_nominal_phase_space(h2_npair_data  , hdata_npair[bin], bin + 1);

                hmc_npair[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_npair[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");

                hmc_npair_whf[bin]   = new TH1F(Form("hmc_npair_whf%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_npair_whf[bin] = new TH1F(Form("hdata_npair_whf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_npair_whf[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_npair_whf[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_npair_truth_whf , hmc_npair_whf[bin]  , bin + 1);
                project_nominal_phase_space(h2_npair_data_whf  , hdata_npair_whf[bin], bin + 1);

                hmc_npair_whf[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_npair_whf[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");

                hmc_npair_wohf[bin]   = new TH1F(Form("hmc_npair_wohf%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_npair_wohf[bin] = new TH1F(Form("hdata_npair_wohf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_npair_wohf[bin]  , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_npair_wohf[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_npair_truth_wohf , hmc_npair_wohf[bin]  , bin + 1);
                project_nominal_phase_space(h2_npair_data_wohf  , hdata_npair_wohf[bin], bin + 1);

                hmc_npair_wohf[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_npair_wohf[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");

                hdata_npair_corr[bin]->Add(hdata_npair_whf[bin], hdata_npair_wohf[bin], 1, 1);

                hmc_eec[bin]->Write();
                hdata_eec[bin]->Write();
                hdata_eec_corr[bin]->Write();
                hmc_npair[bin]->Write();
                hdata_npair[bin]->Write();
                hdata_npair_corr[bin]->Write();

                hmc_eec_whf[bin]->Write();
                hdata_eec_whf[bin]->Write();
                hmc_npair_whf[bin]->Write();
                hdata_npair_whf[bin]->Write();

                hmc_eec_wohf[bin]->Write();
                hdata_eec_wohf[bin]->Write();
                hmc_npair_wohf[bin]->Write();
                hdata_npair_wohf[bin]->Write();
        }

        TH1F* hcorr_tau[nbin_jet_pt];        
        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                double avge_pt2_jet = jet_pt_bjet_data_avge[i];
                
                get_tau_binning_from_eec_binning(tau_binning[i], rl_nominal_binning, avge_pt2_jet);
                
                hcorr_tau[i] = new TH1F(Form("hcorr_tau%i",i),"",nbin_rl_nominal,tau_binning[i]);
                get_tau_from_uoflow_eec(hdata_eec[i + 4], hcorr_tau[i], avge_pt2_jet);
                
                set_histogram_style(hcorr_tau[i], corr_marker_color_jet_pt[i], std_line_width, corr_marker_style_jet_pt[i], std_marker_size + 1);

                file_write->cd();
                hcorr_tau[i]->Write();
                gROOT->cd();
        }

        file_write->Close();
}
