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
#include "../include/names.h"
#include "../include/TBJetsMC.h"
#include "../include/TBJetsMC.C"
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

using namespace std;

void MakeCorrectedObservable(int niter_njet = 4, int niter_npair = 4, std::string variation = "nominal")
{
        if(gSystem->AccessPathName((output_folder + namef_3duncorrecteddistributions_data[variation]).c_str())) {
                std::cout<<"Data file not found for smearing. Check file or variation given as input."<<std::endl;

                return;
        }

        if(gSystem->AccessPathName((output_folder + namef_corrections[variation]).c_str())) {
                std::cout<<"Corrections file not found. Check file or variation given as input."<<std::endl;

                return;
        }

        if (variation == "regpar-upper")
                niter_npair++;
        if (variation == "regpar-lower")
                niter_npair--;

        // Necessary files
        TFile *file_data   = new TFile((output_folder + namef_3duncorrecteddistributions_data[variation]).c_str(), "READ");
        TFile *file_truth  = new TFile((output_folder + "bjets_3duncorrecteddistributions_mc.root").c_str(), "READ"); 
        TFile *file_unfold = new TFile((output_folder + namef_corrections[variation]).c_str(), "READ"); 

        TFile *file_write = new TFile((output_folder + namef_correctedobservable_data[variation]).c_str(), "RECREATE");

        std::cout << "############################## Unfolding 1D Jet distribution ##############################" << std::endl;
        TH1D *h1_jetpt_data = (TH1D *)file_data->Get("Jet_pT");  

        TH1F *h1_purity_jetpt     = (TH1F*) file_unfold->Get("purity_jetpt");
        TH1F *h1_efficiency_jetpt = (TH1F*) file_unfold->Get("efficiency_jetpt");
        TH1F *h1_purity_HFpt      = (TH1F*) file_unfold->Get("purity_HFpt");
        TH1F *h1_efficiency_HFpt  = (TH1F*) file_unfold->Get("efficiency_HFpt");

        TH1D *h1_jetpt_truth      = (TH1D*) file_truth->Get("Jet_pT");
        
        RooUnfoldResponse* response_jetpt = (RooUnfoldResponse*) file_unfold->Get("Roo_response_jetpt");

        TH2 *h2_response_jetpt = (TH2*) response_jetpt->Hresponse();

        h1_jetpt_data->Multiply(h1_purity_jetpt);

        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_data, niter_njet);
        
        h1_jetpt_data = (TH1D *)unfold_jetpt.Hreco();

        h1_jetpt_data->Divide(h1_efficiency_jetpt);

        h1_purity_jetpt->Write("purity_jetpt");
        h1_efficiency_jetpt->Write("efficiency_jetpt");
        h2_response_jetpt->Write("response_jetpt");
        h1_jetpt_truth->Write("h1_jetpt_truth");
        h1_jetpt_data->Write("h1_jetpt_corrected");
        
        std::cout << "############################## Unfolding 3D Jet distribution ##############################" << std::endl;
        TH3D* h3_jet_data                = (TH3D*) file_data->Get("h3_HFptetajetpt");  
        TH3D* h3_jet_truth               = (TH3D*) file_truth->Get("h3_HFptetajetpt");
        TH3F* h3_purity_HFptetajetpt     = (TH3F*) file_unfold->Get("purity_HFptetajetpt");
        TH3F* h3_efficiency_HFptetajetpt = (TH3F*) file_unfold->Get("efficiency_HFptetajetpt");

        RooUnfoldResponse *response_jet = (RooUnfoldResponse *)file_unfold->Get("Roo_response_HFptetajetpt");

        h3_jet_data->Multiply(h3_purity_HFptetajetpt);

        RooUnfoldBayes unfold_jet(response_jet, h3_jet_data, niter_njet);
        
        h3_jet_data = (TH3D *)unfold_jet.Hreco();

        h3_jet_data->Divide(h3_efficiency_HFptetajetpt);

        TH1F* h1_jet_3dcorr = new TH1F("h1_jet_3dcorr", "", ptbinsize, pt_binedges);

        h1_jet_3dcorr = (TH1F*) h3_jet_data->ProjectionZ();

        h1_jet_3dcorr->Write("h1_jet_3dcorr");

        std::cout << "############################## Unfolding 3D Npair distribution ##############################" << std::endl;        
        TH3D* h3_rl_jetpt_weight_truth  = (TH3D*) file_truth->Get("h_npair_mc");
        TH3D* h3_rl_jetpt_weight   = (TH3D*) file_data->Get("h3_rl_jetpt_weight");
        
        TH3D* h3_eff_rl_jetpt_weight    = (TH3D*) file_unfold->Get("efficiency_rl_jetpt_weight");
        TH3D* h3_purity_rl_jetpt_weight = (TH3D*) file_unfold->Get("purity_rl_jetpt_weight");
        
        RooUnfoldResponse *response_rl_jetpt_weight = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair");
        
        h3_rl_jetpt_weight->Multiply(h3_purity_rl_jetpt_weight);
        
        RooUnfoldBayes unfold_rl_jetpt_weight(response_rl_jetpt_weight, h3_rl_jetpt_weight, niter_npair);

        h3_rl_jetpt_weight = (TH3D *)unfold_rl_jetpt_weight.Hreco();
        
        h3_rl_jetpt_weight->Divide(h3_eff_rl_jetpt_weight);

        std::cout << "############################## Unfolding 3D Npair with HF distribution ##############################" << std::endl;
        TH3D *h3_rl_jetpt_weight_truth_whf  = (TH3D*) file_truth->Get("h_npair_whf_mc");
        TH3D *h3_rl_jetpt_weight_whf   = (TH3D*) file_data->Get("h3_rl_jetpt_weight_whf");
        
        TH3D *h3_eff_rl_jetpt_weight_whf    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_weight_whf");
        TH3D *h3_purity_rl_jetpt_weight_whf = (TH3D *)file_unfold->Get("purity_rl_jetpt_weight_whf");
        
        RooUnfoldResponse *response_rl_jetpt_weight_whf = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair_whf");
        
        h3_rl_jetpt_weight_whf->Multiply(h3_purity_rl_jetpt_weight_whf);
        
        RooUnfoldBayes unfold_rl_jetpt_weight_whf(response_rl_jetpt_weight_whf, h3_rl_jetpt_weight_whf, niter_npair);

        h3_rl_jetpt_weight_whf = (TH3D *)unfold_rl_jetpt_weight_whf.Hreco();
        
        h3_rl_jetpt_weight_whf->Divide(h3_eff_rl_jetpt_weight_whf);

        std::cout << "############################## Unfolding 3D Npair without HF distribution ##############################" << std::endl;
        TH3D *h3_rl_jetpt_weight_truth_wohf  = (TH3D*) file_truth->Get("h_npair_wohf_mc");
        TH3D *h3_rl_jetpt_weight_wohf   = (TH3D*) file_data->Get("h3_rl_jetpt_weight_wohf");
        
        TH3D *h3_eff_rl_jetpt_weight_wohf    = (TH3D *)file_unfold->Get("efficiency_rl_jetpt_weight_wohf");
        TH3D *h3_purity_rl_jetpt_weight_wohf = (TH3D *)file_unfold->Get("purity_rl_jetpt_weight_wohf");
        
        RooUnfoldResponse *response_rl_jetpt_weight_wohf = (RooUnfoldResponse *)file_unfold->Get("Roo_response_npair_wohf");
        
        // Correct the distributions
        h3_rl_jetpt_weight_wohf->Multiply(h3_purity_rl_jetpt_weight_wohf);
        
        RooUnfoldBayes unfold_rl_jetpt_weight_wohf(response_rl_jetpt_weight_wohf, h3_rl_jetpt_weight_wohf, niter_npair);

        h3_rl_jetpt_weight_wohf = (TH3D *)unfold_rl_jetpt_weight_wohf.Hreco();
        
        h3_rl_jetpt_weight_wohf->Divide(h3_eff_rl_jetpt_weight_wohf);

        // Calculate EECs from the npair 3D distributions
        TH2D *h2_eec_truth = new TH2D("h2_eec_truth", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec       = new TH2D("h2_eec"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        TH2D *h2_eec_truth_whf = new TH2D("h2_eec_truth_whf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_whf       = new TH2D("h2_eec_whf"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        TH2D *h2_eec_truth_wohf = new TH2D("h2_eec_truth_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D *h2_eec_wohf       = new TH2D("h2_eec_wohf"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        TH2D *h2_npair_jetpt_weight_truth = new TH2D("h2_npair_jetpt_weight_truth", "", ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH2D *h2_npair_jetpt_weight       = new TH2D("h2_npair_jetpt_weight"      , "", ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        apply_unfolded_weights(h3_rl_jetpt_weight_truth     , h2_eec_truth);
        apply_unfolded_weights(h3_rl_jetpt_weight           , h2_eec);
        apply_unfolded_weights(h3_rl_jetpt_weight_truth_whf , h2_eec_truth_whf);
        apply_unfolded_weights(h3_rl_jetpt_weight_whf       , h2_eec_whf);
        apply_unfolded_weights(h3_rl_jetpt_weight_truth_wohf, h2_eec_truth_wohf);
        apply_unfolded_weights(h3_rl_jetpt_weight_wohf      , h2_eec_wohf);

        // h2_eec_truth      = (TH2D*) h3_rl_jetpt_weight_truth->Project3D("yx");
        // h2_eec            = (TH2D*) h3_rl_jetpt_weight->Project3D("yx");
        // h2_eec_truth_whf  = (TH2D*) h3_rl_jetpt_weight_truth_whf->Project3D("yx");
        // h2_eec_whf        = (TH2D*) h3_rl_jetpt_weight_whf->Project3D("yx");
        // h2_eec_truth_wohf = (TH2D*) h3_rl_jetpt_weight_truth_wohf->Project3D("yx");
        // h2_eec_wohf       = (TH2D*) h3_rl_jetpt_weight_wohf->Project3D("yx");

        h2_npair_jetpt_weight_truth      = (TH2D*) h3_rl_jetpt_weight_truth->Project3D("zy");
        h2_npair_jetpt_weight            = (TH2D*) h3_rl_jetpt_weight->Project3D("zy");

        // Project EECs into 1D histograms
        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hdata_eec[ptbinsize]; 
        
        TH1F* hmc_eec_whf[ptbinsize]; 
        TH1F* hdata_eec_whf[ptbinsize]; 
        
        TH1F* hmc_eec_wohf[ptbinsize]; 
        TH1F* hdata_eec_wohf[ptbinsize]; 
        
        TH1F* hdata_eec_jet3dcorr[ptbinsize]; 
        
        TH1F* hmc_npair_weight[ptbinsize]; 
        TH1F* hdata_npair_weight[ptbinsize]; 
        
        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                hmc_eec[bin]             = new TH1F(Form("hmc_eec%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec[bin]           = new TH1F(Form("hdata_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                hmc_eec_whf[bin]         = new TH1F(Form("hmc_eec_whf%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec_whf[bin]       = new TH1F(Form("hdata_eec_whf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                hmc_eec_wohf[bin]        = new TH1F(Form("hmc_eec_wohf%i",bin)  , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hdata_eec_wohf[bin]      = new TH1F(Form("hdata_eec_wohf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                hdata_eec_jet3dcorr[bin] = new TH1F(Form("hdata_jet3dcorr_eec%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                hmc_npair_weight[bin]    = new TH1F(Form("hmc_npair_weight%i",bin)  , "", nbin_weight, weight_binning);
                hdata_npair_weight[bin]  = new TH1F(Form("hdata_npair_weight%i",bin), "", nbin_weight, weight_binning);
                
                set_histogram_style(hmc_eec[bin]            , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec[bin]          , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec_jet3dcorr[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                set_histogram_style(hmc_eec_whf[bin]        , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec_whf[bin]      , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                set_histogram_style(hmc_eec_wohf[bin]        , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hdata_eec_wohf[bin]      , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h2_eec_truth, hmc_eec[bin]            , bin + 1);
                project_nominal_phase_space(h2_eec      , hdata_eec[bin]          , bin + 1);
                project_nominal_phase_space(h2_eec      , hdata_eec_jet3dcorr[bin], bin + 1);

                project_nominal_phase_space(h2_eec_truth_whf, hmc_eec_whf[bin]  , bin + 1);
                project_nominal_phase_space(h2_eec_whf      , hdata_eec_whf[bin], bin + 1);

                project_nominal_phase_space(h2_eec_truth_wohf, hmc_eec_wohf[bin]  , bin + 1);
                project_nominal_phase_space(h2_eec_wohf      , hdata_eec_wohf[bin], bin + 1);

                // Note weight dimension
                for (int i = 1 ; i <= h2_npair_jetpt_weight->GetNbinsY() ; i++) {
                        hmc_npair_weight[bin]->SetBinContent(i, h2_npair_jetpt_weight_truth->GetBinContent(bin + 1, i));
                        hdata_npair_weight[bin]->SetBinContent(i, h2_npair_jetpt_weight->GetBinContent(bin + 1, i));

                        hmc_npair_weight[bin]->SetBinError(i, h2_npair_jetpt_weight_truth->GetBinError(bin + 1, i));
                        hdata_npair_weight[bin]->SetBinError(i, h2_npair_jetpt_weight->GetBinError(bin + 1, i));
                }
                // Note weight dimension

                hmc_eec[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_eec[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");
                hdata_eec_jet3dcorr[bin]->Scale(1./h1_jet_3dcorr->GetBinContent(bin + 1),"width");

                hmc_eec_whf[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_eec_whf[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");
                
                hmc_eec_wohf[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_eec_wohf[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");
                
                hmc_npair_weight[bin]->Scale(1./h1_jetpt_truth->GetBinContent(bin + 1),"width");
                hdata_npair_weight[bin]->Scale(1./h1_jetpt_data->GetBinContent(bin + 1),"width");
                
                hmc_eec[bin]->Write();
                hdata_eec[bin]->Write();
                
                hmc_eec_whf[bin]->Write();
                hdata_eec_whf[bin]->Write();
                
                hmc_eec_wohf[bin]->Write();
                hdata_eec_wohf[bin]->Write();
                
                hdata_eec_jet3dcorr[bin]->Write();
                
                hmc_npair_weight[bin]->Write();
                hdata_npair_weight[bin]->Write();
                
                hdata_eec_whf[bin]->Add(hdata_eec_wohf[bin]);
                hdata_eec_whf[bin]->Divide(hmc_eec[bin]);
                hdata_eec_whf[bin]->Write(Form("hdata_eec_whfwohf_%i", bin));
        }
        
        if (variation != "nominal") {
                if(gSystem->AccessPathName((output_folder + namef_correctedobservable_data["nominal"]).c_str())) {
                        std::cout<<"Nominal observable not found. Cant calculate syst uncertainty."<<std::endl;

                        return;
                }

                TFile* f_nominal = new TFile((output_folder + namef_correctedobservable_data["nominal"]).c_str());

                TH1F* hdata_eec_nominal[ptbinsize]; 
                TH1F* hdata_eec_nominal_whf[ptbinsize]; 
                TH1F* hdata_eec_nominal_wohf[ptbinsize]; 
                TH1F* hdata_eec_nominal_jet3dcorr[ptbinsize]; 

                TH1F* h_syst_variation[ptbinsize];
                TH1F* h_syst_variation_whf[ptbinsize];
                TH1F* h_syst_variation_wohf[ptbinsize];
                TH1F* h_syst_variation_jet3dcorr[ptbinsize];
                
                for (int bin = 0 ; bin < ptbinsize ; bin++) {
                        hdata_eec_nominal[bin]           = (TH1F*) f_nominal->Get(Form("hdata_eec%i", bin));
                        hdata_eec_nominal_whf[bin]       = (TH1F*) f_nominal->Get(Form("hdata_eec_whf%i", bin));
                        hdata_eec_nominal_wohf[bin]      = (TH1F*) f_nominal->Get(Form("hdata_eec_wohf%i", bin));
                        hdata_eec_nominal_jet3dcorr[bin] = (TH1F*) f_nominal->Get(Form("hdata_jet3dcorr_eec%i", bin));

                        h_syst_variation[bin]           = new TH1F(Form("h_syst_variation%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        h_syst_variation_whf[bin]       = new TH1F(Form("h_syst_variation_whf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        h_syst_variation_wohf[bin]      = new TH1F(Form("h_syst_variation_wohf%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                        h_syst_variation_jet3dcorr[bin] = new TH1F(Form("h_syst_variation_jet3dcorr%i",bin), "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);

                        h_syst_variation[bin]->Add(hdata_eec_nominal[bin], hdata_eec[bin], 1., -1.);
                        h_syst_variation_whf[bin]->Add(hdata_eec_nominal_whf[bin], hdata_eec_whf[bin], 1., -1.);
                        h_syst_variation_wohf[bin]->Add(hdata_eec_nominal_wohf[bin], hdata_eec_wohf[bin], 1., -1.);
                        h_syst_variation_jet3dcorr[bin]->Add(hdata_eec_nominal_jet3dcorr[bin], hdata_eec_jet3dcorr[bin], 1., -1.);
                        
                        h_syst_variation[bin]->Divide(hdata_eec_nominal[bin]);
                        h_syst_variation_whf[bin]->Divide(hdata_eec_nominal_whf[bin]);
                        h_syst_variation_wohf[bin]->Divide(hdata_eec_nominal_wohf[bin]);
                        h_syst_variation_jet3dcorr[bin]->Divide(hdata_eec_nominal_jet3dcorr[bin]);

                        h_syst_variation[bin]->Write();
                        h_syst_variation_whf[bin]->Write();
                        h_syst_variation_wohf[bin]->Write();
                        h_syst_variation_jet3dcorr[bin]->Write();
                }
        }

        file_write->Close();
}