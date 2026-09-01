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

// Visual specs
double size_xlabel_text = 33;
double size_ylabel_text = 33;

double size_xtitle_text = 35;
double size_ytitle_text = 35;

double size_legend_text = 52;

double size_latex_tag = 35;

double xtitle_offset = 1.0;
double ytitle_offset = 1.50;

void PrintBreakdownRelSyst2RelStat()
{
        if(gSystem->AccessPathName((output_folder + namef_correctedobservable_data["nominal"]).c_str())) {
                std::cout<<"Data file not found for smearing. Check file or variation given as input."<<std::endl;

                return;
        }

        TFile *f_nominal = new TFile((output_folder + namef_correctedobservable_data["nominal"]).c_str(), "READ");

        // Define visual immediately
        TCanvas* c = new TCanvas("c", "", 1920, 1480);
        c->Draw();

        gStyle->SetLegendFont(133);
        gStyle->SetLegendTextSize(size_legend_text);

        THStack* s[ptbinsize];
        TLegend* l[ptbinsize];

        for (int i = 0 ; i < ptbinsize ; i++) {
                l[i] = new TLegend(1 - 0.35 - gPad->GetRightMargin(), 1 - 0.60 - gPad->GetTopMargin(),1 - gPad->GetRightMargin(), 1 - 0.03 - gPad->GetTopMargin());
                s[i] = new THStack(Form("hs%i",i),Form("hs%i",i));
        }

        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hcorr_eec[ptbinsize]; 
        TH1F* hcorr_eec_syst[ptbinsize]; 
        
        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                hmc_eec[bin]        = (TH1F*) f_nominal->Get(Form("hmc_eec%i",bin));
                hcorr_eec[bin]      = (TH1F*) f_nominal->Get(Form("hdata_eec%i",bin));
                hcorr_eec_syst[bin] = (TH1F*) hcorr_eec[bin]->Clone(Form("hcorr_eec_syst%i",bin));

                set_histogram_style(hmc_eec[bin]       , corr_marker_color_jet_pt[bin], std_line_width - 1, std_marker_style_jet_pt[bin] , std_marker_size);
                set_histogram_style(hcorr_eec[bin]     , corr_marker_color_jet_pt[bin], std_line_width - 1, corr_marker_style_jet_pt[bin], std_marker_size);
                set_histogram_style(hcorr_eec_syst[bin], corr_marker_color_jet_pt[bin], std_line_width - 1, corr_marker_style_jet_pt[bin], std_marker_size);
                
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                
                l[bin]->SetHeader(Form("%.0f < p_{T,jet} < %.0f GeV",pt_binedges[bin],pt_binedges[bin + 1]));
                l[bin]->AddEntry(hcorr_eec_syst[bin],"Statistical Error","f");
        }

        // Include the systematics in the whole deal
        const int nsyst = sizeof(available_systematics)/sizeof(available_systematics[0]);
        TFile* fsyst[nsyst];
        TH1F* hsyst_eec[ptbinsize];

        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                if (gSystem->AccessPathName((output_folder + namef_correctedobservable_data[available_systematics[syst_index]]).c_str()))
                        continue;
                
                fsyst[syst_index] = new TFile((output_folder + namef_correctedobservable_data[available_systematics[syst_index]]).c_str());
                
                for (int bin = 0 ; bin < ptbinsize ; bin++) {
                        hsyst_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_syst_variation%i",bin));

                        set_histo_with_systematics(hsyst_eec[bin], hcorr_eec[bin], hcorr_eec_syst[bin], available_systematics[syst_index]);
                }
        }

        // Extract the statistical error from the total error and get the abs. syst error
        // In order to get the systematic total in the background I have to insert that first and then the other systematics
        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                // substract_stat_error(hcorr_eec[bin], hcorr_eec_syst[bin]);
                set_histoa_relerrors_as_histob_content(hcorr_eec[bin],hcorr_eec_syst[bin]);

                s[bin]->Add(hcorr_eec_syst[bin],"HIST");
        }

        const int nsyst_final = sizeof(final_systematics_names)/sizeof(final_systematics_names[0]);
        
        TH1F* h_final_syst_variation[nsyst_final][ptbinsize];

        for (int final_syst_index = 0 ; final_syst_index < nsyst_final ; final_syst_index++) {
                std::string final_syst_name = final_systematics_names[final_syst_index];

                for (int bin = 0 ; bin < ptbinsize ; bin++) {
                        h_final_syst_variation[final_syst_index][bin] = new TH1F(Form("%i%i",final_syst_index,bin), 
                                                                        "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                        for (int syst_index = 0 ; syst_index < nsyst ; syst_index++) {
                                // Consider systematics that correspond to the group given by final_syst_name
                                if (available_systematics_group[available_systematics[syst_index]] == final_syst_name) {
                                        if (gSystem->AccessPathName((output_folder + namef_correctedobservable_data[available_systematics[syst_index]]).c_str()))
                                                continue;
                                        
                                        fsyst[syst_index] = new TFile((output_folder + namef_correctedobservable_data[available_systematics[syst_index]]).c_str());

                                        hsyst_eec[bin] = (TH1F*) fsyst[syst_index]->Get(Form("h_syst_variation%i",bin));

                                        set_histob_with_grouped_relsyst(hsyst_eec[bin], h_final_syst_variation[final_syst_index][bin], available_systematics[syst_index]);
                                }
                        } // Finished grouping
                        
                        set_histogram_style(h_final_syst_variation[final_syst_index][bin], 
                                            corr_marker_color_jet_pt[final_syst_index], std_line_width - 1, 
                                            std_marker_style_jet_pt[final_syst_index] , std_marker_size + 2);

                        s[bin]->Add(h_final_syst_variation[final_syst_index][bin],"HIST PL");
                        l[bin]->AddEntry(h_final_syst_variation[final_syst_index][bin],final_systematics_legends[final_syst_name].c_str(),"p");
                }
        }

        for (int i = 0 ; i < ptbinsize ; i++) {
                s[i]->Draw("NOSTACK");
                s[i]->SetTitle(";R_{L};Relative uncertainty");
                s[i]->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
                s[i]->SetMinimum(0);
                s[i]->SetMaximum(0.48);
                l[i]->Draw("SAME");
                gPad->SetLogx(1);
                gPad->SetLogy(0);
                
                c->Print(Form("./plots/breakdown_relsyst_2_relstat_jetpt%i.pdf",i));
        }
}