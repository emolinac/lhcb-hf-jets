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

void macro_print_tau_light_heavy()
{
        TFile* f_heavy = new TFile((output_folder + "bjets_nominal_eec.root").c_str());
        TFile* f_light = new TFile((output_folder + "zjets_nominal_eec.root").c_str());
        
        TH1F* hcorr_tau_light[nbin_jet_pt]; 
        TH1F* hcorr_tau[nbin_jet_pt]; 

        double marker_increments = 0.5;

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                if (bin == nbin_jet_pt -1 || bin == 0)
                        marker_increments += 1;

                hcorr_tau[bin] = (TH1F*) f_heavy->Get(Form("hcorr_tau%i",bin));
                hcorr_tau_light[bin] = (TH1F*) f_light->Get(Form("hcorr_tau%i",bin));

                set_histogram_style(hcorr_tau[bin]      , corr_marker_color_jet_pt[bin], std_line_width, std_marker_style_jet_pt[bin], std_marker_size + 1.5 + marker_increments);
                set_histogram_style(hcorr_tau_light[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size + 1.5 + marker_increments);
                
                hcorr_tau[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_tau_light[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);

                marker_increments = 1.5;
        }

        TCanvas* c = new TCanvas("c", "", 1920, 1480);
        c->Draw();

        THStack* s_data = new THStack();
        TLegend* l_data = new TLegend(0.02 + gPad->GetLeftMargin(), 1 - 0.21 - gPad->GetTopMargin(),0.32 + gPad->GetLeftMargin(), 1 - 0.03 - gPad->GetTopMargin());
        
        // Print the EECs as a function of tau
        double tau_binning[nbin_jet_pt][nbin_rl_nominal + 1];
        for (int i = 0 ; i < nbin_jet_pt ; i++) {
                double avge_pt2_jet = jet_pt_bjet_data_avge[i];
                
                get_tau_binning_from_eec_binning(tau_binning[i], rl_nominal_binning, avge_pt2_jet);
        }

        for (int bin = 0 ; bin < nbin_jet_pt ; bin++) {
                s_data->Add(hcorr_tau[bin],"E P X0");
                s_data->Add(hcorr_tau_light[bin],"E P X0");
                l_data->AddEntry(hcorr_tau[bin],Form("B-tagged jet %.0f<#it{p}_{T, jet}<%.0f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
                l_data->AddEntry(hcorr_tau_light[bin],Form("Z-tagged jet : %.0f<#it{p}_{T, jet}<%.0f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]),"p");
        }
        
        TH1F* frame = gPad->DrawFrame(tau_binning[0][0], 0.002, tau_binning[2][nbin_rl_nominal], 0.11);
        s_data->Draw("NOSTACK SAME");
        frame->SetTitle(";#it{R_{L}} #LT #it{p}_{T, jet} #GT (GeV);#Sigma_{EEC}(#it{R_{L}})#times ln(#LT #it{p}_{T, jet} #GT/GeV)/#LT #it{p}_{T, jet} #GT (GeV^{#minus1})");
        l_data->Draw("SAME");
        gPad->SetLogx(1);
        gPad->SetLogy(0);
}