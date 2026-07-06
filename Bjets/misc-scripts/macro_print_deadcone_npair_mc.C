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

void macro_print_deadcone_npair_mc(const int bin = 0)
{
        TFile* f_heavy      = new TFile((output_folder + "bjets_nominal_eec.root").c_str());
        TFile* f_light      = new TFile((output_folder + "zjets_ct_eec.root").c_str());
        TFile* f_light_data = new TFile((output_folder + "zjets_nominal_eec.root").c_str());
        
        THStack* hs = new THStack();
        TLegend* l = new TLegend();
        
        TH1F* h_light[3];
        TH1F* h_heavy[3];
        TH1F* h_deadcone[3];
        
        TH1F* h_light_data[3];
        TH1F* h_heavy_data[3];
        TH1F* h_deadcone_data[3];
        
        for (int bin = 0 ; bin < 3 ; bin++) {
                h_heavy[bin] = (TH1F*) f_heavy->Get(Form("hmc_npair%i", bin + 4));
                h_light[bin] = (TH1F*) f_light->Get(Form("hcorr_npair_truth%i", bin));

                h_deadcone[bin] = new TH1F(Form("h_deadcone%i", bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                h_heavy_data[bin] = (TH1F*) f_heavy->Get(Form("hdata_npair%i", bin + 4));
                h_light_data[bin] = (TH1F*) f_light_data->Get(Form("hcorr_npair%i", bin));
                h_deadcone_data[bin] = new TH1F(Form("h_deadcone_data%i", bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                for (int i = 1 ; i < h_deadcone[bin]->GetNbinsX(); i++) {
                        if (i == 1)
                                continue;

                        double light_npair_interpolation = h_light[bin]->Interpolate(h_heavy[bin]->GetBinCenter(i));

                        h_deadcone[bin]->SetBinContent(i, h_heavy[bin]->GetBinContent(i) / light_npair_interpolation);

                        double light_npair_interpolation_data = h_light_data[bin]->Interpolate(h_heavy_data[bin]->GetBinCenter(i));

                        h_deadcone_data[bin]->SetBinContent(i, h_heavy_data[bin]->GetBinContent(i) / light_npair_interpolation_data);
                }

                h_deadcone[bin]->SetLineStyle(10);
                
                set_histogram_style(h_deadcone[bin], corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin], std_marker_size+2);

                set_histogram_style(h_deadcone_data[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+2);

                hs->Add(h_deadcone[bin], "APC");
                hs->Add(h_deadcone_data[bin], "APC");
                l->AddEntry(h_deadcone[bin],Form("PYTHIA8 : %.1f<p_{T,jet}<%.1f", pt_binedges[bin + 4], pt_binedges[bin + 5]),"p");
                l->AddEntry(h_deadcone_data[bin],Form("Data : %.1f<p_{T,jet}<%.1f", pt_binedges[bin + 4], pt_binedges[bin + 5]),"p");
        }
        
        hs->Draw("NOSTACK");

        hs->SetTitle(";R_{L};#it{N}_{pair}(b-jet) / #it{N}_{pair}(Z+jet)");

        l->Draw("SAME");

        gPad->SetLogx(1);

        hs->GetXaxis()->SetRangeUser(unfolding_rl_nominal_binning[1], 0.5);
}
