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

void macro_print_deadcone_eec(const int bin = 0)
{
        TFile* f_heavy = new TFile((output_folder + "bjets_nominal_eec.root").c_str());
        TFile* f_light = new TFile((output_folder + "zjets_nominal_eec.root").c_str());

        THStack* hs = new THStack();
        TLegend* l = new TLegend();
        TH1F* h_light[3];
        TH1F* h_heavy[3];
        TH1F* h_deadcone[3];
        

        for (int bin = 0 ; bin < 3 ; bin++) {
                h_heavy[bin] = (TH1F*) f_heavy->Get(Form("hdata_eec%i", bin + 4));
                h_light[bin] = (TH1F*) f_light->Get(Form("hcorr_eec%i", bin));

                h_deadcone[bin] = new TH1F(Form("h_deadcone%i", bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);

                for (int i = 1 ; i < h_deadcone[bin]->GetNbinsX(); i++) {
                        if (i == 1)
                                continue;

                        double light_eec_interpolation = h_light[bin]->Interpolate(h_heavy[bin]->GetBinCenter(i));

                        std::cout<<h_heavy[bin]->GetBinContent(i) / light_eec_interpolation<<std::endl;

                        h_deadcone[bin]->SetBinContent(i, h_heavy[bin]->GetBinContent(i) / light_eec_interpolation);
                }

                // h_deadcone[bin]->Divide(h_heavy[bin], h_light[bin], 1, 1);
                
                set_histogram_style(h_deadcone[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size);

                hs->Add(h_deadcone[bin], "");
                l->AddEntry(h_deadcone[bin],Form("%.1f<p_{T,jet}<%.1f", pt_binedges[bin + 4], pt_binedges[bin + 5]),"p");
        }

        hs->Draw("NOSTACK");

        hs->SetTitle(";R_{L};b-jet EEC / Z+jet EEC");

        l->Draw("SAME");

        gPad->SetLogx(1);

        hs->GetXaxis()->SetLimits(unfolding_rl_nominal_binning[1], 0.5);
}
