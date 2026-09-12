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
#include "../../include/names.h"
#include "../../include/TBJetsMC.h"
#include "../../include/TBJetsMC.C"
#include "../../include/utils.cpp"
#include "../../include/utils.h"
#include "../../include/utils-visual.cpp"
#include "../../include/utils-visual.h"

double size_xlabel_text = 33;
double size_ylabel_text = 33;

double size_xtitle_text = 35;
double size_ytitle_text = 35;

double size_legend_text = 30;

double size_latex_tag = 35;

double xtitle_offset = 1.25;
double ytitle_offset = 1.50;

using namespace std;

void macro_print_closure_test()
{
        TCanvas* c = new TCanvas("c","",1080,720);
        c->Draw();

        gStyle->SetPaintTextFormat("4.2f");
        gStyle->SetLegendFont(133);
        gStyle->SetLegendTextSize(size_legend_text);
        gStyle->SetEndErrorSize(10);

        TFile *f = new TFile((output_folder + "bjets_closuretest_eec.root").c_str());

        TH1F* h_ct[ptbinsize];
        THStack* hs = new THStack();
        TLegend* l = new TLegend(0.7,0.7,0.95 - gPad->GetRightMargin(),1- gPad->GetTopMargin());

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                if (bin < 4)
                        continue;

                h_ct[bin] = (TH1F*) f->Get(Form("pseudodata_to_truth_eec%i",bin));

                set_histogram_style(h_ct[bin], corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                hs->Add(h_ct[bin], "E1 P X0");
                l->AddEntry(h_ct[bin],Form("%.1f<p_{T,jet}<%.1f", pt_binedges[bin], pt_binedges[bin + 1]),"p");
        }

        hs->Draw("NOSTACK");
        
        hs->GetXaxis()->SetLabelFont(133);
        hs->GetXaxis()->SetTitleFont(133);
        hs->GetXaxis()->SetLabelSize(size_xlabel_text);
        hs->GetXaxis()->SetLabelOffset(0.005);
        hs->GetXaxis()->SetTitleSize(size_xtitle_text);
        hs->GetXaxis()->SetTitleOffset(xtitle_offset);        
        hs->GetYaxis()->SetLabelFont(133);
        hs->GetYaxis()->SetTitleFont(133);
        hs->GetYaxis()->SetLabelSize(size_ylabel_text);
        hs->GetYaxis()->SetLabelOffset(0.0075);
        hs->GetYaxis()->SetTitleSize(size_ytitle_text);
        hs->GetYaxis()->SetTitleOffset(ytitle_offset);
        
        hs->SetTitle(";#it{R_{L}};Pseudodata/Truth");

        hs->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
        hs->SetMaximum(1.5);
        hs->SetMinimum(0.5);

        l->Draw("same");

        gPad->SetLogx(1);

        TLine* line = new TLine(unfolding_rl_nominal_binning[1], 1, unfolding_rl_nominal_binning[nbin_rl_nominal_unfolding-1], 1);

        line->Draw("same");

        c->Print("./plots/closure_test.pdf");
}
