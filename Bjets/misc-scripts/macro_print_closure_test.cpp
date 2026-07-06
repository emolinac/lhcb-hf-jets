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

void macro_print_closure_test()
{
        gStyle->SetEndErrorSize(10);

        TFile *f = new TFile((output_folder + "bjets_closuretest_eec.root").c_str());

        TH1F* h_ct[ptbinsize];
        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                if (bin < 4)
                        continue;

                h_ct[bin] = (TH1F*) f->Get(Form("pseudodata_to_truth_eec%i",bin));

                set_histogram_style(h_ct[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+2);

                hs->Add(h_ct[bin], "E1 P X0");
                l->AddEntry(h_ct[bin],Form("%.1f<p_{T,jet}<%.1f", pt_binedges[bin], pt_binedges[bin + 1]),"p");
        }

        hs->Draw("NOSTACK");
        hs->SetTitle(";R_{L};Pseudodata/Truth");

        hs->SetMaximum(1.4);
        hs->SetMinimum(0.6);

        l->Draw("same");

        gPad->SetLogx(1);
}
