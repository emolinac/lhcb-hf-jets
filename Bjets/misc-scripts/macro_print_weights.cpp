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

void macro_print_weights()
{
        TFile *f = new TFile((output_folder + "ntuple_bjets_mc.root").c_str());

        TTree* t = (TTree*) f->Get("BTree");

        TH1F* h[ptbinsize];
        THStack* hs = new THStack();
        TLegend* l = new TLegend();

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                h[bin] = new TH1F(Form("weight_%i", bin), "", nbin_weight, weight_binning);

                set_histogram_style(h[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+2);

                t->Project(Form("weight_%i", bin),"pair_weight",Form("jet_pt>%f&&jet_pt<%f", pt_binedges[bin], pt_binedges[bin + 1]));

                h[bin]->Scale(1./h[bin]->Integral());

                hs->Add(h[bin], "HIST P C X0");
                l->AddEntry(h[bin],Form("%.1f<p_{T,jet}<%.1f", pt_binedges[bin], pt_binedges[bin + 1]),"p");
        }

        hs->Draw("NOSTACK");
        hs->SetTitle(";weights;norm.");

        // hs->SetMaximum(1.4);
        // hs->SetMinimum(0.6);

        l->Draw("same");

        gPad->SetLogx(1);
}
