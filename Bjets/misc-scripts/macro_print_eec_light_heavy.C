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

// Visual specs
double size_xlabel_text = 33;
double size_ylabel_text = 33;

double size_xtitle_text = 35;
double size_ytitle_text = 35;

double size_legend_text = 32;

double size_latex_tag = 35;

double xtitle_offset = 1.0;
double ytitle_offset = 1.50;


void draw_three_pads_eec(THStack* h_jetpt_0, THStack* h_jetpt_1, THStack* h_jetpt_2, TLegend* l0, TLegend* l1, TLegend* l2, TLatex* lhcb_print);

void macro_print_eec_light_heavy(const int bin = 0)
{
        TFile* f_heavy      = new TFile((output_folder + "bjets_nominal_eec.root").c_str());
        TFile* f_light      = new TFile((output_folder + "zjets_ct_eec.root").c_str());
        TFile* f_light_data = new TFile((output_folder + "zjets_nominal_eec.root").c_str());
        
        THStack* hs[3];
        TLegend* l[3];
        
        TH1F* h_light[3];
        TH1F* h_heavy[3];
        TH1F* h_deadcone[3];

        TLatex* lhcbprint = new TLatex();
        double l_xcoord[nbin_jet_pt][2] = {{0.52, 0.80}, {0.425, 0.705}, {0.38, 0.66}};
        double l_ycoord[nbin_jet_pt][2] = {{0.68, 0.91}, {0.68, 0.91}, {0.68, 0.91}};
        
        gStyle->SetPaintTextFormat("4.2f");
        gStyle->SetLegendFont(133);
        gStyle->SetLegendTextSize(size_legend_text);

        for (int bin = 0 ; bin < 3 ; bin++) {
                hs[bin] = new THStack();
                l[bin]  = new TLegend(l_xcoord[bin][0], l_ycoord[bin][0], l_xcoord[bin][1], l_ycoord[bin][1]);

                l[bin]->SetHeader(Form("%.0f<#it{p}_{T, jet}<%.0f GeV",jet_pt_binning[bin],jet_pt_binning[bin + 1]));

                h_heavy[bin] = (TH1F*) f_heavy->Get(Form("hdata_eec%i", bin + 4));
                h_light[bin] = (TH1F*) f_light_data->Get(Form("hcorr_eec%i", bin));
                
                h_heavy[bin] = (TH1F*) f_heavy->Get(Form("hdata_eec%i", bin + 4));
                h_light[bin] = (TH1F*) f_light_data->Get(Form("hcorr_eec%i", bin));
                
                set_histogram_style(h_heavy[bin], corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(h_light[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);

                hs[bin]->Add(h_heavy[bin], "APE");
                hs[bin]->Add(h_light[bin], "APE");
                l[bin]->AddEntry(h_heavy[bin],"AK5 B-tagged jet","p");
                l[bin]->AddEntry(h_light[bin],"AK5 Z-tagged jet","p");
        }

        draw_three_pads_eec(hs[0], hs[1], hs[2], l[0], l[1], l[2], lhcbprint);
}

void draw_three_pads_eec(THStack* h_jetpt_0, THStack* h_jetpt_1, THStack* h_jetpt_2, TLegend* l0, TLegend* l1, TLegend* l2, TLatex* lhcb_print)
{
        // Canvas size
        TCanvas *c = new TCanvas("c", "Three Pads with Shared Y Axis", 1500, 600);

        // Overall canvas-level margins in normalized coordinates
        double leftMarginCanvas  = 0.07;
        double rightMarginCanvas = 0.03;

        // Equal plot/frame width for each pad
        double frameWidth = (1.0 - leftMarginCanvas - rightMarginCanvas) / 3.0;

        // Pad x boundaries
        double x0 = 0.0;
        double x1 = leftMarginCanvas + frameWidth;
        double x2 = leftMarginCanvas + 2.0 * frameWidth;
        double x3 = 1.0;

        // Create pads with no gaps between them
        TPad *pad1 = new TPad("pad1", "pad1", x0, 0.0, x1, 1.0);
        TPad *pad2 = new TPad("pad2", "pad2", x1, 0.0, x2, 1.0);
        TPad *pad3 = new TPad("pad3", "pad3", x2, 0.0, x3, 1.0);

        pad1->Draw();
        pad2->Draw();
        pad3->Draw();

        // Margins inside each pad.
        // These are chosen so that the actual plotting areas are equal in size.
        pad1->SetLeftMargin(leftMarginCanvas / x1);
        pad1->SetRightMargin(0.0);

        pad2->SetLeftMargin(0.0);
        pad2->SetRightMargin(0.0);

        pad3->SetLeftMargin(0.0);
        pad3->SetRightMargin(rightMarginCanvas / (x3 - x2));

        // Same vertical margins for all pads
        for (auto pad : {pad1, pad2, pad3}) {
                pad->SetTopMargin(0.08);
                pad->SetBottomMargin(0.14);
                pad->SetTicks(1, 1);
                pad->SetFrameBorderMode(0);
        }

        // Use the same y-axis range for all three plots
        h_jetpt_0->SetMinimum(0.);
        h_jetpt_1->SetMinimum(0.);
        h_jetpt_2->SetMinimum(0.);

        h_jetpt_0->SetMaximum(1.45);
        h_jetpt_1->SetMaximum(1.45);
        h_jetpt_2->SetMaximum(1.45);

        // Draw tag
        lhcb_print->SetTextAlign(22);
        lhcb_print->SetLineWidth(1);
        lhcb_print->SetTextFont(133);
        lhcb_print->SetTextSize(size_latex_tag);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb, #it{pp} collisions, #sqrt{#it{s}} = 13 TeV");

        // Draw
        pad1->cd();
        gPad->SetLogx(1);
        h_jetpt_0->Draw("NOSTACK");
        h_jetpt_0->SetTitle(";#it{R_{L}};#Sigma_{EEC}(#it{R_{L}})");
        h_jetpt_0->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        h_jetpt_0->SetMaximum(1.2);
        
        l0->Draw("SAME");
        
        pad2->cd();
        gPad->SetLogx(1);
        h_jetpt_1->Draw("NOSTACK");
        h_jetpt_1->SetTitle(";#it{R_{L}};#Sigma_{EEC}(#it{R_{L}})");
        h_jetpt_1->GetXaxis()->SetRangeUser(0.02,0.5);
        h_jetpt_1->SetMaximum(1.2);

        l1->Draw("SAME");

        pad3->cd();
        gPad->SetLogx(1);
        h_jetpt_2->Draw("NOSTACK");
        h_jetpt_2->SetTitle(";#it{R_{L}};#Sigma_{EEC}(#it{R_{L}})");
        h_jetpt_2->GetXaxis()->SetRangeUser(0.02,0.5);
        h_jetpt_2->SetMaximum(1.2);

        l2->Draw("SAME");

        // Labels and text sizes!
        for (auto h : {h_jetpt_0, h_jetpt_1, h_jetpt_2}) {
                h->GetXaxis()->SetLabelFont(133);
                h->GetXaxis()->SetTitleFont(133);

                h->GetXaxis()->SetLabelSize(size_xlabel_text);  // same x-axis number size in pixels
                h->GetXaxis()->SetTitleSize(size_xtitle_text);

                h->GetXaxis()->SetTitleOffset(xtitle_offset);
        }

        // Y axis only on the left pad
        h_jetpt_0->GetYaxis()->SetLabelFont(133);
        h_jetpt_0->GetYaxis()->SetTitleFont(133);
        h_jetpt_0->GetYaxis()->SetLabelSize(size_ylabel_text);
        h_jetpt_0->GetYaxis()->SetTitleSize(size_ytitle_text);
        h_jetpt_0->GetYaxis()->SetTitleOffset(ytitle_offset);

        // Hide y-axis labels and titles for middle and right plots
        h_jetpt_1->GetYaxis()->SetLabelSize(0.0);
        h_jetpt_1->GetYaxis()->SetTitleSize(0.0);
        // h_jetpt_1->GetYaxis()->SetTickLength(0.0);

        h_jetpt_2->GetYaxis()->SetLabelSize(0.0);
        h_jetpt_2->GetYaxis()->SetTitleSize(0.0);
        // h_jetpt_2->GetYaxis()->SetTickLength(0.0);

        c->Update();

        c->Print("../plots/eec-light-heavy-data.pdf");
}