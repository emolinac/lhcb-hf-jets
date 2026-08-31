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

double size_legend_text = 32;

double size_latex_tag = 35;

double xtitle_offset = 1.0;
double ytitle_offset = 1.50;

void set_histo_with_systematics(TH1F* hrelerror, TH1F* hnominal, TH1F* hsystematic, std::string variation);
void draw_three_pads_eec(THStack* h_jetpt_0, THStack* h_jetpt_1, THStack* h_jetpt_2, TLegend* l0, TLegend* l1, TLegend* l2, TLatex* lhcb_print);

void MakeFinalObservable()
{
        if(gSystem->AccessPathName((output_folder + namef_correctedobservable_data["nominal"]).c_str())) {
                std::cout<<"Data file not found for smearing. Check file or variation given as input."<<std::endl;

                return;
        }

        TFile *f_nominal = new TFile((output_folder + namef_correctedobservable_data["nominal"]).c_str(), "READ");

        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hcorr_eec[ptbinsize]; 
        TH1F* hcorr_eec_syst[ptbinsize]; 
        
        double marker_increments = 0.5;

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                if (bin == ptbinsize - 1)
                        marker_increments = 1.5;

                if (bin == 0)
                        marker_increments = 0.8;

                hmc_eec[bin]        = (TH1F*) f_nominal->Get(Form("hmc_eec%i",bin));
                hcorr_eec[bin]      = (TH1F*) f_nominal->Get(Form("hdata_eec%i",bin));
                hcorr_eec_syst[bin] = (TH1F*) hcorr_eec[bin]->Clone(Form("hcorr_eec_syst%i",bin));

                set_histogram_style(hmc_eec[bin]       , corr_marker_color_jet_pt[bin], std_line_width - 1, std_marker_style_jet_pt[bin] , std_marker_size + marker_increments);
                set_histogram_style(hcorr_eec[bin]     , corr_marker_color_jet_pt[bin], std_line_width - 1, corr_marker_style_jet_pt[bin], std_marker_size + marker_increments);
                set_histogram_style(hcorr_eec_syst[bin], corr_marker_color_jet_pt[bin], std_line_width - 1, corr_marker_style_jet_pt[bin], std_marker_size + marker_increments);
                
                hcorr_eec[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
                hcorr_eec_syst[bin]->SetFillColorAlpha(corr_marker_color_jet_pt[bin], 0.3);
        
                marker_increments = 0.5;
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

                delete fsyst[syst_index];
        }

        // Print inclusive eecs in three pads
        gStyle->SetPaintTextFormat("4.2f");
        gStyle->SetLegendFont(133);
        gStyle->SetLegendTextSize(size_legend_text);

        double shift = 0.35;
        
        THStack* hs[ptbinsize];
        TLegend* l_data = new TLegend(0.22 + shift, 0.71, 0.50 + shift, 0.90, "Data:","NDC");
        TLegend* l_mc   = new TLegend(0.05 + shift, 0.70, 0.33 + shift, 0.89, "PYTHIA 8:","NDC");
        TLegend* l[ptbinsize];
        TLegend* l_justdata[ptbinsize];
        TLatex* lhcbprint = new TLatex();
        double l_xcoord[ptbinsize][2] = {{0.56, 0.84},{0.56, 0.84}, {0.465, 0.745}, {0.42, 0.70}, {0.56, 0.84}, {0.465, 0.745}, {0.42, 0.70}};
        double l_ycoord[ptbinsize][2] = {{0.68, 0.91},{0.68, 0.91}, {0.68, 0.91}, {0.68, 0.91}, {0.68, 0.91}, {0.68, 0.91}, {0.68, 0.91}};
        
        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                hs[bin] = new THStack();
                l[bin] = new TLegend(l_xcoord[bin][0], l_ycoord[bin][0], l_xcoord[bin][1], l_ycoord[bin][1]);
                l[bin]->SetHeader(Form("%.0f<#it{p}_{T, jet}<%.0f GeV",pt_binedges[bin],pt_binedges[bin + 1]));

                hs[bin]->Add(hcorr_eec[bin],"E X0");
                hs[bin]->Add(hcorr_eec_syst[bin],"E2");
                hs[bin]->Add(hmc_eec[bin],"E");

                l_data->AddEntry(hcorr_eec[bin],Form("%.0f<#it{p}_{T, jet}<%.0f GeV",pt_binedges[bin],pt_binedges[bin + 1]),"p");
                l_mc->AddEntry(hmc_eec[bin],Form("%.0f<#it{p}_{T, jet}<%.0f GeV",pt_binedges[bin],pt_binedges[bin + 1]),"p");

                l[bin]->AddEntry(hcorr_eec[bin],"Data","p");
                l[bin]->AddEntry(hmc_eec[bin],"PYTHIA 8","p");
        }

        draw_three_pads_eec(hs[4], hs[5], hs[6], l[4], l[5], l[6], lhcbprint);
        // draw_three_pads_eec(hs[0], hs[1], hs[2], l[0], l[1], l[2], lhcbprint);
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

        // h_jetpt_0->SetMaximum(1.45);
        // h_jetpt_1->SetMaximum(1.45);
        // h_jetpt_2->SetMaximum(1.45);

        // Draw tag
        lhcb_print->SetTextAlign(22);
        lhcb_print->SetLineWidth(1);
        lhcb_print->SetTextFont(133);
        lhcb_print->SetTextSize(size_latex_tag);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb, #it{pp} collisions, #sqrt{#it{s}} = 13 TeV, AK5 #it{b}-jets");

        // Draw
        pad1->cd();
        gPad->SetLogx(1);
        h_jetpt_0->Draw("NOSTACK");
        h_jetpt_0->SetTitle(";#it{R_{L}};#Sigma_{EEC}(#it{R_{L}})");
        h_jetpt_0->GetXaxis()->SetRangeUser(rl_nominal_binning[0]*1.01,rl_nominal_binning[nbin_rl_nominal]);
        h_jetpt_0->SetMaximum(0.7);
        
        l0->Draw("SAME");
        
        pad2->cd();
        gPad->SetLogx(1);
        h_jetpt_1->Draw("NOSTACK");
        h_jetpt_1->SetTitle(";#it{R_{L}};#Sigma_{EEC}(#it{R_{L}})");
        h_jetpt_1->GetXaxis()->SetRangeUser(0.02,0.5);
        h_jetpt_1->SetMaximum(0.7);

        l1->Draw("SAME");

        pad3->cd();
        gPad->SetLogx(1);
        h_jetpt_2->Draw("NOSTACK");
        h_jetpt_2->SetTitle(";#it{R_{L}};#Sigma_{EEC}(#it{R_{L}})");
        h_jetpt_2->GetXaxis()->SetRangeUser(0.02,0.5);
        h_jetpt_2->SetMaximum(0.7);

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

        c->Print("./plots-analysis/threepad_eec_data_mc.pdf");
}

void set_histo_with_systematics(TH1F* hrelerror, TH1F* hnominal, TH1F* hsystematic, std::string variation)
{
        for (int hbin = 1 ; hbin <= hrelerror->GetNbinsX() ; hbin++)
        {
                double dev     = hrelerror->GetBinContent(hbin);
                double dev_err = hrelerror->GetBinError(hbin);
                
                double syst_error            = std::abs(dev * hnominal->GetBinContent(hbin));
                double syst_error_percentage = std::abs(dev);
                
                hsystematic->SetBinError(hbin, std::sqrt(syst_error * syst_error / available_systematics_scale[variation] + std::pow(hsystematic->GetBinError(hbin),2)));
        }
}