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

using namespace std;

double size_xlabel_text = 33;
double size_ylabel_text = 33;

double size_xtitle_text = 35;
double size_ytitle_text = 35;

double size_legend_text = 25;

double size_latex_tag = 35;

double xtitle_offset = 1.25;
double ytitle_offset = 1.50;

void macro_print_corrections(std::string variation = "nominal")
{
        if(gSystem->AccessPathName((output_folder + "bjets_corrections.root").c_str())) {
                std::cout<<"Corrections file not found. Check file or variation given as input."<<std::endl;

                return;
        }

        TFile* f = new TFile((output_folder + "bjets_corrections.root").c_str());

        TH2D* h2_efficiency_HFptjetpt = (TH2D*) f->Get("efficiency_HFptjetpt");
        TH2D* h2_purity_HFptjetpt     = (TH2D*) f->Get("purity_HFptjetpt");

        TH2D* h2_efficiency_rl_jetpt = (TH2D*) f->Get("efficiency_rl_jetpt");
        TH2D* h2_purity_rl_jetpt     = (TH2D*) f->Get("purity_rl_jetpt");

        TH2D* h2_efficiency_rl_weight = (TH2D*) f->Get("efficiency_rl_weight");
        TH2D* h2_purity_rl_weight     = (TH2D*) f->Get("purity_rl_weight");
        
        TH2D* response_rl     = (TH2D*) f->Get("response_rl");
        TH2D* response_weight = (TH2D*) f->Get("response_weight");
        TH2D* response_jetpt  = (TH2D*) f->Get("h2_response_jetpt");
        
        TH2D* response_rl_detail     = (TH2D*) f->Get("response_rl_detail");
        TH2D* response_weight_detail = (TH2D*) f->Get("response_weight_detail");
        TH2D* response_jetpt_detail  = (TH2D*) f->Get("h2_response_jetpt_detail");
        
        // Print all the relevant corrections
        // TCanvas* c = new TCanvas("c","",800,600);
        TCanvas* c = new TCanvas("c","",1080,720);
        c->Draw();

        // h2_efficiency_HFptjetpt
        // ptbinsize, ptHFbinsize

        gStyle->SetPaintTextFormat("4.2f");
        gStyle->SetLegendFont(133);
        gStyle->SetLegendTextSize(size_legend_text);
        
        THStack* hs = new THStack();
        TLegend* l = new TLegend(0.60,gPad->GetBottomMargin() + 0.03,1 - gPad->GetRightMargin(),gPad->GetBottomMargin()+0.45);
        
        TH1D* h_jet_purity_2d[ptbinsize];
        TH1D* h_jet_efficiency_2d[ptbinsize];
        
        // Jet efficiency
        for (int i = 0 ; i < ptbinsize ; i++) {
                h_jet_efficiency_2d[i] = (TH1D*) h2_efficiency_HFptjetpt->ProjectionX(Form("h_jet_efficiency_2d%i", i), i, i);

                set_histogram_style(h_jet_efficiency_2d[i], corr_marker_color_jet_pt[i], std_line_width-1, corr_marker_style_jet_pt[i], std_marker_size + 0.5);

                hs->Add(h_jet_efficiency_2d[i], "APE");
                l->AddEntry(h_jet_efficiency_2d[i], Form("%.0f<#it{p}_{T, jet}<%.0f GeV",pt_binedges[i],pt_binedges[i + 1]), "P");
        }
        
        hs->SetMaximum(0.34);
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{p}_{T, HF}[GeV];#it{b}-jet reconstruction efficiency");
        l->Draw("SAME");

        c->Print("./plots/bjet-reconstruction-efficiency.pdf");

        // Jet purity
        hs = new THStack();
        l = new TLegend(0.60,gPad->GetBottomMargin() + 0.03,1 - gPad->GetRightMargin(),gPad->GetBottomMargin()+0.45);
        
        for (int i = 0 ; i < ptbinsize ; i++) {
                h_jet_purity_2d[i]     = (TH1D*) h2_purity_HFptjetpt->ProjectionX(Form("h_jet_purity_2d%i", i), i, i);
                
                set_histogram_style(h_jet_purity_2d[i], corr_marker_color_jet_pt[i], std_line_width-1, corr_marker_style_jet_pt[i], std_marker_size + 0.5);

                hs->Add(h_jet_purity_2d[i], "APE");
                l->AddEntry(h_jet_purity_2d[i]    , Form("%.0f<#it{p}_{T, jet}<%.0f GeV",pt_binedges[i],pt_binedges[i + 1]), "P");
        }
        
        hs->SetMinimum(0.1);
        hs->SetMaximum(1.1);
        hs->Draw("NOSTACK");
        hs->SetTitle(";#it{p}_{T, HF}[GeV];#it{b}-jet reconstruction purity");
        l->Draw("SAME");

        c->Print("./plots/bjet-reconstruction-purity.pdf");

        // Pair purity rl jet pt
        TLatex latex;
        latex.SetTextAlign(22); // center alignment
        latex.SetTextSize(0.025);
        latex.SetTextColor(kBlack);
        
        h2_purity_rl_jetpt->Draw("col");
        
        for (int i = 2; i < h2_purity_rl_jetpt->GetNbinsX(); ++i) {
                for (int j = 5; j <= h2_purity_rl_jetpt->GetNbinsY(); ++j) {
                        double x = h2_purity_rl_jetpt->GetXaxis()->GetBinCenter(i);
                        double y = h2_purity_rl_jetpt->GetYaxis()->GetBinCenter(j);
                        double content = h2_purity_rl_jetpt->GetBinContent(i, j);
                        double error = h2_purity_rl_jetpt->GetBinError(i, j);

                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }
        
        h2_purity_rl_jetpt->SetTitle("Purity Correction;#it{R_{L}};#it{p}_{T,jet}[GeV]");
        h2_purity_rl_jetpt->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
        h2_purity_rl_jetpt->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/npairs-reconstruction-purity.pdf");

        // Pair efficiency rl jet pt
        h2_efficiency_rl_jetpt->Draw("col");
        
        for (int i = 2; i < h2_efficiency_rl_jetpt->GetNbinsX(); ++i) {
                for (int j = 5; j <= h2_efficiency_rl_jetpt->GetNbinsY(); ++j) {
                        double x = h2_efficiency_rl_jetpt->GetXaxis()->GetBinCenter(i);
                        double y = h2_efficiency_rl_jetpt->GetYaxis()->GetBinCenter(j);
                        double content = h2_efficiency_rl_jetpt->GetBinContent(i, j);
                        double error = h2_efficiency_rl_jetpt->GetBinError(i, j);

                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }
        
        h2_efficiency_rl_jetpt->SetTitle("Purity Correction;#it{R_{L}};#it{p}_{T,jet}[GeV]");
        h2_efficiency_rl_jetpt->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
        h2_efficiency_rl_jetpt->GetYaxis()->SetRangeUser(jet_pt_binning[0], jet_pt_binning[3]);
        
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/npairs-reconstruction-efficiency.pdf");
        
        // Pair purity rl weight
        latex.SetTextSize(0.02);
        
        h2_purity_rl_weight->Draw("col");
        
        for (int i = 2; i < h2_purity_rl_weight->GetNbinsX(); ++i) {
                for (int j = 1; j <= h2_purity_rl_weight->GetNbinsY(); ++j) {
                        double x = h2_purity_rl_weight->GetXaxis()->GetBinCenter(i);
                        double y = h2_purity_rl_weight->GetYaxis()->GetBinCenter(j);
                        double content = h2_purity_rl_weight->GetBinContent(i, j);
                        double error = h2_purity_rl_weight->GetBinError(i, j);

                        if (content == 0)
                                continue;

                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }
        
        h2_purity_rl_weight->SetTitle("Purity Correction;#it{R_{L}}; Momentum weights");
        h2_purity_rl_weight->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
        h2_purity_rl_weight->GetYaxis()->SetRangeUser(10E-6, 0.3);
        
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/npairs-reconstruction-purity-rl-weight.pdf");

        // Pair efficiency rl weight
        h2_efficiency_rl_weight->Draw("col");
        
        for (int i = 2; i < h2_efficiency_rl_weight->GetNbinsX(); ++i) {
                for (int j = 1; j <= h2_efficiency_rl_weight->GetNbinsY(); ++j) {
                        double x = h2_efficiency_rl_weight->GetXaxis()->GetBinCenter(i);
                        double y = h2_efficiency_rl_weight->GetYaxis()->GetBinCenter(j);
                        double content = h2_efficiency_rl_weight->GetBinContent(i, j);
                        double error = h2_efficiency_rl_weight->GetBinError(i, j);

                        if (content == 0)
                                continue;

                        latex.DrawLatex(x, y, Form("%.2f #pm %.2f", content, error));
                }
        }
        
        h2_efficiency_rl_weight->SetTitle("Efficiency Correction;#it{R_{L}}; Momentum weights");
        h2_efficiency_rl_weight->GetXaxis()->SetRangeUser(rl_nominal_binning[0],rl_nominal_binning[nbin_rl_nominal]);
        h2_efficiency_rl_weight->GetYaxis()->SetRangeUser(10E-6, 0.3);
        
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        
        c->Print("./plots/npairs-reconstruction-efficiency-rl-weight.pdf");

        // RM rl
        TLatex* lhcb_print = new TLatex();
        lhcb_print->SetTextAlign(22);
        lhcb_print->SetLineWidth(1);
        lhcb_print->SetTextFont(133);
        lhcb_print->SetTextSize(size_latex_tag);

        response_rl->GetXaxis()->SetLabelFont(133);
        response_rl->GetXaxis()->SetTitleFont(133);
        response_rl->GetXaxis()->SetLabelSize(size_xlabel_text);
        response_rl->GetXaxis()->SetLabelOffset(0.005);
        response_rl->GetXaxis()->SetTitleSize(size_xtitle_text);
        response_rl->GetXaxis()->SetTitleOffset(xtitle_offset);        
        response_rl->GetYaxis()->SetLabelFont(133);
        response_rl->GetYaxis()->SetTitleFont(133);
        response_rl->GetYaxis()->SetLabelSize(size_ylabel_text);
        response_rl->GetYaxis()->SetLabelOffset(0.0075);
        response_rl->GetYaxis()->SetTitleSize(size_ytitle_text);
        response_rl->GetYaxis()->SetTitleOffset(ytitle_offset);
        response_rl->GetZaxis()->SetLabelFont(133);
        response_rl->GetZaxis()->SetTitleFont(133);
        response_rl->GetZaxis()->SetLabelSize(size_ylabel_text);
        response_rl->GetZaxis()->SetTitleSize(size_ytitle_text);

        response_rl->Draw("colz");
        
        response_rl->Draw("colz");
        response_rl->GetXaxis()->SetRangeUser(0.01,1);
        response_rl->GetYaxis()->SetRangeUser(0.01,1); 
        response_rl->SetTitle("Response matrix of #it{R_{L}};#it{R}_{#it{L}}^{rec};#it{R}_{#it{L}}^{gen}");

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.08);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        gPad->SetLogz(1);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb simulation");
        
        c->Print("./plots/responsematrix_rl.pdf");

        // RM weight
        response_weight->GetXaxis()->SetLabelFont(133);
        response_weight->GetXaxis()->SetTitleFont(133);
        response_weight->GetXaxis()->SetLabelSize(size_xlabel_text);
        response_weight->GetXaxis()->SetLabelOffset(0.005);
        response_weight->GetXaxis()->SetTitleSize(size_xtitle_text);
        response_weight->GetXaxis()->SetTitleOffset(xtitle_offset);        
        response_weight->GetYaxis()->SetLabelFont(133);
        response_weight->GetYaxis()->SetTitleFont(133);
        response_weight->GetYaxis()->SetLabelSize(size_ylabel_text);
        response_weight->GetYaxis()->SetLabelOffset(0.0075);
        response_weight->GetYaxis()->SetTitleSize(size_ytitle_text);
        response_weight->GetYaxis()->SetTitleOffset(ytitle_offset);
        response_weight->GetZaxis()->SetLabelFont(133);
        response_weight->GetZaxis()->SetTitleFont(133);
        response_weight->GetZaxis()->SetLabelSize(size_ylabel_text);
        response_weight->GetZaxis()->SetTitleSize(size_ytitle_text);

        response_weight->Draw("colz");
        response_weight->GetXaxis()->SetRangeUser(10E-6,0.4);
        response_weight->GetYaxis()->SetRangeUser(10E-6,0.4); 
        response_weight->SetTitle("Response matrix of #it{R_{L}};#it{w}^{rec};#it{w}^{gen}");

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.08);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        gPad->SetLogz(1);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb simulation");
        
        c->Print("./plots/responsematrix_weight.pdf");

        // RM weight
        response_jetpt->GetXaxis()->SetLabelFont(133);
        response_jetpt->GetXaxis()->SetTitleFont(133);
        response_jetpt->GetXaxis()->SetLabelSize(size_xlabel_text);
        response_jetpt->GetXaxis()->SetLabelOffset(0.005);
        response_jetpt->GetXaxis()->SetTitleSize(size_xtitle_text);
        response_jetpt->GetXaxis()->SetTitleOffset(xtitle_offset);        
        response_jetpt->GetYaxis()->SetLabelFont(133);
        response_jetpt->GetYaxis()->SetTitleFont(133);
        response_jetpt->GetYaxis()->SetLabelSize(size_ylabel_text);
        response_jetpt->GetYaxis()->SetLabelOffset(0.0075);
        response_jetpt->GetYaxis()->SetTitleSize(size_ytitle_text);
        response_jetpt->GetYaxis()->SetTitleOffset(ytitle_offset);
        response_jetpt->GetZaxis()->SetLabelFont(133);
        response_jetpt->GetZaxis()->SetTitleFont(133);
        response_jetpt->GetZaxis()->SetLabelSize(size_ylabel_text);
        response_jetpt->GetZaxis()->SetTitleSize(size_ytitle_text);

        response_jetpt->Draw("colz");
        response_jetpt->GetXaxis()->SetRangeUser(10E-6,0.4);
        response_jetpt->GetYaxis()->SetRangeUser(10E-6,0.4); 
        response_jetpt->SetTitle("Response matrix of #it{R_{L}};#it{p}_{T, jet}^{rec}[GeV];#it{p}_{T, jet}^{gen}[GeV]");

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.08);
        gPad->SetLogx(1);
        gPad->SetLogy(1);
        gPad->SetLogz(1);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb simulation");
        
        c->Print("./plots/responsematrix_jetpt.pdf");

        // RM rl
        response_rl_detail->GetXaxis()->SetLabelFont(133);
        response_rl_detail->GetXaxis()->SetTitleFont(133);
        response_rl_detail->GetXaxis()->SetLabelSize(size_xlabel_text);
        response_rl_detail->GetXaxis()->SetLabelOffset(0.005);
        response_rl_detail->GetXaxis()->SetTitleSize(size_xtitle_text);
        response_rl_detail->GetXaxis()->SetTitleOffset(xtitle_offset);        
        response_rl_detail->GetYaxis()->SetLabelFont(133);
        response_rl_detail->GetYaxis()->SetTitleFont(133);
        response_rl_detail->GetYaxis()->SetLabelSize(size_ylabel_text);
        response_rl_detail->GetYaxis()->SetLabelOffset(0.0075);
        response_rl_detail->GetYaxis()->SetTitleSize(size_ytitle_text);
        response_rl_detail->GetYaxis()->SetTitleOffset(ytitle_offset);
        response_rl_detail->GetZaxis()->SetLabelFont(133);
        response_rl_detail->GetZaxis()->SetTitleFont(133);
        response_rl_detail->GetZaxis()->SetLabelSize(size_ylabel_text);
        response_rl_detail->GetZaxis()->SetTitleSize(size_ytitle_text);

        response_rl_detail->Draw("colz");
        response_rl_detail->GetXaxis()->SetRangeUser(0.01,1);
        response_rl_detail->GetYaxis()->SetRangeUser(0.01,1); 
        response_rl_detail->SetTitle("Response matrix of #it{R_{L}};#it{R}_{#it{L}}^{rec};#it{R}_{#it{L}}^{gen}");

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.08);
        gPad->SetLogx(0);
        gPad->SetLogy(0);
        gPad->SetLogz(1);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb simulation");
        
        c->Print("./plots/responsematrix_rl_detail.pdf");

        // RM weight
        response_weight_detail->GetXaxis()->SetLabelFont(133);
        response_weight_detail->GetXaxis()->SetTitleFont(133);
        response_weight_detail->GetXaxis()->SetLabelSize(size_xlabel_text);
        response_weight_detail->GetXaxis()->SetLabelOffset(0.005);
        response_weight_detail->GetXaxis()->SetTitleSize(size_xtitle_text);
        response_weight_detail->GetXaxis()->SetTitleOffset(xtitle_offset);        
        response_weight_detail->GetYaxis()->SetLabelFont(133);
        response_weight_detail->GetYaxis()->SetTitleFont(133);
        response_weight_detail->GetYaxis()->SetLabelSize(size_ylabel_text);
        response_weight_detail->GetYaxis()->SetLabelOffset(0.0075);
        response_weight_detail->GetYaxis()->SetTitleSize(size_ytitle_text);
        response_weight_detail->GetYaxis()->SetTitleOffset(ytitle_offset);
        response_weight_detail->GetZaxis()->SetLabelFont(133);
        response_weight_detail->GetZaxis()->SetTitleFont(133);
        response_weight_detail->GetZaxis()->SetLabelSize(size_ylabel_text);
        response_weight_detail->GetZaxis()->SetTitleSize(size_ytitle_text);

        response_weight_detail->Draw("colz");
        response_weight_detail->GetXaxis()->SetRangeUser(10E-6,0.4);
        response_weight_detail->GetYaxis()->SetRangeUser(10E-6,0.4); 
        response_weight_detail->SetTitle("Response matrix of #it{R_{L}};#it{w}^{rec};#it{w}^{gen}");

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.08);
        gPad->SetLogx(0);
        gPad->SetLogy(0);
        gPad->SetLogz(1);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb simulation");
        
        c->Print("./plots/responsematrix_weight_detail.pdf");

        // RM weight
        response_jetpt_detail->GetXaxis()->SetLabelFont(133);
        response_jetpt_detail->GetXaxis()->SetTitleFont(133);
        response_jetpt_detail->GetXaxis()->SetLabelSize(size_xlabel_text);
        response_jetpt_detail->GetXaxis()->SetLabelOffset(0.005);
        response_jetpt_detail->GetXaxis()->SetTitleSize(size_xtitle_text);
        response_jetpt_detail->GetXaxis()->SetTitleOffset(xtitle_offset);        
        response_jetpt_detail->GetYaxis()->SetLabelFont(133);
        response_jetpt_detail->GetYaxis()->SetTitleFont(133);
        response_jetpt_detail->GetYaxis()->SetLabelSize(size_ylabel_text);
        response_jetpt_detail->GetYaxis()->SetLabelOffset(0.0075);
        response_jetpt_detail->GetYaxis()->SetTitleSize(size_ytitle_text);
        response_jetpt_detail->GetYaxis()->SetTitleOffset(ytitle_offset);
        response_jetpt_detail->GetZaxis()->SetLabelFont(133);
        response_jetpt_detail->GetZaxis()->SetTitleFont(133);
        response_jetpt_detail->GetZaxis()->SetLabelSize(size_ylabel_text);
        response_jetpt_detail->GetZaxis()->SetTitleSize(size_ytitle_text);

        response_jetpt_detail->Draw("colz");
        response_jetpt_detail->GetXaxis()->SetRangeUser(10E-6,0.4);
        response_jetpt_detail->GetYaxis()->SetRangeUser(10E-6,0.4); 
        response_jetpt_detail->SetTitle("Response matrix of #it{R_{L}};#it{p}_{T, jet}^{rec}[GeV];#it{p}_{T, jet}^{gen}[GeV]");

        gPad->SetRightMargin(0.15);
        gPad->SetTopMargin(0.08);
        gPad->SetLogx(0);
        gPad->SetLogy(0);
        gPad->SetLogz(1);
        
        lhcb_print->DrawLatexNDC(0.52, 0.96, "LHCb simulation");
        
        c->Print("./plots/responsematrix_jetpt_detail.pdf");
}