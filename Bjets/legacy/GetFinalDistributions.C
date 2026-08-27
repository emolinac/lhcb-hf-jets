#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "Settings.h"
#include "../Helpers.h"
using namespace std;

void GetFinalDistributions(int NumEvts = -1,
                           int dataset = 91599, bool onlysim9 = true)
{ 
  gStyle->SetImageScaling(3.);
  gStyle->SetPaintTextFormat("3.3f");

  const int sigfigs = 3;
  TString string_final, string_truth, string_eff, string_unfold, string_sys, extension;
  TString str_followHard, str_ghost, str_charged, str_Mag, str_flavor, str_PID(""), str_simversion("");
  TString loc_rootfiles_data("../../root_files/Bjets/");
  TString loc_rootfiles_MC("../../root_files/BjetsMC/");
  TString loc_plots("../../plots/Bjets/results/");

  int JetMeth = (dataset / 1000) % 10;
  int flavor = (dataset / 100) % 10;
  int ptRange = (dataset / 10) % 10;
  int Mag = (dataset / 1) % 10;
  
  int NumIters = 4;

  // int NumIters = 2;
  //int NumIters_z = 1;
  //int NumIters_kt = 1;

  if (Mag == 0)
    str_Mag = "_MD";
  else if (Mag == 1)
    str_Mag = "_MU";
  if (flavor == 1)
    str_flavor = "_udsg";
  else if (flavor == 4)
    str_flavor = "_c";
  else if (flavor == 5)
    str_flavor = "_b";

  if (PID_cut)
    str_PID = "_PID";
  if (flavor == 1)
    str_followHard = "_hard";
  else
  {
    //   if(followHardest) str_followHard = "_hard";
    // else
    str_followHard = "_HF";
  }
  if (ghostCut)
    str_ghost = Form("_ghost_%.1f", ghostProb);

  TString str_pythia = "Pythia 8";
  if (onlysim9)
  {
    str_simversion = "_sim9";
    str_pythia += ".186";
  }


  
  string_truth = TString("truth") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_Mag + str_flavor + str_simversion +  Form("_%d", dataset);

  string_final = TString("corrected_data") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_PID + Form("_iters_%d", NumIters) + Form("_%d", dataset);

  string_sys = TString("total_sys") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost  + str_Mag + str_flavor + str_PID + Form("_iters_%d", NumIters) + Form("_%d", dataset);
  
  extension = TString("publish") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_PID + Form("_iters_%d", NumIters) + str_simversion + Form("_%d", dataset);
  

  cout << string_final << endl;
  cout << string_truth << endl;

  /////////////////////   Get Files /////////////////////////////////

  TFile *file_final = new TFile(loc_rootfiles_data + string_final + TString(".root"), "READ");
  TFile *file_truth = new TFile(loc_rootfiles_MC + string_truth + TString(".root"), "READ");
  TFile *file_sys = new TFile(loc_rootfiles_data + "Systematics/" + string_sys + TString(".root"), "READ");
  //TFile *file_write = new TFile(loc_hists + extension + TString(".root"), "RECREATE");
  //TFile *file_write = TFile::Open(loc_rootfiles_data + extension + TString(".root"), "RECREATE");    

  TH2D *h2_ptz_data = (TH2D *)file_final->Get("ptz_final");
  TH2D *h2_ptz_truth = (TH2D *)file_truth->Get("ptz");
  
  TH2D *h2_ptjt_data = (TH2D *)file_final->Get("ptjt_final");
  TH2D *h2_ptjt_truth = (TH2D *)file_truth->Get("ptjt");
  
  TH2D *h2_ptr_data = (TH2D *)file_final->Get("ptr_final");
  TH2D *h2_ptr_truth = (TH2D *)file_truth->Get("ptr");    

  TH3D *h3_ptzjt_data = (TH3D *)file_final->Get("ptzjt_final");
  TH3D *h3_ptzjt_truth = (TH3D *)file_truth->Get("ptzjt");
  
  TH3D *h3_ptzr_data = (TH3D *)file_final->Get("ptzr_final");
  TH3D *h3_ptzr_truth = (TH3D *)file_truth->Get("ptzr");
  
  TH3D *h3_ptjtr_data = (TH3D *)file_final->Get("ptjtr_final");
  TH3D *h3_ptjtr_truth = (TH3D *)file_truth->Get("ptjtr");    
       
  TH1D *h1_z_data[ptbinsize-1], *h1_z_truth[ptbinsize-1];
  TH1D *h1_jt_data[ptbinsize-1], *h1_jt_truth[ptbinsize-1];
  TH1D *h1_r_data[ptbinsize-1], *h1_r_truth[ptbinsize-1]; 
  TH1D *h1_z_ratio_MCData[ptbinsize-1], *h1_jt_ratio_MCData[ptbinsize-1], *h1_r_ratio_MCData[ptbinsize-1];     
  TH1D *h1_z_sys[ptbinsize-1], *h1_jt_sys[ptbinsize-1], *h1_r_sys[ptbinsize-1];  
  
  TH2D *h2_zjt_data[ptbinsize-1], *h2_zjt_truth[ptbinsize-1];
  TH2D *h2_zr_data[ptbinsize-1], *h2_zr_truth[ptbinsize-1];
  TH2D *h2_jtr_data[ptbinsize-1], *h2_jtr_truth[ptbinsize-1];  
  TH2D *h2_zjt_ratio_MCData[ptbinsize-1], *h2_zr_ratio_MCData[ptbinsize-1], *h2_jtr_ratio_MCData[ptbinsize-1];     
  TH2D *h2_zjt_sys[ptbinsize-1], *h2_zr_sys[ptbinsize-1], *h2_jtr_sys[ptbinsize-1];    
  TH2D *h2_zjt_tot_unc[ptbinsize-1], *h2_zr_tot_unc[ptbinsize-1], *h2_jtr_tot_unc[ptbinsize-1];     


  /////////////////////  Combine last 2 z bins to include z=1.0 entries in the 0.95-1.0 bin /////////////////////////////////

  for (int iy = 1; iy <= h3_ptzjt_truth->GetNbinsY(); iy++)
  {
    for (int iz = 1; iz <= h3_ptzjt_truth->GetNbinsZ(); iz++)
    {
      h3_ptzjt_truth->SetBinContent(zbinsize2D-1, iy, iz, h3_ptzjt_truth->GetBinContent(zbinsize2D-1, iy, iz) + h3_ptzjt_truth->GetBinContent(zbinsize2D, iy, iz));
      h3_ptzjt_truth->SetBinError(zbinsize2D-1, iy, iz, sqrt(pow(h3_ptzjt_truth->GetBinError(zbinsize2D-1, iy, iz),2) + pow(h3_ptzjt_truth->GetBinError(zbinsize2D, iy, iz),2)));
      h3_ptzjt_truth->SetBinContent(zbinsize2D, iy, iz, 0.);
      h3_ptzjt_truth->SetBinError(zbinsize2D, iy, iz, 0.);
    }
  }
  
  for (int iy = 1; iy <= h3_ptzr_truth->GetNbinsY(); iy++)
  {
    for (int iz = 1; iz <= h3_ptzr_truth->GetNbinsZ(); iz++)
    {
      h3_ptzr_truth->SetBinContent(zbinsize2D-1, iy, iz, h3_ptzr_truth->GetBinContent(zbinsize2D-1, iy, iz) + h3_ptzr_truth->GetBinContent(zbinsize2D, iy, iz));
      h3_ptzr_truth->SetBinError(zbinsize2D-1, iy, iz, sqrt(pow(h3_ptzr_truth->GetBinError(zbinsize2D-1, iy, iz),2) + pow(h3_ptzr_truth->GetBinError(zbinsize2D, iy, iz),2)));
      h3_ptzr_truth->SetBinContent(zbinsize2D, iy, iz, 0.);
      h3_ptzr_truth->SetBinError(zbinsize2D, iy, iz, 0.);
    }
  }  
  
  for (int iy = 1; iy <= h2_ptz_truth->GetNbinsY(); iy++)
  {
    std::cout << " Bin content for 0.95 - 1.0 : " << h2_ptz_truth->GetBinContent(zbinsize-1, iy) << " Bin error for 0.95 - 1.0 : " << h2_ptz_truth->GetBinError(zbinsize-1, iy) << std::endl;
    std::cout << " Bin content for 1.0 - 1.1 : " << h2_ptz_truth->GetBinContent(zbinsize, iy) << " Bin error for 1.0 - 1.1 : " << h2_ptz_truth->GetBinError(zbinsize, iy) << std::endl;
    h2_ptz_truth->SetBinContent(zbinsize-1, iy, h2_ptz_truth->GetBinContent(zbinsize-1, iy) + h2_ptz_truth->GetBinContent(zbinsize, iy));
    h2_ptz_truth->SetBinError(zbinsize-1, iy, sqrt(pow(h2_ptz_truth->GetBinError(zbinsize-1, iy),2) + pow(h2_ptz_truth->GetBinError(zbinsize, iy),2)));
    //std::cout << " Bin error for 0.95 - 1.0 : " << h2_ptz_truth->GetBinError(zbinsize-1) << " 1.0 - 1.1 : " << h2_ptz_truth->GetBinError(zbinsize) << std::endl;
    h2_ptz_truth->SetBinContent(zbinsize, iy, 0.);
    h2_ptz_truth->SetBinError(zbinsize, iy, 0.); 
  }

  ////

  //---- paint setup...
  //
  int ican = -1, iframe = -1, itext = -1;
  TCanvas *ccan[1000];
  TH1F *frame[1000];
  TLatex *text[1000];
  for (int i = 0; i < 1000; i++)
  {
    text[i] = new TLatex();
    text[i]->SetNDC(kTRUE);
    text[i]->SetTextSize(0.06);
  }
  TLatex Tl;
  Tl.SetNDC(kTRUE);
  Tl.SetTextSize(0.06);
  Tl.SetTextFont(1);
  Int_t lhcbFont = 132;
  Double_t lhcbTSize = 0.05;
  Double_t lhcbWidth = 2.00;

  Tl.SetTextFont(lhcbFont);
  Tl.SetTextColor(1);
  Tl.SetTextSize(lhcbTSize);
  Tl.SetTextAlign(12);
  //

  // add LHCb label
  auto lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
                                1.01 - gStyle->GetPadTopMargin(),
                                gStyle->GetPadLeftMargin() + 0.30,
                                1.07 - gStyle->GetPadTopMargin(),
                                "BRNDC");
  lhcbName->AddText("LHCb");
  lhcbName->SetFillColor(0);
  lhcbName->SetTextAlign(12);
  lhcbName->SetBorderSize(0);
  // lhcbName->SetTextFont(lhcbFont);

  gStyle->SetOptStat(0);
  // gROOT->SetStyle("Plain");
  // gStyle->SetCanvasDefH(1414); gStyle->SetCanvasDefW(1000);
  // set the paper & margin sizes
  gStyle->SetFillColor(1);
  gStyle->SetFillStyle(1001); // solid
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPaperSize(20, 26);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetTitleStyle(0);
  gStyle->SetStatStyle(0);
  // use large fonts
  gStyle->SetTextFont(lhcbFont);
  gStyle->SetTextSize(lhcbTSize);
  gStyle->SetLabelFont(lhcbFont, "x");
  gStyle->SetLabelFont(lhcbFont, "y");
  gStyle->SetLabelFont(lhcbFont, "z");
  gStyle->SetLabelSize(lhcbTSize, "x");
  gStyle->SetLabelSize(lhcbTSize, "y");
  gStyle->SetLabelSize(lhcbTSize, "z");
  gStyle->SetTitleFont(lhcbFont);
  gStyle->SetTitleFont(lhcbFont, "x");
  gStyle->SetTitleFont(lhcbFont, "y");
  gStyle->SetTitleFont(lhcbFont, "z");
  gStyle->SetTitleSize(1.2 * lhcbTSize, "x");
  gStyle->SetTitleSize(1.2 * lhcbTSize, "y");
  gStyle->SetTitleSize(1.2 * lhcbTSize, "z");

  // use medium bold lines and thick markers
  // gStyle->SetLineWidth(lhcbWidth);
  // gStyle->SetFrameLineWidth(lhcbWidth);
  // gStyle->SetHistLineWidth(lhcbWidth);
  // gStyle->SetFuncWidth(lhcbWidth);
  // gStyle->SetGridWidth(lhcbWidth);
  // gStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes
  // gStyle->SetMarkerStyle(20);
  // gStyle->SetMarkerSize(1.0);

  // label offsets
  gStyle->SetLabelOffset(0.010, "X");
  gStyle->SetLabelOffset(0.010, "Y");

  // by default, do not display histogram decorations:
  // gStyle->SetOptStat(0);
  gStyle->SetOptStat("emr"); // show only nent -e , mean - m , rms -r
  // full opts at http://root.cern.ch/root/html/TStyle.html#TStyle:SetOptStat
  gStyle->SetStatFormat("6.3g");      // specified as c printf options
  gStyle->SetPaintTextFormat("6.3g"); // specified as c printf options

  gStyle->SetOptTitle(0);
  gStyle->SetOptFit(0);
  // gStyle->SetOptFit(1011); // order is probability, Chi2, errors, parameters
  // titles
  gStyle->SetTitleOffset(0.95, "X");
  gStyle->SetTitleOffset(0.95, "Y");
  gStyle->SetTitleOffset(1.2, "Z");
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFont(lhcbFont, "title");
  gStyle->SetTitleX(0.0);
  gStyle->SetTitleY(1.0);
  gStyle->SetTitleW(1.0);
  gStyle->SetTitleH(0.05);

  // look of the statistics box:
  gStyle->SetStatBorderSize(0);
  gStyle->SetStatFont(lhcbFont);
  gStyle->SetStatFontSize(0.05);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.15);

  // put tick marks on top and RHS of plots
  // gStyle->SetPadTickX(1);
  // gStyle->SetPadTickY(1);

  // histogram divisions: only 5 in x to avoid label overlaps
  gStyle->SetNdivisions(505, "x");
  gStyle->SetNdivisions(510, "y");

  // gStyle->SetLineWidth(3);
  gStyle->SetOptStat(0);

  //---- paint...
  char buf[100];
  char bufb[100];
  TString rootfile;
  TString plotfile;
  TString plotfilePDF;
  TString plotfileO;
  TString plotfileC;
  // TString OutputFileBase	= outbase+outinfo;
  rootfile = loc_rootfiles_data + extension + TString(".root");
  plotfile = loc_plots + extension + TString(".ps");
  plotfilePDF = loc_plots + extension + TString(".pdf");
  plotfileO = plotfilePDF + TString("(");
  plotfileC = plotfilePDF + TString("]");
  // c->SaveAs("plots/"+extension+".pdf");

  // Begin plotting
  //
  auto legend_stack_z = new TLegend(0.7, 0.67, 0.9, 0.87);
  legend_stack_z->SetTextSize(0.07);
  
  auto legend_stack_jt = new TLegend(0.7, 0.67, 0.9, 0.87);
  legend_stack_jt->SetTextSize(0.07);
  
  auto legend_stack_r = new TLegend(0.7, 0.67, 0.9, 0.87);
  legend_stack_r->SetTextSize(0.07);  
  
  double left = 0.2;
  double step = 0.1;
  double top = 0.8 - step;
  double top2 = 0.5 - step;
  double top2d = 0.9 - step;
  double step2d = 0.06;

  ofstream textables;
  textables.open("textables.txt");
  // z histograms // 
  textables << "\\begin{longtable}{ccccc}\n";
  textables << "    \\caption{Final results of $\\frac{1}{N_{j}}\\frac{dN}{dz}$ from Figure~\\ref{fig:z_final_pt1}-~\\ref{fig:z_final_pt6}}\\\\\n";
  textables << "    \\label{tab:z_final}\n";
  textables << "    \\ptj (GeV/c) & $z$ & $\\frac{1}{N_{j}}\\frac{dN}{dz}$ & $\\sigma_{stat}$ & $\\sigma_{sys}$\\\\\n";
  textables << "    \\hline\n";
  textables << "    \\hline\n";
  for (int i = 1; i < ptbinsize; i++)
  {
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
    ccan[ican]->SetFillColor(10);


    ccan[ican]->cd();
    ccan[ican]->Divide(1, 1, 0.0001, 0.0001);
    ccan[ican]->cd(1);

    TPad* top_pad = new TPad(Form("top_pad_%d_zpt", i), "Top pad", 0, 0.45, 1, 1);
    top_pad->SetBottomMargin(0);
    top_pad->SetLeftMargin(0.16);
    top_pad->SetRightMargin(0.05);
    top_pad->SetTicks(1, 1);
    top_pad->Draw();

    TPad* bot_pad = new TPad(Form("bot_pad_%d_zpt", i), "Bottom pad", 0, 0, 1, 0.45);
    gStyle->cd();
    gROOT->SetStyle(gStyle->GetName());
    bot_pad->UseCurrentStyle();
    bot_pad->SetTopMargin(0);
    bot_pad->SetBottomMargin(0.25);
    bot_pad->SetLeftMargin(0.16);
    bot_pad->SetRightMargin(0.05);
    bot_pad->SetTicks(1, 1);
    bot_pad->Draw();

    h2_ptz_truth->GetNbinsY();    
    h1_z_data[i-1] = (TH1D *)h2_ptz_data->ProjectionX(Form("z_data_pt%d",i), i+1, i+1);
    h1_z_truth[i-1] = (TH1D *)h2_ptz_truth->ProjectionX(Form("z_truth_pt%d",i), i+1, i+1);
    
    NormalizeHist(h1_z_data[i-1]);
    NormalizeHist(h1_z_truth[i-1]);   
       
    h1_z_sys[i-1] = (TH1D*)file_sys->Get(Form("z_sys_total_pt%d",i));   
    TH1D *h1_z_data_w_sys = (TH1D*)h1_z_data[i-1]->Clone("h1_z_data_w_sys");
    h1_z_sys[i-1]->Multiply(h1_z_sys[i-1], h1_z_data[i-1]);
    if (i!=1)
    {
      textables << "    \\hline\n";
    }
    textables << Form("    $%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]);
    for (int ibin=0; ibin<h1_z_data_w_sys->GetNbinsX(); ++ibin)
    { 
      h1_z_data_w_sys->SetBinError(ibin+1, h1_z_sys[i-1]->GetBinContent(ibin+1));
      if (ibin==0)
        textables << " & " << Form("$%.2f<z<%.2f$", z_binedges[ibin], z_binedges[ibin+1]) << " & " << std::fixed << std::setprecision(sigfigs) << h1_z_data[i-1]->GetBinContent(ibin+1) << " & " << h1_z_data[i-1]->GetBinError(ibin+1) << " & " <<  h1_z_data_w_sys->GetBinError(ibin+1)  << " \\\\\n"; 
      else
        textables << "                     & " << Form("$%.2f<z<%.2f$", z_binedges[ibin], z_binedges[ibin+1]) << " & " << std::fixed << std::setprecision(sigfigs) << h1_z_data[i-1]->GetBinContent(ibin+1) << " & " << h1_z_data[i-1]->GetBinError(ibin+1) << " & " <<  h1_z_data_w_sys->GetBinError(ibin+1)  << " \\\\\n"; 
      //std::cout << "stat error pt: " << i << " : " << h1_z_data[i-1]->GetBinError(ibin+1) << std::endl;
      //std::cout << "sys error pt: " << i << " : " << h1_z_data_w_sys->GetBinError(ibin+1) << std::endl;      
    } 
    
     // Calcualte ratios to truth
     TH1D *h1_z_data_stat_ratio = (TH1D *)h1_z_truth[i-1]->Clone("h1_z_data_stat_ratio");
     h1_z_data_stat_ratio->Divide(h1_z_data[i-1]);
     TH1D *h1_z_data_sys_ratio = (TH1D *)h1_z_truth[i-1]->Clone("h1_z_data_sys_ratio");
     SetHistErrorZero(h1_z_data_sys_ratio);  // No sys. uncertainties on Pythia
     h1_z_data_sys_ratio->Divide(h1_z_data_w_sys);
    
    SetTruthStyle(h1_z_truth[i-1]);
    SetDataStyle(h1_z_data[i-1], 0.7);
    SetDataStyle(h1_z_data_w_sys, 0.7);
    SetDataStyle(h1_z_data_stat_ratio, 0.7);
    SetDataStyle(h1_z_data_sys_ratio, 0.7);
    h1_z_data_w_sys->SetFillColorAlpha(16, 0.5);
    h1_z_data_w_sys->SetLineColorAlpha(0, 0.5);
    //h1_z_truth[i-1]->GetXaxis()->SetRange(1, h1_z_truth[i-1]->GetNbinsX() - 1);

    top_pad->cd();
    h1_z_truth[i-1]->GetYaxis()->SetMaxDigits(2);
    h1_z_truth[i-1]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{dz}");
    h1_z_truth[i-1]->GetYaxis()->SetTitleOffset(0.65);
    h1_z_truth[i-1]->GetYaxis()->SetRangeUser(0, 4.);    
    h1_z_truth[i-1]->GetYaxis()->SetTitleSize(0.1);
    h1_z_truth[i-1]->GetYaxis()->SetTitleFont(lhcbFont);
    h1_z_truth[i-1]->GetYaxis()->SetLabelSize(0.08);
    h1_z_truth[i-1]->GetYaxis()->SetLabelFont(lhcbFont);
    h1_z_truth[i-1]->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_z_truth[i-1]->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_z_truth[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h1_z_truth[i-1]->Draw("HIST SAME");
    h1_z_truth[i-1]->Draw("PE SAME");
    h1_z_data_w_sys->Draw("PE2 SAME");
    h1_z_data[i-1]->Draw("PE SAME");    
 
    lhcbName->Draw("SAME");
    Tl.DrawLatex(left, top + step, "#scale[1.2]{AK5 B^{#pm}-tagged jets}");

    Tl.DrawLatex(left, top, "#scale[1.2]{pp #sqrt{s} = 13 TeV}");
    Tl.DrawLatex(left, top - step, Form("#scale[1.2]{%.1f<p_{T,jet}<%.1f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    Tl.DrawLatex(left, top - 2 * step, "#scale[1.2]{2.5 < #eta_{jet} < 4}");

    if (i == 1)
    {
      legend_stack_z->SetBorderSize(0);
      legend_stack_z->SetFillStyle(0);
      legend_stack_z->SetFillColor(3);
      legend_stack_z->AddEntry(h1_z_data[i-1], "Data");
      legend_stack_z->AddEntry(h1_z_truth[i-1], str_pythia);
      legend_stack_z->Draw("SAME");
    }
    else
    {
      legend_stack_z->Draw("SAME");
    }

    bot_pad->cd();
    h1_z_data_sys_ratio->GetYaxis()->SetTitle("Sim/Data");
    h1_z_data_sys_ratio->GetYaxis()->SetTitleOffset(0.5);
    h1_z_data_sys_ratio->GetYaxis()->SetTitleSize(0.1);
    h1_z_data_sys_ratio->GetYaxis()->SetTitleFont(lhcbFont);
    h1_z_data_sys_ratio->GetYaxis()->SetLabelSize(0.09);
    h1_z_data_sys_ratio->GetYaxis()->SetLabelFont(lhcbFont);
    h1_z_data_sys_ratio->GetXaxis()->SetTitleSize(0.12);
    h1_z_data_sys_ratio->GetXaxis()->SetTitleOffset(0.85);
    h1_z_data_sys_ratio->GetXaxis()->SetTitleFont(lhcbFont);
    h1_z_data_sys_ratio->GetXaxis()->SetLabelSize(0.09);
    h1_z_data_sys_ratio->GetXaxis()->SetLabelFont(lhcbFont);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_z_data_sys_ratio->SetMinimum(0.0);
    h1_z_data_sys_ratio->SetMaximum(2.0);
    h1_z_data_sys_ratio->SetFillColorAlpha(16, 0.5);
    h1_z_data_sys_ratio->SetLineColorAlpha(0, 0.5);
    h1_z_data_sys_ratio->GetXaxis()->SetRangeUser(0.1, 1.0);
    h1_z_data_sys_ratio->Draw("PE2 SAME");
    h1_z_data_stat_ratio->Draw("PE SAME");
    TLine *line = new TLine(0.1, 1.0, 1.0, 1.0);
    line->SetLineStyle(2);
    line->Draw("SAME");


    ccan[ican]->cd();
    ccan[ican]->Update();
    if (ican == 0)
    {
      ccan[ican]->Print(plotfileO.Data());
    }
    else
    {
      ccan[ican]->Print(plotfilePDF.Data());
    }
    ccan[ican]->SaveAs(Form(loc_plots + "z_final_w_sys_pt%d.png",i));
        
  }     
  textables << "\\end{longtable}\n";


  // jT historgrams //
  textables << "\\begin{longtable}{ccccc}\n";
  textables << "    \\caption{Final results of $\\frac{1}{N_{j}}\\frac{dN}{d\\jt}$ from Figure~\\ref{fig:jt_final_pt1}-~\\ref{fig:jt_final_pt6}}\\\\\n";
  textables << "    \\label{tab:jt_final}\n";
  textables << "    \\ptj (GeV/c) & \\jt (GeV/c) & $\\frac{1}{N_{j}}\\frac{dN}{d\\jt}$ & $\\sigma_{stat}$ & $\\sigma_{sys}$\\\\\n";
  textables << "    \\hline\n";
  textables << "    \\hline\n";
  for (int i = 1; i < ptbinsize; i++)
  {    
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
    ccan[ican]->SetFillColor(10);
   
    ccan[ican]->cd();
    ccan[ican]->Divide(1, 1, 0.0001, 0.0001);
    ccan[ican]->cd(1);

    TPad* top_pad = new TPad(Form("top_pad_%d_jtpt", i), "Top pad", 0, 0.45, 1, 1);
    top_pad->SetBottomMargin(0);
    top_pad->SetLeftMargin(0.16);
    top_pad->SetRightMargin(0.05);
    top_pad->SetTicks(1, 1);
    top_pad->Draw();

    TPad* bot_pad = new TPad(Form("bot_pad_%d_jtpt", i), "Bottom pad", 0, 0, 1, 0.45);
    gStyle->cd();
    gROOT->SetStyle(gStyle->GetName());
    bot_pad->UseCurrentStyle();
    bot_pad->SetTopMargin(0);
    bot_pad->SetBottomMargin(0.25);
    bot_pad->SetLeftMargin(0.16);
    bot_pad->SetRightMargin(0.05);
    bot_pad->SetTicks(1, 1);
    bot_pad->Draw();
    
    h1_jt_data[i-1] = (TH1D *)h2_ptjt_data->ProjectionX(Form("jt_data_pt%d",i), i+1, i+1);
    h1_jt_truth[i-1] = (TH1D *)h2_ptjt_truth->ProjectionX(Form("jt_truth_pt%d",i), i+1, i+1);
    
    NormalizeHist(h1_jt_data[i-1]);
    NormalizeHist(h1_jt_truth[i-1]); 
        
    h1_jt_sys[i-1] = (TH1D*)file_sys->Get(Form("jt_sys_total_pt%d",i));   
   /* 
    for (int ibin=0; ibin<h1_jt_sys[i-1]->GetNbinsX(); ++ibin)
    { 
      
      std::cout << "jt relative sys error pt: " << i << " : " << h1_jt_sys[i-1]->GetBinContent(ibin+1) << std::endl;    
      std::cout << "jt data point pt " << i << " : " << h1_jt_data[i-1]->GetBinContent(ibin+1) << std::endl;
      
    }
   */        
    TH1D *h1_jt_data_w_sys = (TH1D*)h1_jt_data[i-1]->Clone("h1_jt_data_w_sys");
    h1_jt_sys[i-1]->Multiply(h1_jt_sys[i-1], h1_jt_data[i-1]);
    if (i!=1)
    {
      textables << "    \\hline\n";
    }
    textables << Form("    $%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]);
    for (int ibin=0; ibin<h1_jt_data_w_sys->GetNbinsX(); ++ibin)
    { 
      h1_jt_data_w_sys->SetBinError(ibin+1, h1_jt_sys[i-1]->GetBinContent(ibin+1));
      if (h1_jt_data[i-1]->GetBinContent(ibin+1) == 0)
      {
        continue;
      }
      if (ibin==0)
        textables << " & " << Form("$%.2f<\\jt<%.2f$", jt_binedges[ibin], jt_binedges[ibin+1]) << " & "  << std::fixed << std::setprecision(sigfigs) << h1_jt_data[i-1]->GetBinContent(ibin+1) << " & " << h1_jt_data[i-1]->GetBinError(ibin+1) << " & " << h1_jt_data_w_sys->GetBinError(ibin+1)  << " \\\\\n"; 
      else
        textables << "                     & " << Form("$%.2f<\\jt<%.2f$", jt_binedges[ibin], jt_binedges[ibin+1]) << " & "  << std::fixed << std::setprecision(sigfigs) << h1_jt_data[i-1]->GetBinContent(ibin+1) << " & " << h1_jt_data[i-1]->GetBinError(ibin+1) << " & " << h1_jt_data_w_sys->GetBinError(ibin+1) << " \\\\\n";       
      //std::cout << "stat error pt: " << i << " : " << h1_jt_data[i-1]->GetBinError(ibin+1) << std::endl;
      //std::cout << "sys error pt: " << i << " : " << h1_jt_data_w_sys->GetBinError(ibin+1) << std::endl;      
      
    } 
    
     // Calcualte ratios to truth
     TH1D *h1_jt_data_stat_ratio = (TH1D *)h1_jt_truth[i-1]->Clone("h1_jt_data_stat_ratio");
     h1_jt_data_stat_ratio->Divide(h1_jt_data[i-1]);
     TH1D *h1_jt_data_sys_ratio = (TH1D *)h1_jt_truth[i-1]->Clone("h1_jt_data_sys_ratio");
     SetHistErrorZero(h1_jt_data_sys_ratio);  // No sys. uncertainties on Pythia
     h1_jt_data_sys_ratio->Divide(h1_jt_data_w_sys);
       

    SetTruthStyle(h1_jt_truth[i-1]);
    SetDataStyle(h1_jt_data[i-1], 0.7);
    SetDataStyle(h1_jt_data_w_sys, 0.7);
    SetDataStyle(h1_jt_data_stat_ratio, 0.7);
    SetDataStyle(h1_jt_data_sys_ratio, 0.7); 
    //h1_jt_data_w_sys->
    h1_jt_data_w_sys->SetFillColorAlpha(16, 0.5);
    h1_jt_data_w_sys->SetLineColorAlpha(0, 0.5);
    //h1_jt_truth[i-1]->GetXaxis()->SetRange(1, h1_jt_truth[i-1]->GetNbinsX() - 1);

    top_pad->cd();
    top_pad->SetLogy();
    h1_jt_truth[i-1]->GetYaxis()->SetMaxDigits(2);
    h1_jt_truth[i-1]->GetXaxis()->SetTitle("j_{T} [GeV/c]");
    h1_jt_truth[i-1]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{dj_{T}}");
    h1_jt_truth[i-1]->GetYaxis()->SetTitleOffset(0.65);
    h1_jt_truth[i-1]->GetYaxis()->SetRangeUser(0, 100.);    
    h1_jt_truth[i-1]->GetYaxis()->SetTitleSize(0.1);
    h1_jt_truth[i-1]->GetYaxis()->SetTitleFont(lhcbFont);
    h1_jt_truth[i-1]->GetYaxis()->SetLabelSize(0.08);
    h1_jt_truth[i-1]->GetYaxis()->SetLabelFont(lhcbFont);
    h1_jt_truth[i-1]->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_jt_truth[i-1]->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_jt_truth[i-1]->SetMinimum(0.00000001);
    h1_jt_truth[i-1]->Draw("HIST SAME");
    h1_jt_truth[i-1]->Draw("PE SAME");
    //h1_ptktdR_data_stat->Draw("HIST SAME");
    h1_jt_data_w_sys->SetMinimum(0.00000001);    
    h1_jt_data[i-1]->SetMinimum(0.00000001);
    h1_jt_data_w_sys->Draw("PE2 SAME");    
    h1_jt_data[i-1]->Draw("PE SAME");

    //gPad->SetLogy();
    // TGaxis *axis_dR = new TGaxis(gPad->GetUxmax(),gPad->GetUymax(),gPad->GetUxmin(),gPad->GetUymax(),jetradius/exp(dR_max), jetradius/exp(dR_min),510,"+L");
    // TGaxis *axis_dR = new TGaxis(gPad->GetUxmax(),gPad->GetUymax(),gPad->GetUxmin(),gPad->GetUymax(),jetradius/exp(dR_max), jetradius/exp(dR_min) - 0.0001,510,"+L");
    // axis_dR->SetLabelOffset(-0.01);
    // Tl.DrawLatex(0.85, 0.95, Form("#DeltaR"));
    // axis_dR->Draw("SAME");
    // axis_dR->SetNdivisions(-510);
      

    // Tl.DrawLatex(left, top + step, "#scale[1.4]{LHCb Unofficial}");
    lhcbName->Draw("SAME");
    Tl.DrawLatex(left, top2 + step, "#scale[1.2]{AK5 B^{#pm}-tagged jets}");

    Tl.DrawLatex(left, top2, "#scale[1.2]{pp #sqrt{s} = 13 TeV}");
    Tl.DrawLatex(left, top2 - step, Form("#scale[1.2]{%.1f<p_{T,jet}<%.1f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    Tl.DrawLatex(left, top2 - 2 * step, "#scale[1.2]{2.5 < #eta_{jet} < 4}");

    if (i == 1)
    {
      legend_stack_jt->SetBorderSize(0);
      legend_stack_jt->SetFillStyle(0);
      legend_stack_jt->SetFillColor(3);
      legend_stack_jt->AddEntry(h1_jt_data[i-1], "Data");
      legend_stack_jt->AddEntry(h1_jt_truth[i-1], str_pythia);
      legend_stack_jt->Draw("SAME");
    }
    else
    {
      legend_stack_jt->Draw("SAME");
    }
 

    bot_pad->cd();
    h1_jt_data_sys_ratio->GetYaxis()->SetTitle("Sim/Data");
    h1_jt_data_sys_ratio->GetYaxis()->SetTitleOffset(0.5);
    h1_jt_data_sys_ratio->GetYaxis()->SetTitleSize(0.1);
    h1_jt_data_sys_ratio->GetYaxis()->SetTitleFont(lhcbFont);
    h1_jt_data_sys_ratio->GetYaxis()->SetLabelSize(0.09);
    h1_jt_data_sys_ratio->GetYaxis()->SetLabelFont(lhcbFont);
    h1_jt_data_sys_ratio->GetXaxis()->SetTitleSize(0.12);
    h1_jt_data_sys_ratio->GetXaxis()->SetTitleOffset(0.85);
    h1_jt_data_sys_ratio->GetXaxis()->SetTitleFont(lhcbFont);
    h1_jt_data_sys_ratio->GetXaxis()->SetLabelSize(0.09);
    h1_jt_data_sys_ratio->GetXaxis()->SetLabelFont(lhcbFont);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_jt_data_sys_ratio->SetMinimum(0.0);
    h1_jt_data_sys_ratio->SetMaximum(2.0);
    h1_jt_data_sys_ratio->SetFillColorAlpha(16, 0.5);
    h1_jt_data_sys_ratio->SetLineColorAlpha(0, 0.5);
    h1_jt_data_sys_ratio->Draw("PE2 SAME");
    h1_jt_data_stat_ratio->Draw("PE SAME");
    TLine *line_jt = new TLine(0.0, 1.0, 10.0, 1.0);
    line_jt->SetLineStyle(2);
    line_jt->Draw("SAME");

    ccan[ican]->cd();
    ccan[ican]->Update();
    if (ican == 0)
    {
      ccan[ican]->Print(plotfileO.Data());
    }
    else
    {
      ccan[ican]->Print(plotfilePDF.Data());
    }      
    ccan[ican]->SaveAs(Form(loc_plots + "jt_final_w_sys_pt%d.png",i));   
        
  }
  textables << "\\end{longtable}\n";

  // r histograms //
  textables << "\\begin{longtable}{ccccc}\n";
  textables << "    \\caption{Final results of $\\frac{1}{N_{j}}\\frac{dN}{dr}$ from Figure~\\ref{fig:r_final_pt1}-~\\ref{fig:r_final_pt6}}\\\\\n";
  textables << "    \\label{tab:r_final}\n";
  textables << "    \\ptj (GeV/c) & $r$ & $\\frac{1}{N_{j}}\\frac{dN}{dr}$ & $\\sigma_{stat}$ & $\\sigma_{sys}$\\\\\n";
  textables << "    \\hline\n";
  textables << "    \\hline\n"; 
  for (int i = 1; i < ptbinsize; i++)
  {     
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
    ccan[ican]->SetFillColor(10);

    ccan[ican]->cd();
    ccan[ican]->Divide(1, 1, 0.0001, 0.0001);
    ccan[ican]->cd(1);
    
    TPad* top_pad = new TPad(Form("top_pad_%d_rpt", i), "Top pad", 0, 0.45, 1, 1);
    top_pad->SetBottomMargin(0);
    top_pad->SetLeftMargin(0.16);
    top_pad->SetRightMargin(0.05);
    top_pad->SetTicks(1, 1);
    top_pad->Draw();

    TPad* bot_pad = new TPad(Form("bot_pad_%d_rpt", i), "Bottom pad", 0, 0, 1, 0.45);
    gStyle->cd();
    gROOT->SetStyle(gStyle->GetName());
    bot_pad->UseCurrentStyle();
    bot_pad->SetTopMargin(0);
    bot_pad->SetBottomMargin(0.25);
    bot_pad->SetLeftMargin(0.16);
    bot_pad->SetRightMargin(0.05);
    bot_pad->SetTicks(1, 1);
    bot_pad->Draw();

    h1_r_data[i-1] = (TH1D *)h2_ptr_data->ProjectionX(Form("r_data_pt%d",i), i+1, i+1);
    h1_r_truth[i-1] = (TH1D *)h2_ptr_truth->ProjectionX(Form("r_truth_pt%d",i), i+1, i+1);
    
    NormalizeHist(h1_r_data[i-1]);
    NormalizeHist(h1_r_truth[i-1]); 
        
    h1_r_sys[i-1] = (TH1D*)file_sys->Get(Form("r_sys_total_pt%d",i));   
    TH1D *h1_r_data_w_sys = (TH1D*)h1_r_data[i-1]->Clone("h1_r_data_w_sys");
    h1_r_sys[i-1]->Multiply(h1_r_sys[i-1], h1_r_data[i-1]);
    if (i!=1)
    {
      textables << "    \\hline\n";
    }
    textables << Form("    $%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]);   
    for (int ibin=0; ibin<h1_r_data_w_sys->GetNbinsX(); ++ibin)
    { 
      h1_r_data_w_sys->SetBinError(ibin+1, h1_r_sys[i-1]->GetBinContent(ibin+1));
      if (h1_r_data[i-1]->GetBinContent(ibin+1) == 0)
      {
        continue;
      }
      if (ibin==0)
        textables << " & " << Form("$%.2f<r<%.2f$", r_binedges[ibin], r_binedges[ibin+1]) << " & "  << std::fixed << std::setprecision(sigfigs) << h1_r_data[i-1]->GetBinContent(ibin+1) << " & " << h1_r_data[i-1]->GetBinError(ibin+1) << " & " << h1_r_data_w_sys->GetBinError(ibin+1)  << " \\\\\n"; 
      else
        textables << "                     & " << Form("$%.2f<r<%.2f$", r_binedges[ibin], r_binedges[ibin+1]) << " & "  << std::fixed << std::setprecision(sigfigs) << h1_r_data[i-1]->GetBinContent(ibin+1) << " & " << h1_r_data[i-1]->GetBinError(ibin+1) << " & " << h1_r_data_w_sys->GetBinError(ibin+1) << " \\\\\n";       
      //std::cout << "stat error pt: " << i << " : " << h1_r_data[i-1]->GetBinError(ibin+1) << std::endl;
      //std::cout << "sys error pt: " << i << " : " << h1_r_data_w_sys->GetBinError(ibin+1) << std::endl;      
      
    } 
    
    // Calcualte ratios to truth
    TH1D *h1_r_data_stat_ratio = (TH1D *)h1_r_truth[i-1]->Clone("h1_r_data_stat_ratio");
    h1_r_data_stat_ratio->Divide(h1_r_data[i-1]);
    TH1D *h1_r_data_sys_ratio = (TH1D *)h1_r_truth[i-1]->Clone("h1_r_data_sys_ratio");
    SetHistErrorZero(h1_r_data_sys_ratio);  // No sys. uncertainties on Pythia
    h1_r_data_sys_ratio->Divide(h1_r_data_w_sys);
        
    SetTruthStyle(h1_r_truth[i-1]);
    SetDataStyle(h1_r_data[i-1], 0.7);
    SetDataStyle(h1_r_data_w_sys, 0.7);
    SetDataStyle(h1_r_data_stat_ratio, 0.7);
    SetDataStyle(h1_r_data_sys_ratio, 0.7);    

    h1_r_data_w_sys->SetFillColorAlpha(16, 0.5);
    h1_r_data_w_sys->SetLineColorAlpha(0, 0.5);

    top_pad->cd();
    top_pad->SetLogy();

    h1_r_truth[i-1]->GetYaxis()->SetMaxDigits(2);
    h1_r_truth[i-1]->GetXaxis()->SetTitle("r");
    h1_r_truth[i-1]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{dr}");
    h1_r_truth[i-1]->GetYaxis()->SetTitleOffset(0.65);
    h1_r_truth[i-1]->GetYaxis()->SetRangeUser(0, 100.);    
    h1_r_truth[i-1]->GetYaxis()->SetTitleSize(0.1);
    h1_r_truth[i-1]->GetYaxis()->SetTitleFont(lhcbFont);
    h1_r_truth[i-1]->GetYaxis()->SetLabelSize(0.08);
    h1_r_truth[i-1]->GetYaxis()->SetLabelFont(lhcbFont);
    h1_r_truth[i-1]->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_r_truth[i-1]->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_r_truth[i-1]->SetMinimum(0.00001);
    h1_r_truth[i-1]->GetXaxis()->SetRangeUser(0.0, 0.4);
    h1_r_truth[i-1]->Draw("HIST SAME");
    h1_r_truth[i-1]->Draw("PE SAME");
    h1_r_data_w_sys->SetMinimum(0.00001);
    h1_r_data[i-1]->SetMinimum(0.00001);
    h1_r_data_w_sys->Draw("PE2 SAME");    
    h1_r_data[i-1]->Draw("PE SAME");

    // Tl.DrawLatex(left, top + step, "#scale[1.4]{LHCb Unofficial}");
    lhcbName->Draw("SAME");
    Tl.DrawLatex(left, top2 + step, "#scale[1.2]{AK5 B^{#pm}-tagged jets}");

    Tl.DrawLatex(left, top2, "#scale[1.2]{pp #sqrt{s} = 13 TeV}");
    Tl.DrawLatex(left, top2 - step, Form("#scale[1.2]{%.1f<p_{T,jet}<%.1f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    Tl.DrawLatex(left, top2 - 2 * step, "#scale[1.2]{2.5 < #eta_{jet} < 4}");

    if (i == 1)
    {
      legend_stack_r->SetBorderSize(0);
      legend_stack_r->SetFillStyle(0);
      legend_stack_r->SetFillColor(3);
      legend_stack_r->AddEntry(h1_r_data[i-1], "Data");
      legend_stack_r->AddEntry(h1_r_truth[i-1], str_pythia);
      legend_stack_r->Draw("SAME");
    }
    else
    {
      legend_stack_r->Draw("SAME");
    }


    bot_pad->cd();
    h1_r_data_sys_ratio->GetYaxis()->SetTitle("Sim/Data");
    h1_r_data_sys_ratio->GetYaxis()->SetTitleOffset(0.5);
    h1_r_data_sys_ratio->GetYaxis()->SetTitleSize(0.1);
    h1_r_data_sys_ratio->GetYaxis()->SetTitleFont(lhcbFont);
    h1_r_data_sys_ratio->GetYaxis()->SetLabelSize(0.09);
    h1_r_data_sys_ratio->GetYaxis()->SetLabelFont(lhcbFont);
    h1_r_data_sys_ratio->GetXaxis()->SetTitleSize(0.12);
    h1_r_data_sys_ratio->GetXaxis()->SetTitleOffset(0.85);
    h1_r_data_sys_ratio->GetXaxis()->SetTitleFont(lhcbFont);
    h1_r_data_sys_ratio->GetXaxis()->SetLabelSize(0.09);
    h1_r_data_sys_ratio->GetXaxis()->SetLabelFont(lhcbFont);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(-1, -1, 0);    
    h1_r_data_sys_ratio->SetMinimum(0.0);
    h1_r_data_sys_ratio->SetMaximum(2.0);
    h1_r_data_sys_ratio->SetFillColorAlpha(16, 0.5);
    h1_r_data_sys_ratio->SetLineColorAlpha(0, 0.5);
    h1_r_data_sys_ratio->GetXaxis()->SetRangeUser(0.0, 0.4);
    h1_r_data_sys_ratio->Draw("PE2 SAME");
    h1_r_data_stat_ratio->Draw("PE SAME");
    TLine *line_r = new TLine(0.0, 1.0, 0.4, 1.0);
    line_r->SetLineStyle(2);
    line_r->Draw("SAME");

    ccan[ican]->cd();
    ccan[ican]->Update();
    if (ican == 0)
    {
      ccan[ican]->Print(plotfileO.Data());
    }
    else
    {
      ccan[ican]->Print(plotfilePDF.Data());
    }
    ccan[ican]->SaveAs(Form(loc_plots + "r_final_w_sys_pt%d.png",i));     
  }  
  textables << "\\end{longtable}\n";

  // zjt histos // 
  for (int i = 1; i < ptbinsize; i++)
  {
    textables << "\\begin{sidewaystable}\n";
    textables << "    \\centering\n";
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
    ccan[ican]->SetFillColor(10);

    ccan[ican]->cd();
    ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
    ccan[ican]->cd(1);
    
    h3_ptzjt_data->GetZaxis()->SetRange(i+1, i+1);
    h3_ptzjt_truth->GetZaxis()->SetRange(i+1, i+1);   
    
    h2_zjt_data[i-1] = (TH2D *)h3_ptzjt_data->Project3D("yx");
    h2_zjt_truth[i-1] = (TH2D *)h3_ptzjt_truth->Project3D("yx");     
    h2_zjt_data[i-1]->SetName(Form("zjt_pt%d",i));     
    h2_zjt_truth[i-1]->SetName(Form("zjt_truth_pt%d",i));     
    
    NormalizeHist(h2_zjt_data[i-1]);
    NormalizeHist(h2_zjt_truth[i-1]); 
        
    h2_zjt_sys[i-1] = (TH2D*)file_sys->Get(Form("zjt_sys_total_pt%d",i));   
    h2_zjt_sys[i-1]->Multiply(h2_zjt_sys[i-1], h2_zjt_data[i-1]);
    
    h2_zjt_tot_unc[i-1] = (TH2D*)h2_zjt_data[i-1]->Clone(Form("zjt_tot_unc_pt%d", i));
    for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
      {
        h2_zjt_tot_unc[i-1]->SetBinContent( xbin+1, ybin+1,  std::sqrt(std::pow(h2_zjt_data[i-1]->GetBinError(xbin+1, ybin+1),2) + std::pow(h2_zjt_sys[i-1]->GetBinContent(xbin+1, ybin+1), 2) ) );
        //std::cout << "stat error pt: " << i << " : " << h1_r_data[i-1]->GetBinError(xbin+1) << std::endl;
        //std::cout << "sys error pt: " << i << " : " << h1_r_data_w_sys->GetBinError(xbin+1) << std::endl; 
      }     
      
    } 

    // ( z , j t )   d a t a   t a b l e //
    textables << "    \\caption{Final results for $\\frac{1}{N_{j}}\\frac{dN}{dzd\\jt}$ from " << Form("Figure~\\ref{fig:zjt_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:zjt_final_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.875\\paperheight}{c|cccccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{14}{c}{\\normalsize{$\\frac{1}{N_{j}}\\frac{dN}{dzd\\jt}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $\\jt \\backslash z$ & "; 
    for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
    { 
      if(ybin!=0)
      {
        if (ybin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[ybin], jt_binedges[ybin+1]) << " & "; 
      }
      for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
      {
        if (xbin==0) {continue;}
        if (ybin==0)
        {
          if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << " & "; 
        }
        else if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
        {
          if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }       
    textables << "    \\end{tabularx}\n";

    // ( z , j t )   s t a t   u n c   t a b l e //
    textables << "    \\caption{Final results for $\\sigma_{stat}$ corresponding to " << Form("Table~\\ref{tab:zjt_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:zjt_final_sigma_stat_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.875\\paperheight}{c|cccccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{14}{c}{\\normalsize{$\\sigma_{stat}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $\\jt \\backslash z$ & "; 
    for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
    { 
      if(ybin!=0)
      {
        if (ybin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[ybin], jt_binedges[ybin+1]) << " & "; 
      }
      for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
      {
        if (xbin==0) {continue;}
        if (ybin==0)
        {
          if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << " & "; 
        }
        else if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
        {
          if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zjt_data[i-1]->GetBinError(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zjt_data[i-1]->GetBinError(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }   
    textables << "    \\end{tabularx}\n";

    // ( z , j t )   s y s   u n c   t a b l e //
    textables << "    \\caption{Final results for $\\sigma_{sys}$ corresponding to " << Form("Table~\\ref{tab:zjt_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:zjt_final_sigma_sys_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.875\\paperheight}{c|cccccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{14}{c}{\\normalsize{$\\sigma_{sys}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $\\jt \\backslash z$ & "; 
    for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
    { 
      if(ybin!=0)
      {
        if (ybin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[ybin], jt_binedges[ybin+1]) << " & "; 
      }
      for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
      {
        if (xbin==0) {continue;}
        if (ybin==0)
        {
          if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << " & "; 
        }
        else if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
        {
          if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zjt_sys[i-1]->GetBinContent(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zjt_sys[i-1]->GetBinContent(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }   
    textables << "    \\end{tabularx}\n";

    gPad->SetLogz();
    h2_zjt_data[i-1]->GetXaxis()->SetTitle("z");
    h2_zjt_data[i-1]->GetYaxis()->SetTitle("j_{T} [GeV/c]");
    h2_zjt_data[i-1]->SetMinimum(0.00001);  
    h2_zjt_data[i-1]->SetMaximum(10.);      
    //h2_zjt_data[i-1]->GetYaxis()->SetLimits(0.0, 20.0);
    h2_zjt_data[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zjt_data[i-1]->Draw("COLZ");

    lhcbName->Draw("SAME");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{AK5 B^{#pm}-tagged jets              pp #sqrt{s} = 13 TeV}");    
    //Tl.DrawLatex(left, top2d, "#scale[0.8]{}");
    Tl.DrawLatex(left, top2d, Form("#scale[0.8]{%.1f<p_{T,jet}<%.1f GeV/c        2.5 < #eta_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top2d - 2 * step2d, "#scale[0.8]{}");    
    
    ccan[ican]->cd(2); 
    h2_zjt_tot_unc[i-1]->SetMinimum(0.00001);  
    h2_zjt_tot_unc[i-1]->SetMaximum(10.);        
    gPad->SetLogz();
    h2_zjt_tot_unc[i-1]->GetXaxis()->SetTitle("z");
    h2_zjt_tot_unc[i-1]->GetYaxis()->SetTitle("j_{T} [GeV/c]");  
    h2_zjt_tot_unc[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);    
    h2_zjt_tot_unc[i-1]->Draw("COLZ");  
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{Total uncertainty (statistical + systematic)}");     
    
    ccan[ican]->cd(3);   
    h2_zjt_truth[i-1]->SetMinimum(0.00001);  
    h2_zjt_truth[i-1]->SetMaximum(10.);                  
    gPad->SetLogz();
    h2_zjt_truth[i-1]->GetXaxis()->SetTitle("z");
    h2_zjt_truth[i-1]->GetYaxis()->SetTitle("j_{T} [GeV/c]");
    h2_zjt_truth[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zjt_truth[i-1]->Draw("COLZ");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{" + str_pythia + "}");
    
    ccan[ican]->cd(4);
    h2_zjt_ratio_MCData[i-1] = (TH2D*)h2_zjt_truth[i-1]->Clone(Form("h2_zjt_ratio_pt%d",i));
    h2_zjt_ratio_MCData[i-1]->Divide(h2_zjt_truth[i-1], h2_zjt_data[i-1]);
    h2_zjt_ratio_MCData[i-1]->SetMaximum(3.);
    //h2_zjt_ratio_MCData[i-1]->SetMinimum(0.0001);        
    h2_zjt_ratio_MCData[i-1]->GetXaxis()->SetTitle("z");
    h2_zjt_ratio_MCData[i-1]->GetYaxis()->SetTitle("j_{T} [GeV/c]");
    h2_zjt_ratio_MCData[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zjt_ratio_MCData[i-1]->Draw("COLZ");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{MC/Data Ratio}");     
    
    ccan[ican]->cd();
    //gPad->SetLogz();
    ccan[ican]->Update();
    if (ican == 0)
    {
      ccan[ican]->Print(plotfileO.Data());
    }
    else
    {
      ccan[ican]->Print(plotfilePDF.Data());
    }
    ccan[ican]->SaveAs(Form(loc_plots + "zjt_final_w_sys_pt%d.png",i));   

    textables << "\\end{sidewaystable}\n";  
  }    
 
  // zr plots //
  for (int i = 1; i < ptbinsize; i++)
  {
    textables << "\\begin{sidewaystable}\n";
    textables << "    \\centering\n";
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
    ccan[ican]->SetFillColor(10);

    ccan[ican]->cd();
    ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
    ccan[ican]->cd(1);
    
    h3_ptzr_data->GetZaxis()->SetRange(i+1, i+1);
    h3_ptzr_truth->GetZaxis()->SetRange(i+1, i+1);   
    
    h2_zr_data[i-1] = (TH2D *)h3_ptzr_data->Project3D("yx");
    h2_zr_truth[i-1] = (TH2D *)h3_ptzr_truth->Project3D("yx");     
    h2_zr_data[i-1]->SetName(Form("zr_pt%d",i));     
    h2_zr_truth[i-1]->SetName(Form("zr_truth_pt%d",i));     
    
    NormalizeHist(h2_zr_data[i-1]);
    NormalizeHist(h2_zr_truth[i-1]); 
        
    h2_zr_sys[i-1] = (TH2D*)file_sys->Get(Form("zr_sys_total_pt%d",i));
    h2_zr_sys[i-1]->Multiply(h2_zr_sys[i-1], h2_zr_data[i-1]);
    
    h2_zr_tot_unc[i-1] = (TH2D*)h2_zr_data[i-1]->Clone(Form("zr_tot_unc_pt%d", i));
    for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
      {
        h2_zr_tot_unc[i-1]->SetBinContent( xbin+1, ybin+1,  std::sqrt(std::pow(h2_zr_data[i-1]->GetBinError(xbin+1, ybin+1),2) + std::pow(h2_zr_sys[i-1]->GetBinContent(xbin+1, ybin+1), 2) ) );
        //std::cout << "stat error pt: " << i << " : " << h1_r_data[i-1]->GetBinError(xbin+1) << std::endl;
        //std::cout << "sys error pt: " << i << " : " << h1_r_data_w_sys->GetBinError(xbin+1) << std::endl; 
      }     
      
    } 

    // ( z , r )   d a t a   t a b l e //
    textables << "    \\caption{Final results for $\\frac{1}{N_{j}}\\frac{dN}{dzdr}$ from " << Form("Figure~\\ref{fig:zr_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:zr_final_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.875\\paperheight}{c|cccccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{14}{c}{\\normalsize{$\\frac{1}{N_{j}}\\frac{dN}{dzdr}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $r \\backslash z$ & "; 
    for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
    { 
      if(ybin!=0)
      {
        if (ybin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << " & "; 
      }
      for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
      {
        if (xbin==0) {continue;}
        if (ybin==0)
        {
          if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << " & "; 
        }
        else if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
        {
          if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }       
    textables << "    \\end{tabularx}\n";

    // ( z , r )   s t a t   u n c   t a b l e //
    textables << "    \\caption{Final results for $\\sigma_{stat}$ corresponding to " << Form("Table~\\ref{tab:zr_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:zr_final_sigma_stat_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.875\\paperheight}{c|cccccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{14}{c}{\\normalsize{$\\sigma_{stat}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $r \\backslash z$ & "; 
    for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
    { 
      if(ybin!=0)
      {
        if (ybin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[ybin], jt_binedges[ybin+1]) << " & "; 
      }
      for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
      {
        if (xbin==0) {continue;}
        if (ybin==0)
        {
          if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << " & "; 
        }
        else if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
        {
          if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zr_data[i-1]->GetBinError(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zr_data[i-1]->GetBinError(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }   
    textables << "    \\end{tabularx}\n";

    // ( z , r )   s y s   u n c   t a b l e //
    textables << "    \\caption{Final results for $\\sigma_{sys}$ corresponding to " << Form("Table~\\ref{tab:zr_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:zr_final_sigma_sys_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.875\\paperheight}{c|cccccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{14}{c}{\\normalsize{$\\sigma_{sys}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $r \\backslash z$ & "; 
    for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
    { 
      if(ybin!=0)
      {
        if (ybin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << " & "; 
      }
      for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
      {
        if (xbin==0) {continue;}
        if (ybin==0)
        {
          if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", z_binedges[xbin], z_binedges[xbin+1]) << " & "; 
        }
        else if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
        {
          if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zr_sys[i-1]->GetBinContent(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_zr_sys[i-1]->GetBinContent(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }   
    textables << "    \\end{tabularx}\n";
        
    gPad->SetLogz();
    h2_zr_data[i-1]->SetMinimum(0.00001);  
    h2_zr_data[i-1]->SetMaximum(200.);      
    h2_zr_data[i-1]->GetXaxis()->SetTitle("z");
    h2_zr_data[i-1]->GetYaxis()->SetTitle("r");
    h2_zr_data[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zr_data[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_zr_data[i-1]->Draw("COLZ");

    lhcbName->Draw("SAME");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{AK5 B^{#pm}-tagged jets              pp #sqrt{s} = 13 TeV}");    
    //Tl.DrawLatex(left, top2d, "#scale[0.8]{}");
    Tl.DrawLatex(left, top2d, Form("#scale[0.8]{%.1f<p_{T,jet}<%.1f GeV/c        2.5 < #eta_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top2d - 2 * step2d, "#scale[0.8]{}");  
    
    ccan[ican]->cd(2);   
    h2_zr_tot_unc[i-1]->SetMinimum(0.00001);  
    h2_zr_tot_unc[i-1]->SetMaximum(200.);               
    gPad->SetLogz();
    h2_zr_tot_unc[i-1]->GetXaxis()->SetTitle("z");
    h2_zr_tot_unc[i-1]->GetYaxis()->SetTitle("r");
    h2_zr_tot_unc[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zr_tot_unc[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_zr_tot_unc[i-1]->Draw("COLZ");  
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{Total uncertainty (statistical + systematic)}");     
    
    ccan[ican]->cd(3);  
    h2_zr_truth[i-1]->SetMinimum(0.00001);  
    h2_zr_truth[i-1]->SetMaximum(200.);                  
    gPad->SetLogz();
    h2_zr_truth[i-1]->GetXaxis()->SetTitle("z");
    h2_zr_truth[i-1]->GetYaxis()->SetTitle("r");
    h2_zr_truth[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zr_truth[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_zr_truth[i-1]->Draw("COLZ");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{" + str_pythia + "}");
    
    ccan[ican]->cd(4);
    h2_zr_ratio_MCData[i-1] = (TH2D*)h2_zr_truth[i-1]->Clone(Form("h2_zr_ratio_pt%d",i));
    h2_zr_ratio_MCData[i-1]->Divide(h2_zr_truth[i-1], h2_zr_data[i-1]);    
    //h2_zr_ratio_MCData[i-1]->SetMinimum(0.0001);
    h2_zr_ratio_MCData[i-1]->SetMaximum(3.); 
    h2_zr_ratio_MCData[i-1]->GetXaxis()->SetTitle("z");
    h2_zr_ratio_MCData[i-1]->GetYaxis()->SetTitle("r");
    h2_zr_ratio_MCData[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zr_ratio_MCData[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_zr_ratio_MCData[i-1]->Draw("COLZ");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{MC/Data Ratio}");     
    
    ccan[ican]->cd();
    //gPad->SetLogz();
    ccan[ican]->Update();
    if (ican == 0)
    {
      ccan[ican]->Print(plotfileO.Data());
    }
    else
    {
      ccan[ican]->Print(plotfilePDF.Data());
    }
    ccan[ican]->SaveAs(Form(loc_plots + "zr_final_w_sys_pt%d.png",i));
    
    textables << "\\end{sidewaystable}\n";  
  }      
  
  // jtr plots //
  for (int i = 1; i < ptbinsize; i++)
  {
    textables << "\\begin{sidewaystable}\n";
    textables << "    \\centering\n";
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
    ccan[ican]->SetFillColor(10);

    ccan[ican]->cd();
    ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
    ccan[ican]->cd(1);
    
    h3_ptjtr_data->GetZaxis()->SetRange(i+1, i+1);
    h3_ptjtr_truth->GetZaxis()->SetRange(i+1, i+1);   
    
    h2_jtr_data[i-1] = (TH2D *)h3_ptjtr_data->Project3D("yx");
    h2_jtr_truth[i-1] = (TH2D *)h3_ptjtr_truth->Project3D("yx");     
    h2_jtr_data[i-1]->SetName(Form("jtr_pt%d",i));     
    h2_jtr_truth[i-1]->SetName(Form("jtr_truth_pt%d",i));     
    
    NormalizeHist(h2_jtr_data[i-1]);
    NormalizeHist(h2_jtr_truth[i-1]); 
        
    h2_jtr_sys[i-1] = (TH2D*)file_sys->Get(Form("jtr_sys_total_pt%d",i));   
    h2_jtr_sys[i-1]->Multiply(h2_jtr_sys[i-1], h2_jtr_data[i-1]);
    
    h2_jtr_tot_unc[i-1] = (TH2D*)h2_jtr_data[i-1]->Clone(Form("jtr_tot_unc_pt%d", i));
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      {
        h2_jtr_tot_unc[i-1]->SetBinContent( xbin+1, ybin+1,  std::sqrt(std::pow(h2_jtr_data[i-1]->GetBinError(xbin+1, ybin+1),2) + std::pow(h2_jtr_sys[i-1]->GetBinContent(xbin+1, ybin+1), 2) ) );
        //std::cout << "stat error pt: " << i << " : " << h1_r_data[i-1]->GetBinError(xbin+1) << std::endl;
        //std::cout << "sys error pt: " << i << " : " << h1_r_data_w_sys->GetBinError(xbin+1) << std::endl; 
      }     
      
    }
    
    // ( j t , r  )   d a t a   t a b l e //
    textables << "    \\caption{Final results for $\\frac{1}{N_{j}}\\frac{dN}{d\\jt dr}$ from " << Form("Figure~\\ref{fig:jtr_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:jtr_final_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.7\\paperheight}{c|ccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{11}{c}{\\normalsize{$\\frac{1}{N_{j}}\\frac{dN}{d\\jt dr}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $\\jt \\backslash r$ & "; 
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      if(xbin!=0)
      {
        if (xbin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[xbin], jt_binedges[xbin+1]) << " & "; 
      }
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      {
        if (xbin==0)
        {
          if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
            textables << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << " & "; 
        }
        else if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
        {
          if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }       
    textables << "    \\end{tabularx}\n";

    // ( j t , r )   s t a t   u n c   t a b l e //
    textables << "    \\caption{Final results for $\\sigma_{stat}$ corresponding to " << Form("Table~\\ref{tab:jtr_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:jtr_final_sigma_stat_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.7\\paperheight}{c|ccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{11}{c}{\\normalsize{$\\sigma_{stat}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $\\jt \\backslash r$ & "; 
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      if(xbin!=0)
      {
        if (xbin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[xbin], jt_binedges[xbin+1]) << " & "; 
      }
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      {
        if (xbin==0)
        {
          if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
            textables << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << " & "; 
        }
        else if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
        {
          if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_jtr_data[i-1]->GetBinError(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_jtr_data[i-1]->GetBinError(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }   
    textables << "    \\end{tabularx}\n";

    // ( j t , r )   s y s   u n c   t a b l e //
    textables << "    \\caption{Final results for $\\sigma_{sys}$ corresponding to " << Form("Table~\\ref{tab:jtr_final_pt%d}}",i)<<"\n";
    textables << "    " << Form("\\label{tab:jtr_final_sigma_sys_pt%d}",i) << "\n";
    textables << "    \\tiny\n";
    textables << "    \\begin{tabularx}{0.7\\paperheight}{c|ccccccccccc}\n";
    textables << "    " << Form("$%.1f<\\ptj<%.1f$", pt_binedges[i], pt_binedges[i+1]) << "& \\multicolumn{11}{c}{\\normalsize{$\\sigma_{sys}$}}\\\\\n";
    textables << "        \\hline\n";
    textables << "        \\hline\n"; 
    textables << "        $\\jt \\backslash r$ & "; 
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      if(xbin!=0)
      {
        if (xbin==1)
        {
          textables << "        \\hline\n";
        }
        textables << "        " << Form("$%.2f-%.2f$", jt_binedges[xbin], jt_binedges[xbin+1]) << " & "; 
      }
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      {
        if (xbin==0)
        {
          if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
            textables << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << "\\\\\n"; 
          else
            textables << Form("$%.2f-%.2f$", r_binedges[ybin], r_binedges[ybin+1]) << " & "; 
        }
        else if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
        {
          if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - \\\\\n";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_jtr_sys[i-1]->GetBinContent(xbin+1, ybin+1) << "\\\\\n";
          }
        }
        else
        {
          if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
          {
            textables << " - &";
          }
          else
          {
            textables << std::defaultfloat << std::setprecision(sigfigs) << h2_jtr_sys[i-1]->GetBinContent(xbin+1, ybin+1) << " & ";
          }
        }
      }           
    }   
    textables << "    \\end{tabularx}\n";
        
    h2_jtr_data[i-1]->SetMinimum(0.00001);  
    h2_jtr_data[i-1]->SetMaximum(50.);        
    gPad->SetLogz();
    h2_jtr_data[i-1]->GetXaxis()->SetTitle("j_{T} [GeV/c]");
    h2_jtr_data[i-1]->GetYaxis()->SetTitle("r");
    h2_jtr_data[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);  
    h2_jtr_data[i-1]->Draw("COLZ");

    lhcbName->Draw("SAME");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{AK5 B^{#pm}-tagged jets              pp #sqrt{s} = 13 TeV}");    
    //Tl.DrawLatex(left, top2d, "#scale[0.8]{}");
    Tl.DrawLatex(left, top2d, Form("#scale[0.8]{%.1f<p_{T,jet}<%.1f GeV/c        2.5 < #eta_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top2d - 2 * step2d, "#scale[0.8]{}");  
    
    ccan[ican]->cd(2); 
    h2_jtr_tot_unc[i-1]->SetMinimum(0.00001);
    h2_jtr_tot_unc[i-1]->SetMaximum(50.);
    gPad->SetLogz();
    h2_jtr_tot_unc[i-1]->GetXaxis()->SetTitle("j_{T} [GeV/c]");
    h2_jtr_tot_unc[i-1]->GetYaxis()->SetTitle("r");
    h2_jtr_tot_unc[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_jtr_tot_unc[i-1]->Draw("COLZ");  
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{Total uncertainty (statistical + systematic)}");     
    
    ccan[ican]->cd(3);
    h2_jtr_truth[i-1]->SetMinimum(0.00001);
    h2_jtr_truth[i-1]->SetMaximum(50.);
    gPad->SetLogz();
    h2_jtr_truth[i-1]->GetXaxis()->SetTitle("j_{T} [GeV/c]");
    h2_jtr_truth[i-1]->GetYaxis()->SetTitle("r");
    h2_jtr_truth[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_jtr_truth[i-1]->Draw("COLZ");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{" + str_pythia + "}");
    
    ccan[ican]->cd(4);
    h2_jtr_ratio_MCData[i-1] = (TH2D*)h2_jtr_truth[i-1]->Clone(Form("h2_jtr_ratio_pt%d",i));
    h2_jtr_ratio_MCData[i-1]->Divide(h2_jtr_truth[i-1], h2_jtr_data[i-1]);    
    //h2_jtr_ratio_MCData[i-1]->SetMinimum(0.0001);
    h2_jtr_ratio_MCData[i-1]->SetMaximum(3.);
    h2_jtr_ratio_MCData[i-1]->GetXaxis()->SetTitle("j_{T} [GeV/c]");
    h2_jtr_ratio_MCData[i-1]->GetYaxis()->SetTitle("r");
    h2_jtr_ratio_MCData[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_jtr_ratio_MCData[i-1]->Draw("COLZ");
    Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{MC/Data Ratio}");
    
    ccan[ican]->cd();
    //gPad->SetLogz();
    ccan[ican]->Update();
    if (ican == 0)
    {
      ccan[ican]->Print(plotfileO.Data());
    }
    else
    {
      ccan[ican]->Print(plotfilePDF.Data());
    }
    ccan[ican]->SaveAs(Form(loc_plots + "jtr_final_w_sys_pt%d.png",i));    

    textables << "\\end{sidewaystable}\n";  
  }    
  
  if (ican > -1)
  {
    cout << " You plotted " << ican + 1 << " canvasses......." << endl;
    ccan[ican]->Print(plotfileC.Data());
  }



  //file_write->Write();
  //file_write->Close();
  //delete file_write;
  //
}


