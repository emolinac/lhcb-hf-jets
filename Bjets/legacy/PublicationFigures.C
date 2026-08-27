#include <TCanvas.h>
#include <TImage.h>
#include <TPad.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "Settings.h"
#include "../Helpers.h"
using namespace std;

void PublicationFigures(int NumEvts = -1,
                        int dataset = 91599,
                        bool onlysim9 = true)
{ 
  gROOT->ForceStyle(1);
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

  string_truth = TString("truth") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_Mag + str_flavor + str_simversion + Form("_%d", dataset);

  string_final = TString("corrected_data") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_PID + Form("_iters_%d", NumIters) + Form("_%d", dataset);

  string_sys = TString("total_sys") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost  + str_Mag + str_flavor + str_PID + Form("_iters_%d", NumIters) + Form("_%d", dataset);  
  
  extension = TString("publish_use") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_PID + Form("_iters_%d", NumIters) + str_simversion + Form("_%d", dataset);
  

  cout << string_final << endl;
  cout << string_truth << endl;
  
  /////////////////////   Get Files /////////////////////////////////

  std::cout << "made it to get files" << std::endl;

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

  /////////////////////  Combine last 2 z bins to include z=1.0 entries in the previous bin /////////////////////////////////

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
    h2_ptz_truth->SetBinContent(zbinsize-1, iy, h2_ptz_truth->GetBinContent(zbinsize-1, iy) + h2_ptz_truth->GetBinContent(zbinsize, iy));
    h2_ptz_truth->SetBinError(zbinsize-1, iy, sqrt(pow(h2_ptz_truth->GetBinError(zbinsize-1, iy),2) + pow(h2_ptz_truth->GetBinError(zbinsize, iy),2)));
    h2_ptz_truth->SetBinContent(zbinsize, iy, 0.);
    h2_ptz_truth->SetBinError(zbinsize, iy, 0.); 
  }

  // h2_ptktdR_truth->GetYaxis()->SetRange(1, ktbinsize);

  // NormalizeHist(h2_ptktdR_truth);
  // NormalizeHist(h2_ptzdR_truth);
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
  auto legend_stack_z = new TLegend(0.6, 0.65, 0.8, 0.85);
  legend_stack_z->SetTextSize(0.08);
  
  auto legend_stack_jt = new TLegend(0.2, 0.65, 0.4, 0.85);
  legend_stack_jt->SetTextSize(0.08);
  
  auto legend_stack_r = new TLegend(0.6, 0.6, 0.8, 0.8);
  legend_stack_r->SetTextSize(0.08);  
  
  double left = 0.2;
  double step = 0.1;
  double top = 0.8 - step;
  double top2 = 0.5 - step;
  double top2d = 0.9 - step;
  double step2d = 0.06;

  std::cout << "made it before 1D z plots" << std::endl;

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 800, 1100);
  ccan[ican]->Divide(2, 3, 0.0001, 0.0001);

  std::cout << "made it after canvas divide" << std::endl;
  //ccan[ican]->cd();
  // z histograms // 
  top = 0.8;
  step = 0.12;
  ofstream hepdata;
  for (int i = 1; i < ptbinsize; i++)
  {
    sprintf(buf, "ccan%d", ican);
    ccan[ican]->SetFillColor(10);

    ccan[ican]->cd(i);

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
    for (int ibin=0; ibin<h1_z_data_w_sys->GetNbinsX(); ++ibin)
    { 
      h1_z_data_w_sys->SetBinError(ibin+1, h1_z_sys[i-1]->GetBinContent(ibin+1));
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
    h1_z_truth[i-1]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#it{z}}");
    h1_z_truth[i-1]->GetYaxis()->SetTitleOffset(0.65);
    h1_z_truth[i-1]->GetYaxis()->SetRangeUser(0, 4.);    
    h1_z_truth[i-1]->GetYaxis()->SetTitleSize(0.1);
    h1_z_truth[i-1]->GetYaxis()->SetTitleFont(lhcbFont);
    h1_z_truth[i-1]->GetYaxis()->SetLabelSize(0.08);
    h1_z_truth[i-1]->GetYaxis()->SetLabelFont(lhcbFont);
    h1_z_truth[i-1]->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_z_truth[i-1]->GetYaxis()->ChangeLabel(-1, -1, 0); 

    h1_z_truth[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0); 
    //h1_z_data[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0); 
    //h1_z_data_w_sys->GetXaxis()->SetRangeUser(0.1, 1.0); 


    h1_z_truth[i-1]->Draw("HIST SAME");
    h1_z_truth[i-1]->Draw("PE SAME");
    h1_z_data_w_sys->Draw("PE2 SAME");
    h1_z_data[i-1]->Draw("PE SAME");    
 
    Tl.DrawLatex(left, top, "#scale[1.5]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
    Tl.DrawLatex(left, top - step, "#scale[1.5]{AK5 #it{B}^{#pm}-tagged jets}");
    Tl.DrawLatex(left, top - 2 * step, Form("#scale[1.5]{%.0f < #it{p}_{T,jet} < %.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    Tl.DrawLatex(left, top - 3 * step, "#scale[1.5]{2.5 < #it{y}_{jet} < 4}");
    if (i==1)
    {
      //lhcbName->Draw("SAME");

      legend_stack_z->SetBorderSize(0);
      legend_stack_z->SetFillStyle(0);
      legend_stack_z->SetFillColor(3);
      legend_stack_z->AddEntry(h1_z_data[i-1], "Data");
      legend_stack_z->AddEntry(h1_z_truth[i-1], str_pythia);
    }
    legend_stack_z->Draw("SAME");

    //else if (i == 2)
    //{

    //}
    //else
    //{
    //  legend_stack_z->Draw("SAME");
    //}

    bot_pad->cd();
    h1_z_data_sys_ratio->GetYaxis()->SetTitle("Sim/Data");
    h1_z_data_sys_ratio->GetYaxis()->SetTitleOffset(0.4);
    h1_z_data_sys_ratio->GetYaxis()->SetTitleSize(0.12);
    h1_z_data_sys_ratio->GetYaxis()->SetTitleFont(lhcbFont);
    h1_z_data_sys_ratio->GetYaxis()->SetLabelSize(0.1);
    h1_z_data_sys_ratio->GetYaxis()->SetLabelFont(lhcbFont);
    h1_z_data_sys_ratio->GetXaxis()->SetTitle("#it{z}");
    h1_z_data_sys_ratio->GetXaxis()->SetTitleSize(0.15);
    h1_z_data_sys_ratio->GetXaxis()->SetTitleOffset(0.65);
    h1_z_data_sys_ratio->GetXaxis()->SetTitleFont(lhcbFont);
    h1_z_data_sys_ratio->GetXaxis()->SetLabelSize(0.1);
    h1_z_data_sys_ratio->GetXaxis()->SetLabelFont(lhcbFont);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(3, -1, 0);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(5, -1, 0);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(-3, -1, 0);
    h1_z_data_sys_ratio->GetYaxis()->ChangeLabel(-5, -1, 0);
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


    ccan[ican]->cd(i);
    ccan[ican]->Update();

    //pad = (TPad*)ccan[ican]->cd(i);
    //pad->Draw();
    //pad->SaveAs(Form("hepdata/Figure1.%d.png", i));

    hepdata.open(Form("hepdata/fig1-%d.yaml",i));
    hepdata << "independent_variables:\n";
    hepdata << "- header: {name: '$z$'}\n";
    hepdata << "  values:\n";
    for (int ibin=0; ibin<h1_z_data[i-1]->GetNbinsX(); ++ibin)
    { 
      if (ibin == 0 || ibin == h1_z_data[i-1]->GetNbinsX()-1)
        continue;
      hepdata << "  - { low: " << z_binedges[ibin] << ", high: " << z_binedges[ibin+1]<< "}\n";
    }
    hepdata << "dependent_variables:\n";
    hepdata << "- header: {name: '$\\frac{1}{N_{\\rm{jets}}}\\frac{dN_{B}}{dz}$'}\n";
    hepdata << "  qualifiers:\n";
    hepdata << "  - {name: 'region', value: '$" << pt_binedges[i] << " < p_{T,\\rm{jet}} < " <<  pt_binedges[i+1] <<"$, $2.5 < y_{jet} < 4.0$'}\n";
    hepdata << "  values:\n";
    for (int ibin=0; ibin<h1_z_data[i-1]->GetNbinsX(); ++ibin)
    { 
      if (ibin == 0 || ibin == h1_z_data[i-1]->GetNbinsX()-1)
        continue;
      hepdata << "  - value: " << h1_z_data[i-1]->GetBinContent(ibin+1) << "\n";
      hepdata << "    errors:\n";
      hepdata << "    - {symerror: " << h1_z_data[i-1]->GetBinError(ibin+1) << ", label: 'Statistical'}\n"; 
      hepdata << "    - {symerror: " << h1_z_data_w_sys->GetBinError(ibin+1) << ", label: 'Systematic'}\n"; 
    }
    hepdata.close();
  }     

  if (ican == 0)
  {
    ccan[ican]->Print(plotfileO.Data());
  }
  else
  {
    ccan[ican]->Print(plotfilePDF.Data());
  }
  ccan[ican]->SaveAs(TString(loc_plots + "z_final_w_sys.png"));
  for (int i = 1; i < ptbinsize; ++i) {
      TPad *pad = (TPad*)ccan[ican]->GetPad(i);
      if (!pad) continue;

      TCanvas *tmp = new TCanvas(Form("pad_canvas_%d", i), "", 800, 600);

      tmp->cd();

      // Copy all primitives from the pad
      TIter next(pad->GetListOfPrimitives());
      TObject *obj;
      while ((obj = next())) {
          // Draw a clone of the object
          obj->DrawClone();
      }

      // Optional: copy margins
      tmp->SetLeftMargin(pad->GetLeftMargin());
      tmp->SetRightMargin(pad->GetRightMargin());
      tmp->SetTopMargin(pad->GetTopMargin());
      tmp->SetBottomMargin(pad->GetBottomMargin());

      tmp->Update();
      tmp->SaveAs(Form("hepdata/fig1-%d.png", i));

      delete tmp;
  }

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 800, 1100);
  ccan[ican]->Divide(2, 3, 0.0001, 0.0001);
  ccan[ican]->cd();


  // jT historgrams //
  top = 0.825;
  step = 0.12;
  left = 0.55;
  for (int i = 1; i < ptbinsize; i++)
  {    
    sprintf(buf, "ccan%d", ican);
    ccan[ican]->SetFillColor(10);
   
    ccan[ican]->cd(i);

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
    for (int ibin=0; ibin<h1_jt_data_w_sys->GetNbinsX(); ++ibin)
    { 
      h1_jt_data_w_sys->SetBinError(ibin+1, h1_jt_sys[i-1]->GetBinContent(ibin+1));
      if (h1_jt_data[i-1]->GetBinContent(ibin+1) == 0)
      {
        continue;
      }      
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
    h1_jt_truth[i-1]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#it{j}_{T}}");
    h1_jt_truth[i-1]->GetYaxis()->SetTitleOffset(0.65);
    h1_jt_truth[i-1]->GetYaxis()->SetRangeUser(0, 100.);    
    h1_jt_truth[i-1]->GetYaxis()->SetTitleSize(0.1);
    h1_jt_truth[i-1]->GetYaxis()->SetTitleFont(lhcbFont);
    h1_jt_truth[i-1]->GetYaxis()->SetLabelSize(0.08);
    h1_jt_truth[i-1]->GetYaxis()->SetLabelFont(lhcbFont);
    h1_jt_truth[i-1]->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_jt_truth[i-1]->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_jt_truth[i-1]->SetMinimum(0.00001);
    h1_jt_truth[i-1]->Draw("HIST SAME");
    h1_jt_truth[i-1]->Draw("PE SAME");
    //h1_ptktdR_data_stat->Draw("HIST SAME");
    h1_jt_data_w_sys->SetMinimum(0.00001);    
    h1_jt_data[i-1]->SetMinimum(0.00001);
    h1_jt_data_w_sys->Draw("PE2 SAME");    
    h1_jt_data[i-1]->Draw("PE SAME");
      
    //Tl.DrawLatex(left, top, Form("#scale[1.5]{%.0f<p_{T,jet}<%.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    //if (i == 1)
    //{
      //lhcbName->Draw("SAME");
      //Tl.DrawLatex(left, top - step, "#scale[1.5]{2.5 < y_{jet} < 4}");
      //Tl.DrawLatex(left, top - 2 * step, "#scale[1.5]{AK5 B^{#pm}-tagged jets}");
      //Tl.DrawLatex(left, top - 3 * step, "#scale[1.5]{pp #sqrt{s} = 13 TeV}");
    //}
    Tl.DrawLatex(left, top, "#scale[1.5]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
    Tl.DrawLatex(left, top - step, "#scale[1.5]{AK5 #it{B}^{#pm}-tagged jets}");
    Tl.DrawLatex(left, top - 2 * step, Form("#scale[1.5]{%.0f < #it{p}_{T,jet} < %.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    if (i != 6 )
    { 
      Tl.DrawLatex(left, top - 3 * step, "#scale[1.5]{2.5 < #it{y}_{jet} < 4}");
    }
    else
    {
      Tl.DrawLatex(0.725, top - 3 * step, "#scale[1.5]{2.5 < #it{y}_{jet} < 4}");
    }
    if (i == 1)
    {
      legend_stack_jt->SetBorderSize(0);
      legend_stack_jt->SetFillStyle(0);
      legend_stack_jt->SetFillColor(3);
      legend_stack_jt->AddEntry(h1_jt_data[i-1], "Data");
      legend_stack_jt->AddEntry(h1_jt_truth[i-1], str_pythia);
    }
    legend_stack_jt->Draw("SAME");


    bot_pad->cd();
    h1_jt_data_sys_ratio->GetYaxis()->SetTitle("Sim/Data");
    h1_jt_data_sys_ratio->GetYaxis()->SetTitleOffset(0.4);
    h1_jt_data_sys_ratio->GetYaxis()->SetTitleSize(0.12);
    h1_jt_data_sys_ratio->GetYaxis()->SetTitleFont(lhcbFont);
    h1_jt_data_sys_ratio->GetYaxis()->SetLabelSize(0.1);
    h1_jt_data_sys_ratio->GetYaxis()->SetLabelFont(lhcbFont);
    h1_jt_data_sys_ratio->GetXaxis()->SetTitle("#it{j}_{T} [GeV/c]");
    h1_jt_data_sys_ratio->GetXaxis()->SetTitleSize(0.12);
    h1_jt_data_sys_ratio->GetXaxis()->SetTitleOffset(0.9);
    h1_jt_data_sys_ratio->GetXaxis()->SetTitleFont(lhcbFont);
    h1_jt_data_sys_ratio->GetXaxis()->SetLabelSize(0.1);
    h1_jt_data_sys_ratio->GetXaxis()->SetLabelFont(lhcbFont);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(3, -1, 0);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(5, -1, 0);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(-3, -1, 0);
    h1_jt_data_sys_ratio->GetYaxis()->ChangeLabel(-5, -1, 0);
    h1_jt_data_sys_ratio->SetMinimum(0.0);
    h1_jt_data_sys_ratio->SetMaximum(2.0);
    h1_jt_data_sys_ratio->SetFillColorAlpha(16, 0.5);
    h1_jt_data_sys_ratio->SetLineColorAlpha(0, 0.5);
    h1_jt_data_sys_ratio->Draw("PE2 SAME");
    h1_jt_data_stat_ratio->Draw("PE SAME");
    TLine *line_jt = new TLine(0.0, 1.0, 10.0, 1.0);
    line_jt->SetLineStyle(2);
    line_jt->Draw("SAME");

    ccan[ican]->cd(i);
    ccan[ican]->Update(); 

    hepdata.open(Form("hepdata/fig2-%d.yaml",i));
    hepdata << "independent_variables:\n";
    hepdata << "- header: {name: '$j_{T}$ [GeV/c]'}\n";
    hepdata << "  values:\n";
    for (int ibin=0; ibin<h1_jt_data[i-1]->GetNbinsX(); ++ibin)
    { 
      if (h1_jt_data[i-1]->GetBinContent(ibin+1) !=0 )
      {
        hepdata << "  - { low: " << jt_binedges[ibin] << ", high: " << jt_binedges[ibin+1]<< "}\n";
      }
    }
    hepdata << "dependent_variables:\n";
    hepdata << "- header: {name: '$\\frac{1}{N_{\\rm{jets}}}\\frac{dN_{B}}{dj_{T}}$'}\n";
    hepdata << "  qualifiers:\n";
    hepdata << "  - {name: 'region', value: '$" << pt_binedges[i] << " < p_{T,\\rm{jet}} < " <<  pt_binedges[i+1] <<"$, $2.5 < y_{jet} < 4.0$'}\n";
    hepdata << "  values:\n";
    for (int ibin=0; ibin<h1_jt_data[i-1]->GetNbinsX(); ++ibin)
    { 
      if (h1_jt_data[i-1]->GetBinContent(ibin+1) !=0 )
      {
        hepdata << "  - value: " << h1_jt_data[i-1]->GetBinContent(ibin+1) << "\n";
        hepdata << "    errors:\n";
        hepdata << "    - {symerror: " << h1_jt_data[i-1]->GetBinError(ibin+1) << ", label: 'Statistical'}\n"; 
        hepdata << "    - {symerror: " << h1_jt_data_w_sys->GetBinError(ibin+1) << ", label: 'Systematic'}\n"; 
      }
    }
    hepdata.close();
  }

  if (ican == 0)
  {
    ccan[ican]->Print(plotfileO.Data());
  }
  else
  {
    ccan[ican]->Print(plotfilePDF.Data());
  }      
  ccan[ican]->SaveAs(TString(loc_plots + "jt_final_w_sys.png"));  
  for (int i = 1; i < ptbinsize; ++i) {
    TPad *pad = (TPad*)ccan[ican]->GetPad(i);
    if (!pad) continue;

    TCanvas *tmp = new TCanvas(Form("pad_canvas_%d", i), "", 800, 600);

    tmp->cd();

    // Copy all primitives from the pad
    TIter next(pad->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        // Draw a clone of the object
        obj->DrawClone();
    }

    // Optional: copy margins
    tmp->SetLeftMargin(pad->GetLeftMargin());
    tmp->SetRightMargin(pad->GetRightMargin());
    tmp->SetTopMargin(pad->GetTopMargin());
    tmp->SetBottomMargin(pad->GetBottomMargin());

    tmp->Update();
    tmp->SaveAs(Form("hepdata/fig2-%d.png", i));

    delete tmp;
}


  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 800, 1100);
  ccan[ican]->Divide(2, 3, 0.0001, 0.0001);
  ccan[ican]->cd();

  // r histograms //
  //top = 0.96;
  //step = 0.06;
  //left = 0.175;
  top = 0.35;
  step = 0.12;
  left = 0.2;
  for (int i = 1; i < ptbinsize; i++)
  {     
    ccan[ican]->SetFillColor(10);

    ccan[ican]->cd(i);
    
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
    for (int ibin=0; ibin<h1_r_data_w_sys->GetNbinsX(); ++ibin)
    { 
      h1_r_data_w_sys->SetBinError(ibin+1, h1_r_sys[i-1]->GetBinContent(ibin+1));
      if (h1_r_data[i-1]->GetBinContent(ibin+1) == 0)
      {
        continue;
      }   
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
    h1_r_truth[i-1]->GetYaxis()->SetTitle("#frac{1}{N_{jets}}#frac{dN}{d#it{r}}");
    h1_r_truth[i-1]->GetYaxis()->SetTitleOffset(0.65);
    h1_r_truth[i-1]->GetYaxis()->SetRangeUser(0, 100.);    
    h1_r_truth[i-1]->GetYaxis()->SetTitleSize(0.1);
    h1_r_truth[i-1]->GetYaxis()->SetTitleFont(lhcbFont);
    h1_r_truth[i-1]->GetYaxis()->SetLabelSize(0.08);
    h1_r_truth[i-1]->GetYaxis()->SetLabelFont(lhcbFont);
    h1_r_truth[i-1]->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_r_truth[i-1]->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_r_truth[i-1]->SetMinimum(0.01);
    h1_r_truth[i-1]->GetXaxis()->SetRangeUser(0.0, 0.4);
    h1_r_truth[i-1]->Draw("HIST SAME");
    h1_r_truth[i-1]->Draw("PE SAME");
    h1_r_data_w_sys->SetMinimum(0.01);
    h1_r_data[i-1]->SetMinimum(0.01);
    h1_r_data_w_sys->Draw("PE2 SAME");    
    h1_r_data[i-1]->Draw("PE SAME");

    /*
    Tl.DrawLatex(left, top, Form("#scale[1.5]{%.0f<p_{T,jet}<%.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    if (i==1)
    {
      //lhcbName->Draw("SAME");
      Tl.DrawLatex(left, top - step, "#scale[1.5]{2.5 < y_{jet} < 4}");
      Tl.DrawLatex(left, top - 2 * step, "#scale[1.5]{AK5 B^{#pm}-tagged jets}");
      Tl.DrawLatex(0.525, top - step, "#scale[1.5]{pp #sqrt{s} = 13 TeV}");
    }
    */
    Tl.DrawLatex(left , top, "#scale[1.5]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
    Tl.DrawLatex(left, top - step, "#scale[1.5]{AK5 #it{B}^{#pm}-tagged jets}");
    Tl.DrawLatex(left, top - 2 * step, Form("#scale[1.5]{%.0f < #it{p}_{T,jet} < %.0f GeV/c,   2.5 < #it{y}_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top - 3 * step, "#scale[1.5]{2.5 < #it{y}_{jet} < 4}");
    if (i == 1)
    {
      legend_stack_r->SetBorderSize(0);
      legend_stack_r->SetFillStyle(0);
      legend_stack_r->SetFillColor(3);
      legend_stack_r->AddEntry(h1_r_data[i-1], "Data");
      legend_stack_r->AddEntry(h1_r_truth[i-1], str_pythia);
    }
    legend_stack_r->Draw("SAME");


    bot_pad->cd();
    h1_r_data_sys_ratio->GetYaxis()->SetTitle("Sim/Data");
    h1_r_data_sys_ratio->GetYaxis()->SetTitleOffset(0.4);
    h1_r_data_sys_ratio->GetYaxis()->SetTitleSize(0.12);
    h1_r_data_sys_ratio->GetYaxis()->SetTitleFont(lhcbFont);
    h1_r_data_sys_ratio->GetYaxis()->SetLabelSize(0.1);
    h1_r_data_sys_ratio->GetYaxis()->SetLabelFont(lhcbFont);
    h1_r_data_sys_ratio->GetXaxis()->SetTitle("#it{r}");
    h1_r_data_sys_ratio->GetXaxis()->SetTitleSize(0.15);
    h1_r_data_sys_ratio->GetXaxis()->SetTitleOffset(0.65);
    h1_r_data_sys_ratio->GetXaxis()->SetTitleFont(lhcbFont);
    h1_r_data_sys_ratio->GetXaxis()->SetLabelSize(0.1);
    h1_r_data_sys_ratio->GetXaxis()->SetLabelFont(lhcbFont);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(1, -1, 0);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(3, -1, 0);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(5, -1, 0);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(-1, -1, 0);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(-3, -1, 0);
    h1_r_data_sys_ratio->GetYaxis()->ChangeLabel(-5, -1, 0); 
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

    ccan[ican]->cd(i);
    ccan[ican]->Update();

    hepdata.open(Form("hepdata/fig3-%d.yaml",i));
    hepdata << "independent_variables:\n";
    hepdata << "- header: {name: '$r$'}\n";
    hepdata << "  values:\n";
    for (int ibin=0; ibin<h1_r_data[i-1]->GetNbinsX(); ++ibin)
    { 
      if (ibin == h1_r_data[i-1]->GetNbinsX()-1)
        continue;
      hepdata << "  - { low: " << r_binedges[ibin] << ", high: " << r_binedges[ibin+1]<< "}\n";
    }
    hepdata << "dependent_variables:\n";
    hepdata << "- header: {name: '$\\frac{1}{N_{\\rm{jets}}}\\frac{dN_{B}}{dr}$'}\n";
    hepdata << "  qualifiers:\n";
    hepdata << "  - {name: 'region', value: '$" << pt_binedges[i] << " < p_{T,\\rm{jet}} < " <<  pt_binedges[i+1] <<"$, $2.5 < y_{jet} < 4.0$'}\n";
    hepdata << "  values:\n";
    for (int ibin=0; ibin<h1_r_data[i-1]->GetNbinsX(); ++ibin)
    { 
      if (ibin == h1_r_data[i-1]->GetNbinsX()-1)
        continue;
      hepdata << "  - value: " << h1_r_data[i-1]->GetBinContent(ibin+1) << "\n";
      hepdata << "    errors:\n";
      hepdata << "    - {symerror: " << h1_r_data[i-1]->GetBinError(ibin+1) << ", label: 'Statistical'}\n"; 
      hepdata << "    - {symerror: " << h1_r_data_w_sys->GetBinError(ibin+1) << ", label: 'Systematic'}\n"; 
    }
    hepdata.close();
  
  }  

  if (ican == 0)
  {
    ccan[ican]->Print(plotfileO.Data());
  }
  else
  {
    ccan[ican]->Print(plotfilePDF.Data());
  }
  ccan[ican]->SaveAs(TString(loc_plots + "r_final_w_sys.png"));   
  for (int i = 1; i < ptbinsize; ++i) {
    TPad *pad = (TPad*)ccan[ican]->GetPad(i);
    if (!pad) continue;

    TCanvas *tmp = new TCanvas(Form("pad_canvas_%d", i), "", 800, 600);

    tmp->cd();

    // Copy all primitives from the pad
    TIter next(pad->GetListOfPrimitives());
    TObject *obj;
    while ((obj = next())) {
        // Draw a clone of the object
        obj->DrawClone();
    }

    // Optional: copy margins
    tmp->SetLeftMargin(pad->GetLeftMargin());
    tmp->SetRightMargin(pad->GetRightMargin());
    tmp->SetTopMargin(pad->GetTopMargin());
    tmp->SetBottomMargin(pad->GetBottomMargin());

    tmp->Update();
    tmp->SaveAs(Form("hepdata/fig3-%d.png", i));

    delete tmp;
}


  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 800, 1100);
  ccan[ican]->SetFillColor(10);

  ccan[ican]->cd();
  ccan[ican]->Divide(2, 3, 0.0001, 0.0001);
  // zjt histos // 
  top = 0.96;
  step = 0.06;
  left = 0.175;
  for (int i = 1; i < ptbinsize; i++)
  {

    //gStyle->cd();
    //gROOT->SetStyle(gStyle->GetName());

    ccan[ican]->cd(i);
    
    h3_ptzjt_data->GetZaxis()->SetRange(i+1, i+1);
    //h3_ptzjt_truth->GetZaxis()->SetRange(i+1, i+1);   
    
    h2_zjt_data[i-1] = (TH2D *)h3_ptzjt_data->Project3D("yx");
    //h2_zjt_truth[i-1] = (TH2D *)h3_ptzjt_truth->Project3D("yx");     
    h2_zjt_data[i-1]->SetName(Form("zjt_pt%d",i));     
    //h2_zjt_truth[i-1]->SetName(Form("zjt_truth_pt%d",i));     
    
    NormalizeHist(h2_zjt_data[i-1]);
    //NormalizeHist(h2_zjt_truth[i-1]); 

    gPad->SetTopMargin(0.2);
    gPad->SetBottomMargin(0.125);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    h2_zjt_data[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zjt_data[i-1]->GetXaxis()->SetTitle("#it{z}");
    //h2_zjt_data[i-1]->GetXaxis()->SetTitleSize(0.05);
    //h2_zjt_data[i-1]->GetXaxis()->SetTitleOffset(0.85);
    //h2_zjt_data[i-1]->GetXaxis()->SetLabelSize(0.04);
    h2_zjt_data[i-1]->GetYaxis()->SetTitle("#it{j}_{T} [GeV/c]");
    //h2_zjt_data[i-1]->GetYaxis()->SetTitleSize(0.05);
    //h2_zjt_data[i-1]->GetYaxis()->SetTitleOffset(0.85);
    //h2_zjt_data[i-1]->GetYaxis()->SetLabelSize(0.04);
    //h2_zjt_data[i-1]->GetZaxis()->SetLabelSize(0.04);
    h2_zjt_data[i-1]->SetMinimum(0.00001);  
    h2_zjt_data[i-1]->SetMaximum(10.);      
    //h2_zjt_data[i-1]->GetYaxis()->SetLimits(0.0, 20.0);
    h2_zjt_data[i-1]->UseCurrentStyle();
    h2_zjt_data[i-1]->Draw("COLZ");


    Tl.DrawLatex(left, top, "#scale[1.0]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
    Tl.DrawLatex(left, top - step, "#scale[1.0]{AK5 #it{B}^{#pm}-tagged jets}");
    Tl.DrawLatex(left, top - 2.1 * step, Form("#scale[0.85]{%.0f < #it{p}_{T,jet} < %.0f GeV/c,   2.5 < #it{y}_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top - 3 * step, "#scale[0.75]{2.5 < #it{y}_{jet} < 4}");
    //lhcbName->Draw("SAME");
    /*
    Tl.DrawLatex(left, top, Form("#scale[1.0]{%.0f<p_{T,jet}<%.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    if (i==1)
    {
      //lhcbName->Draw("SAME");
      Tl.DrawLatex(left, top - step, "#scale[1.0]{2.5 < y_{jet} < 4}");
      Tl.DrawLatex(left, top - 2*step, "#scale[1.0]{AK5 B^{#pm}-tagged jets}");
      Tl.DrawLatex(left, top - 3*step, "#scale[1.0]{pp #sqrt{s} = 13 TeV}");

    }
    */
    ccan[ican]->cd();
    //gPad->SetLogz();
    ccan[ican]->Update(); 

    h2_zjt_sys[i-1] = (TH2D*)file_sys->Get(Form("zjt_sys_total_pt%d",i));   
    TH2D *h2_zjt_data_w_sys = (TH2D*)h2_zjt_data[i-1]->Clone("h2_zjt_data_w_sys");
    h2_zjt_sys[i-1]->Multiply(h2_zjt_sys[i-1], h2_zjt_data[i-1]);
    for (int xbin=0; xbin<h2_zjt_data_w_sys->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_zjt_data_w_sys->GetNbinsY(); ++ybin)
      { 
        h2_zjt_data_w_sys->SetBinError(xbin+1, ybin+1, h2_zjt_sys[i-1]->GetBinContent(xbin+1, ybin+1));
        if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
        {
          continue;
        }   
      }    
    } 

    hepdata.open(Form("hepdata/fig4-%d.yaml",i));
    hepdata << "independent_variables:\n";
    hepdata << "- header: {name: '$z$'}\n";
    //hepdata << "- header: {name: '$j_{T} [GeV/c]$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
    { 
      if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
        continue;
      for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - { low: " << z_binedges2D[xbin] << ", high: " << z_binedges2D[xbin+1]<< "}\n";
        }
      }
    }
    hepdata << "- header: {name: '$j_{T}$ [GeV/c]'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
    {
      if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
        continue;
      for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - { low: " << jt_binedges2D[ybin] << ", high: " << jt_binedges2D[ybin+1]<< "}\n";
        }
      }
    }
    hepdata << "dependent_variables:\n";
    hepdata << "- header: {name: '$\\frac{1}{N_{\\rm{jets}}}\\frac{dN_{B}}{dzdj_{T}}$'}\n";
    hepdata << "  qualifiers:\n";
    hepdata << "  - {name: 'region', value: '$" << pt_binedges[i] << " < p_{T,\\rm{jet}} < " <<  pt_binedges[i+1] <<"$, $2.5 < y_{jet} < 4.0$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_zjt_data[i-1]->GetNbinsX(); ++xbin)
    {
      if (xbin == h2_zjt_data[i-1]->GetNbinsX()-1)
        continue; 
      for (int ybin=0; ybin<h2_zjt_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - value: " << h2_zjt_data[i-1]->GetBinContent(xbin+1, ybin+1) << "\n";
          hepdata << "    errors:\n";
          hepdata << "    - {symerror: " << h2_zjt_data[i-1]->GetBinError(xbin+1, ybin+1) << ", label: 'Statistical'}\n"; 
          hepdata << "    - {symerror: " << h2_zjt_data_w_sys->GetBinError(xbin+1, ybin+1) << ", label: 'Systematic'}\n"; 
        }
      }
    }
    hepdata.close();
  }    

  if (ican == 0)
  {
    ccan[ican]->Print(plotfileO.Data());
  }
  else
  {
    ccan[ican]->Print(plotfilePDF.Data());
  }
  ccan[ican]->SaveAs(TString(loc_plots + "zjt_final_w_sys.png"));  
  for (int i = 1; i < ptbinsize; ++i) {
    TPad *pad = (TPad*)ccan[ican]->GetPad(i);
    if (!pad) continue;
    TCanvas *tmp = new TCanvas(Form("pad_canvas_%d", i), "", 800, 600);

    // ---- copy pad state ----
    tmp->SetLogx(pad->GetLogx());
    tmp->SetLogy(pad->GetLogy());
    tmp->SetLogz(pad->GetLogz());

    tmp->SetLeftMargin(pad->GetLeftMargin());
    tmp->SetRightMargin(pad->GetRightMargin());
    tmp->SetTopMargin(pad->GetTopMargin());
    tmp->SetBottomMargin(pad->GetBottomMargin());

    tmp->cd();

    // ---- find the main TH2 ----
    TH2 *h2 = nullptr;
    TIter it(pad->GetListOfPrimitives());
    TObject *obj;

    while ((obj = it())) {
        if ((h2 = dynamic_cast<TH2*>(obj))) break;
    }

    // ---- redraw TH2 properly ----
    if (h2) {
        TH2 *h2c = (TH2*)h2->Clone();
        h2c->Draw("COLZ");   // force palette regeneration
    }

    // ---- redraw everything else on top ----
    it.Reset();
    while ((obj = it())) {
      if (obj == h2) continue;
      if (obj->InheritsFrom("TPaletteAxis")) continue;
      if (obj->InheritsFrom("TH1")) continue;    // includes TH2
      if (obj->InheritsFrom("TFrame")) continue; // very important

      // Draw annotations without "same"
      if (obj->InheritsFrom("TLatex") ||
          obj->InheritsFrom("TText")  ||
          obj->InheritsFrom("TLegend")||
          obj->InheritsFrom("TLine")  ||
          obj->InheritsFrom("TArrow") ||
          obj->InheritsFrom("TBox")) {

            obj->DrawClone();
      }
      // Graph-like overlays
      else if (obj->InheritsFrom("TGraph")) {
        obj->DrawClone("same");
      }
    }

    tmp->SaveAs(Form("hepdata/fig4-%d.png", i));

    delete tmp; 
  }
 

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 800, 1100);
  ccan[ican]->SetFillColor(10);

  ccan[ican]->cd();
  ccan[ican]->Divide(2, 3, 0.0001, 0.0001);
  // zr plots //
  top = 0.96;
  step = 0.06;
  left = 0.175;
  for (int i = 1; i < ptbinsize; i++)
  {
    ccan[ican]->cd(i);
    
    h3_ptzr_data->GetZaxis()->SetRange(i+1, i+1);
    //h3_ptzr_truth->GetZaxis()->SetRange(i+1, i+1);   
    
    h2_zr_data[i-1] = (TH2D *)h3_ptzr_data->Project3D("yx");
    //h2_zr_truth[i-1] = (TH2D *)h3_ptzr_truth->Project3D("yx");     
    h2_zr_data[i-1]->SetName(Form("zr_pt%d",i));     
    //h2_zr_truth[i-1]->SetName(Form("zr_truth_pt%d",i));     
    
    NormalizeHist(h2_zr_data[i-1]);
    //NormalizeHist(h2_zr_truth[i-1]); 
    

    gPad->SetTopMargin(0.2);
    gPad->SetBottomMargin(0.125);
    gPad->SetRightMargin(0.15);
    gPad->SetLogz();
    h2_zr_data[i-1]->SetMinimum(0.00001);  
    h2_zr_data[i-1]->SetMaximum(200.);      
    h2_zr_data[i-1]->GetXaxis()->SetRangeUser(0.1, 1.0);
    h2_zr_data[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_zr_data[i-1]->GetXaxis()->SetTitle("#it{z}");
    //h2_zr_data[i-1]->GetXaxis()->SetTitleSize(0.05);
    //h2_zr_data[i-1]->GetXaxis()->SetTitleOffset(0.65);
    //h2_zr_data[i-1]->GetXaxis()->SetLabelSize(0.04);
    h2_zr_data[i-1]->GetYaxis()->SetTitle("#it{r}");
    //h2_zr_data[i-1]->GetYaxis()->SetTitleSize(0.05);
    //h2_zr_data[i-1]->GetYaxis()->SetTitleOffset(1.0);
    //h2_zr_data[i-1]->GetYaxis()->SetLabelSize(0.04);
    //h2_zr_data[i-1]->GetZaxis()->SetLabelSize(0.04);
    h2_zr_data[i-1]->UseCurrentStyle();

    h2_zr_data[i-1]->Draw("COLZ");


    Tl.DrawLatex(left, top, "#scale[1.0]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
    Tl.DrawLatex(left, top - step, "#scale[1.0]{AK5 #it{B}^{#pm}-tagged jets}");
    Tl.DrawLatex(left, top - 2.1 * step, Form("#scale[0.85]{%.0f < #it{p}_{T,jet} < %.0f GeV/c,   2.5 < #it{y}_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top - 3 * step, "#scale[0.75]{2.5 < #it{y}_{jet} < 4}");
    //lhcbName->Draw("SAME");
    /*
    Tl.DrawLatex(left, top, Form("#scale[1.0]{%.0f<p_{T,jet}<%.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    if (i==1)
    {
      //lhcbName->Draw("SAME");
      Tl.DrawLatex(left, top - step, "#scale[1.0]{2.5 < y_{jet} < 4}");
      Tl.DrawLatex(left, top - 2*step, "#scale[1.0]{AK5 B^{#pm}-tagged jets}");
      Tl.DrawLatex(left, top - 3*step, "#scale[1.0]{pp #sqrt{s} = 13 TeV}");

    }
    */
    //Tl.DrawLatex(left, top2d + step2d, "#scale[0.8]{AK5 B^{#pm}-tagged jets              pp #sqrt{s} = 13 TeV}");    
    //Tl.DrawLatex(left, top2d, "#scale[0.8]{}");
    //Tl.DrawLatex(left, top2d, Form("#scale[0.8]{%.0f<p_{T,jet}<%.0f GeV/c        2.5 < y_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top2d - 2 * step2d, "#scale[0.8]{}");  
    
    ccan[ican]->cd();
    //gPad->SetLogz();

    h2_zr_sys[i-1] = (TH2D*)file_sys->Get(Form("zr_sys_total_pt%d",i));   
    TH2D *h2_zr_data_w_sys = (TH2D*)h2_zr_data[i-1]->Clone("h2_zr_data_w_sys");
    h2_zr_sys[i-1]->Multiply(h2_zr_sys[i-1], h2_zr_data[i-1]);
    for (int xbin=0; xbin<h2_zr_data_w_sys->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_zr_data_w_sys->GetNbinsY(); ++ybin)
      { 
        h2_zr_data_w_sys->SetBinError(xbin+1, ybin+1, h2_zr_sys[i-1]->GetBinContent(xbin+1, ybin+1));
        if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
        {
          continue;
        }   
      }    
    } 

    hepdata.open(Form("hepdata/fig5-%d.yaml",i));
    hepdata << "independent_variables:\n";
    hepdata << "- header: {name: '$z$'}\n";
    //hepdata << "- header: {name: '$j_{T} [GeV/c]$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
        continue;
      for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - { low: " << z_binedges2D[xbin] << ", high: " << z_binedges2D[xbin+1]<< "}\n";
        }
      }
    }
    hepdata << "- header: {name: '$r$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
    {
      if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
        continue;
      for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - { low: " << r_binedges2D[ybin] << ", high: " << r_binedges2D[ybin+1]<< "}\n";
        }
      }
    }
    hepdata << "dependent_variables:\n";
    hepdata << "- header: {name: '$\\frac{1}{N_{\\rm{jets}}}\\frac{dN_{B}}{dzdr}$'}\n";
    hepdata << "  qualifiers:\n";
    hepdata << "  - {name: 'region', value: '$" << pt_binedges[i] << " < p_{T,\\rm{jet}} < " <<  pt_binedges[i+1] <<"$, $2.5 < y_{jet} < 4.0$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_zr_data[i-1]->GetNbinsX(); ++xbin)
    {
      if (xbin == h2_zr_data[i-1]->GetNbinsX()-1)
        continue; 
      for (int ybin=0; ybin<h2_zr_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - value: " << h2_zr_data[i-1]->GetBinContent(xbin+1, ybin+1) << "\n";
          hepdata << "    errors:\n";
          hepdata << "    - {symerror: " << h2_zr_data[i-1]->GetBinError(xbin+1, ybin+1) << ", label: 'Statistical'}\n"; 
          hepdata << "    - {symerror: " << h2_zr_data_w_sys->GetBinError(xbin+1, ybin+1) << ", label: 'Systematic'}\n"; 
        }
      }
    }
    hepdata.close();
  }       
  ccan[ican]->Update();
  if (ican == 0)
  {
    ccan[ican]->Print(plotfileO.Data());
  }
  else
  {
    ccan[ican]->Print(plotfilePDF.Data());
  }
  ccan[ican]->SaveAs(TString(loc_plots + "zr_final_w_sys.png"));
  for (int i = 1; i < ptbinsize; ++i) {
    TPad *pad = (TPad*)ccan[ican]->GetPad(i);
    if (!pad) continue;
    TCanvas *tmp = new TCanvas(Form("pad_canvas_%d", i), "", 800, 600);

    // ---- copy pad state ----
    tmp->SetLogx(pad->GetLogx());
    tmp->SetLogy(pad->GetLogy());
    tmp->SetLogz(pad->GetLogz());

    tmp->SetLeftMargin(pad->GetLeftMargin());
    tmp->SetRightMargin(pad->GetRightMargin());
    tmp->SetTopMargin(pad->GetTopMargin());
    tmp->SetBottomMargin(pad->GetBottomMargin());

    tmp->cd();

    // ---- find the main TH2 ----
    TH2 *h2 = nullptr;
    TIter it(pad->GetListOfPrimitives());
    TObject *obj;

    while ((obj = it())) {
        if ((h2 = dynamic_cast<TH2*>(obj))) break;
    }

    // ---- redraw TH2 properly ----
    if (h2) {
        TH2 *h2c = (TH2*)h2->Clone();
        h2c->Draw("COLZ");   // force palette regeneration
    }

    // ---- redraw everything else on top ----
    it.Reset();
    while ((obj = it())) {
      if (obj == h2) continue;
      if (obj->InheritsFrom("TPaletteAxis")) continue;
      if (obj->InheritsFrom("TH1")) continue;    // includes TH2
      if (obj->InheritsFrom("TFrame")) continue; // very important

      // Draw annotations without "same"
      if (obj->InheritsFrom("TLatex") ||
          obj->InheritsFrom("TText")  ||
          obj->InheritsFrom("TLegend")||
          obj->InheritsFrom("TLine")  ||
          obj->InheritsFrom("TArrow") ||
          obj->InheritsFrom("TBox")) {

            obj->DrawClone();
      }
      // Graph-like overlays
      else if (obj->InheritsFrom("TGraph")) {
        obj->DrawClone("same");
      }
    }

    tmp->SaveAs(Form("hepdata/fig5-%d.png", i));

    delete tmp; 
  }

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 800, 1100);
  ccan[ican]->SetFillColor(10);

  ccan[ican]->cd();
  ccan[ican]->Divide(2, 3, 0.0001, 0.0001);

  // jtr plots //
  top = 0.96;
  step = 0.06;
  left = 0.175;
  for (int i = 1; i < ptbinsize; i++)
  {
    ccan[ican]->cd(i);
    gPad->SetTopMargin(0.1);
    gPad->SetBottomMargin(0.12);

    h3_ptjtr_data->GetZaxis()->SetRange(i+1, i+1);
    //h3_ptjtr_truth->GetZaxis()->SetRange(i+1, i+1);   
    
    h2_jtr_data[i-1] = (TH2D *)h3_ptjtr_data->Project3D("yx");
    //h2_jtr_truth[i-1] = (TH2D *)h3_ptjtr_truth->Project3D("yx");     
    h2_jtr_data[i-1]->SetName(Form("jtr_pt%d",i));     
    //h2_jtr_truth[i-1]->SetName(Form("jtr_truth_pt%d",i));     
    
    NormalizeHist(h2_jtr_data[i-1]);
    //NormalizeHist(h2_jtr_truth[i-1]); 
    
        
    h2_jtr_data[i-1]->SetMinimum(0.00001);  
    h2_jtr_data[i-1]->SetMaximum(50.);
    h2_jtr_data[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    

    gPad->SetTopMargin(0.2);
    gPad->SetBottomMargin(0.125);
    gPad->SetRightMargin(0.15);   
    gPad->SetLogz();

    //h2_jtr_data[i-1]->GetYaxis()->SetRangeUser(0.0, 0.4);
    h2_jtr_data[i-1]->GetXaxis()->SetTitle("#it{j}_{T} [GeV/c]");
    //h2_jtr_data[i-1]->GetXaxis()->SetTitleSize(0.01);
    //h2_jtr_data[i-1]->GetXaxis()->SetTitleOffset(0.1);
    //h2_jtr_data[i-1]->GetXaxis()->SetLabelSize(0.04);
    h2_jtr_data[i-1]->GetYaxis()->SetTitle("#it{r}");
    //h2_jtr_data[i-1]->GetYaxis()->SetTitleSize(0.05);
    //h2_jtr_data[i-1]->GetYaxis()->SetTitleOffset(1.0);
    //h2_jtr_data[i-1]->GetYaxis()->SetLabelSize(0.04);
    //h2_jtr_data[i-1]->GetZaxis()->SetLabelSize(0.04);
    h2_jtr_data[i-1]->UseCurrentStyle();
    h2_jtr_data[i-1]->GetXaxis()->SetTitleOffset(0.9);


    h2_jtr_data[i-1]->Draw("COLZ");

    Tl.DrawLatex(left, top, "#scale[1.0]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
    Tl.DrawLatex(left, top - step, "#scale[1.0]{AK5 #it{B}^{#pm}-tagged jets}");
    Tl.DrawLatex(left, top - 2.1 * step, Form("#scale[0.85]{%.0f < #it{p}_{T,jet} < %.0f GeV/c,   2.5 < #it{y}_{jet} < 4}", pt_binedges[i], pt_binedges[i+1]));
    //Tl.DrawLatex(left, top - 3 * step, "#scale[0.75]{2.5 < #it{y}_{jet} < 4}");
    //lhcbName->Draw("SAME");
    /*
    Tl.DrawLatex(left, top, Form("#scale[1.0]{%.0f<p_{T,jet}<%.0f GeV/c}", pt_binedges[i], pt_binedges[i+1]));
    if (i==1)
    {
      //lhcbName->Draw("SAME");
      Tl.DrawLatex(left, top -  step, "#scale[1.0]{2.5 < y_{jet} < 4}");
      Tl.DrawLatex(left, top - 2 * step, "#scale[1.0]{AK5 B^{#pm}-tagged jets}");
      Tl.DrawLatex(left, top - 3 * step, "#scale[1.0]{pp #sqrt{s} = 13 TeV}");

    }
    */

    ccan[ican]->cd();
    //gPad->SetLogz();
    ccan[ican]->Update();

    h2_jtr_sys[i-1] = (TH2D*)file_sys->Get(Form("jtr_sys_total_pt%d",i));   
    TH2D *h2_jtr_data_w_sys = (TH2D*)h2_jtr_data[i-1]->Clone("h2_jtr_data_w_sys");
    h2_jtr_sys[i-1]->Multiply(h2_jtr_sys[i-1], h2_jtr_data[i-1]);
    for (int xbin=0; xbin<h2_jtr_data_w_sys->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_jtr_data_w_sys->GetNbinsY(); ++ybin)
      { 
        h2_jtr_data_w_sys->SetBinError(xbin+1, ybin+1, h2_jtr_sys[i-1]->GetBinContent(xbin+1, ybin+1));
        if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) == 0)
        {
          continue;
        }   
      }    
    } 

    hepdata.open(Form("hepdata/fig6-%d.yaml",i));
    hepdata << "independent_variables:\n";
    hepdata << "- header: {name: '$j_{T}$ [GeV/c]'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    { 
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
          continue;
        if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - { low: " << jt_binedges2D[xbin] << ", high: " << jt_binedges2D[xbin+1]<< "}\n";
        }
      }
    }
    hepdata << "- header: {name: '$r$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    {
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      {
        if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
          continue; 
        if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - { low: " << r_binedges2D[ybin] << ", high: " << r_binedges2D[ybin+1]<< "}\n";
        }
      }
    }
    hepdata << "dependent_variables:\n";
    hepdata << "- header: {name: '$\\frac{1}{N_{\\rm{jets}}}\\frac{dN_{B}}{dj_{T}dr}$'}\n";
    hepdata << "  qualifiers:\n";
    hepdata << "  - {name: 'region', value: '$" << pt_binedges[i] << " < p_{T,\\rm{jet}} < " <<  pt_binedges[i+1] <<"$, $2.5 < y_{jet} < 4.0$'}\n";
    hepdata << "  values:\n";
    for (int xbin=0; xbin<h2_jtr_data[i-1]->GetNbinsX(); ++xbin)
    {
      for (int ybin=0; ybin<h2_jtr_data[i-1]->GetNbinsY(); ++ybin)
      { 
        if (ybin == h2_jtr_data[i-1]->GetNbinsY()-1)
          continue;
        if (h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) !=0 )
        {
          hepdata << "  - value: " << h2_jtr_data[i-1]->GetBinContent(xbin+1, ybin+1) << "\n";
          hepdata << "    errors:\n";
          hepdata << "    - {symerror: " << h2_jtr_data[i-1]->GetBinError(xbin+1, ybin+1) << ", label: 'Statistical'}\n"; 
          hepdata << "    - {symerror: " << h2_jtr_data_w_sys->GetBinError(xbin+1, ybin+1) << ", label: 'Systematic'}\n"; 
        }
      }
    }
    hepdata.close();
  }    
  
  if (ican == 0)
  {
    ccan[ican]->Print(plotfileO.Data());
  }
  else
  {
    ccan[ican]->Print(plotfilePDF.Data());
  }
  ccan[ican]->SaveAs(TString(loc_plots + "jtr_final_w_sys.png"));
  for (int i = 1; i < ptbinsize; ++i) {
    TPad *pad = (TPad*)ccan[ican]->GetPad(i);
    if (!pad) continue;
    TCanvas *tmp = new TCanvas(Form("pad_canvas_%d", i), "", 800, 600);

    // ---- copy pad state ----
    tmp->SetLogx(pad->GetLogx());
    tmp->SetLogy(pad->GetLogy());
    tmp->SetLogz(pad->GetLogz());

    tmp->SetLeftMargin(pad->GetLeftMargin());
    tmp->SetRightMargin(pad->GetRightMargin());
    tmp->SetTopMargin(pad->GetTopMargin());
    tmp->SetBottomMargin(pad->GetBottomMargin());

    tmp->cd();

    // ---- find the main TH2 ----
    TH2 *h2 = nullptr;
    TIter it(pad->GetListOfPrimitives());
    TObject *obj;

    while ((obj = it())) {
        if ((h2 = dynamic_cast<TH2*>(obj))) break;
    }

    // ---- redraw TH2 properly ----
    if (h2) {
        TH2 *h2c = (TH2*)h2->Clone();
        h2c->Draw("COLZ");   // force palette regeneration
    }

    // ---- redraw everything else on top ----
    it.Reset();
    while ((obj = it())) {
      if (obj == h2) continue;
      if (obj->InheritsFrom("TPaletteAxis")) continue;
      if (obj->InheritsFrom("TH1")) continue;    // includes TH2
      if (obj->InheritsFrom("TFrame")) continue; // very important

      // Draw annotations without "same"
      if (obj->InheritsFrom("TLatex") ||
          obj->InheritsFrom("TText")  ||
          obj->InheritsFrom("TLegend")||
          obj->InheritsFrom("TLine")  ||
          obj->InheritsFrom("TArrow") ||
          obj->InheritsFrom("TBox")) {

            obj->DrawClone();
      }
      // Graph-like overlays
      else if (obj->InheritsFrom("TGraph")) {
        obj->DrawClone("same");
      }
    }

    tmp->SaveAs(Form("hepdata/fig6-%d.png", i));

    delete tmp; 
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


