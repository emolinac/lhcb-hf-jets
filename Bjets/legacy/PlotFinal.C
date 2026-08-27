#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"
#include <TF1.h>
#include <TLatex.h>
#include <THStack.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include <TRatioPlot.h>
#include "../Helpers.h"
using namespace std;

void PlotFinal(int NumEvts = -1, int dataset = 91599, bool L0MuonDiMuon = false)
{
  bool WTA_cut = false;
  TString string_data, string_data_test, string_unfold, string_unfold_test, extension;
  TString str_followHard, str_ghost, str_Mag, str_flavor;
  TString loc_rootfiles_data("../../root_files/Bjets/");
  TString loc_rootfiles_MC("../../root_files/BjetsMC/");
  TString loc_plots("../../plots/Bjets/");

  bool ghostCut1(true), ghostCut2(true);
  int chargedJetCut1 = chargedJetCut;
  int chargedJetCut2 = chargedJetCut;
  bool followHardest1 = false;
  bool followHardest2 = false;

  int JetMeth1 = (dataset / 1000) % 10;
  int flavor1 = (dataset / 100) % 10;
  int ptRange1 = (dataset / 10) % 10;
  int Mag1 = (dataset / 1) % 10;


  double ghostProb_test = ghostProb;
  int NumIters = 4;


  TString str_Mag1 = "";
  TString str_pT1 = "";
  TString str_level1 = "";
  TString str_followHard1 = "";
  TString str_flavor1 = "";
  TString str_ghost1 = "";


  TString str_PID("");
  TString str_WTA("");

  if (PID_cut)
    str_PID = "_PID";
  if (WTA_cut)
    str_WTA = "_WTA";

  if (Mag1 == 0)
    str_Mag1 = "_MD";
  else if (Mag1 == 1)
    str_Mag1 = "_MU";

  if (flavor1 == 1)
    str_flavor1 = "_udsg";
  else if (flavor1 == 4)
    str_flavor1 = "_c";
  else if (flavor1 == 5)
    str_flavor1 = "_b";

  if (followHardest1)
    str_followHard1 = "_hard";
  else
    str_followHard1 = "_HF";

  if (ghostCut1)
    str_ghost1 = Form("_ghost_%.1f", ghostProb);


  TString str_L0 = "";
  if (L0MuonDiMuon)
    str_L0 = "_L0MuonDiMuon";

  string_data = loc_rootfiles_data + TString("data") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard1 + str_ghost1 + str_Mag1 + str_flavor1 + str_PID + str_WTA + str_L0 +  Form("_%d", dataset);
  string_unfold = loc_rootfiles_MC + TString("unfold_reco") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard1 + str_ghost1 + str_flavor1 + str_PID + str_WTA + str_L0 + Form("_%d", dataset);

  cout << string_data << endl;
  cout << string_unfold << endl;

  extension = TString("corrected_data") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard1 + str_ghost1 + str_Mag1 + str_flavor1 + str_PID + str_WTA + Form("_iters_%d", NumIters) + str_L0 + Form("_%d", dataset);

  /////////////////////   Get Files /////////////////////////////////

  TFile *file_data = new TFile(string_data + TString(".root"), "READ");
  TFile *file_unfold = new TFile(string_unfold + TString(".root"), "READ");
  TFile *file_save = new TFile(loc_rootfiles_data + extension + ".root", "RECREATE");


  /////////////////////   Create observable histograms /////////////////////////////////

  TH3D *h3_ptzjt = (TH3D *)file_data->Get("ptzjt");
  TH3D *h3_ptzjt_final = (TH3D *)h3_ptzjt->Clone("ptzjt_final");
  
  TH3D *h3_ptzr = (TH3D *)file_data->Get("ptzr"); 
  TH3D *h3_ptzr_final = (TH3D *)h3_ptzr->Clone("ptzr_final");
  
  TH3D *h3_ptjtr = (TH3D *)file_data->Get("ptjtr");
  TH3D *h3_ptjtr_final = (TH3D *)h3_ptjtr->Clone("ptjtr_final");   
  
  TH2D *h2_ptz = (TH2D *)file_data->Get("ptz");
  TH2D *h2_ptz_final = (TH2D *)h2_ptz->Clone("ptz_final");
  
  TH2D *h2_ptjt = (TH2D *)file_data->Get("ptjt");  
  TH2D *h2_ptjt_final = (TH2D *)h2_ptjt->Clone("ptjt_final");  
  
  TH2D *h2_ptr = (TH2D *)file_data->Get("ptr");
  TH2D *h2_ptr_final = (TH2D *)h2_ptr->Clone("ptr_final"); 


  /////////////////////   Get Efficiency Hists /////////////////////////////////

  TH3D *h3_eff_ptzjt = (TH3D *)file_unfold->Get("efficiency_ptzjt");
  TH3D *h3_purity_ptzjt = (TH3D *)file_unfold->Get("purity_ptzjt");
  
  TH3D *h3_eff_ptzr = (TH3D *)file_unfold->Get("efficiency_ptzr");
  TH3D *h3_purity_ptzr = (TH3D *)file_unfold->Get("purity_ptzr");     
  
  TH3D *h3_eff_ptjtr = (TH3D *)file_unfold->Get("efficiency_ptjtr");
  TH3D *h3_purity_ptjtr = (TH3D *)file_unfold->Get("purity_ptjtr");  

  TH2D *h2_eff_ptz = (TH2D *)file_unfold->Get("efficiency_ptz");
  TH2D *h2_purity_ptz = (TH2D *)file_unfold->Get("purity_ptz");
  
  TH2D *h2_eff_ptjt = (TH2D *)file_unfold->Get("efficiency_ptjt");
  TH2D *h2_purity_ptjt = (TH2D *)file_unfold->Get("purity_ptjt");     
  
  TH2D *h2_eff_ptr = (TH2D *)file_unfold->Get("efficiency_ptr");
  TH2D *h2_purity_ptr = (TH2D *)file_unfold->Get("purity_ptr");  

  /////////////////////   Get Response Matrices /////////////////////////////////

  RooUnfoldResponse *response_ptzjt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptzjt");
  RooUnfoldResponse *response_ptzr = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptzr");
  RooUnfoldResponse *response_ptjtr = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptjtr");    
  RooUnfoldResponse *response_ptz = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptz"); 
  RooUnfoldResponse *response_ptjt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptjt");  
  RooUnfoldResponse *response_ptr = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptr");  
 
  /////////////////////  Correct Purities /////////////////////////////////

  h3_ptzjt_final->Multiply(h3_ptzjt_final, h3_purity_ptzjt); 
  h3_ptzr_final->Multiply(h3_ptzr_final, h3_purity_ptzr); 
  h3_ptjtr_final->Multiply(h3_ptjtr_final, h3_purity_ptjtr);   
  h2_ptz_final->Multiply(h2_ptz_final, h2_purity_ptz); 
  h2_ptjt_final->Multiply(h2_ptjt_final, h2_purity_ptjt);  
  h2_ptr_final->Multiply(h2_ptr_final, h2_purity_ptr);     
  
  cout << "About to Unfold" << endl;

  // RooUnfoldids    unfold1D(response_jetpt, h1_jetpt_reco_purity, NumIters);

  /////////////////////   Create Unfold Objects and unfold /////////////////////////////////

  RooUnfoldBayes unfold_ptzjt(response_ptzjt, h3_ptzjt_final, NumIters);
  RooUnfoldBayes unfold_ptzr(response_ptzr, h3_ptzr_final, NumIters);
  RooUnfoldBayes unfold_ptjtr(response_ptjtr, h3_ptjtr_final, NumIters);  
  RooUnfoldBayes unfold_ptz(response_ptz, h2_ptz_final, NumIters); 
  RooUnfoldBayes unfold_ptjt(response_ptjt, h2_ptjt_final, NumIters);
  RooUnfoldBayes unfold_ptr(response_ptr, h2_ptr_final, NumIters);

  cout << "Unfold successful" << endl;
  /////////////////////   Get Unfolded Distribution /////////////////////////////////
  // h1_jetpt_reco_purity_unfold = (TH1D*) unfold1D.Hreco();

  h3_ptzjt_final = (TH3D *)unfold_ptzjt.Hreco();
  h3_ptzr_final = (TH3D *)unfold_ptzr.Hreco();
  h3_ptjtr_final = (TH3D *)unfold_ptjtr.Hreco();   
  h2_ptz_final = (TH2D *)unfold_ptz.Hreco();
  h2_ptjt_final = (TH2D *)unfold_ptjt.Hreco();
  h2_ptr_final = (TH2D *)unfold_ptr.Hreco();
 
  cout << "Got Unfolded Dist" << endl;
  /////////////////////  Correct for efficiency if needed /////////////////////////////////

  h3_ptzjt_final->Divide(h3_ptzjt_final, h3_eff_ptzjt);
  h3_ptzr_final->Divide(h3_ptzr_final, h3_eff_ptzr);
  h3_ptjtr_final->Divide(h3_ptjtr_final, h3_eff_ptjtr);
  h2_ptz_final->Divide(h2_ptz_final, h2_eff_ptz);
  h2_ptjt_final->Divide(h2_ptjt_final, h2_eff_ptjt);
  h2_ptr_final->Divide(h2_ptr_final, h2_eff_ptr);

  
  cout << "Corrected for Eff" << endl;

  /////////////////////  Combine last 2 z bins to include z=1.0 entries in the 0.95-1.0 bin /////////////////////////////////

  for (int iy = 1; iy <= h3_ptzjt_final->GetNbinsY(); iy++)
  {
    for (int iz = 1; iz <= h3_ptzjt_final->GetNbinsZ(); iz++)
    {
      h3_ptzjt_final->SetBinContent(zbinsize2D-1, iy, iz, h3_ptzjt_final->GetBinContent(zbinsize2D-1, iy, iz) + h3_ptzjt_final->GetBinContent(zbinsize2D, iy, iz));
      h3_ptzjt_final->SetBinError(zbinsize2D-1, iy, iz, sqrt(pow(h3_ptzjt_final->GetBinError(zbinsize2D-1, iy, iz),2) + pow(h3_ptzjt_final->GetBinError(zbinsize2D, iy, iz),2)));
      h3_ptzjt_final->SetBinContent(zbinsize2D, iy, iz, 0.);
      h3_ptzjt_final->SetBinError(zbinsize2D, iy, iz, 0.);
    }
  }
  
  for (int iy = 1; iy <= h3_ptzr_final->GetNbinsY(); iy++)
  {
    for (int iz = 1; iz <= h3_ptzr_final->GetNbinsZ(); iz++)
    {
      h3_ptzr_final->SetBinContent(zbinsize2D-1, iy, iz, h3_ptzr_final->GetBinContent(zbinsize2D-1, iy, iz) + h3_ptzr_final->GetBinContent(zbinsize2D, iy, iz));
      h3_ptzr_final->SetBinError(zbinsize2D-1, iy, iz, sqrt(pow(h3_ptzr_final->GetBinError(zbinsize2D-1, iy, iz),2) + pow(h3_ptzr_final->GetBinError(zbinsize2D, iy, iz),2)));
      h3_ptzr_final->SetBinContent(zbinsize2D, iy, iz, 0.);
      h3_ptzr_final->SetBinError(zbinsize2D, iy, iz, 0.);
    }
  }  
  
  for (int iy = 1; iy <= h2_ptz_final->GetNbinsY(); iy++)
  {
    h2_ptz_final->SetBinContent(zbinsize-1, iy, h2_ptz_final->GetBinContent(zbinsize-1, iy) + h2_ptz_final->GetBinContent(zbinsize, iy));
    h2_ptz_final->SetBinError(zbinsize-1, iy, sqrt(pow(h2_ptz_final->GetBinError(zbinsize-1, iy),2) + pow(h2_ptz_final->GetBinError(zbinsize, iy),2)));
    h2_ptz_final->SetBinContent(zbinsize, iy, 0.);
    h2_ptz_final->SetBinError(zbinsize, iy, 0.); 
  }
  
  /////////////////////   Write corrected histograms out for later processing ////////////////////
  
  int binlow = h2_ptz->GetYaxis()->FindBin(ptMin + 0.1);
  int binhigh = h2_ptz->GetYaxis()->FindBin(ptMax - 0.1); 

  h3_ptzjt_final->GetZaxis()->SetRange(binlow, binhigh); 
  h3_ptzr_final->GetZaxis()->SetRange(binlow, binhigh); 
  h3_ptjtr_final->GetZaxis()->SetRange(binlow, binhigh); 
  h2_ptz_final->GetYaxis()->SetRange(binlow, binhigh);  
  h2_ptjt_final->GetYaxis()->SetRange(binlow, binhigh);  
  h2_ptr_final->GetYaxis()->SetRange(binlow, binhigh);  
    
  h3_ptzjt_final->Write("ptzjt_final");
  h3_ptzr_final->Write("ptzr_final");
  h3_ptjtr_final->Write("ptjtr_final");
  h2_ptz_final->Write("ptz_final");
  h2_ptjt_final->Write("ptjt_final");
  h2_ptr_final->Write("ptr_final");
  
  /////////////////////   Normalize histograms /////////////////////////////////

  
  TH2D *h2_zjt_ptbinned_final[ptbinsize-1], *h2_zjt_ptbinned[ptbinsize-1];
  TH2D *h2_zr_ptbinned_final[ptbinsize-1], *h2_zr_ptbinned[ptbinsize-1];
  TH2D *h2_jtr_ptbinned_final[ptbinsize-1], *h2_jtr_ptbinned[ptbinsize-1];     
  TH1D *h1_z_ptbinned_final[ptbinsize-1], *h1_z_ptbinned[ptbinsize-1];
  TH1D *h1_jt_ptbinned_final[ptbinsize-1], *h1_jt_ptbinned[ptbinsize-1];
  TH1D *h1_r_ptbinned_final[ptbinsize-1], *h1_r_ptbinned[ptbinsize-1];    
    
  TCanvas *can_z[ptbinsize-1], *can_jt[ptbinsize-1], *can_r[ptbinsize-1];
  TCanvas *can_zjt[ptbinsize-1], *can_zr[ptbinsize-1], *can_jtr[ptbinsize-1];  
  TString plotdir = "../../plots/Bjets/data_unfold_compare/";
  for (int j=1; j < ptbinsize; ++j)
  {

    h3_ptzjt_final->GetZaxis()->SetRange(j+1, j+1);
    h3_ptzjt->GetZaxis()->SetRange(j+1, j+1); 

    h3_ptzr_final->GetZaxis()->SetRange(j+1, j+1);
    h3_ptzr->GetZaxis()->SetRange(j+1, j+1); 
    
    h3_ptjtr_final->GetZaxis()->SetRange(j+1, j+1);
    h3_ptjtr->GetZaxis()->SetRange(j+1, j+1); 
      
    h2_zjt_ptbinned_final[j-1] = (TH2D *)h3_ptzjt_final->Project3D("yx");
    h2_zjt_ptbinned[j-1] = (TH2D *)h3_ptzjt->Project3D("yx");     
    h2_zjt_ptbinned_final[j-1]->SetName(Form("zjt_final_pt%d",j));     
    h2_zjt_ptbinned[j-1]->SetName(Form("zjt_pt%d",j));
        
    h2_zr_ptbinned_final[j-1] = (TH2D *)h3_ptzr_final->Project3D("yx");
    h2_zr_ptbinned[j-1] = (TH2D *)h3_ptzr->Project3D("yx");  
    h2_zr_ptbinned_final[j-1]->SetName(Form("zr_final_pt%d",j));     
    h2_zr_ptbinned[j-1]->SetName(Form("zr_pt%d",j));         
    
    h2_jtr_ptbinned_final[j-1] = (TH2D *)h3_ptjtr_final->Project3D("yx");
    h2_jtr_ptbinned[j-1] = (TH2D *)h3_ptjtr->Project3D("yx"); 
    h2_jtr_ptbinned_final[j-1]->SetName(Form("jtr_final_pt%d",j));     
    h2_jtr_ptbinned[j-1]->SetName(Form("jtr_pt%d",j));       
    
    h1_z_ptbinned_final[j-1] = (TH1D *)h2_ptz_final->ProjectionX(Form("z_final_pt%d",j), j+1, j+1);
    h1_z_ptbinned[j-1] = (TH1D *)h2_ptz->ProjectionX(Form("z_pt%d",j), j+1, j+1);    
    
    h1_jt_ptbinned_final[j-1] = (TH1D *)h2_ptjt_final->ProjectionX(Form("jt_final_pt%d",j), j+1, j+1);
    h1_jt_ptbinned[j-1] = (TH1D *)h2_ptjt->ProjectionX(Form("jt_pt%d",j), j+1, j+1);       
    
    h1_r_ptbinned_final[j-1] = (TH1D *)h2_ptr_final->ProjectionX(Form("r_final_pt%d",j), j+1, j+1);
    h1_r_ptbinned[j-1] = (TH1D *)h2_ptr->ProjectionX(Form("r_pt%d",j), j+1, j+1);                  

    //h1_z_ptbinned_final[j-1]->GetXaxis()->SetRange(binlow, binhigh);
    //h1_z_ptbinned_final_test[j-1]->GetXaxis()->SetRange(binlow, binhigh);    

    //NormalizeHist(h1_z_ptbinned[j-1]_final);
    //NormalizeHist(h1_z_ptbinned[j-1]_final_test);

    NormalizeHist(h2_zjt_ptbinned_final[j-1]);    
    NormalizeHist(h2_zjt_ptbinned[j-1]);      
    
    NormalizeHist(h2_zr_ptbinned_final[j-1]);
    NormalizeHist(h2_zr_ptbinned[j-1]);    
    
    NormalizeHist(h2_jtr_ptbinned_final[j-1]);
    NormalizeHist(h2_jtr_ptbinned[j-1]);  
    
    NormalizeHist(h1_z_ptbinned_final[j-1]);    
    NormalizeHist(h1_z_ptbinned[j-1]);      
    
    NormalizeHist(h1_jt_ptbinned_final[j-1]);
    NormalizeHist(h1_jt_ptbinned[j-1]);    
    
    NormalizeHist(h1_r_ptbinned_final[j-1]);
    NormalizeHist(h1_r_ptbinned[j-1]);            
    
    gStyle->SetOptStat(0);
    
    TLegend *zleg = new TLegend(0.2, 0.5, 0.5, 0.8);
    zleg->SetLineColor(0);
    can_z[j-1] = new TCanvas(Form("can_z_pt%d", j));
    can_z[j-1]->cd();
    h1_z_ptbinned_final[j-1]->SetLineColor(kBlack);
    h1_z_ptbinned_final[j-1]->SetMarkerColor(kBlack);
    h1_z_ptbinned_final[j-1]->SetMarkerStyle(kFullCircle);
    h1_z_ptbinned[j-1]->SetLineColor(kBlue);
    h1_z_ptbinned[j-1]->SetMarkerColor(kBlue);
    h1_z_ptbinned[j-1]->SetMarkerStyle(kFullSquare);
    h1_z_ptbinned_final[j-1]->GetXaxis()->SetTitle("z");
    h1_z_ptbinned_final[j-1]->Draw();
    h1_z_ptbinned[j-1]->Draw("same");
    h1_z_ptbinned_final[j-1]->GetYaxis()->SetRangeUser(0.0, 3.5);
    zleg->AddEntry(h1_z_ptbinned_final[j-1], "corrected data", "lp");
    zleg->AddEntry(h1_z_ptbinned[j-1], "uncorrected data", "lp");
    zleg->Draw("same");
    TString plotname = Form("z_compare_pt%d", j);
    plotname = plotdir + plotname;
    plotname += str_L0;
    plotname += ".png";
    can_z[j-1]->SaveAs(plotname);
    
    TLegend *jtleg = new TLegend(0.6, 0.6, 0.85, 0.85);
    jtleg->SetLineColor(0);
    can_jt[j-1] = new TCanvas(Form("can_jt_pt%d", j));
    can_jt[j-1]->cd();
    can_jt[j-1]->SetLogy();
    h1_jt_ptbinned_final[j-1]->SetLineColor(kBlack);
    h1_jt_ptbinned_final[j-1]->SetMarkerColor(kBlack);
    h1_jt_ptbinned_final[j-1]->SetMarkerStyle(kFullCircle);
    h1_jt_ptbinned[j-1]->SetLineColor(kBlue);
    h1_jt_ptbinned[j-1]->SetMarkerColor(kBlue);
    h1_jt_ptbinned[j-1]->SetMarkerStyle(kFullSquare);
    h1_jt_ptbinned_final[j-1]->GetXaxis()->SetTitle("j_{T} [GeV/c]");
    h1_jt_ptbinned_final[j-1]->Draw();
    h1_jt_ptbinned[j-1]->Draw("same");
    jtleg->AddEntry(h1_z_ptbinned_final[j-1], "corrected data", "lp");
    jtleg->AddEntry(h1_z_ptbinned[j-1], "uncorrected data", "lp");
    jtleg->Draw("same");
    plotname = Form("jt_compare_pt%d", j);
    plotname = plotdir + plotname;
    plotname += str_L0;
    plotname += ".png";
    can_jt[j-1]->SaveAs(plotname);
    
    TLegend *rleg = new TLegend(0.6, 0.6, 0.85, 0.85);
    rleg->SetLineColor(0);
    can_r[j-1] = new TCanvas(Form("can_r_pt%d", j));
    can_r[j-1]->cd();
    can_r[j-1]->SetLogy();
    h1_r_ptbinned_final[j-1]->SetLineColor(kBlack);
    h1_r_ptbinned_final[j-1]->SetMarkerColor(kBlack);
    h1_r_ptbinned_final[j-1]->SetMarkerStyle(kFullCircle);
    h1_r_ptbinned[j-1]->SetLineColor(kBlue);
    h1_r_ptbinned[j-1]->SetMarkerColor(kBlue);
    h1_r_ptbinned[j-1]->SetMarkerStyle(kFullSquare);
    h1_r_ptbinned_final[j-1]->GetXaxis()->SetTitle("r");
    h1_r_ptbinned_final[j-1]->Draw();
    h1_r_ptbinned[j-1]->Draw("same");
    rleg->AddEntry(h1_z_ptbinned_final[j-1], "corrected data", "lp");
    rleg->AddEntry(h1_z_ptbinned[j-1], "uncorrected data", "lp");
    rleg->Draw("same");
    plotname = Form("r_compare_pt%d", j);
    plotname = plotdir + plotname;
    plotname += str_L0;
    plotname += ".png";
    can_r[j-1]->SaveAs(plotname);
                     
  }
  
/*
  /////////////////////   Create pull hists /////////////////////////////////

  TH2D *h2_ptktdR_pull = (TH2D *)h2_ptktdR_final->Clone("ptktdR_pull");
  TH1D *h1_ptktdR_pulldist = new TH1D("ptktdR_pulldist", "", 16, -4.5, 4.5);

  TH2D *h2_ptzdR_pull = (TH2D *)h2_ptzdR_final->Clone("ptzdR_pull");
  TH1D *h1_ptzdR_pulldist = new TH1D("ptzdR_pulldist", "", 16, -4.5, 4.5);

  TH2D *h2_ptthetaErad_pull = (TH2D *)h2_ptthetaErad_final->Clone("ptthetaErad_pull");
  TH1D *h1_ptthetaErad_pulldist = new TH1D("ptthetaErad_pulldist", "", 16, -4.5, 4.5);
*/


  //GetPulls(h2_ptktdR_diff, h2_ptktdR_pull, h1_ptktdR_pulldist);
/*
  //---- paint setup...
  //
  int ican = -1, iframe = -1, itext = -1;
  TCanvas *ccan[1000];
  TH1F *frame[1000];
  TLatex *text[1000];
  TRatioPlot *rp[1000];
  for (int i = 0; i < 1000; i++)
  {
    text[i] = new TLatex();
    text[i]->SetNDC(kTRUE);
    text[i]->SetTextSize(0.06);
  }
  TLatex Tl;
  Tl.SetNDC(kTRUE);
  Tl.SetTextSize(0.04);
  //
  gStyle->SetOptStat(0);
  // gStyle->SetPaperSize(TStyle::kUSLetter);
  // gStyle->SetPadBottomMargin(0.08);
  // gStyle->SetPadTopMargin(0.005);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.13);
  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetLabelSize(0.05, "Y");
  gStyle->SetTitleXSize(0.055);
  gStyle->SetTitleYSize(0.055);
  gStyle->SetTitleOffset(0.55, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetStatW(0.2);
  // gStyle->SetPalette(kBird);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(100);
  gStyle->SetPaintTextFormat("4.5g");
  // gStyle->SetErrorX(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetStatStyle(0);
  // gStyle->SetLineWidth(3);

  //---- paint...
  char buf[100];
  char bufb[100];
  TString rootfile;
  TString plotfile;
  TString plotfilePDF;
  TString plotfileO;
  TString plotfileC;
  // TString OutputFileBase	= outbase+outinfo;
  rootfile = loc_hists + extension + TString(".root");
  plotfile = loc_plots + extension + TString(".ps");
  plotfilePDF = loc_plots + extension + TString(".pdf");
  plotfileO = plotfilePDF + TString("(");
  plotfileC = plotfilePDF + TString("]");
  // c->SaveAs("plots/"+extension+".pdf");

  // Begin plotting
  //

  //

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
  ccan[ican]->SetFillColor(10);
  // gPad->SetLeftMargin(0.16);
  // gPad->SetBottomMargin(0.06);
  ccan[ican]->cd();
  ccan[ican]->Divide(2, 2, 0.0001, 0.0001);

  ccan[ican]->cd(1);
  h2_ptktdR_final_test->Draw("COL0Z");

  ccan[ican]->cd(2);
  h1_ptktdR_pulldist->Draw("E SAME");
  h1_ptktdR_pulldist->Draw("HIST SAME");
  ccan[ican]->cd(3);
  // h2_ktdR_diff->Set
  h2_ptktdR_ratio->Draw("COL0Z TEXTE MIN0");
  h2_ptktdR_ratio->SetMinimum(-0.1);
  h2_ptktdR_ratio->SetMaximum(0.1);

  ccan[ican]->cd(4);
  h2_ptktdR_pull->Draw("COL0Z TEXT MIN0");
  h2_ptktdR_pull->SetMaximum(4.);
  h2_ptktdR_pull->SetMinimum(-4.);

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
  //

  //

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
  ccan[ican]->SetFillColor(10);
  // gPad->SetLeftMargin(0.16);
  // gPad->SetBottomMargin(0.06);
  ccan[ican]->cd();
  ccan[ican]->Divide(2, 2, 0.0001, 0.0001);

  ccan[ican]->cd(1);
  h2_ptzdR_final_test->Draw("COL0Z");

  ccan[ican]->cd(2);
  h1_ptzdR_pulldist->Draw("E SAME");
  h1_ptzdR_pulldist->Draw("HIST SAME");
  // h2_ktdR_->SetTitle("Data+Purity+Unfold; ;");
  ccan[ican]->cd(3);
  // h2_ktdR_diff->Set
  h2_ptzdR_ratio->Draw("COL0Z TEXTE MIN0");
  h2_ptzdR_ratio->SetMinimum(-0.1);
  h2_ptzdR_ratio->SetMaximum(0.1);

  ccan[ican]->cd(4);
  h2_ptzdR_pull->Draw("COL0Z TEXT MIN0");
  h2_ptzdR_pull->SetMaximum(4.);
  h2_ptzdR_pull->SetMinimum(-4.);

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
  //

  if (ican > -1)
  {
    cout << " You plotted " << ican + 1 << " canvasses......." << endl;
    ccan[ican]->Print(plotfileC.Data());
  }
  //
*/
  //file_save->Write();
  //file_save->Close();
}
