#include <TFile.h>
#include <TH2.h>
#include "Settings.h"
#include "../Helpers.h"

void SPlotCompare()
{


    gStyle->SetPaintTextFormat("3.3f");
    gStyle->SetOptStat(0);
    gROOT->ForceStyle();
    
    TFile infile = TFile("../../root_files/Bjets/data_ev_-1_ptj_5100_eta_2.54.0_HF_ghost_0.4_b_PID_splotfit_91599.root");
    //TFile infile = TFile("../../root_files/Bjets/data_ev_-1_ptj_10100_eta_2.54.0_HF_ghost_0.4_b_PID_91599.root");
    TFile outfile = TFile("../../root_files/Bjets/splot_compare.root", "RECREATE");
    

    TH2D *ptz = (TH2D*)infile.Get("ptz");
    NormalizeHist(ptz);
    TH2D *ptz_sweight = (TH2D*)infile.Get("ptz_sweight");
    NormalizeHist(ptz_sweight);
    TH2D *ptz_uncorrected = (TH2D*)infile.Get("ptz_uncorrected");
    NormalizeHist(ptz_uncorrected);
    TH2D *ptz_uncorrected_nomasscond = (TH2D*)infile.Get("ptz_uncorrected_nomasscond");
    NormalizeHist(ptz_uncorrected_nomasscond);
    TH2D *ptz_nobgsub = (TH2D*)infile.Get("ptz_nobgsub");
    NormalizeHist(ptz_nobgsub);

    TH2D *ptz_ratio = (TH2D*)ptz->Clone("ptz_ratio");
    TH2D *ptz_ratio_uncorr_sweight = (TH2D*)ptz_uncorrected_nomasscond->Clone("ptz_ratio_uncorr_sweight");
    TH2D *ptz_ratio_uncorr_sbsub = (TH2D*)ptz_uncorrected->Clone("ptz_ratio_uncorr_sbsub");
    TH2D *ptz_ratio_evtweight = (TH2D*)ptz_uncorrected->Clone("ptz_ratio_evtweight");

    ptz_ratio->Divide(ptz, ptz_sweight);
    ptz_ratio_evtweight->Divide(ptz_uncorrected, ptz_nobgsub);
    ptz_ratio_uncorr_sweight->Divide(ptz_uncorrected_nomasscond, ptz_sweight);
    ptz_ratio_uncorr_sbsub->Divide(ptz_uncorrected, ptz);

    ptz_ratio->SetOption("COLZ, text");
    ptz_ratio_evtweight->SetOption("COLZ, text");
    ptz_ratio_uncorr_sweight->SetOption("COLZ, text");
    ptz_ratio_uncorr_sbsub->SetOption("COLZ, text");
    ptz_ratio->Write();
    ptz_ratio_evtweight->Write();
    ptz_ratio_uncorr_sweight->Write();
    ptz_ratio_uncorr_sbsub->Write();

    THStack *hs_ptz_ratio = new THStack("z_sweight_ratio", ";z;#frac{SB subtraction}{sPlot}");      
    THStack *hs_ptz_ratio_evtweight = new THStack("z_ratio_evtweight", ";z;#frac{Uncorrected}{event weighted}");      
    THStack *hs_ptz_ratio_uncorr_sweight = new THStack("z_ratio_uncorr_sweight", ";z;#frac{Uncorrected}{sPlot}");   
    THStack *hs_ptz_ratio_uncorr_sbsub = new THStack("z_ratio_uncorr_sbsub", ";z;#frac{Uncorrected}{SB subtraction}");      
   

    for (int j = 1; j < ptbinsize; j++)
    {   
        
      TH1D *h1_ptz_temp = (TH1D *)ptz_ratio->ProjectionX(Form("z_ratio_pt%d", j), j + 1, j + 1); 
      TH1D *h1_ptz_temp_evtweight = (TH1D *)ptz_ratio_evtweight->ProjectionX(Form("z_ratio_evtweight_pt%d", j), j + 1, j + 1);
      TH1D *h1_ptz_temp_uncorr_sweight = (TH1D *)ptz_ratio_uncorr_sweight->ProjectionX(Form("z_ratio_uncorr_sweight_pt%d", j), j + 1, j + 1); 
      TH1D *h1_ptz_temp_uncorr_sbsub = (TH1D *)ptz_ratio_uncorr_sbsub->ProjectionX(Form("z_ratio_uncorr_sbsub_pt%d", j), j + 1, j + 1); 

      h1_ptz_temp->SetStats(0);
      h1_ptz_temp_evtweight->SetStats(0);
      h1_ptz_temp_uncorr_sweight->SetStats(0);
      h1_ptz_temp_uncorr_sbsub->SetStats(0);

      h1_ptz_temp->SetMarkerStyle(j + 20);
      h1_ptz_temp_evtweight->SetMarkerStyle(j + 20);
      h1_ptz_temp_uncorr_sweight->SetMarkerStyle(j + 20);
      h1_ptz_temp_uncorr_sbsub->SetMarkerStyle(j + 20);

      if (j!=5) 
      {
        h1_ptz_temp->SetMarkerColor(j);
        h1_ptz_temp->SetLineColor(j);
        h1_ptz_temp_evtweight->SetMarkerColor(j);
        h1_ptz_temp_evtweight->SetLineColor(j);
        h1_ptz_temp_uncorr_sweight->SetMarkerColor(j);
        h1_ptz_temp_uncorr_sweight->SetLineColor(j);
        h1_ptz_temp_uncorr_sbsub->SetMarkerColor(j);
        h1_ptz_temp_uncorr_sbsub->SetLineColor(j);
      }
      else
      {
        h1_ptz_temp->SetMarkerColor(j*j+3);
        h1_ptz_temp->SetLineColor(j*j+3);
        h1_ptz_temp_evtweight->SetMarkerColor(j*j+3);
        h1_ptz_temp_evtweight->SetLineColor(j*j+3);
        h1_ptz_temp_uncorr_sweight->SetMarkerColor(j*j+3);
        h1_ptz_temp_uncorr_sweight->SetLineColor(j*j+3);
        h1_ptz_temp_uncorr_sbsub->SetMarkerColor(j*j+3);
        h1_ptz_temp_uncorr_sbsub->SetLineColor(j*j+3);
      }
      h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptz_temp_evtweight->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptz_temp_uncorr_sweight->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptz_temp_uncorr_sbsub->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));

      //h1_ptz_temp->SetOption("PE HIST");
      hs_ptz_ratio->Add(h1_ptz_temp);
      hs_ptz_ratio_evtweight->Add(h1_ptz_temp_evtweight);
      hs_ptz_ratio_uncorr_sweight->Add(h1_ptz_temp_uncorr_sweight);
      hs_ptz_ratio_uncorr_sbsub->Add(h1_ptz_temp_uncorr_sbsub);
    }
    hs_ptz_ratio->SetMaximum(2.0);
    hs_ptz_ratio->SetMinimum(0.);
    hs_ptz_ratio_evtweight->SetMaximum(2.0);
    hs_ptz_ratio_evtweight->SetMinimum(0.);
    hs_ptz_ratio_uncorr_sweight->SetMaximum(2.0);
    hs_ptz_ratio_uncorr_sweight->SetMinimum(0.);
    hs_ptz_ratio_uncorr_sbsub->SetMaximum(2.0);
    hs_ptz_ratio_uncorr_sbsub->SetMinimum(0.);
    hs_ptz_ratio->Write();
    hs_ptz_ratio_evtweight->Write();
    hs_ptz_ratio_uncorr_sweight->Write();
    hs_ptz_ratio_uncorr_sbsub->Write();

    // jT hists //
    TH2D *ptjt = (TH2D*)infile.Get("ptjt");
    NormalizeHist(ptjt);
    TH2D *ptjt_sweight = (TH2D*)infile.Get("ptjt_sweight");
    NormalizeHist(ptjt_sweight);
    TH2D *ptjt_uncorrected = (TH2D*)infile.Get("ptjt_uncorrected");
    NormalizeHist(ptjt_uncorrected);
    TH2D *ptjt_uncorrected_nomasscond = (TH2D*)infile.Get("ptjt_uncorrected_nomasscond");
    NormalizeHist(ptjt_uncorrected_nomasscond);
    TH2D *ptjt_nobgsub = (TH2D*)infile.Get("ptjt_nobgsub");
    NormalizeHist(ptjt_nobgsub);

    TH2D *ptjt_ratio = (TH2D*)ptjt->Clone("ptjt_ratio");
    TH2D *ptjt_ratio_uncorr_sweight = (TH2D*)ptjt_uncorrected_nomasscond->Clone("ptjt_ratio_uncorr_sweight");
    TH2D *ptjt_ratio_uncorr_sbsub = (TH2D*)ptjt_uncorrected->Clone("ptjt_ratio_uncorr_sbsub");
    TH2D *ptjt_ratio_evtweight = (TH2D*)ptjt_uncorrected->Clone("ptjt_ratio_evtweight");

    ptjt_ratio->Divide(ptjt, ptjt_sweight);
    ptjt_ratio_evtweight->Divide(ptjt_uncorrected, ptjt_nobgsub);
    ptjt_ratio_uncorr_sweight->Divide(ptjt_uncorrected_nomasscond, ptjt_sweight);
    ptjt_ratio_uncorr_sbsub->Divide(ptjt_uncorrected, ptjt);

    ptjt_ratio->SetOption("COLZ, text");
    ptjt_ratio_evtweight->SetOption("COLZ, text");
    ptjt_ratio_uncorr_sweight->SetOption("COLZ, text");
    ptjt_ratio_uncorr_sbsub->SetOption("COLZ, text");
    ptjt_ratio->Write();
    ptjt_ratio_evtweight->Write();
    ptjt_ratio_uncorr_sweight->Write();
    ptjt_ratio_uncorr_sbsub->Write();

    THStack *hs_ptjt_ratio = new THStack("jt_sweight_ratio", ";j_{T};#frac{SB subtraction}{sPlot}");      
    THStack *hs_ptjt_ratio_evtweight = new THStack("jt_ratio_evtweight", ";j_{T};#frac{Uncorrected}{event weighted}");      
    THStack *hs_ptjt_ratio_uncorr_sweight = new THStack("jt_ratio_uncorr_sweight", ";j_{T};#frac{Uncorrected}{sPlot}");   
    THStack *hs_ptjt_ratio_uncorr_sbsub = new THStack("jt_ratio_uncorr_sbsub", ";j_{T};#frac{Uncorrected}{SB subtraction}");      
   

    for (int j = 1; j < ptbinsize; j++)
    {   
        
      TH1D *h1_ptjt_temp = (TH1D *)ptjt_ratio->ProjectionX(Form("jt_ratio_pt%d", j), j + 1, j + 1); 
      TH1D *h1_ptjt_temp_evtweight = (TH1D *)ptjt_ratio_evtweight->ProjectionX(Form("jt_ratio_evtweight_pt%d", j), j + 1, j + 1);
      TH1D *h1_ptjt_temp_uncorr_sweight = (TH1D *)ptjt_ratio_uncorr_sweight->ProjectionX(Form("jt_ratio_uncorr_sweight_pt%d", j), j + 1, j + 1); 
      TH1D *h1_ptjt_temp_uncorr_sbsub = (TH1D *)ptjt_ratio_uncorr_sbsub->ProjectionX(Form("jt_ratio_uncorr_sbsub_pt%d", j), j + 1, j + 1); 

      h1_ptjt_temp->SetStats(0);
      h1_ptjt_temp_evtweight->SetStats(0);
      h1_ptjt_temp_uncorr_sweight->SetStats(0);
      h1_ptjt_temp_uncorr_sbsub->SetStats(0);

      h1_ptjt_temp->SetMarkerStyle(j + 20);
      h1_ptjt_temp_evtweight->SetMarkerStyle(j + 20);
      h1_ptjt_temp_uncorr_sweight->SetMarkerStyle(j + 20);
      h1_ptjt_temp_uncorr_sbsub->SetMarkerStyle(j + 20);

      if (j!=5) 
      {
        h1_ptjt_temp->SetMarkerColor(j);
        h1_ptjt_temp->SetLineColor(j);
        h1_ptjt_temp_evtweight->SetMarkerColor(j);
        h1_ptjt_temp_evtweight->SetLineColor(j);
        h1_ptjt_temp_uncorr_sweight->SetMarkerColor(j);
        h1_ptjt_temp_uncorr_sweight->SetLineColor(j);
        h1_ptjt_temp_uncorr_sbsub->SetMarkerColor(j);
        h1_ptjt_temp_uncorr_sbsub->SetLineColor(j);
      }
      else
      {
        h1_ptjt_temp->SetMarkerColor(j*j+3);
        h1_ptjt_temp->SetLineColor(j*j+3);
        h1_ptjt_temp_evtweight->SetMarkerColor(j*j+3);
        h1_ptjt_temp_evtweight->SetLineColor(j*j+3);
        h1_ptjt_temp_uncorr_sweight->SetMarkerColor(j*j+3);
        h1_ptjt_temp_uncorr_sweight->SetLineColor(j*j+3);
        h1_ptjt_temp_uncorr_sbsub->SetMarkerColor(j*j+3);
        h1_ptjt_temp_uncorr_sbsub->SetLineColor(j*j+3);
      }
      h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptjt_temp_evtweight->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptjt_temp_uncorr_sweight->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptjt_temp_uncorr_sbsub->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));

      //h1_ptjt_temp->SetOption("PE HIST");
      hs_ptjt_ratio->Add(h1_ptjt_temp);
      hs_ptjt_ratio_evtweight->Add(h1_ptjt_temp_evtweight);
      hs_ptjt_ratio_uncorr_sweight->Add(h1_ptjt_temp_uncorr_sweight);
      hs_ptjt_ratio_uncorr_sbsub->Add(h1_ptjt_temp_uncorr_sbsub);
    }
    hs_ptjt_ratio->SetMaximum(2.0);
    hs_ptjt_ratio->SetMinimum(0.);
    hs_ptjt_ratio_evtweight->SetMaximum(2.0);
    hs_ptjt_ratio_evtweight->SetMinimum(0.);
    hs_ptjt_ratio_uncorr_sweight->SetMaximum(2.0);
    hs_ptjt_ratio_uncorr_sweight->SetMinimum(0.);
    hs_ptjt_ratio_uncorr_sbsub->SetMaximum(2.0);
    hs_ptjt_ratio_uncorr_sbsub->SetMinimum(0.);
    hs_ptjt_ratio->Write();
    hs_ptjt_ratio_evtweight->Write();
    hs_ptjt_ratio_uncorr_sweight->Write();
    hs_ptjt_ratio_uncorr_sbsub->Write();


    // r hists //
    TH2D *ptr = (TH2D*)infile.Get("ptr");
    NormalizeHist(ptr);
    TH2D *ptr_sweight = (TH2D*)infile.Get("ptr_sweight");
    NormalizeHist(ptr_sweight);
    TH2D *ptr_uncorrected = (TH2D*)infile.Get("ptr_uncorrected");
    NormalizeHist(ptr_uncorrected);
    TH2D *ptr_uncorrected_nomasscond = (TH2D*)infile.Get("ptr_uncorrected_nomasscond");
    NormalizeHist(ptr_uncorrected_nomasscond);
    TH2D *ptr_nobgsub = (TH2D*)infile.Get("ptr_nobgsub");
    NormalizeHist(ptr_nobgsub);

    TH2D *ptr_ratio = (TH2D*)ptr->Clone("ptr_ratio");
    TH2D *ptr_ratio_uncorr_sweight = (TH2D*)ptr_uncorrected_nomasscond->Clone("ptr_ratio_uncorr_sweight");
    TH2D *ptr_ratio_uncorr_sbsub = (TH2D*)ptr_uncorrected->Clone("ptr_ratio_uncorr_sbsub");
    TH2D *ptr_ratio_evtweight = (TH2D*)ptr_uncorrected->Clone("ptr_ratio_evtweight");

    ptr_ratio->Divide(ptr, ptr_sweight);
    ptr_ratio_evtweight->Divide(ptr_uncorrected, ptr_nobgsub);
    ptr_ratio_uncorr_sweight->Divide(ptr_uncorrected_nomasscond, ptr_sweight);
    ptr_ratio_uncorr_sbsub->Divide(ptr_uncorrected, ptr);

    ptr_ratio->SetOption("COLZ, text");
    ptr_ratio_evtweight->SetOption("COLZ, text");
    ptr_ratio_uncorr_sweight->SetOption("COLZ, text");
    ptr_ratio_uncorr_sbsub->SetOption("COLZ, text");
    ptr_ratio->Write();
    ptr_ratio_evtweight->Write();
    ptr_ratio_uncorr_sweight->Write();
    ptr_ratio_uncorr_sbsub->Write();

    THStack *hs_ptr_ratio = new THStack("r_sweight_ratio", ";r;#frac{SB subtraction}{sPlot}");      
    THStack *hs_ptr_ratio_evtweight = new THStack("r_ratio_evtweight", ";r;#frac{Uncorrected}{event weighted}");      
    THStack *hs_ptr_ratio_uncorr_sweight = new THStack("r_ratio_uncorr_sweight", ";r;#frac{Uncorrected}{sPlot}");   
    THStack *hs_ptr_ratio_uncorr_sbsub = new THStack("r_ratio_uncorr_sbsub", ";r;#frac{Uncorrected}{SB subtraction}");      
   

    for (int j = 1; j < ptbinsize; j++)
    {   
        
      TH1D *h1_ptr_temp = (TH1D *)ptr_ratio->ProjectionX(Form("r_ratio_pt%d", j), j + 1, j + 1); 
      TH1D *h1_ptr_temp_evtweight = (TH1D *)ptr_ratio_evtweight->ProjectionX(Form("r_ratio_evtweight_pt%d", j), j + 1, j + 1);
      TH1D *h1_ptr_temp_uncorr_sweight = (TH1D *)ptr_ratio_uncorr_sweight->ProjectionX(Form("r_ratio_uncorr_sweight_pt%d", j), j + 1, j + 1); 
      TH1D *h1_ptr_temp_uncorr_sbsub = (TH1D *)ptr_ratio_uncorr_sbsub->ProjectionX(Form("r_ratio_uncorr_sbsub_pt%d", j), j + 1, j + 1); 

      h1_ptr_temp->SetStats(0);
      h1_ptr_temp_evtweight->SetStats(0);
      h1_ptr_temp_uncorr_sweight->SetStats(0);
      h1_ptr_temp_uncorr_sbsub->SetStats(0);

      h1_ptr_temp->SetMarkerStyle(j + 20);
      h1_ptr_temp_evtweight->SetMarkerStyle(j + 20);
      h1_ptr_temp_uncorr_sweight->SetMarkerStyle(j + 20);
      h1_ptr_temp_uncorr_sbsub->SetMarkerStyle(j + 20);

      if (j!=5) 
      {
        h1_ptr_temp->SetMarkerColor(j);
        h1_ptr_temp->SetLineColor(j);
        h1_ptr_temp_evtweight->SetMarkerColor(j);
        h1_ptr_temp_evtweight->SetLineColor(j);
        h1_ptr_temp_uncorr_sweight->SetMarkerColor(j);
        h1_ptr_temp_uncorr_sweight->SetLineColor(j);
        h1_ptr_temp_uncorr_sbsub->SetMarkerColor(j);
        h1_ptr_temp_uncorr_sbsub->SetLineColor(j);
      }
      else
      {
        h1_ptr_temp->SetMarkerColor(j*j+3);
        h1_ptr_temp->SetLineColor(j*j+3);
        h1_ptr_temp_evtweight->SetMarkerColor(j*j+3);
        h1_ptr_temp_evtweight->SetLineColor(j*j+3);
        h1_ptr_temp_uncorr_sweight->SetMarkerColor(j*j+3);
        h1_ptr_temp_uncorr_sweight->SetLineColor(j*j+3);
        h1_ptr_temp_uncorr_sbsub->SetMarkerColor(j*j+3);
        h1_ptr_temp_uncorr_sbsub->SetLineColor(j*j+3);
      }
      h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptr_temp_evtweight->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptr_temp_uncorr_sweight->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
      h1_ptr_temp_uncorr_sbsub->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));

      //h1_ptr_temp->SetOption("PE HIST");
      hs_ptr_ratio->Add(h1_ptr_temp);
      hs_ptr_ratio_evtweight->Add(h1_ptr_temp_evtweight);
      hs_ptr_ratio_uncorr_sweight->Add(h1_ptr_temp_uncorr_sweight);
      hs_ptr_ratio_uncorr_sbsub->Add(h1_ptr_temp_uncorr_sbsub);
    }
    hs_ptr_ratio->SetMaximum(2.0);
    hs_ptr_ratio->SetMinimum(0.);
    hs_ptr_ratio_evtweight->SetMaximum(2.0);
    hs_ptr_ratio_evtweight->SetMinimum(0.);
    hs_ptr_ratio_uncorr_sweight->SetMaximum(2.0);
    hs_ptr_ratio_uncorr_sweight->SetMinimum(0.);
    hs_ptr_ratio_uncorr_sbsub->SetMaximum(2.0);
    hs_ptr_ratio_uncorr_sbsub->SetMinimum(0.);
    hs_ptr_ratio->Write();
    hs_ptr_ratio_evtweight->Write();
    hs_ptr_ratio_uncorr_sweight->Write();
    hs_ptr_ratio_uncorr_sbsub->Write();

    outfile.Close();
}
