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
#include "../include/TBJetsMC.h"
#include "../include/TBJetsMC.C"

using namespace std;

void ClosureTest(int NumEvts = -1,
                 int dataset1 = 91599,
                 int dataset2 = 91599,
                 int NumIters = 4,  
                 bool DoShapeClosure = false)
{
        gStyle->SetPaintTextFormat("3.3f");
        gROOT->ForceStyle();
        gStyle->SetOptStat(0);

        TRandom3 *myRNG = new TRandom3(); // Set initializing value to 0 in order to set a random seed each time this is executed
        const int nRuns =  DoShapeClosure ? 1 : 10;

        bool smear_by_data = (!DoShapeClosure);  
        
        const int fixSmear = 42;

        // int NumIters = 1;
        // int RegIDS = 5;

        /////////////////////   Get Files /////////////////////////////////

        TFile *file_reco = new TFile((output_folder + "bjets_simpleobservable_mcreco.root").c_str(), "READ"); 
        TFile *file_data = new TFile((output_folder + "bjets_simpleobservable_data.root").c_str(), "READ");
        TFile *file_truth = new TFile((output_folder + "bjets_simpleobservable_mc.root").c_str(), "READ"); 
        TFile *file_unfold = new TFile((output_folder + "bjets_corrections.root").c_str(), "READ"); 

        TFile *file_write = new TFile((output_folder + "bjets_closuretest.root").c_str(), "RECREATE");

        /////////////////////   Get histograms /////////////////////////////////

        TH1D *h1_jetpt_reco = (TH1D *)file_reco->Get("Jet_pT");
        TH1D *h1_jetpt_data = (TH1D *)file_data->Get("Jet_pT");  
        h1_jetpt_reco->Draw();
        TH1D *h1_jetpt_final = (TH1D *)h1_jetpt_reco->Clone("jetpt_final");

        /////////////////////   Get Truth histograms /////////////////////////////////
        TH1D *h1_jetpt_truth = (TH1D *)file_truth->Get("Jet_pT");

        /////////////////////   Get Purity & Efficiency Hists /////////////////////////////////
        TH1D *h1_purity_jetpt = (TH1D *)file_unfold->Get("purity_jetpt");
        TH1D *h1_efficiency_jetpt = (TH1D *)file_unfold->Get("efficiency_jetpt");

        /////////////////////   Get Response Matrices /////////////////////////////////
        RooUnfoldResponse *response_jetpt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_jetpt");

        TH2 *h2_response_jetpt = (TH2 *)response_jetpt->Hresponse();

        response_jetpt->UseOverflow();

        // Multiply by purity
        h1_jetpt_final->Multiply(h1_purity_jetpt);

        // Unfold
        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_final, NumIters);
        
        h1_jetpt_final = (TH1D *)unfold_jetpt.Hreco();
        
        // Divide by efficiency
        h1_jetpt_final->Divide(h1_efficiency_jetpt);

        TH1D *h1_jetpt_final_ratio = (TH1D *)h1_jetpt_final->Clone("h1_jetpt_final_ratio");
        
        h1_jetpt_final_ratio->Divide(h1_jetpt_truth);

        h1_jetpt_reco->Write("h1_jetpt_reco");
        h1_jetpt_truth->Write("h1_jetpt_truth");
        h1_jetpt_final->Write("h1_jetpt_final");  

        int binlow_jet  = h1_jetpt_truth->FindBin(ptMin);
        int binhigh_jet = h1_jetpt_truth->FindBin(ptMax);

        double Njets_reco  = h1_jetpt_reco->Integral(binlow_jet, binhigh_jet);
        double Njets_truth = h1_jetpt_truth->Integral(binlow_jet, binhigh_jet);
        double Njets_final = h1_jetpt_final->Integral(binlow_jet, binhigh_jet);
                
        ////////////////////////////////////
        // Smearing the jet pt distribution
        ///////////////////////////////////
        TH1D *h1_jetpt_closure_error;
        for (int i = 0; i < nRuns; i++) {
                TH1D *h1_jetpt_smear = (TH1D *)h1_jetpt_reco->Clone(Form("jetpt_smear%d", i));
                //myRNG->SetSeed();
                //myRNG->SetSeed(fixSmear);
                if (smear_by_data)
                        SmearObservables(h1_jetpt_smear, h1_jetpt_data, myRNG);
                
                // Multiply by purity
                h1_jetpt_smear->Multiply(h1_purity_jetpt);
                
                RooUnfoldBayes unfold_jetpt_smear(response_jetpt, h1_jetpt_smear, NumIters);      
                
                h1_jetpt_smear = (TH1D *)unfold_jetpt_smear.Hreco();

                // Divide by efficiency
                h1_jetpt_smear->Divide(h1_efficiency_jetpt);
                
                TH1D *h1_jetpt_ratio_smear = (TH1D *)h1_jetpt_reco->Clone(Form("h1_jetpt_ratio_smear%d", i));
                //TH1D *h1_jetpt_pull_smear = (TH1D *)h1_jetpt_reco->Clone(Form("h1_jetpt_pull_smear%d", i));
                h1_jetpt_ratio_smear->Divide(h1_jetpt_smear, h1_jetpt_truth);
                h1_jetpt_ratio_smear->Write();
                //h1_jetpt_pull_smear->Write(); 
        }


        /////////////////////   Create z histograms /////////////////////////////////
        
        std::cout << "############################## Unfolding for ptz ##############################" << std::endl;

        TH2D *h2_ptz_truth  = (TH2D *)file_truth->Get("ptz");
        TH2D *h2_ptz        = (TH2D *)file_reco->Get("ptz");
        TH2D *h2_ptz_data   = (TH2D *)file_data->Get("ptz");  
        TH2D *h2_ptz_final  = (TH2D *)h2_ptz->Clone("ptz_final");
        TH2D *h2_eff_ptz    = (TH2D *)file_unfold->Get("efficiency_ptz");
        TH2D *h2_purity_ptz = (TH2D *)file_unfold->Get("purity_ptz");
        
        RooUnfoldResponse *response_ptz = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptz");
        
        // Multiply by purity  
        h2_ptz_final->Multiply(h2_ptz_final, h2_purity_ptz);
        RooUnfoldBayes unfold_ptz(response_ptz, h2_ptz_final, NumIters);

        h2_ptz_final = (TH2D *)unfold_ptz.Hreco();
        // Divide by efficiency  
        h2_ptz_final->Divide(h2_ptz_final, h2_eff_ptz);
        
        int binlow = h2_ptz->GetYaxis()->FindBin(ptMin + 0.1);
        int binhigh = h2_ptz->GetYaxis()->FindBin(ptMax - 0.1);
        int binlowz = h2_ptz->GetXaxis()->FindBin(z_binedges[0]);
        int binhighz = h2_ptz->GetXaxis()->FindBin(z_binedges[zbinsize-1]);
        //int binhighz = h2_ptz->GetXaxis()->FindBin(0.95);  

        h2_ptz_final->GetYaxis()->SetRangeUser(ptMin, ptMax);
        h2_ptz_truth->GetYaxis()->SetRangeUser(ptMin, ptMax);

        h2_ptz->GetYaxis()->SetRangeUser(ptMin, ptMax);

        //
        TH1D *h1_z_final = (TH1D *)h2_ptz_final->ProjectionX("z_final");
        TH1D *h1_z = (TH1D *)h2_ptz->ProjectionX("z");
        TH1D *h1_z_truth = (TH1D *)h2_ptz_truth->ProjectionX("z_truth");
        //
        
        h1_z_final->GetXaxis()->SetRange(binlowz, binhighz);
        h1_z_truth->GetXaxis()->SetRange(binlowz, binhighz);
        
        // This will depend on results of 1D jet pt unfolding
        //h1_z_final->Scale(1. / Njets_final, "width");
        //h1_z_truth->Scale(1. / Njets_truth, "width");
        //h1_z->Scale(1. / Njets_reco, "width");  
        
        NormalizeHist(h1_z_final);
        NormalizeHist(h1_z_truth);
        NormalizeHist(h1_z);
        

        TH1D *h1_z_ratio = (TH1D *)h1_z->Clone("h1_z_ratio");
        TH1D *h1_z_ratio_final = (TH1D *)h1_z_final->Clone("h1_z_ratio_final"); 


        h1_z_ratio->Divide(h1_z, h1_z_truth);
        h1_z_ratio_final->Divide(h1_z_final, h1_z_truth);

        h1_z_ratio_final->GetXaxis()->SetRange(binlowz, binhighz);
        
        
        ////////////////////////////////////
        // Smearing the Observables
        ///////////////////////////////////
        TH1D *h1_z_ptbinned_closure_error[ptbinsize-1];
        
        for (int i = 0; i < nRuns; i++) {
                TH2D *h2_ptz_smear = (TH2D *)h2_ptz->Clone(Form("ptz_smear%d", i));
                //myRNG->SetSeed();
                //myRNG->SetSeed(fixSmear);
                if (smear_by_data) 
                        SmearObservables(h2_ptz_smear, h2_ptz_data, myRNG);
                
                // Multiply by purity
                h2_ptz_smear->Multiply(h2_ptz_smear, h2_purity_ptz);

                RooUnfoldBayes unfold_ptz_smear(response_ptz, h2_ptz_smear, NumIters);
                h2_ptz_smear = (TH2D *)unfold_ptz_smear.Hreco();
                
                // Divide by efficiency    
                h2_ptz_smear->Divide(h2_ptz_smear, h2_eff_ptz);    
                
                h2_ptz_smear->GetYaxis()->SetRangeUser(ptMin, ptMax);
                
                TH1D *h1_z_smear = (TH1D *)h2_ptz_smear->ProjectionX(Form("z_smear%d", i));    
                //TH1D *h1_z_smear_20_100 = (TH1D *)h2_ptz_smear->ProjectionX(Form("z_smear_20_100_%d", i), 3, 5);  

                h1_z_smear->GetXaxis()->SetRange(binlowz, binhighz);
                h1_z_smear->GetXaxis()->SetTitle("z"); 
                        
                //h1_z_smear->Scale(1. / Njets_final, "width");
                NormalizeHist(h1_z_smear);

                std::cout << endl;    
                std::cout << "h1_z after norm : " << h1_z_smear->GetEntries() << std::endl;
                std::cout << endl;    
                        
                TH1D *h1_z_ratio = (TH1D*)h1_z_smear->Clone(Form("z_ratio_smear%d", i));
                h1_z_ratio->Divide(h1_z_smear, h1_z_truth);
                h1_z_ratio->Write();
                
                // pT bins [5-10, 10-12, 12-15, 15-20, 20-30, 30-50, 50-100]
                TH1D *h1_z_ptbinned[ptbinsize-1], *h1_z_truth_ptbinned[ptbinsize-1];
                TH1D *h1_z_eff_ptbinned[ptbinsize-1], *h1_z_purity_ptbinned[ptbinsize-1];   
                TH1D *h1_z_ptbinned_ratio[ptbinsize-1];    

                for (int j=1; j < ptbinsize; ++j)
                {
                        h1_z_ptbinned[j-1] = (TH1D *)h2_ptz_smear->ProjectionX(Form("z_smear%d_pt%d", i,j), j+1, j+1);     
                        h1_z_truth_ptbinned[j-1] = (TH1D *)h2_ptz_truth->ProjectionX(Form("z_truth_smear%d_pt%d", i,j), j+1, j+1);
                        h1_z_eff_ptbinned[j-1] = (TH1D *)h2_eff_ptz->ProjectionX(Form("z_eff_smear%d_pt%d", i,j), j+1, j+1);   
                        h1_z_purity_ptbinned[j-1] = (TH1D *)h2_purity_ptz->ProjectionX(Form("z_purty_smear%d_pt%d", i,j), j+1, j+1);   
                        h1_z_ptbinned[j-1]->GetXaxis()->SetRange(binlowz, binhighz);
                        h1_z_truth_ptbinned[j-1]->GetXaxis()->SetRange(binlowz, binhighz);
                        h1_z_ptbinned[j-1]->GetXaxis()->SetTitle("z");
                        h1_z_truth_ptbinned[j-1]->GetXaxis()->SetTitle("z");
                                
                        //h1_z_ptbinned[j-1]->Scale(1./h1_jetpt_final->Integral(j+1, j+1), "width");
                        //h1_z_truth_ptbinned[j-1]->Scale(1./h1_jetpt_truth->Integral(j+1, j+1), "width");       
                        NormalizeHist(h1_z_ptbinned[j-1]);
                        NormalizeHist(h1_z_truth_ptbinned[j-1]);                 
                        h1_z_ptbinned_ratio[j-1] = (TH1D*)h1_z_ptbinned[j-1]->Clone();
                        h1_z_ptbinned_ratio[j-1]->Divide(h1_z_ptbinned[j-1], h1_z_truth_ptbinned[j-1]);
                        h1_z_ptbinned_ratio[j-1]->Write(Form("z_ratio_smear%d_pt%d", i,j));           
                        
                        TString sys_title = DoShapeClosure ? "shape nonclosure" : "nonclosure";

                        sys_title += Form(" %.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]);

                        if (i == 0) {
                                h1_z_ptbinned_closure_error[j-1] = (TH1D *)h1_z_ptbinned_ratio[j-1]->Clone(Form("z_pt%d_closure_error", j));
                                
                                if (!DoShapeClosure)
                                        h1_z_ptbinned_closure_error[j-1]->Reset(); // Reset the histogram to zero
                        }

                        h1_z_ptbinned_closure_error[j-1]->GetYaxis()->SetTitle("#sigma_{sys}");
                        h1_z_ptbinned_closure_error[j-1]->SetTitle(sys_title);
                        
                        if (DoShapeClosure) {
                                SubtractUnity(h1_z_ptbinned_closure_error[j-1]);
                                h1_z_ptbinned_closure_error[j-1]->Write();
                        } else {
                                // Loop over all bins and add in quadrature
                                for (int x = 1; x <= h1_z_ptbinned_ratio[j-1]->GetNbinsX(); ++x) {
                                        double binContent = 0;
                                        if ( h1_z_ptbinned_ratio[j-1]->GetBinContent(x) != 0 )
                                                binContent = 1. - h1_z_ptbinned_ratio[j-1]->GetBinContent(x);
                                        
                                        double sumContent = h1_z_ptbinned_closure_error[j-1]->GetBinContent(x);

                                        h1_z_ptbinned_closure_error[j-1]->SetBinContent(x, sumContent + binContent * binContent);

                                        double binError = h1_z_ptbinned_ratio[j-1]->GetBinError(x);
                                        double sumError = h1_z_ptbinned_closure_error[j-1]->GetBinError(x);
                                        h1_z_ptbinned_closure_error[j-1]->SetBinError(x, sumError + binError * binError);
                                }                      
                        }   
                }
        }

        if (nRuns > 1) {
                for (int j=1; j < ptbinsize; ++j) {
                        // Take the square root of the final sum histogram to get the quadrature sum
                        for (int x = 1; x <= h1_z_ptbinned_closure_error[j-1]->GetNbinsX(); ++x) {
                                double sumContent = h1_z_ptbinned_closure_error[j-1]->GetBinContent(x);
                                h1_z_ptbinned_closure_error[j-1]->SetBinContent(x, TMath::Sqrt(sumContent) / TMath::Sqrt(nRuns));

                                double sumError = h1_z_ptbinned_closure_error[j-1]->GetBinError(x);
                                h1_z_ptbinned_closure_error[j-1]->SetBinError(x, TMath::Sqrt(sumError) / TMath::Sqrt(nRuns));
                        }

                        h1_z_ptbinned_closure_error[j-1]->Write();                      
                }
        }  
        
        
        /////////////////////   Create jt histograms /////////////////////////////////

        
        std::cout << "############################## Unfolding for ptjt ##############################" << std::endl;

        //TH1D *h1_jt_truth = (TH1D *)file_truth->Get("jt");
        TH2D *h2_ptjt_truth = (TH2D *)file_truth->Get("ptjt");
        TH2D *h2_ptjt = (TH2D *)file_reco->Get("ptjt");
        TH2D *h2_ptjt_data = (TH2D *)file_data->Get("ptjt");  
        TH2D *h2_ptjt_final = (TH2D *)h2_ptjt->Clone("ptjt_final");
        TH2D *h2_eff_ptjt = (TH2D *)file_unfold->Get("efficiency_ptjt");
        TH2D *h2_purity_ptjt = (TH2D *)file_unfold->Get("purity_ptjt");
        RooUnfoldResponse *response_ptjt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptjt");
        
        // Multiply by purity  
        h2_ptjt_final->Multiply(h2_ptjt_final, h2_purity_ptjt);
        RooUnfoldBayes unfold_ptjt(response_ptjt, h2_ptjt_final, NumIters);

        h2_ptjt_final = (TH2D *)unfold_ptjt.Hreco();
        // Divide by efficiency  
        h2_ptjt_final->Divide(h2_ptjt_final, h2_eff_ptjt);
        
        int binlowjt = h2_ptjt->GetXaxis()->FindBin(jt_binedges[0]);
        int binhighjt = h2_ptjt->GetXaxis()->FindBin(jt_binedges[jtbinsize-1]);      
        h2_ptjt_final->GetYaxis()->SetRangeUser(ptMin, ptMax);
        h2_ptjt_truth->GetYaxis()->SetRangeUser(ptMin, ptMax);

        h2_ptjt->GetYaxis()->SetRangeUser(ptMin, ptMax);

        //
        TH1D *h1_jt_final = (TH1D *)h2_ptjt_final->ProjectionX("jt_final");
        TH1D *h1_jt = (TH1D *)h2_ptjt->ProjectionX("jt");
        TH1D *h1_jt_truth = (TH1D *)h2_ptjt_truth->ProjectionX("jt_truth");
        //
        
        h1_jt_final->GetXaxis()->SetRange(binlowjt, binhighjt);
        h1_jt_truth->GetXaxis()->SetRange(binlowjt, binhighjt);
        
        // This will depend on results of 1D jet pt unfolding
        //h1_jt_final->Scale(1. / Njets_final, "width");
        //h1_jt_truth->Scale(1. / Njets_truth, "width");
        //h1_jt->Scale(1. / Njets_reco, "width");  
        
        NormalizeHist(h1_jt_final);
        NormalizeHist(h1_jt_truth);
        NormalizeHist(h1_jt);

        TH1D *h1_jt_ratio = (TH1D *)h1_jt->Clone("h1_jt_ratio");
        TH1D *h1_jt_ratio_final = (TH1D *)h1_jt_final->Clone("h1_jt_ratio_final"); 


        h1_jt_ratio->Divide(h1_jt, h1_jt_truth);
        h1_jt_ratio_final->Divide(h1_jt_final, h1_jt_truth);

        h1_jt_ratio_final->GetXaxis()->SetRange(binlowjt, binhighjt);
        
        
        ////////////////////////////////////
        // Smearing the Observables
        ///////////////////////////////////
        TH1D *h1_jt_ptbinned_closure_error[ptbinsize-1];
        for (int i = 0; i < nRuns; i++)
        {
        TH2D *h2_ptjt_smear = (TH2D *)h2_ptjt->Clone(Form("ptjt_smear%d", i));
        //myRNG->SetSeed();
        //myRNG->SetSeed(fixSmear);
        if (smear_by_data)
        {
                SmearObservables(h2_ptjt_smear, h2_ptjt_data, myRNG);
        }  

        // Multiply by purity
        h2_ptjt_smear->Multiply(h2_ptjt_smear, h2_purity_ptjt);
        RooUnfoldBayes unfold_ptjt_smear(response_ptjt, h2_ptjt_smear, NumIters);
        h2_ptjt_smear = (TH2D *)unfold_ptjt_smear.Hreco();
        // Divide by efficiency    
        h2_ptjt_smear->Divide(h2_ptjt_smear, h2_eff_ptjt);
        h2_ptjt_smear->GetYaxis()->SetRangeUser(ptMin, ptMax);  
        
        TH1D *h1_jt_smear = (TH1D *)h2_ptjt_smear->ProjectionX(Form("jt_smear%d", i));                        
        h1_jt_smear->GetXaxis()->SetRange(binlowjt, binhighjt);   
        h1_jt_smear->GetXaxis()->SetTitle("j_{T} (GeV/c)");
        //h1_jt_smear->Scale(1. / Njets_final, "width");     
        NormalizeHist(h1_jt_smear);
        TH1D *h1_jt_ratio = (TH1D*)h1_jt_smear->Clone(Form("jt_ratio_smear%d", i));
        h1_jt_ratio->Divide(h1_jt_smear, h1_jt_truth);
        h1_jt_ratio->Write();     

        // pT bins [5-10, 10-12, 12-15, 15-20, 20-30, 30-50, 50-100]
        TH1D *h1_jt_ptbinned[ptbinsize-1], *h1_jt_truth_ptbinned[ptbinsize-1];
        TH1D *h1_jt_ptbinned_ratio[ptbinsize-1];
        for (int j=1; j < ptbinsize; ++j)
        {
        h1_jt_ptbinned[j-1] = (TH1D *)h2_ptjt_smear->ProjectionX(Form("jt_smear%d_pt%d", i,j), j+1, j+1);    
        h1_jt_truth_ptbinned[j-1] = (TH1D *)h2_ptjt_truth->ProjectionX(Form("jt_truth_smear%d_pt%d", i,j), j+1, j+1);
        h1_jt_ptbinned[j-1]->GetXaxis()->SetRange(binlowjt, binhighjt);
        h1_jt_truth_ptbinned[j-1]->GetXaxis()->SetRange(binlowjt, binhighjt);
        h1_jt_ptbinned[j-1]->GetXaxis()->SetTitle("j_{T} (GeV/c)");
        h1_jt_truth_ptbinned[j-1]->GetXaxis()->SetTitle("j_{T} (GeV/c)");
                
        //h1_jt_ptbinned[j-1]->Scale(1./h1_jetpt_final->Integral(j+1, j+1), "width");
        //h1_jt_truth_ptbinned[j-1]->Scale(1./h1_jetpt_truth->Integral(j+1, j+1), "width");
        NormalizeHist(h1_jt_ptbinned[j-1]);
        NormalizeHist(h1_jt_truth_ptbinned[j-1]);
        h1_jt_ptbinned_ratio[j-1] = (TH1D*)h1_jt_ptbinned[j-1]->Clone();
        h1_jt_ptbinned_ratio[j-1]->Divide(h1_jt_ptbinned[j-1], h1_jt_truth_ptbinned[j-1]);       
        h1_jt_ptbinned_ratio[j-1]->Write(Form("jt_ratio_smear%d_pt%d", i,j));       
        
        sys_title = DoShapeClosure ? "shape nonclosure" : "nonclosure";
        sys_title += Form(" %.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]);      
        if (i == 0)
        {
                h1_jt_ptbinned_closure_error[j-1] = (TH1D *)h1_jt_ptbinned_ratio[j-1]->Clone(Form("jt_pt%d_closure_error", j));
                if (!DoShapeClosure)
                {
                h1_jt_ptbinned_closure_error[j-1]->Reset(); // Reset the histogram to zero
                }
        }
        h1_jt_ptbinned_closure_error[j-1]->GetYaxis()->SetTitle("#sigma_{sys}");
        h1_jt_ptbinned_closure_error[j-1]->SetTitle(sys_title);
        
        if (DoShapeClosure)
        {
                SubtractUnity(h1_jt_ptbinned_closure_error[j-1]);
                h1_jt_ptbinned_closure_error[j-1]->Write();
        }
        else 
        {
        // Loop over all bins and add in quadrature
                for (int x = 1; x <= h1_jt_ptbinned_ratio[j-1]->GetNbinsX(); ++x)
                {       
                double binContent = 0;
                if ( h1_jt_ptbinned_ratio[j-1]->GetBinContent(x) != 0 )
                { 
                binContent = 1. - h1_jt_ptbinned_ratio[j-1]->GetBinContent(x);
                }
                double sumContent = h1_jt_ptbinned_closure_error[j-1]->GetBinContent(x);
                h1_jt_ptbinned_closure_error[j-1]->SetBinContent(x, sumContent + binContent * binContent);

                double binError = h1_jt_ptbinned_ratio[j-1]->GetBinError(x);
                double sumError = h1_jt_ptbinned_closure_error[j-1]->GetBinError(x);
                h1_jt_ptbinned_closure_error[j-1]->SetBinError(x, sumError + binError * binError);
                }                    
        }  
        }        
        }  
        if (nRuns > 1)
        {
        for (int j=1; j < ptbinsize; ++j)
        {
        // Take the square root of the final sum histogram to get the quadrature sum
        for (int x = 1; x <= h1_jt_ptbinned_closure_error[j-1]->GetNbinsX(); ++x)
        {
                double sumContent = h1_jt_ptbinned_closure_error[j-1]->GetBinContent(x);
                h1_jt_ptbinned_closure_error[j-1]->SetBinContent(x, TMath::Sqrt(sumContent) / TMath::Sqrt(nRuns));

                double sumError = h1_jt_ptbinned_closure_error[j-1]->GetBinError(x);
                h1_jt_ptbinned_closure_error[j-1]->SetBinError(x, TMath::Sqrt(sumError) / TMath::Sqrt(nRuns));
        }
        h1_jt_ptbinned_closure_error[j-1]->Write();                
        }
        } 
        
        
        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        /////////////////////   Create r histograms /////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        std::cout << "############################## Unfolding for ptr ##############################" << std::endl;
        
        //TH1D *h1_r_truth = (TH1D *)file_truth->Get("r");
        TH2D *h2_ptr_truth = (TH2D *)file_truth->Get("ptr");
        TH2D *h2_ptr = (TH2D *)file_reco->Get("ptr");
        TH2D *h2_ptr_data = (TH2D *)file_data->Get("ptr");  
        TH2D *h2_ptr_final = (TH2D *)h2_ptr->Clone("ptr_final");
        TH2D *h2_eff_ptr = (TH2D *)file_unfold->Get("efficiency_ptr");
        TH2D *h2_purity_ptr = (TH2D *)file_unfold->Get("purity_ptr");
        RooUnfoldResponse *response_ptr = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptr");
        
        // Multiply by purity  
        h2_ptr_final->Multiply(h2_ptr_final, h2_purity_ptr);
        RooUnfoldBayes unfold_ptr(response_ptr, h2_ptr_final, NumIters);

        h2_ptr_final = (TH2D *)unfold_ptr.Hreco();
        // Divide by efficiency  
        h2_ptr_final->Divide(h2_ptr_final, h2_eff_ptr);
        
        int binlowr = h2_ptr->GetXaxis()->FindBin(r_binedges[0]);
        int binhighr = h2_ptr->GetXaxis()->FindBin(r_binedges[rbinsize-1]);  

        h2_ptr_final->GetYaxis()->SetRangeUser(ptMin, ptMax);
        h2_ptr_truth->GetYaxis()->SetRangeUser(ptMin, ptMax);

        h2_ptr->GetYaxis()->SetRangeUser(ptMin, ptMax);

        //
        TH1D *h1_r_final = (TH1D *)h2_ptr_final->ProjectionX("r_final");
        TH1D *h1_r = (TH1D *)h2_ptr->ProjectionX("r");
        TH1D *h1_r_truth = (TH1D *)h2_ptr_truth->ProjectionX("r_truth");
        //
        
        h1_r_final->GetXaxis()->SetRange(binlowr, binhighr);
        h1_r_truth->GetXaxis()->SetRange(binlowr, binhighr);
        
        // This will depend on results of 1D jet pt unfolding
        //h1_r_final->Scale(1. / Njets_final, "width");
        //h1_r_truth->Scale(1. / Njets_truth, "width");
        //h1_r->Scale(1. / Njets_reco, "width");  
        
        NormalizeHist(h1_r_final);
        NormalizeHist(h1_r_truth);
        NormalizeHist(h1_r);

        TH1D *h1_r_ratio = (TH1D *)h1_r->Clone("h1_r_ratio");
        TH1D *h1_r_ratio_final = (TH1D *)h1_r_final->Clone("h1_r_ratio_final"); 


        h1_r_ratio->Divide(h1_r, h1_r_truth);
        h1_r_ratio_final->Divide(h1_r_final, h1_r_truth);

        h1_r_ratio_final->GetXaxis()->SetRange(binlowr, binhighr);
        
        
        ////////////////////////////////////
        // Smearing the Observables
        ///////////////////////////////////
        TH1D *h1_r_ptbinned_closure_error[ptbinsize-1];
        for (int i = 0; i < nRuns; i++)
        {
        TH2D *h2_ptr_smear = (TH2D *)h2_ptr->Clone(Form("ptr_smear%d", i));
        //myRNG->SetSeed();
        //myRNG->SetSeed(fixSmear);
        if (smear_by_data)
        {
                SmearObservables(h2_ptr_smear, h2_ptr_data, myRNG);
        }
        
        // Multiply by purity
        h2_ptr_smear->Multiply(h2_ptr_smear, h2_purity_ptr);
        
        RooUnfoldBayes unfold_ptr_smear(response_ptr, h2_ptr_smear, NumIters);
        h2_ptr_smear = (TH2D *)unfold_ptr_smear.Hreco();
        
        // Divide by efficiency    
        h2_ptr_smear->Divide(h2_ptr_smear, h2_eff_ptr);
        
        h2_ptr_smear->GetYaxis()->SetRangeUser(ptMin, ptMax);
                
        TH1D *h1_r_smear = (TH1D *)h2_ptr_smear->ProjectionX(Form("r_smear%d", i));                     
        h1_r_smear->GetXaxis()->SetRange(binlowr, binhighr);    
        h1_r_smear->GetXaxis()->SetTitle("r");    
        //h1_r_smear->Scale(1. / Njets_final, "width");
        NormalizeHist(h1_r_smear); 
        
        TH1D *h1_r_ratio = (TH1D*)h1_r_smear->Clone(Form("r_ratio_smear%d", i));
        h1_r_ratio->Divide(h1_r_smear, h1_r_truth);
        h1_r_ratio->Write();
        
        // pT bins [5-10, 10-12, 12-15, 15-20, 20-30, 30-50, 50-100]
        TH1D *h1_r_ptbinned[ptbinsize-1], *h1_r_truth_ptbinned[ptbinsize-1];
        TH1D *h1_r_ptbinned_ratio[ptbinsize-1];
        for (int j=1; j < ptbinsize; ++j)
        {
        h1_r_ptbinned[j-1] = (TH1D *)h2_ptr_smear->ProjectionX(Form("r_smear%d_pt%d", i,j), j+1, j+1);    
        h1_r_truth_ptbinned[j-1] = (TH1D *)h2_ptr_truth->ProjectionX(Form("r_truth_smear%d_pt%d", i,j), j+1, j+1);
        h1_r_ptbinned[j-1]->GetXaxis()->SetRange(binlowr, binhighr);
        h1_r_truth_ptbinned[j-1]->GetXaxis()->SetRange(binlowr, binhighr);
        h1_r_ptbinned[j-1]->GetXaxis()->SetTitle("r");
        h1_r_truth_ptbinned[j-1]->GetXaxis()->SetTitle("r");
                
        //h1_r_ptbinned[j-1]->Scale(1./h1_jetpt_final->Integral(j+1, j+1), "width");
        //h1_r_truth_ptbinned[j-1]->Scale(1./h1_jetpt_truth->Integral(j+1, j+1), "width");
        NormalizeHist(h1_r_ptbinned[j-1]);
        NormalizeHist(h1_r_truth_ptbinned[j-1]);
        h1_r_ptbinned_ratio[j-1] = (TH1D*)h1_r_ptbinned[j-1]->Clone();
        h1_r_ptbinned_ratio[j-1]->Divide(h1_r_ptbinned[j-1], h1_r_truth_ptbinned[j-1]);       
        h1_r_ptbinned_ratio[j-1]->Write(Form("r_ratio_smear%d_pt%d", i,j));      
        
        sys_title = DoShapeClosure ? "shape nonclosure" : "nonclosure";
        sys_title += Form(" %.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]);        
        if (i == 0)
        {
                h1_r_ptbinned_closure_error[j-1] = (TH1D *)h1_r_ptbinned_ratio[j-1]->Clone(Form("r_pt%d_closure_error", j));
                if (!DoShapeClosure)
                {
                h1_r_ptbinned_closure_error[j-1]->Reset(); // Reset the histogram to zero
                }
        }     
        h1_r_ptbinned_closure_error[j-1]->GetYaxis()->SetTitle("#sigma_{sys}");
        h1_r_ptbinned_closure_error[j-1]->SetTitle(sys_title);
                
        if (DoShapeClosure)
        {
                SubtractUnity(h1_r_ptbinned_closure_error[j-1]);
                h1_r_ptbinned_closure_error[j-1]->Write();
        }
        else 
        {
                // Loop over all bins and add in quadrature
                for (int x = 1; x <= h1_r_ptbinned_ratio[j-1]->GetNbinsX(); ++x)
                {
                double binContent = 0;
                if ( h1_r_ptbinned_ratio[j-1]->GetBinContent(x) != 0 )
                { 
                binContent = 1. - h1_r_ptbinned_ratio[j-1]->GetBinContent(x);
                }
                double sumContent = h1_r_ptbinned_closure_error[j-1]->GetBinContent(x);
                h1_r_ptbinned_closure_error[j-1]->SetBinContent(x, sumContent + binContent * binContent);

                double binError = h1_r_ptbinned_ratio[j-1]->GetBinError(x);
                double sumError = h1_r_ptbinned_closure_error[j-1]->GetBinError(x);
                h1_r_ptbinned_closure_error[j-1]->SetBinError(x, sumError + binError * binError);
                }                   
        }
        }         
        }  
        if (nRuns > 1)
        {
        for (int j=1; j < ptbinsize; ++j)
        {
        // Take the square root of the final sum histogram to get the quadrature sum
        for (int x = 1; x <= h1_r_ptbinned_closure_error[j-1]->GetNbinsX(); ++x)
        {
                double sumContent = h1_r_ptbinned_closure_error[j-1]->GetBinContent(x);
                h1_r_ptbinned_closure_error[j-1]->SetBinContent(x, TMath::Sqrt(sumContent) / TMath::Sqrt(nRuns));

                double sumError = h1_r_ptbinned_closure_error[j-1]->GetBinError(x);
                h1_r_ptbinned_closure_error[j-1]->SetBinError(x, TMath::Sqrt(sumError) / TMath::Sqrt(nRuns));
        }   
        h1_r_ptbinned_closure_error[j-1]->Write();                
        }
        }    
        
        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        /////////////////////   Create zjt histograms /////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        std::cout << "############################## Unfolding for ptzjt ##############################" << std::endl;
        
        //TH2D *h2_zjt_truth = (TH2D *)file_truth->Get("zjt");
        TH3D *h3_ptzjt_truth = (TH3D *)file_truth->Get("ptzjt");
        TH3D *h3_ptzjt = (TH3D *)file_reco->Get("ptzjt");
        TH3D *h3_ptzjt_data = (TH3D *)file_data->Get("ptzjt");  
        TH3D *h3_ptzjt_final = (TH3D *)h3_ptzjt->Clone("ptzjt_final");
        TH3D *h3_eff_ptzjt = (TH3D *)file_unfold->Get("efficiency_ptzjt");
        TH3D *h3_purity_ptzjt = (TH3D *)file_unfold->Get("purity_ptzjt");
        RooUnfoldResponse *response_ptzjt = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptzjt");
        
        // Multiply by purity  
        h3_ptzjt_final->Multiply(h3_ptzjt_final, h3_purity_ptzjt);
        RooUnfoldBayes unfold_ptzjt(response_ptzjt, h3_ptzjt_final, NumIters);
        h3_ptzjt_final = (TH3D *)unfold_ptzjt.Hreco();
        // Divide by efficiency  
        h3_ptzjt_final->Divide(h3_ptzjt_final, h3_eff_ptzjt);

        int binlowz2D = h3_ptzjt->GetXaxis()->FindBin(z_binedges2D[0]);
        int binhighz2D = h3_ptzjt->GetXaxis()->FindBin(z_binedges2D[zbinsize2D-1]);
        int binlowjt2D = h3_ptzjt->GetYaxis()->FindBin(jt_binedges2D[0]);
        int binhighjt2D = h3_ptzjt->GetYaxis()->FindBin(jt_binedges2D[jtbinsize2D-1]);

        h3_ptzjt_final->GetZaxis()->SetRangeUser(ptMin, ptMax);
        h3_ptzjt_truth->GetZaxis()->SetRangeUser(ptMin, ptMax);

        h3_ptzjt->GetZaxis()->SetRangeUser(ptMin, ptMax);

        TH2D *h2_zjt_final = (TH2D *)h3_ptzjt_final->Project3D("yx");
        h2_zjt_final->SetName("zjt_final");
        TH2D *h2_zjt = (TH2D *)h3_ptzjt->Project3D("yx");
        h2_zjt->SetName("zjt");
        TH2D *h2_zjt_truth = (TH2D *)h3_ptzjt_truth->Project3D("yx");
        h2_zjt_truth->SetName("zjt_truth");
        //
        h2_zjt_final->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_final->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zjt_truth->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_truth->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        
        
        // This will depend on results of 1D jet unfolding
        //h2_zjt_final->Scale(1. / Njets_final, "width");
        //h2_zjt_truth->Scale(1. / Njets_truth, "width");
        NormalizeHist(h2_zjt_final);
        NormalizeHist(h2_zjt_truth);

        TH2D *h2_zjt_ratio = (TH2D *)h2_zjt->Clone("h2_zjt_ratio");
        TH2D *h2_zjt_ratio_final = (TH2D *)h2_zjt->Clone("h2_zjt_ratio_final");
        TH2D *h2_zjt_pullratio_final = (TH2D *)h2_zjt->Clone("h2_zjt_pullratio_final");

        //TH2D *h2_zjt_pull = (TH2D *)h2_zjt_truth->Clone("zjt_pull");

        //TH1D *h1_zjt_pulldist = new TH1D("zjt_pulldist", "", 20, -4.5, 4.5);
        TH1D *h1_zjt_pulldist_final = new TH1D("zjt_pulldist_final", "", 20, -4.5, 4.5);
        //TH1D *h1_ptzjt_pulldist_vecfinal = new TH1D("ptzjt_pulldist_vecfinal", "", 20, -4.5, 4.5);

        h2_zjt_ratio->Divide(h2_zjt, h2_zjt_truth);
        h2_zjt_ratio_final->Divide(h2_zjt_final, h2_zjt_truth);

        h2_zjt_ratio_final->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_ratio_final->GetXaxis()->SetRange(binlowz2D, binhighz2D);

        h2_zjt_pullratio_final->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);      
        h2_zjt_pullratio_final->GetXaxis()->SetRange(binlowz2D, binhighz2D);

        //GetPullsRatio(h2_zjt_ratio, h2_zjt_pull, h1_zjt_pulldist);
        GetPullsRatio(h2_zjt_ratio_final, h2_zjt_pullratio_final, h1_zjt_pulldist_final);

        //double zjt_pullperf1 = GetPullPerformance(h2_zjt_pull, 1);
        //double zjt_pullperf2 = GetPullPerformance(h2_zjt_pull, 2);

        //double zjt_pullperf1_final = GetPullPerformance(h2_zjt_pull_final, 1);
        //double zjt_pullperf2_final = GetPullPerformance(h2_zjt_pull_final, 2);
        
        h2_zjt_ratio_final->SetOption("COLZ, text");  
        h2_zjt_ratio_final->Write();

        h2_zjt_pullratio_final->SetOption("COLZ, text");
        h2_zjt_pullratio_final->Write();
        h1_zjt_pulldist_final->Write();

        
        ////////////////////////////////////
        // Smearing the Observables
        ///////////////////////////////////
        TH2D *h2_zjt_ptbinned_closure_error[ptbinsize-1];
        for (int i = 0; i < nRuns; i++)
        {
        TH3D *h3_ptzjt_smear = (TH3D *)h3_ptzjt->Clone(Form("ptzjt_smear%d", i));
        //myRNG->SetSeed();
        //myRNG->SetSeed(fixSmear);
        if (smear_by_data)
        {
                SmearObservables(h3_ptzjt_smear, h3_ptzjt_data, myRNG);
        }
        // Multiply by purity
        h3_ptzjt_smear->Multiply(h3_ptzjt_smear, h3_purity_ptzjt);
        RooUnfoldBayes unfold_ptzjt_smear(response_ptzjt, h3_ptzjt_smear, NumIters);

        h3_ptzjt_smear = (TH3D *)unfold_ptzjt_smear.Hreco();
        // Divide by efficiency
        h3_ptzjt_smear->Divide(h3_ptzjt_smear, h3_eff_ptzjt);
        h3_ptzjt_smear->GetZaxis()->SetRangeUser(ptMin, ptMax);

        TH2D *h2_zjt_smear = (TH2D *)h3_ptzjt_smear->Project3D("yx");
        h2_zjt_smear->SetName(Form("zjt_smear%d", i));

        h2_zjt_smear->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_smear->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        //h2_zjt_smear->Scale(1. / Njets_final, "width");
        NormalizeHist(h2_zjt_smear);
        TH2D *h2_zjt_ratio_smear = (TH2D *)h2_zjt_smear->Clone(Form("zjt_ratio_smear%d", i));
        //TH2D *h2_zjt_pull_smear = (TH2D *)h2_zjt_smear->Clone(Form("h2_zjt_pull_smear%d", i));
        h2_zjt_ratio_smear->Divide(h2_zjt_smear, h2_zjt_truth); 

        h2_zjt_ratio_smear->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_ratio_smear->GetXaxis()->SetRange(binlowz2D, binhighz2D);

        //h2_zjt_pull_smear->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        //h2_zjt_pull_smear->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zjt_ratio_smear->SetOption("COLZ, text");    
        h2_zjt_ratio_smear->Write();      
        //h2_zjt_pull_smear->Write(); 
        
        // pT bins [5-10, 10-12, 12-15, 15-20, 20-30, 30-50, 50-100]
        TH2D *h2_zjt_ptbinned[ptbinsize-1], *h2_zjt_truth_ptbinned[ptbinsize-1], *h2_zjt_ptbinned_ratio[ptbinsize-1], *h2_zjt_ptbinned_pullratio[ptbinsize-1];
        TH1D *h1_zjt_ptbinned_pulldist[ptbinsize-1];
        for (int j=1; j < ptbinsize; ++j)
        {
        h3_ptzjt_smear->GetZaxis()->SetRange(j+1, j+1);
        h3_ptzjt_truth->GetZaxis()->SetRange(j+1, j+1);        
        h2_zjt_ptbinned[j-1] = (TH2D *)h3_ptzjt_smear->Project3D("yx");  
        h2_zjt_ptbinned[j-1]->SetName(Form("zjt_smear%d_pt%d", i,j)); 
        h2_zjt_truth_ptbinned[j-1] = (TH2D *)h3_ptzjt_truth->Project3D("yx");
        h2_zjt_truth_ptbinned[j-1]->SetName(Form("zjt_truth_smear%d_pt%d", i,j));
        h2_zjt_ptbinned_pullratio[j-1] = (TH2D *)h2_zjt_ptbinned[j-1]->Clone(Form("zjt_pullratio_smear%d_pt%d", i, j));
        h1_zjt_ptbinned_pulldist[j-1] = new TH1D(Form("zjt_pulldist_smear%d_pt%d", i, j), ";pulls", 20, -4.5, 4.5);


        h2_zjt_ptbinned[j-1]->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zjt_truth_ptbinned[j-1]->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zjt_ptbinned[j-1]->GetXaxis()->SetTitle("z");
        h2_zjt_truth_ptbinned[j-1]->GetXaxis()->SetTitle("z");
        h2_zjt_ptbinned[j-1]->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_truth_ptbinned[j-1]->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_zjt_ptbinned[j-1]->GetYaxis()->SetTitle("j_{T}");
        h2_zjt_truth_ptbinned[j-1]->GetYaxis()->SetTitle("j_{T}");    
        h2_zjt_ptbinned_pullratio[j-1]->GetYaxis()->SetRange(binlowjt2D, binhighjt2D);      
        h2_zjt_ptbinned_pullratio[j-1]->GetXaxis()->SetRange(binlowz2D, binhighz2D);    
                
        //h2_zjt_ptbinned[j-1]->Scale(1./h1_jetpt_final->Integral(j+1, j+1));
        //h2_zjt_truth_ptbinned[j-1]->Scale(1./h1_jetpt_truth->Integral(j+1, j+1));
        NormalizeHist(h2_zjt_ptbinned[j-1]);
        NormalizeHist(h2_zjt_truth_ptbinned[j-1]);
        h2_zjt_ptbinned_ratio[j-1] = (TH2D *)h2_zjt_ptbinned[j-1]->Clone(Form("zjt_ratio_smear%d_pt%d", i,j));
        h2_zjt_ptbinned_ratio[j-1]->Divide(h2_zjt_ptbinned[j-1], h2_zjt_truth_ptbinned[j-1]);     
        h2_zjt_ptbinned_ratio[j-1]->SetOption("COLZ, text");         
        h2_zjt_ptbinned_ratio[j-1]->Write();
        
        GetPullsRatio(h2_zjt_ptbinned_ratio[j-1], h2_zjt_ptbinned_pullratio[j-1], h1_zjt_ptbinned_pulldist[j-1]);
        h2_zjt_ptbinned_pullratio[j-1]->SetOption("COLZ, text");
        h2_zjt_ptbinned_pullratio[j-1]->Write();
        h1_zjt_ptbinned_pulldist[j-1]->Write();

                
        sys_title = DoShapeClosure ? "shape nonclosure" : "nonclosure";
        sys_title += Form(" %.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]);        
        if (i == 0)
        {
                h2_zjt_ptbinned_closure_error[j-1] = (TH2D *)h2_zjt_ptbinned_ratio[j-1]->Clone(Form("zjt_pt%d_closure_error", j));
                if (!DoShapeClosure)
                {
                h2_zjt_ptbinned_closure_error[j-1]->Reset(); // Reset the histogram to zero
                }
        }
        h2_zjt_ptbinned_closure_error[j-1]->GetZaxis()->SetTitle("#sigma_{sys}");
        h2_zjt_ptbinned_closure_error[j-1]->SetTitle(sys_title);      
        
        if (DoShapeClosure)
        {
                SubtractUnity(h2_zjt_ptbinned_closure_error[j-1]);
                h2_zjt_ptbinned_closure_error[j-1]->Write();
        }
        else 
        {
                // Loop over all bins and add in quadrature
                for (int x = 1; x <= h2_zjt_ptbinned_ratio[j-1]->GetNbinsX(); ++x)
                {
                for (int y = 1; y <= h2_zjt_ptbinned_ratio[j-1]->GetNbinsY(); ++y)
                {
        
                double binContent = 0;
                if ( h2_zjt_ptbinned_ratio[j-1]->GetBinContent(x,y) != 0 )
                { 
                binContent = 1. - h2_zjt_ptbinned_ratio[j-1]->GetBinContent(x, y);
                }
                double sumContent = h2_zjt_ptbinned_closure_error[j-1]->GetBinContent(x, y);
                h2_zjt_ptbinned_closure_error[j-1]->SetBinContent(x, y, sumContent + binContent * binContent);
        
                double binError = h2_zjt_ptbinned_ratio[j-1]->GetBinError(x, y);
                double sumError = h2_zjt_ptbinned_closure_error[j-1]->GetBinError(x, y);
                h2_zjt_ptbinned_closure_error[j-1]->SetBinError(x, y, sumError + binError * binError);
                }
                }                     
        }
        }    
        }
        if (nRuns > 1)
        {
        for (int j=1; j < ptbinsize; ++j)
        {
        // Take the square root of the final sum histogram to get the quadrature sum
        for (int x = 1; x <= h2_zjt_ptbinned_closure_error[j-1]->GetNbinsX(); ++x)
        {
                for (int y = 1; y <= h2_zjt_ptbinned_closure_error[j-1]->GetNbinsY(); ++y)
                {      
                double sumContent = h2_zjt_ptbinned_closure_error[j-1]->GetBinContent(x, y);
                h2_zjt_ptbinned_closure_error[j-1]->SetBinContent(x, y, TMath::Sqrt(sumContent) / TMath::Sqrt(nRuns));

                double sumError = h2_zjt_ptbinned_closure_error[j-1]->GetBinError(x, y);
                h2_zjt_ptbinned_closure_error[j-1]->SetBinError(x, y, TMath::Sqrt(sumError) / TMath::Sqrt(nRuns));
                }
        }
        h2_zjt_ptbinned_closure_error[j-1]->Write();                
        }  
        }                          


        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        /////////////////////   Create zr histograms /////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        std::cout << "############################## Unfolding for ptzr ##############################" << std::endl;
        
        //TH2D *h2_zr_truth = (TH2D *)file_truth->Get("zr");
        TH3D *h3_ptzr_truth = (TH3D *)file_truth->Get("ptzr");
        TH3D *h3_ptzr = (TH3D *)file_reco->Get("ptzr");
        TH3D *h3_ptzr_data = (TH3D *)file_data->Get("ptzr");  
        TH3D *h3_ptzr_final = (TH3D *)h3_ptzr->Clone("ptzr_final");
        TH3D *h3_eff_ptzr = (TH3D *)file_unfold->Get("efficiency_ptzr");
        TH3D *h3_purity_ptzr = (TH3D *)file_unfold->Get("purity_ptzr");
        RooUnfoldResponse *response_ptzr = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptzr");
        
        // Multiply by purity  
        h3_ptzr_final->Multiply(h3_ptzr_final, h3_purity_ptzr);
        RooUnfoldBayes unfold_ptzr(response_ptzr, h3_ptzr_final, NumIters);
        h3_ptzr_final = (TH3D *)unfold_ptzr.Hreco();
        // Divide by efficiency  
        h3_ptzr_final->Divide(h3_ptzr_final, h3_eff_ptzr);

        int binlowr2D = h3_ptzr->GetYaxis()->FindBin(r_binedges2D[0]);
        int binhighr2D = h3_ptzr->GetYaxis()->FindBin(r_binedges2D[rbinsize2D-1]);

        h3_ptzr_final->GetZaxis()->SetRangeUser(ptMin, ptMax);
        h3_ptzr_truth->GetZaxis()->SetRangeUser(ptMin, ptMax);

        h3_ptzr->GetZaxis()->SetRangeUser(ptMin, ptMax);

        TH2D *h2_zr_final = (TH2D *)h3_ptzr_final->Project3D("yx");
        h2_zr_final->SetName("zr_final");
        TH2D *h2_zr = (TH2D *)h3_ptzr->Project3D("yx");
        h2_zr->SetName("zr");
        TH2D *h2_zr_truth = (TH2D *)h3_ptzr_truth->Project3D("yx");
        h2_zr_truth->SetName("zr_truth");
        //
        h2_zr_final->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_final->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zr_truth->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_truth->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        
        
        // This will depend on results of 1D jet unfolding
        //h2_zr_final->Scale(1. / Njets_final, "width");  
        //h2_zr_truth->Scale(1. / Njets_truth, "width");
        NormalizeHist(h2_zr_final);
        NormalizeHist(h2_zr_truth);
        
        TH2D *h2_zr_ratio = (TH2D *)h2_zr->Clone("h2_zr_ratio");
        TH2D *h2_zr_ratio_final = (TH2D *)h2_zr->Clone("h2_zr_ratio_final");
        //TH2D *h2_zr_pull_final = (TH2D *)h2_zr->Clone("h2_zr_pull_final");

        //TH2D *h2_zr_pull = (TH2D *)h2_zr_truth->Clone("zr_pull");

        //TH1D *h1_zr_pulldist = new TH1D("zr_pulldist", "", 20, -4.5, 4.5);
        //TH1D *h1_zr_pulldist_final = new TH1D("zr_pulldist_final", "", 20, -4.5, 4.5);
        //TH1D *h1_ptzr_pulldist_vecfinal = new TH1D("ptzr_pulldist_vecfinal", "", 20, -4.5, 4.5);

        h2_zr_ratio->Divide(h2_zr, h2_zr_truth);
        h2_zr_ratio_final->Divide(h2_zr_final, h2_zr_truth);

        h2_zr_ratio_final->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_ratio_final->GetXaxis()->SetRange(binlowz2D, binhighz2D);

        //h2_zr_pull_final->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        //h2_zr_pull_final->GetXaxis()->SetRange(binlowz2D, binhighz2D);

        //GetPullsRatio(h2_zr_ratio, h2_zr_pull, h1_zr_pulldist);
        //GetPullsRatio(h2_zr_ratio_final, h2_zr_pull_final, h1_zr_pulldist_final);

        //double zr_pullperf1 = GetPullPerformance(h2_zr_pull, 1);
        //double zr_pullperf2 = GetPullPerformance(h2_zr_pull, 2);

        //double zr_pullperf1_final = GetPullPerformance(h2_zr_pull_final, 1);
        //double zr_pullperf2_final = GetPullPerformance(h2_zr_pull_final, 2);
        h2_zr_ratio_final->SetOption("COLZ, text");  
        h2_zr_ratio_final->Write();

        
        ////////////////////////////////////
        // Smearing the Observables
        ///////////////////////////////////
        TH2D *h2_zr_ptbinned_closure_error[ptbinsize-1];
        for (int i = 0; i < nRuns; i++)
        {
        TH3D *h3_ptzr_smear = (TH3D *)h3_ptzr->Clone(Form("ptzr_smear%d", i));
        //myRNG->SetSeed();
        //myRNG->SetSeed(fixSmear);
        if (smear_by_data)
        {
                SmearObservables(h3_ptzr_smear, h3_ptzr_data, myRNG);
        }
        // Multiply by purity
        h3_ptzr_smear->Multiply(h3_ptzr_smear, h3_purity_ptzr);
        RooUnfoldBayes unfold_ptzr_smear(response_ptzr, h3_ptzr_smear, NumIters);

        h3_ptzr_smear = (TH3D *)unfold_ptzr_smear.Hreco();
        // Divide by efficiency
        h3_ptzr_smear->Divide(h3_ptzr_smear, h3_eff_ptzr);
        h3_ptzr_smear->GetZaxis()->SetRangeUser(ptMin, ptMax);

        TH2D *h2_zr_smear = (TH2D *)h3_ptzr_smear->Project3D("yx");
        h2_zr_smear->SetName(Form("zr_smear%d", i));

        h2_zr_smear->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_smear->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        //h2_zr_smear->Scale(1. / Njets_final, "width");
        NormalizeHist(h2_zr_smear);
        TH2D *h2_zr_ratio_smear = (TH2D *)h2_zr_smear->Clone(Form("zr_ratio_smear%d", i));
        //TH2D *h2_zr_pull_smear = (TH2D *)h2_zr_smear->Clone(Form("h2_zr_pull_smear%d", i));
        h2_zr_ratio_smear->Divide(h2_zr_smear, h2_zr_truth); 

        h2_zr_ratio_smear->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_ratio_smear->GetXaxis()->SetRange(binlowz2D, binhighz2D);

        //h2_zr_pull_smear->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        //h2_zr_pull_smear->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zr_ratio_smear->SetOption("COLZ, text");
        h2_zr_ratio_smear->Write();      
        //h2_zr_pull_smear->Write();   

        // pT bins [5-10, 10-12, 12-15, 15-20, 20-30, 30-50, 50-100]
        TH2D *h2_zr_ptbinned[ptbinsize-1], *h2_zr_truth_ptbinned[ptbinsize-1], *h2_zr_ptbinned_ratio[ptbinsize-1], *h2_zr_ptbinned_pullratio[ptbinsize-1];
        TH1D *h1_zr_ptbinned_pulldist[ptbinsize-1];
        for (int j=1; j < ptbinsize; ++j)
        {
        h3_ptzr_smear->GetZaxis()->SetRange(j+1, j+1);
        h3_ptzr_truth->GetZaxis()->SetRange(j+1, j+1);        
        h2_zr_ptbinned[j-1] = (TH2D *)h3_ptzr_smear->Project3D("yx");    
        h2_zr_ptbinned[j-1]->SetName(Form("zr_smear%d_pt%d", i,j));
        h2_zr_truth_ptbinned[j-1] = (TH2D *)h3_ptzr_truth->Project3D("yx");
        h2_zr_truth_ptbinned[j-1]->SetName(Form("zr_truth_smear%d_pt%d", i,j));
        h2_zr_ptbinned_pullratio[j-1] = (TH2D *)h2_zr_ptbinned[j-1]->Clone(Form("zr_pullratio_smear%d_pt%d", i, j));
        h1_zr_ptbinned_pulldist[j-1] = new TH1D(Form("zr_pulldist_smear%d_pt%d", i, j), ";pulls", 20, -4.5, 4.5);

        h2_zr_ptbinned[j-1]->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zr_truth_ptbinned[j-1]->GetXaxis()->SetRange(binlowz2D, binhighz2D);
        h2_zr_ptbinned[j-1]->GetXaxis()->SetTitle("z");
        h2_zr_truth_ptbinned[j-1]->GetXaxis()->SetTitle("z");
        h2_zr_ptbinned[j-1]->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_truth_ptbinned[j-1]->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_zr_ptbinned[j-1]->GetYaxis()->SetTitle("r");
        h2_zr_truth_ptbinned[j-1]->GetYaxis()->SetTitle("r");
        h2_zr_ptbinned_pullratio[j-1]->GetYaxis()->SetRange(binlowr2D, binhighr2D);      
        h2_zr_ptbinned_pullratio[j-1]->GetXaxis()->SetRange(binlowz2D, binhighz2D);           
                
        //h2_zr_ptbinned[j-1]->Scale(1./h1_jetpt_final->Integral(j+1, j+1));
        //h2_zr_truth_ptbinned[j-1]->Scale(1./h1_jetpt_truth->Integral(j+1, j+1));
        NormalizeHist(h2_zr_ptbinned[j-1]);
        NormalizeHist(h2_zr_truth_ptbinned[j-1]);
        h2_zr_ptbinned_ratio[j-1] = (TH2D *)h2_zr_ptbinned[j-1]->Clone(Form("zr_ratio_smear%d_pt%d", i,j));
        h2_zr_ptbinned_ratio[j-1]->Divide(h2_zr_ptbinned[j-1], h2_zr_truth_ptbinned[j-1]);     
        h2_zr_ptbinned_ratio[j-1]->SetOption("COLZ, text");         
        h2_zr_ptbinned_ratio[j-1]->Write();   
        
        GetPullsRatio(h2_zr_ptbinned_ratio[j-1], h2_zr_ptbinned_pullratio[j-1], h1_zr_ptbinned_pulldist[j-1]);
        h2_zr_ptbinned_pullratio[j-1]->SetOption("COLZ, text");
        h2_zr_ptbinned_pullratio[j-1]->Write();
        h1_zr_ptbinned_pulldist[j-1]->Write();

        sys_title = DoShapeClosure ? "shape nonclosure" : "nonclosure";
        sys_title += Form(" %.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]);          
        if (i == 0)
        {
                h2_zr_ptbinned_closure_error[j-1] = (TH2D *)h2_zr_ptbinned_ratio[j-1]->Clone(Form("zr_pt%d_closure_error", j));
                if (!DoShapeClosure)
                {
                h2_zr_ptbinned_closure_error[j-1]->Reset(); // Reset the histogram to zero
                }
        }
        h2_zr_ptbinned_closure_error[j-1]->GetZaxis()->SetTitle("#sigma_{sys}");
        h2_zr_ptbinned_closure_error[j-1]->SetTitle(sys_title);      
                
        if (DoShapeClosure)
        {
                SubtractUnity(h2_zr_ptbinned_closure_error[j-1]);
                h2_zr_ptbinned_closure_error[j-1]->Write();
        }
        else 
        {
                // Loop over all bins and add in quadrature   
                for (int x = 1; x <= h2_zr_ptbinned_ratio[j-1]->GetNbinsX(); ++x)
                {
                for (int y = 1; y <= h2_zr_ptbinned_ratio[j-1]->GetNbinsY(); ++y)
                {
                double binContent = 0;
                if ( h2_zr_ptbinned_ratio[j-1]->GetBinContent(x,y) != 0 )
                { 
                binContent = 1. - h2_zr_ptbinned_ratio[j-1]->GetBinContent(x, y);
                }
                double sumContent = h2_zr_ptbinned_closure_error[j-1]->GetBinContent(x, y);         
                h2_zr_ptbinned_closure_error[j-1]->SetBinContent(x, y, sumContent + binContent * binContent);          
        
                double binError = h2_zr_ptbinned_ratio[j-1]->GetBinError(x, y);
                double sumError = h2_zr_ptbinned_closure_error[j-1]->GetBinError(x, y);
                h2_zr_ptbinned_closure_error[j-1]->SetBinError(x, y, sumError + binError * binError);
                }
                }                 
        }
        }         
        }

        if (nRuns > 1)
        {
        for (int j=1; j < ptbinsize; ++j)
        {
        // Take the square root of the final sum histogram to get the quadrature sum
        for (int x = 1; x <= h2_zr_ptbinned_closure_error[j-1]->GetNbinsX(); ++x)
        {
                for (int y = 1; y <= h2_zr_ptbinned_closure_error[j-1]->GetNbinsY(); ++y)
                {      
                double sumContent = h2_zr_ptbinned_closure_error[j-1]->GetBinContent(x, y);
                h2_zr_ptbinned_closure_error[j-1]->SetBinContent(x, y, TMath::Sqrt(sumContent) / TMath::Sqrt(nRuns));

                double sumError = h2_zr_ptbinned_closure_error[j-1]->GetBinError(x, y);
                h2_zr_ptbinned_closure_error[j-1]->SetBinError(x, y, TMath::Sqrt(sumError) / TMath::Sqrt(nRuns));
                }
        }
        h2_zr_ptbinned_closure_error[j-1]->Write();                
        }  
        }     
                        
        
        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        /////////////////////   Create jtr histograms /////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////

        std::cout << "############################## Unfolding for ptjtr ##############################" << std::endl;

        //TH2D *h2_jtr_truth = (TH2D *)file_truth->Get("jtr");
        TH3D *h3_ptjtr_truth = (TH3D *)file_truth->Get("ptjtr");
        TH3D *h3_ptjtr = (TH3D *)file_reco->Get("ptjtr");
        TH3D *h3_ptjtr_data = (TH3D *)file_data->Get("ptjtr");  
        TH3D *h3_ptjtr_final = (TH3D *)h3_ptjtr->Clone("ptjtr_final");
        TH3D *h3_eff_ptjtr = (TH3D *)file_unfold->Get("efficiency_ptjtr");
        TH3D *h3_purity_ptjtr = (TH3D *)file_unfold->Get("purity_ptjtr");
        RooUnfoldResponse *response_ptjtr = (RooUnfoldResponse *)file_unfold->Get("Roo_response_ptjtr");
        
        // Multiply by purity  
        h3_ptjtr_final->Multiply(h3_ptjtr_final, h3_purity_ptjtr);
        RooUnfoldBayes unfold_ptjtr(response_ptjtr, h3_ptjtr_final, NumIters);

        h3_ptjtr_final = (TH3D *)unfold_ptjtr.Hreco();
        // Divide by efficiency  
        h3_ptjtr_final->Divide(h3_ptjtr_final, h3_eff_ptjtr);

        h3_ptjtr_final->GetZaxis()->SetRangeUser(ptMin, ptMax);
        h3_ptjtr_truth->GetZaxis()->SetRangeUser(ptMin, ptMax);

        h3_ptjtr->GetZaxis()->SetRangeUser(ptMin, ptMax);

        TH2D *h2_jtr_final = (TH2D *)h3_ptjtr_final->Project3D("yx");
        h2_jtr_final->SetName("jtr_final");
        TH2D *h2_jtr = (TH2D *)h3_ptjtr->Project3D("yx");
        h2_jtr->SetName("jtr");
        TH2D *h2_jtr_truth = (TH2D *)h3_ptjtr_truth->Project3D("yx");
        h2_jtr_truth->SetName("jtr_truth");
        //
        h2_jtr_final->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_final->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_jtr_truth->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_truth->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        
        
        // This will depend on results of 1D jet unfolding
        //h2_jtr_final->Scale(1. / Njets_final, "width");
        //h2_jtr_truth->Scale(1. / Njets_truth, "width");
        NormalizeHist(h2_jtr_final);
        NormalizeHist(h2_jtr_truth);

        TH2D *h2_jtr_ratio = (TH2D *)h2_jtr->Clone("h2_jtr_ratio");
        TH2D *h2_jtr_ratio_final = (TH2D *)h2_jtr->Clone("h2_jtr_ratio_final");

        h2_jtr_ratio->Divide(h2_jtr, h2_jtr_truth);
        h2_jtr_ratio_final->Divide(h2_jtr_final, h2_jtr_truth);

        h2_jtr_ratio_final->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_ratio_final->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        
        h2_jtr_ratio_final->SetOption("COLZ, text");
        h2_jtr_ratio_final->Write();   

        ////////////////////////////////////
        // Smearing the Observables
        ///////////////////////////////////
        TH2D *h2_jtr_ptbinned_closure_error[ptbinsize-1];
        for (int i = 0; i < nRuns; i++)
        {
        TH3D *h3_ptjtr_smear = (TH3D *)h3_ptjtr->Clone(Form("ptjtr_smear%d", i));
        //myRNG->SetSeed();
        //myRNG->SetSeed(fixSmear);
        if (smear_by_data)
        {
                SmearObservables(h3_ptjtr_smear, h3_ptjtr_data, myRNG);
        }
        // Multiply by purity
        h3_ptjtr_smear->Multiply(h3_ptjtr_smear, h3_purity_ptjtr);
        RooUnfoldBayes unfold_ptjtr_smear(response_ptjtr, h3_ptjtr_smear, NumIters);

        h3_ptjtr_smear = (TH3D *)unfold_ptjtr_smear.Hreco();
        // Divide by efficiency
        h3_ptjtr_smear->Divide(h3_ptjtr_smear, h3_eff_ptjtr);
        h3_ptjtr_smear->GetZaxis()->SetRangeUser(ptMin, ptMax);

        TH2D *h2_jtr_smear = (TH2D *)h3_ptjtr_smear->Project3D("yx");
        h2_jtr_smear->SetName(Form("jtr_smear%d", i));

        h2_jtr_smear->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_smear->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        //h2_jtr_smear->Scale(1. / Njets_final, "width");
        NormalizeHist(h2_jtr_smear);
        TH2D *h2_jtr_ratio_smear = (TH2D *)h2_jtr_smear->Clone(Form("jtr_ratio_smear%d", i));
        //TH2D *h2_jtr_pull_smear = (TH2D *)h2_jtr_smear->Clone(Form("h2_jtr_pull_smear%d", i));
        h2_jtr_ratio_smear->Divide(h2_jtr_smear, h2_jtr_truth); 

        h2_jtr_ratio_smear->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_ratio_smear->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);

        //h2_jtr_pull_smear->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        //h2_jtr_pull_smear->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_jtr_ratio_smear->SetOption("COLZ, text");
        h2_jtr_ratio_smear->Write();      
        //h2_jtr_pull_smear->Write();     
        
        TH2D *h2_jtr_ptbinned[ptbinsize-1], *h2_jtr_truth_ptbinned[ptbinsize-1], *h2_jtr_ptbinned_ratio[ptbinsize-1], *h2_jtr_ptbinned_pullratio[ptbinsize-1];
        TH1D *h1_jtr_ptbinned_pulldist[ptbinsize-1];
        for (int j=1; j < ptbinsize; ++j)
        {
        h3_ptjtr_smear->GetZaxis()->SetRange(j+1, j+1);
        h3_ptjtr_truth->GetZaxis()->SetRange(j+1, j+1);        
        h2_jtr_ptbinned[j-1] = (TH2D *)h3_ptjtr_smear->Project3D("yx");    
        h2_jtr_ptbinned[j-1]->SetName(Form("jtr_smear%d_pt%d", i,j));
        h2_jtr_truth_ptbinned[j-1] = (TH2D *)h3_ptjtr_truth->Project3D("yx");
        h2_jtr_truth_ptbinned[j-1]->SetName(Form("jtr_truth_smear%d_pt%d", i,j));
        h2_jtr_ptbinned_pullratio[j-1] = (TH2D *)h2_jtr_ptbinned[j-1]->Clone(Form("jtr_pullratio_smear%d_pt%d", i, j));
        h1_jtr_ptbinned_pulldist[j-1] = new TH1D(Form("jtr_pulldist_smear%d_pt%d", i, j), ";pulls", 20, -4.5, 4.5);

        h2_jtr_ptbinned[j-1]->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_jtr_truth_ptbinned[j-1]->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);
        h2_jtr_ptbinned[j-1]->GetXaxis()->SetTitle("j_{T}");
        h2_jtr_truth_ptbinned[j-1]->GetXaxis()->SetTitle("j_{T}");
        h2_jtr_ptbinned[j-1]->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_truth_ptbinned[j-1]->GetYaxis()->SetRange(binlowr2D, binhighr2D);
        h2_jtr_ptbinned[j-1]->GetYaxis()->SetTitle("r");
        h2_jtr_truth_ptbinned[j-1]->GetYaxis()->SetTitle("r");
        h2_jtr_ptbinned_pullratio[j-1]->GetYaxis()->SetRange(binlowr2D, binhighr2D);      
        h2_jtr_ptbinned_pullratio[j-1]->GetXaxis()->SetRange(binlowjt2D, binhighjt2D);   
                
                
        //h2_jtr_ptbinned[j-1]->Scale(1./h1_jetpt_final->Integral(j+1, j+1));
        //h2_jtr_truth_ptbinned[j-1]->Scale(1./h1_jetpt_truth->Integral(j+1, j+1));
        NormalizeHist(h2_jtr_ptbinned[j-1]);
        NormalizeHist(h2_jtr_truth_ptbinned[j-1]);
        h2_jtr_ptbinned_ratio[j-1] = (TH2D *)h2_jtr_ptbinned[j-1]->Clone(Form("jtr_ratio_smear%d_pt%d", i,j));
        h2_jtr_ptbinned_ratio[j-1]->Divide(h2_jtr_ptbinned[j-1], h2_jtr_truth_ptbinned[j-1]);     
        h2_jtr_ptbinned_ratio[j-1]->SetOption("COLZ, text");        
        //h2_jtr_ptbinned_ratio[j-1]->Write(Form("zjt_ratio_smear%d_pt%d", i,j));
        h2_jtr_ptbinned_ratio[j-1]->Write();     
        
        GetPullsRatio(h2_jtr_ptbinned_ratio[j-1], h2_jtr_ptbinned_pullratio[j-1], h1_jtr_ptbinned_pulldist[j-1]);
        h2_jtr_ptbinned_pullratio[j-1]->SetOption("COLZ, text");
        h2_jtr_ptbinned_pullratio[j-1]->Write();
        h1_jtr_ptbinned_pulldist[j-1]->Write();
        
        sys_title = DoShapeClosure ? "shape nonclosure" : "nonclosure";
        sys_title += Form(" %.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]);          
        if (i == 0)
        {
                h2_jtr_ptbinned_closure_error[j-1] = (TH2D *)h2_jtr_ptbinned_ratio[j-1]->Clone(Form("jtr_pt%d_closure_error", j));
                if (!DoShapeClosure)
                {
                h2_jtr_ptbinned_closure_error[j-1]->Reset(); // Reset the histogram to zero
                }
        }
        h2_jtr_ptbinned_closure_error[j-1]->GetZaxis()->SetTitle("#sigma_{sys}");
        h2_jtr_ptbinned_closure_error[j-1]->SetTitle(sys_title);      
                
        if (DoShapeClosure)
        {
                SubtractUnity(h2_jtr_ptbinned_closure_error[j-1]);
                h2_jtr_ptbinned_closure_error[j-1]->Write();
        }
        else 
        {
                // Loop over all bins and add in quadrature
                for (int x = 1; x <= h2_jtr_ptbinned_ratio[j-1]->GetNbinsX(); ++x)
                {
                for (int y = 1; y <= h2_jtr_ptbinned_ratio[j-1]->GetNbinsY(); ++y)
                {    
                double binContent = 0;
                if ( h2_jtr_ptbinned_ratio[j-1]->GetBinContent(x,y) != 0 )
                { 
                binContent = 1. - h2_jtr_ptbinned_ratio[j-1]->GetBinContent(x, y);
                }
                double sumContent = h2_jtr_ptbinned_closure_error[j-1]->GetBinContent(x, y);
                h2_jtr_ptbinned_closure_error[j-1]->SetBinContent(x, y, sumContent + binContent * binContent);

                double binError = h2_jtr_ptbinned_ratio[j-1]->GetBinError(x, y);
                double sumError = h2_jtr_ptbinned_closure_error[j-1]->GetBinError(x, y);
                h2_jtr_ptbinned_closure_error[j-1]->SetBinError(x, y, sumError + binError * binError);
                }
                }  
        }            
        }                        
        }  
        if (nRuns > 1)
        {
        for (int j=1; j < ptbinsize; ++j)
        {
        // Take the square root of the final sum histogram to get the quadrature sum
        for (int x = 1; x <= h2_jtr_ptbinned_closure_error[j-1]->GetNbinsX(); ++x)
        {
                for (int y = 1; y <= h2_jtr_ptbinned_closure_error[j-1]->GetNbinsY(); ++y)
                {      
                double sumContent = h2_jtr_ptbinned_closure_error[j-1]->GetBinContent(x, y);
                h2_jtr_ptbinned_closure_error[j-1]->SetBinContent(x, y, TMath::Sqrt(sumContent) / TMath::Sqrt(nRuns));

                double sumError = h2_jtr_ptbinned_closure_error[j-1]->GetBinError(x, y);
                h2_jtr_ptbinned_closure_error[j-1]->SetBinError(x, y, TMath::Sqrt(sumError) / TMath::Sqrt(nRuns));
                }
        }  
        h2_jtr_ptbinned_closure_error[j-1]->Write();                
        }  
        }   
        
        //file_write->Write();
        file_write->Close();
        //
}
