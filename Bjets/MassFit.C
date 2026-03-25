#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"

#include <sstream>
#include <fstream>
#include <iostream>

#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/TBJetsMCReco.h"
#include "../include/TBJetsMCReco.C"

using namespace RooFit;

void MassFit(int NumEvts = -1, int dataset = 91599, bool isData = true,
             bool UseDTF = true,
             bool DoRecSelEff = 0,
             bool DoSystematic = 0,             
             float ptmin_user = pTLow,
             float ptmax_user = 250.)
{
        ptMin = ptmin_user;
        ptMax = ptmax_user;
        
        int flavor = 5;
        
        TString str_Mag = "";
        TString str_pT = "";
        TString str_level = "";
        TString str_followHard = "";
        TString str_flavor = "";
        TString str_ghost = "";
        TString str_charged = "";
        TString str_PID = "";
        TString str_year = "";

        if (PID_cut)
                str_PID = "_PID";

        if (flavor == 1)
                str_flavor = "_udsg";
        else if (flavor == 4)
                str_flavor = "_c";
        else if (flavor == 5)
                str_flavor = "_b";
        
                if (isData) {
                str_level = "data";
        } else {
                if (truthLevel)
                        str_level = "truth";
                else
                        str_level = "reco";
        }

        if (flavor == 1) {
                str_followHard = "_hard";
        } else {
                if (followHardest)
                        str_followHard = "_hard";
                else
                        str_followHard = "_HF";
        }

        if (ghostCut)
                str_ghost = Form("_ghost_%.1f", ghostProb);

        if (chargedJetCut)
                str_charged = "_charge";


        TString extension = TString("massfit_") + str_level + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_ghost + str_charged + str_Mag + str_flavor + str_PID + str_L0 + Form("_%d", dataset);
        
        if (DoRecSelEff)
                extension = "recselsys_" + extension;
        
        if (DoSystematic)
                extension = "sys_" + extension;
        
        cout << extension << endl;
        
        // Setup Tree
        TChain *BTree = new TChain("BTree", "B-jets Tree Variables");

        if (isData)
                BTree->Add((output_folder + "ntuple_bjets_data.root/BTree").c_str());
        else
                BTree->Add((output_folder + "ntuple_bjets_mcreco.root/BTree").c_str());
        
        if (NumEvts > BTree->GetEntries())
                NumEvts = BTree->GetEntries();
        
        if (NumEvts == -1)
                NumEvts = BTree->GetEntries();

        cout << "Number of entries : " << BTree->GetEntries() << endl;
        
        float mup_px, mup_py, mup_pz, mup_e;
        float mum_px, mum_py, mum_pz, mum_e;
        float K_px, K_py, K_pz, K_e, K_PIDK;
        float jet_pt, bmass_dtf;
        float HF_px, HF_py, HF_pz, HF_e;
        float tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e;
        float Jpsi_CHI2NDOF, Jpsi_BPVDLS, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_IPCHI2;
        float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
        float jet_eta, jet_rap;
        double WTA_dist;
        bool isTrueBjet, TOS;
        int nSV;
        
        BTree->SetBranchAddress("HF_px", &HF_px);
        BTree->SetBranchAddress("HF_py", &HF_py);
        BTree->SetBranchAddress("HF_pz", &HF_pz);
        BTree->SetBranchAddress("HF_e", &HF_e);

        BTree->SetBranchAddress("Bu_IPCHI2", &Bu_IPCHI2);
        BTree->SetBranchAddress("Bu_CHI2", &Bu_CHI2);
        BTree->SetBranchAddress("Bu_CHI2NDOF", &Bu_CHI2NDOF);

        BTree->SetBranchAddress("Jpsi_FDCHI2", &Jpsi_FDCHI2);
        BTree->SetBranchAddress("Jpsi_CHI2", &Jpsi_CHI2);
        BTree->SetBranchAddress("Jpsi_CHI2NDOF", &Jpsi_CHI2NDOF);
        BTree->SetBranchAddress("Jpsi_BPVDLS", &Jpsi_BPVDLS);    

        BTree->SetBranchAddress("mum_px", &mum_px);
        BTree->SetBranchAddress("mum_py", &mum_py);
        BTree->SetBranchAddress("mum_pz", &mum_pz);
        BTree->SetBranchAddress("mum_e", &mum_e);

        BTree->SetBranchAddress("mup_px", &mup_px);
        BTree->SetBranchAddress("mup_py", &mup_py);
        BTree->SetBranchAddress("mup_pz", &mup_pz);
        BTree->SetBranchAddress("mup_e", &mup_e);

        BTree->SetBranchAddress("K_px", &K_px);
        BTree->SetBranchAddress("K_py", &K_py);
        BTree->SetBranchAddress("K_pz", &K_pz);
        BTree->SetBranchAddress("K_e", &K_e);
        BTree->SetBranchAddress("K_PIDK", &K_PIDK);

        BTree->SetBranchAddress("jet_pt", &jet_pt);
        BTree->SetBranchAddress("jet_eta", &jet_eta);
        BTree->SetBranchAddress("jet_rap", &jet_rap);    

        BTree->SetBranchAddress("nSV", &nSV);
        BTree->SetBranchAddress("bmass_dtf", &bmass_dtf);
        BTree->SetBranchAddress("isTrueBjet", &isTrueBjet);

        BTree->SetBranchAddress("TOS", &TOS);

        BTree->SetBranchAddress("WTA_reco_dist", &WTA_dist);

        if (!isData) {
                BTree->SetBranchAddress("tr_HF_px", &tr_HF_px);
                BTree->SetBranchAddress("tr_HF_py", &tr_HF_py);
                BTree->SetBranchAddress("tr_HF_pz", &tr_HF_pz);
                BTree->SetBranchAddress("tr_HF_e", &tr_HF_e);
        }
        
        float mass_low = 5.15;
        float mass_high = 5.55;

        vector<TH1D *> h1_mass_HFpt;
        vector<float> vec_bkg_frac, vec_res_frac, vec_bkg_yield, vec_sig_yield;

        std::string type_of_event = (isData) ? "data" : "mcreco";
        TFile f((output_folder + "mass-fits/results_mass_fit_" + type_of_event + ".root").c_str(), "RECREATE");

        int nBins = 80;
        float binsize;
        // int nBins[ptHFbinsize] = {100, 100, 100, 100, 100, 100};
        for (int i = 0; i < ptHFbinsize; i++) {
                if (i >= ptHFbinsize - 2)
                        nBins = 60;
                else
                        nBins = 80;
                
                TH1D *h1_temp = new TH1D(Form("h1_mass%d", i), "", nBins, mass_low, mass_high);
                
                h1_mass_HFpt.push_back(h1_temp);
                
                binsize = (mass_high - mass_low) / (float)nBins;
        }

        TH1D *h1_mass     = new TH1D("h1_mass", "", nBins, mass_low, mass_high);
        TH1D *h1_mass_all = new TH1D("h1_mass_all", "", nBins + 17, 4.8, mass_high);
        
        TH1D *h1_MassMin  = new TH1D("h1_MassMin", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_MassMax  = new TH1D("h1_MassMax", "", ptHFbinsize, ptHF_binedges);
        
        TH1D *h1_Sideband1Min = new TH1D("h1_Sideband1Min", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband1Max = new TH1D("h1_Sideband1Max", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband2Min = new TH1D("h1_Sideband2Min", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband2Max = new TH1D("h1_Sideband2Max", "", ptHFbinsize, ptHF_binedges);
        
        TH1D *h1_Sideband1Min_forSysNear = new TH1D("h1_Sideband1Min_forSysNear", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband1Max_forSysNear = new TH1D("h1_Sideband1Max_forSysNear", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband2Min_forSysNear = new TH1D("h1_Sideband2Min_forSysNear", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband2Max_forSysNear = new TH1D("h1_Sideband2Max_forSysNear", "", ptHFbinsize, ptHF_binedges);
        
        TH1D *h1_Sideband1Min_forSysFar = new TH1D("h1_Sideband1Min_forSysFar", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband1Max_forSysFar = new TH1D("h1_Sideband1Max_forSysFar", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband2Min_forSysFar = new TH1D("h1_Sideband2Min_forSysFar", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_Sideband2Max_forSysFar = new TH1D("h1_Sideband2Max_forSysFar", "", ptHFbinsize, ptHF_binedges);    
        
        TH1D *h1_BkgScale = new TH1D("h1_BkgScale", "", ptHFbinsize, ptHF_binedges);
        
        TH1D *h1_BkgScale_forSysNear  = new TH1D("h1_BkgScale_forSysNear", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_BkgScale_forSysFar   = new TH1D("h1_BkgScale_forSysFar", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_BkgScale_forSysLower = new TH1D("h1_BkgScale_forSysLower", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_BkgScale_forSysUpper = new TH1D("h1_BkgScale_forSysUpper", "", ptHFbinsize, ptHF_binedges);
        
        TLorentzVector HFmeson, HFjet, mum, mup, Kmeson, Jpsi, tr_HFmeson;

        // Fill the mass plots
        for (int ev = 0; ev < NumEvts; ev++) {
                BTree->GetEntry(ev);

                if (ev%10000 == 0) {
                        double percentage = 100.*ev/NumEvts;
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (jet_rap < etaMin || jet_rap > etaMax)
                        continue;

                if (jet_pt > ptMax)
                        continue;

                if (jet_pt < ptMin)
                        continue;
                
                if (!TOS)
                        continue;
                
                if (WTA_dist > WTA_dist_max)
                        continue;

                if (DoRecSelEff) {
                        // cout << Bu_IPCHI2 << ", " << Bu_CHI2 << ", " << Jpsi_CHI2 << ", " << Jpsi_CHI2NDOF << ", " << sqrt(Jpsi_FDCHI2) << endl;
                        if (Bu_IPCHI2 > 22)
                                continue;
                        if (Bu_CHI2 > 42)
                                continue;
                        if (Jpsi_CHI2 > 22)
                                continue;
                        if (Jpsi_CHI2NDOF > 18)
                                continue;
                        if (fabs(Jpsi_BPVDLS) < 3.2)
                                continue;
                }

                mup.SetPxPyPzE(mup_px, 
                               mup_py, 
                               mup_pz, 
                               mup_e);

                mum.SetPxPyPzE(mum_px, 
                               mum_py, 
                               mum_pz, 
                               mum_e);

                Kmeson.SetPxPyPzE(K_px, 
                                  K_py, 
                                  K_pz, 
                                  K_e);
                
                // HFmeson = mup + mum + Kmeson;
                
                HFmeson.SetPxPyPzE(HF_px, 
                                   HF_py, 
                                   HF_pz, 
                                   HF_e);

                tr_HFmeson.SetPxPyPzE(tr_HF_px, 
                                      tr_HF_py, 
                                      tr_HF_pz, 
                                      tr_HF_e);

                // if(dtf_chi2ndf > 5) continue;
                // cout<<bmass_dtf<<endl;
                float bmass = HFmeson.M();
                
                if (UseDTF)
                        bmass = bmass_dtf;
                
                bool PID_cond = (K_PIDK > 0);
                
                if (!PID_cond && PID_cut)
                        continue;
                
                if (!isData)
                        if (!isTrueBjet)
                                continue;

                h1_mass->Fill(bmass);
                h1_mass_all->Fill(bmass);

                for (int i = 0; i < ptHFbinsize; i++) {
                        if (HFmeson.Pt() > ptHF_binedges[i] && HFmeson.Pt() < ptHF_binedges[i + 1]) {
                                h1_mass_HFpt[i]->Fill(bmass);

                                break;
                        }
                }
        }

        //---- paint setup...
        //
        int ican = -1, iframe = -1, itext = -1;
        TCanvas *ccan[1000];
        TH1F *frame[1000];
        TLatex *text[1000];
        TRatioPlot *rp[1000];

        for (int i = 0; i < 1000; i++) {
                text[i] = new TLatex();
                text[i]->SetNDC(kTRUE);
                text[i]->SetTextSize(0.06);
        }

        TLatex Tl;
        Tl.SetNDC(kTRUE);
        Tl.SetTextSize(0.03);
        //
        gStyle->SetCanvasDefH(2000);
        gStyle->SetCanvasDefW(3508);
        //
        gStyle->SetOptStat(0);
        // gStyle->SetPaperSize(TStyle::kUSLetter);
        gStyle->SetPadBottomMargin(0.2);
        // gStyle->SetPadTopMargin(0.005);
        gStyle->SetPadLeftMargin(0.15);
        gStyle->SetPadRightMargin(0.15);
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
        gStyle->SetPaintTextFormat("4.3g");
        // gStyle->SetErrorX(0);
        gStyle->SetTitleStyle(0);
        gStyle->SetStatStyle(0);
        // gStyle->SetLineWidth(3);

        //---- paint...
        char buf[100];
        char bufb[100];
        
        TString plotfilePDF;
        TString plotfileO;
        TString plotfileC;
        
        plotfilePDF = Form("./plots/massfit_%s.pdf",(isData)?"data":"mcreco");

        plotfileO = plotfilePDF + TString("(");
        plotfileC = plotfilePDF + TString("]");
        
        // Plot data and PDF overlaid
        ican = 0;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 1);
        ccan[ican]->SetFillColor(10);

        // gPad->SetBottomMargin(0.06);
        ccan[ican]->cd();
        ccan[ican]->Divide(3, 5, 0.0001, 0.0001);

        ican = 1;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 1);
        ccan[ican]->SetFillColor(10);

        // gPad->SetBottomMargin(0.06);
        ccan[ican]->cd();
        ccan[ican]->Divide(3, 5, 0.0001, 0.0001);
        // c->SaveAs("plots/"+extension+".pdf");
        // h1_mass->Scale(1./h1_mass->GetEntries());

        std::ofstream SigParams("csv/signal_" + extension + ".csv");
        std::ofstream BkgParams("csv/bkg_" + extension + ".csv");
        std::ofstream ResonantParams("csv/res_" + extension + ".csv");

        if (DoSystematic)
                SigParams << "m,t,n,nsig,fraction" << endl;
        else
                SigParams << "mu,sigma1,sigma2,alpha1,n1,n1_2,alpha2,n2,n2_2,nsig" << endl;

        BkgParams << "exp,nbkg" << endl;
        ResonantParams << "mu,sigma,alpha1,n1,alpha2,n2,nres" << endl;

        // S e t u p   c o m p o n e n t   p d f s
        // ---------------------------------------

        // Declare observable x
        std::string latexFileName = "fit_results.tex";
        std::ofstream latexFile(latexFileName);

        TFile *file_workspace;
        TFile *file_workspace_misid;
        
        for (int i = 0; i < ptHFbinsize; i++) {
                if (i >= ptHFbinsize - 2)
                        nBins = 60;
                else
                        nBins = 80;
                
                binsize = (mass_high - mass_low) / (float)nBins;

                // file_workspace works when the mass fits are executed for the MCReco first!

                file_workspace       = new TFile((output_folder + Form("mass-fits/workspace%d_", i) + "testing_massfitsmcreco.root").c_str(), "READ");
                file_workspace_misid = new TFile((output_folder + Form("mass-fits/workspace%d_", i) + "testing_massfits.root").c_str(), "READ");
                
                if (file_workspace == NULL) {
                        std::cout<<"First execute the massfits for mcreco!"<<std::endl;

                        return;
                }
                
                if (file_workspace_misid == NULL) {
                        std::cout<<"First execute the massfits for MisID!"<<std::endl;

                        return;
                }

                RooWorkspace *w_read = (RooWorkspace *)file_workspace->Get(Form("w%d", i));
                RooWorkspace *w_read_misid = (RooWorkspace *)file_workspace_misid->Get(Form("w%d", i));
                RooRealVar *sigma_ratio, *mean, *sigma1, *sigma2;
                RooRealVar *mu_sig, *alpha1_sig, *alpha2_sig, *p1_sig, *p2_sig;
                RooRealVar *mu_sig2, *alpha1_sig2, *alpha2_sig2, *p1_sig2, *p2_sig2;        
                RooRealVar *mu, *width, *alpha1, *alpha2, *p1, *p2;

                RooRealVar HFMass("HFMass", "HFMass", mass_low, mass_high);
                RooDataHist B_mass("B_mass", "B_mass", RooArgList(HFMass), h1_mass_HFpt[i], 1.);

                HFMass.setBins(nBins);

                // Create two Gaussian PDFs g1(HFMass,mean1,sigma) anf g2(HFMass,mean2,sigma) and their parameters
                ////////////////////////////////////////////////////
                // build siganl
                ///////////////////////////////////////////////////

                if (UseDTF) {
                        sigma_ratio = new RooRealVar("sigma_ratio", "sigma_ratio", 1.52571 / 0.856549);
                        mean        = new RooRealVar("mean", "mean of gaussians", 5.27966, 5.27, 5.282);
                        
                        if (isData) {
                                sigma1 = new RooRealVar("sigma1", "width of gaussians", 0.011, 0.001, 0.03);
                                sigma2 = new RooRealVar("sigma2", "width of gaussians", 0.007, 0.001, 0.03); // CHANGE
                        } else {
                                sigma1 = new RooRealVar("sigma1", "width of gaussians", 0.011, 0.001, 0.03);
                                sigma2 = new RooRealVar("sigma2", "width of gaussians", 0.007, 0.001, 0.03); // CHANGE
                        }     
                } else {
                        sigma_ratio = new RooRealVar("sigma_ratio", "sigma_ratio", 1.78528 / 3.15958);
                        mean        = new RooRealVar("mean", "mean of gaussians", 5.279, 5.27, 5.282);
                        sigma1      = new RooRealVar("sigma1", "width of gaussians", 0.002, 0.001, 0.05);
                        sigma2      = new RooRealVar("sigma2", "width of gaussians", 0.003, 0.001, 0.05); // CHANGE
                }

                RooGaussian gauss1("gauss_sig1", "Signal component 1", HFMass, *mean, *sigma1);
                RooGaussian gauss2("gauss_sig2", "Signal component 2", HFMass, *mean, *sigma2);

                mu_sig     = new RooRealVar("mu_sig", "mu", 5.279, 5.27, 5.282);
                alpha1_sig = new RooRealVar("alpha1_sig", "alpha1", 2., 0.01, 10.);
                alpha2_sig = new RooRealVar("alpha2_sig", "alpha2", 2., 0.01, 10.);

                alpha1_sig2 = new RooRealVar("alpha1_sig2", "alpha1", 2., 0.01, 10.);
                alpha2_sig2 = new RooRealVar("alpha2_sig2", "alpha2", 2., 0.01, 10.);
                
                // p1_sig = RooRealVar("p1_sig", "p1", 1., 0.6, 4.);
                // p2_sig = RooRealVar("p2_sig", "p2", 1., 0.6, 4.);
                
                p1_sig = new RooRealVar("p1_sig", "p1", 2., 1., 6.);
                p2_sig = new RooRealVar("p2_sig", "p2", 3., 1., 6.);

                // Or, Create signal from two Crystal Ball functions
                // These parameters have been derived from simulation


                if (isData && w_read != NULL) {
                        mu_sig = (RooRealVar *)w_read->obj("mu_sig");
                        alpha1_sig = (RooRealVar *)w_read->obj("alpha1_sig");
                        alpha2_sig = (RooRealVar *)w_read->obj("alpha2_sig");
                        p1_sig = (RooRealVar *)w_read->obj("p1_sig");
                        p2_sig = (RooRealVar *)w_read->obj("p2_sig");

                        // p1_sig2 = (RooRealVar *)w_read->obj("p1_sig2");
                        // p2_sig2 = (RooRealVar *)w_read->obj("p2_sig2");
                        // mu_sig.setConstant(kTRUE);
                        alpha1_sig->setConstant(kTRUE);
                        alpha2_sig->setConstant(kTRUE);
                        p1_sig->setConstant(kTRUE);
                        p2_sig->setConstant(kTRUE);
                        // p1_sig2->setConstant(kTRUE);
                        // p2_sig2->setConstant(kTRUE);
                }

                // Define the CB functions
                RooCrystalBall dcbPdf_sig1("dcbPdf_sig1", "floatSidedCB_sig1", HFMass, *mu_sig, *sigma1, *alpha1_sig, *p1_sig, *alpha2_sig, *p2_sig);
                RooCrystalBall dcbPdf_sig2("dcbPdf_sig2", "floatSidedCB_sig2", HFMass, *mu_sig, *sigma2, *alpha1_sig, *p1_sig, *alpha2_sig, *p2_sig);

                RooRealVar *nsig1 = new RooRealVar("nsig1", "fraction of component 1 in signal", 0.3, 0., 1.);

                // Choose if you want two CB functions as signal, or two Gaussians.

                //RooAddPdf sig("sig", "Signal", RooArgList(gauss1, gauss2), RooArgList(nsig1));
                RooAddPdf sig("sig", "Signal", RooArgList(dcbPdf_sig1, dcbPdf_sig2), RooArgList(*nsig1));
                
                RooRealVar *sig2_m = new RooRealVar("sig2_m", "sig2_m", 5.279, 5.27, 5.282);
                RooRealVar *sig2_n = new RooRealVar("sig2_n", "sig2_n", 5.);
                sig2_n->setConstant(kTRUE);
                RooRealVar *sig2_t = new RooRealVar("sig2_t", "sig2_t", 0.02, 0.001, 0.3);

                RooGenericPdf sig2("student-t", "", "tgamma((sig2_n+1)/2)/(tgamma(sig2_n/2)*sqrt(TMath::Pi()*sig2_n)*sig2_t) * pow(1+(1/sig2_n * pow((HFMass-sig2_m)/sig2_t, 2)), (-1-sig2_n)/2)", RooArgList(HFMass, *sig2_m, *sig2_n, *sig2_t));

                if (UseDTF && w_read_misid != NULL) {
                        mu = (RooRealVar *)w_read_misid->obj("mu");
                        width = (RooRealVar *)w_read_misid->obj("width");
                        alpha1 = (RooRealVar *)w_read_misid->obj("alpha1");
                        alpha2 = (RooRealVar *)w_read_misid->obj("alpha2");
                        p1 = (RooRealVar *)w_read_misid->obj("p1");
                        p2 = (RooRealVar *)w_read_misid->obj("p2");
                        mu->setConstant(kTRUE);
                        width->setConstant(kTRUE);
                        alpha1->setConstant(kTRUE);
                        alpha2->setConstant(kTRUE);
                        p1->setConstant(kTRUE);
                        p2->setConstant(kTRUE);
                } else {
                        return; // EFMC: sorry, but wtf
                        mu = new RooRealVar("mu", "mu", 5.33183e+00);
                        width = new RooRealVar("width", "width", 2.62897e-02);
                        alpha1 = new RooRealVar("alpha1", "alpha1", 1.72233e+00);
                        alpha2 = new RooRealVar("alpha2", "alpha2", 7.67677e-01);
                        p1 = new RooRealVar("p1", "p1", 1.73799e+00);
                        p2 = new RooRealVar("p2", "p2", 2.09840e+00);
                }

                RooCrystalBall dcbPdf("dcbPdf", "floatSidedCB", HFMass, *mu, *width, *alpha1, *p1, *alpha2, *p2);

                /////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Build linear background PDF, for combinatoric background
                ///////////////////////////////////////////////////

                RooRealVar *a2 = new RooRealVar("c2", "coefficient #2", -0.1, -1., -0.001);
                RooRealVar *a1 = new RooRealVar("exp_c", "exp_c", -2., -4, -0.1);
                RooExponential bkg("bkg", "bkg", HFMass, *a1);

                //  RooRealVar * a2("a2", "a2", -1, -2, -0.1);
                // RooRealVar * a1("a1", "a1", 2, 0., 9.);
                // RooPolynomial bkg("bkg", "Background", HFMass, RooArgList(a1, a2));

                RooRealVar *a0_cheb = new RooRealVar("a0_cheb", "a0_cheb", -0.5, -3, -0.0001);
                RooRealVar *a1_cheb = new RooRealVar("a1_cheb", "a1_cheb", -0.2, -3, -0.0001);
                RooRealVar *a2_cheb = new RooRealVar("a2_cheb", "a2_cheb", -0.2, -3, -0.0001);
                // RooChebychev bkg_cheb("bkg_cheb", "Background_cheb", HFMass, RooArgSet(*a0_cheb, *a1_cheb, *a2_cheb));
                RooChebychev bkg_cheb("bkg_cheb", "Background_cheb", HFMass, RooArgSet(*a0_cheb, *a1_cheb));

                /////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Fractions
                ///////////////////////////////////////////////////

                RooRealVar *nsig = new RooRealVar("nsig", "normalization for signal", 10000, 0., 10000000.);
                RooRealVar *nbkg = new RooRealVar("nbkg", "normalization for background", 1000, 0., 1000000);
                RooRealVar *nres = new RooRealVar("nres", "normalization for  background", 0.00384 * nsig->getVal(), 0., 0.0384 * nsig->getVal());

                // S a m p l e ,   f i t   a n d   p l o t   m o d e l
                // ---------------------------------------------------
                HFMass.setRange("comb1", 5.55, 5.78);
                HFMass.setRange("comb2", RooConst(mass_low), RooConst(mass_low + 0.1));
                HFMass.setRange("Full", mass_low, mass_high);

                // Fit model to data
                RooAddPdf *model;
                if (isData)
                {
                        if (UseDTF){
                                if (DoSystematic)
                                        model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg_cheb, sig2, dcbPdf), RooArgList(*nbkg, *nsig, *nres));
                                else
                                        model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg_cheb, sig, dcbPdf), RooArgList(*nbkg, *nsig, *nres));
                        } else {
                                model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg, sig, dcbPdf), RooArgList(*nbkg, *nsig, *nres));
                        }
                } else {
                model = new RooAddPdf("model", "g1+g2+a", RooArgList(sig), RooArgList(*nsig));
                // model = new RooAddPdf("model", "g1+g2+a", RooArgList(sig2), RooArgList(*nsig));
                }

                // Capture LaTeX output
                RooFitResult* result =  model->fitTo(B_mass, PrintEvalErrors(-1), Save(true));
                RooArgSet params_final = result->floatParsFinal();
                std::streambuf* oldCoutBuffer = std::cout.rdbuf();
                std::cout.rdbuf(latexFile.rdbuf());  // Redirect cout to outFile
                latexFile << "\\begin{table}[H]\n";
                latexFile << "    \\centering\n";
                latexFile << "    \\caption{Fit result for " << ptHF_binedges[i] << " $ < \\pthf < $ " << ptHF_binedges[i+1] << "}\n";
                latexFile << "    \\label{tab:fit_results_pthf" << i+1 << "}\n";
                params_final.printLatex();
                latexFile << "\\end{table}\n";
                std::cout.rdbuf(oldCoutBuffer);


                // Begin plotting
                //
                ican = 0;
                ccan[ican]->cd(i + 1);
                gPad->SetRightMargin(0.12);
                gPad->SetLeftMargin(0.12);
                ccan[ican]->SetFillStyle(4070);
                RooPlot *xframe = HFMass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
                xframe->SetTitleOffset(0.7, "Y");
                xframe->SetTitleOffset(0.9, "X");

                B_mass.plotOn(xframe, Name("B_mass"));
                // model->plotOn(xframe);

                // Overlay the background component of model with a dashed line
                model->plotOn(xframe, Name("sig"), Components(sig), LineStyle(5), LineColor(kGreen + 2), LineWidth(1.));
                model->plotOn(xframe, Name("bkg_cheb"), Components(bkg_cheb), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
                model->plotOn(xframe, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));
                model->plotOn(xframe, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1.));
                RooHist *hpull = xframe->pullHist();
                float chi2;
                if (isData)
                chi2 = xframe->chiSquare(10);
                else
                chi2 = xframe->chiSquare(10);


                // Overlay the background+sig2 components of model with a dotted line
                xframe->Draw("SAME");
                gPad->SetLogy();
                xframe->SetMinimum(1.);

                TLegend *leg1 = new TLegend(0.6, 0.5, 0.7, 0.87);
                leg1->SetTextSize(0.05);
                leg1->SetBorderSize(0);
                leg1->AddEntry("B_mass", "Data", "PE");
                leg1->AddEntry("model", "Fit model", "L");
                leg1->AddEntry("sig", "Signal", "L");
                leg1->AddEntry("bkg", "Comb bkg", "L");
                
                if (isData)
                        leg1->AddEntry("dcbPdf", "misID bkg", "F");
                
                leg1->AddEntry((TObject *)0, Form("#chi^{2}/ndf = %.2f", chi2), "");
                leg1->AddEntry((TObject *)0, Form("%.1f < p_{T}^{HF} < %.1f GeV", ptHF_binedges[i], ptHF_binedges[i + 1]), "");

                leg1->Draw("SAME");

                ican = 1;
                ccan[ican]->cd(i + 1);
                gPad->SetRightMargin(0.12);
                gPad->SetLeftMargin(0.12);
                ccan[ican]->SetFillStyle(4070);
                RooPlot *xframe2 = HFMass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
                xframe2->SetTitleOffset(0.7, "Y");
                xframe2->SetTitleOffset(0.9, "X");
                B_mass.plotOn(xframe2, Name("B_mass"));
                // model->plotOn(xframe2);

                // Overlay the background component of model with a dashed line
                model->plotOn(xframe2, Name("sig"), Components(sig), LineStyle(5), LineColor(kGreen + 2), LineWidth(1.));
                model->plotOn(xframe2, Name("bkg_cheb"), Components(bkg_cheb), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
                model->plotOn(xframe2, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));
                model->plotOn(xframe2, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1.));

                xframe2->Draw("SAME");
                xframe2->SetMinimum(0.);
                leg1->Draw("SAME");

                RooWorkspace *w = new RooWorkspace(Form("w%d", i), Form("workspace%d", i));
                w->import(*model);
                w->import(B_mass);
                w->Print();

                if (isData)
                        w->writeToFile((output_folder + Form("mass-fits/workspace%d_", i) + "testing_massfitsdata.root").c_str());
                else
                        w->writeToFile((output_folder + Form("mass-fits/workspace%d_", i) + "testing_massfitsmcreco.root").c_str());

                cout << "Chi2/dof = " << chi2 << endl;
                cout << "a1 = " << a1->getVal() << endl;
                cout << "a2 = " << a2->getVal() << endl;
                cout << "sigma1 = " << sigma1->getVal() << endl;
                cout << "sigma2 = " << sigma2->getVal() << endl;

                float MassMu = mu_sig->getVal();
                float MassSigma;
                
                if (nsig1->getVal() > 0.5)
                        MassSigma = sigma1->getVal();
                else
                        MassSigma = sigma2->getVal();     
                
                if (DoSystematic) {
                        MassMu = sig2_m->getVal();
                        MassSigma = sig2_t->getVal();
                        cout << "MassMu = " << sig2_m->getVal() << endl;
                        cout << "MassSigma = " << sig2_t->getVal() << endl;
                }
                
                /*
                float MassUp = MassMu + 3 * MassSigma;
                float MassDown = MassMu - 3 * MassSigma;
                Sideband1_Min = MassMu - 11 * MassSigma;
                Sideband1_Max = MassMu - 7 * MassSigma;
                Sideband2_Min = MassMu + 10 * MassSigma;
                Sideband2_Max = MassMu + 14 * MassSigma;
                */
                float MassUp, MassDown;
                if (i < 2) {
                        MassDown = MassMu - 2 * MassSigma;
                        MassUp = MassMu + 2 * MassSigma;
                }
                
                if (i < 5) {
                        MassDown = MassMu - 3 * MassSigma;
                        MassUp = MassMu + 3 * MassSigma;
                } else {
                        MassDown = MassMu - 4 * MassSigma;
                        MassUp = MassMu + 4 * MassSigma;
                }
                
                Sideband1_Min = (MassMu - 14 * MassSigma) > 5.15 ? (MassMu - 14 * MassSigma) : 5.15;
                Sideband1_Max = (MassMu - 7 * MassSigma) > Sideband1_Min ? (MassMu - 7 * MassSigma) : 5.17;

                Sideband2_Min = MassMu + 10 * MassSigma;
                Sideband2_Max = (MassMu + 25 * MassSigma) < 5.55 ? (MassMu + 25 * MassSigma) : 5.55; // Biggest difference wrt IC code
                
                float Sideband1_Min_forSysNear = Sideband1_Min + MassSigma;
                float Sideband1_Max_forSysNear = Sideband1_Max;
                float Sideband2_Min_forSysNear = Sideband2_Min;
                float Sideband2_Max_forSysNear = Sideband2_Max - MassSigma;

                float Sideband1_Min_forSysFar = Sideband1_Min;
                float Sideband1_Max_forSysFar = Sideband1_Max - MassSigma;
                float Sideband2_Min_forSysFar = Sideband2_Min + MassSigma;
                float Sideband2_Max_forSysFar = Sideband2_Max;
                
                cout << "Signal: " << MassDown << ", " << MassUp << endl;
                cout << "Sideband1: " << Sideband1_Min << ", " << Sideband1_Max << endl;
                cout << "Sideband2: " << Sideband2_Min << ", " << Sideband2_Max << endl;
                cout << "Sideband1_forSysNear: " << Sideband1_Min_forSysNear << ", " << Sideband1_Max_forSysNear << endl;
                cout << "Sideband2_forSysNear: " << Sideband2_Min_forSysNear << ", " << Sideband2_Max_forSysNear << endl;
                cout << "Sideband1_forSysFar: " << Sideband1_Min_forSysFar << ", " << Sideband1_Max_forSysFar << endl;
                cout << "Sideband2_forSysFar: " << Sideband2_Min_forSysFar << ", " << Sideband2_Max_forSysFar << endl;
                
                HFMass.setRange("signal", MassDown, MassUp);
                HFMass.setRange("sideband1", Sideband1_Min, Sideband1_Max);
                HFMass.setRange("sideband2", Sideband2_Min, Sideband2_Max);

                HFMass.setRange("sideband1_forSysNear", Sideband1_Min_forSysNear, Sideband1_Max_forSysNear);
                HFMass.setRange("sideband2_forSysNear", Sideband2_Min_forSysNear, Sideband2_Max_forSysNear);
                HFMass.setRange("sideband1_forSysFar", Sideband1_Min_forSysFar, Sideband1_Max_forSysFar);
                HFMass.setRange("sideband2_forSysFar", Sideband2_Min_forSysFar, Sideband2_Max_forSysFar);
                
                if (Sideband2_Min_forSysNear < Sideband1_Max_forSysNear)
                        std::cout << "Overlapping sideband ranges!!" << std::endl;
                
                RooAbsReal *intsigX;
                RooAbsReal *intbkgX, *intresX, *intbkgSide, *intbkgSideLeft, *intbkgSideRight;
                intsigX = sig.createIntegral(HFMass, NormSet(HFMass), Range("signal"));
                intbkgX = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("signal"));

                //intbkgSide = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband1,sideband2"));
                intbkgSideLeft = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband1"));
                intbkgSideRight = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband2"));
                //RooAbsReal *intbkgSide_forSysNear = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband1_forSysNear, sideband2_forSysNear"));
                RooAbsReal *intbkgSideLeft_forSysNear = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband1_forSysNear"));
                RooAbsReal *intbkgSideRight_forSysNear = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband2_forSysNear"));
                //RooAbsReal *intbkgSide_forSysFar = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband1_forSysFar, sideband2_forSysFar"));
                RooAbsReal *intbkgSideLeft_forSysFar = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband1_forSysFar"));
                RooAbsReal *intbkgSideRight_forSysFar = bkg_cheb.createIntegral(HFMass, NormSet(HFMass), Range("sideband2_forSysFar"));
                
                intresX = dcbPdf.createIntegral(HFMass, NormSet(HFMass), Range("signal"));
                float sig_yield = intsigX->getVal() * nsig->getVal();
                float bkg_yield = intbkgX->getVal() * nbkg->getVal();
                float res_yield = intresX->getVal() * nres->getVal();
                float bkg_frac = (bkg_yield) / (sig_yield + bkg_yield + res_yield);
                float res_frac = (res_yield) / (sig_yield + bkg_yield + res_yield);
                float sig_frac = (sig_yield) / (sig_yield + bkg_yield + res_yield);        
                float bkg_scale = intbkgX->getVal() / (intbkgSideLeft->getVal() + intbkgSideRight->getVal());
                float bkg_scale_forSysNear = intbkgX->getVal() / (intbkgSideLeft_forSysNear->getVal() + intbkgSideRight_forSysNear->getVal());
                float bkg_scale_forSysFar = intbkgX->getVal() / (intbkgSideLeft_forSysFar->getVal() + intbkgSideRight_forSysFar->getVal());
                float bkg_scale_forSysLower = intbkgX->getVal() / intbkgSideLeft->getVal();      
                float bkg_scale_forSysUpper = intbkgX->getVal() / intbkgSideRight->getVal();          
                cout << "Signal Yield = " << sig_yield << endl;
                cout << "Bkg Fraction = " << bkg_frac << endl;
                cout << "Bkg Fraction = " << res_frac << endl;

                vec_bkg_frac.push_back(bkg_frac);
                vec_res_frac.push_back(res_frac);
                vec_bkg_yield.push_back((res_yield + bkg_yield));
                vec_sig_yield.push_back(sig_yield);

                h1_MassMax->SetBinContent(i + 1, MassUp);
                h1_MassMin->SetBinContent(i + 1, MassDown);
                
                h1_Sideband1Min->SetBinContent(i + 1, Sideband1_Min);
                h1_Sideband1Max->SetBinContent(i + 1, Sideband1_Max);
                h1_Sideband2Min->SetBinContent(i + 1, Sideband2_Min);
                h1_Sideband2Max->SetBinContent(i + 1, Sideband2_Max);

                h1_Sideband1Min_forSysNear->SetBinContent(i + 1, Sideband1_Min_forSysNear);
                h1_Sideband1Max_forSysNear->SetBinContent(i + 1, Sideband1_Max_forSysNear);
                h1_Sideband2Min_forSysNear->SetBinContent(i + 1, Sideband2_Min_forSysNear);
                h1_Sideband2Max_forSysNear->SetBinContent(i + 1, Sideband2_Max_forSysNear);

                h1_Sideband1Min_forSysFar->SetBinContent(i + 1, Sideband1_Min_forSysFar);
                h1_Sideband1Max_forSysFar->SetBinContent(i + 1, Sideband1_Max_forSysFar);
                h1_Sideband2Min_forSysFar->SetBinContent(i + 1, Sideband2_Min_forSysFar);
                h1_Sideband2Max_forSysFar->SetBinContent(i + 1, Sideband2_Max_forSysFar);
                
                h1_BkgScale->SetBinContent(i + 1, bkg_scale);        

                h1_BkgScale_forSysNear->SetBinContent(i + 1, bkg_scale_forSysNear);
                h1_BkgScale_forSysFar->SetBinContent(i + 1, bkg_scale_forSysFar);

                h1_BkgScale_forSysLower->SetBinContent(i+1, bkg_scale_forSysLower);
                h1_BkgScale_forSysUpper->SetBinContent(i+1, bkg_scale_forSysUpper);

                std::cout << " // // // // // // Upper vs Lower Sideband Crosscheck // // // // // // " << std::endl;
                std::cout << " intbkgSideLeft : " << intbkgSideLeft->getVal() << std::endl;
                std::cout << " intbkgSideRight : " << intbkgSideRight->getVal() << std::endl;
                std::cout << " intbkgSide : " << intbkgSideLeft->getVal() + intbkgSideRight->getVal() << std::endl;
                std::cout << " intbkgSig : " << intbkgX->getVal() << std::endl;
                std::cout << " bkg_scale_forSysLower : " << bkg_scale_forSysLower << std::endl;
                std::cout << " bkg_scale_forSysUpper : " << bkg_scale_forSysUpper << std::endl;
                std::cout << "Area_Sig / Area_Sideband (bkg_scale) = " << bkg_scale << std::endl;
                std::cout << "Area left side / Area right side = " << intbkgSideLeft->getVal() / intbkgSideRight->getVal() << std::endl;



                // Add a vertical dashed line at x = 2

                if (isData)
                {
                        double x_linelow = MassDown;
                        TLine *linelow = new TLine(x_linelow, xframe->GetMinimum(), x_linelow, xframe->GetMaximum());
                        linelow->SetLineColorAlpha(kGreen + 1, 1.0);
                        linelow->SetLineStyle(3); // Dashed line
                        linelow->SetLineWidth(2);

                        double x_linehigh = MassUp;
                        TLine *linehigh = new TLine(x_linehigh, xframe->GetMinimum(), x_linehigh, xframe->GetMaximum());
                        linehigh->SetLineColorAlpha(kGreen + 1, 1.0);
                        linehigh->SetLineStyle(3); // Dashed line
                        linehigh->SetLineWidth(2);
                        
                        double x_sideline1low = Sideband1_Min;
                        TLine *sideline1low = new TLine(x_sideline1low, xframe->GetMinimum(), x_sideline1low, xframe->GetMaximum());
                        sideline1low->SetLineColorAlpha(kBlue + 1, 1.0);
                        sideline1low->SetLineStyle(3); // Dashed line
                        sideline1low->SetLineWidth(2);
                        
                        double x_sideline1high = Sideband1_Max;
                        TLine *sideline1high = new TLine(x_sideline1high, xframe->GetMinimum(), x_sideline1high, xframe->GetMaximum());
                        sideline1high->SetLineColorAlpha(kBlue + 1, 1.0);
                        sideline1high->SetLineStyle(3); // Dashed line
                        sideline1high->SetLineWidth(2);
                        
                        double x_sideline2low = Sideband2_Min;
                        TLine *sideline2low = new TLine(x_sideline2low, xframe->GetMinimum(), x_sideline2low, xframe->GetMaximum());
                        sideline2low->SetLineColorAlpha(kBlue + 1, 1.0);
                        sideline2low->SetLineStyle(3); // Dashed line
                        sideline2low->SetLineWidth(2);
                        
                        double x_sideline2high = Sideband2_Max;
                        TLine *sideline2high = new TLine(x_sideline2high, xframe->GetMinimum(), x_sideline2high, xframe->GetMaximum());
                        sideline2high->SetLineColorAlpha(kBlue + 1, 1.0);
                        sideline2high->SetLineStyle(3); // Dashed line
                        sideline2high->SetLineWidth(2);
                        
                        ican = 0;
                        ccan[ican]->cd(i + 1);
                        linelow->Draw("SAME");
                        linehigh->Draw("SAME");
                        sideline1low->Draw("SAME");
                        sideline1high->Draw("SAME");
                        sideline2low->Draw("SAME");
                        sideline2high->Draw("SAME");
                        
                        linelow->SetY1(xframe2->GetMinimum());
                        linehigh->SetY1(xframe2->GetMinimum());
                        sideline1low->SetY1(xframe2->GetMinimum());
                        sideline1high->SetY1(xframe2->GetMinimum());
                        sideline2low->SetY1(xframe2->GetMinimum());
                        sideline2high->SetY1(xframe2->GetMinimum());
                        
                        linelow->SetY2(xframe2->GetMaximum());
                        linehigh->SetY2(xframe2->GetMaximum());
                        sideline1low->SetY2(xframe2->GetMaximum());
                        sideline1high->SetY2(xframe2->GetMaximum());
                        sideline2low->SetY2(xframe2->GetMaximum());
                        sideline2high->SetY2(xframe2->GetMaximum());
                        
                        ican = 1;
                        ccan[ican]->cd(i + 1);
                        linelow->Draw("SAME");
                        linehigh->Draw("SAME");
                        sideline1low->Draw("SAME");
                        sideline1high->Draw("SAME");
                        sideline2low->Draw("SAME");
                        sideline2high->Draw("SAME");

                        if (i==5) {
                                Int_t lhcbFont = 132;
                                Double_t lhcbTSize = 0.05;
                                Double_t lhcbWidth = 2.00;
                                
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
                                
                                TGaxis::SetMaxDigits(3);
                                
                                TCanvas *c0 = new TCanvas("c0");
                                c0->cd();
                                c0->SetBottomMargin(0.15);
                                c0->SetTopMargin(0.1);
                                c0->SetLeftMargin(0.15);
                                c0->SetRightMargin(0.1);
                                
                                RooPlot *xframe_pub = HFMass.frame(Title(Form(";#it{m}_{#mu^{+}#mu^{-}K^{#pm}} [GeV/c^{2}];Candidates/(%.1f MeV/c^{2}) ", binsize * 1000.)));
                                xframe_pub->SetTitleOffset(1.1, "X");
                                xframe_pub->SetTitleOffset(0.8, "Y");
                        
                                B_mass.plotOn(xframe_pub, Name("B_mass"), MarkerStyle(kFullCircle), MarkerSize(4.), LineWidth(2));
                                // model->plotOn(xframe2);
                        
                                // Overlay the background component of model with a dashed line
                                model->plotOn(xframe_pub, Name("sig"), Components(sig), LineStyle(5), LineColor(kGreen + 2), LineWidth(2.));
                                model->plotOn(xframe_pub, Name("bkg_cheb"), Components(bkg_cheb), LineStyle(kDashed), LineColor(kBlue), LineWidth(2.));
                                model->plotOn(xframe_pub, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));
                                model->plotOn(xframe_pub, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(2.));

                                xframe_pub->Draw("SAME");

                                TLatex Tl;
                                Tl.SetNDC(kTRUE); // Use normalized coordinates (0-1)
                                
                                Tl.SetTextFont(lhcbFont);
                                Tl.SetTextColor(1);
                                Tl.SetTextSize(lhcbTSize);
                                Tl.SetTextAlign(12);
                                Tl.DrawLatex(0.575, 0.85, "#scale[1.0]{LHCb pp #sqrt{#it{s}} = 13 TeV}");
                                Tl.DrawLatex(0.575, 0.75, "#scale[1.0]{AK5 #it{B}^{#pm}-tagged jets}");
                                Tl.DrawLatex(0.575, 0.65, Form("#scale[1.0]{%.1f < #it{p}_{T,B^{#pm}} < %.1f GeV/c}", ptHF_binedges[i], ptHF_binedges[i+1]));
                                Tl.DrawLatex(0.575, 0.55, Form("#scale[1.0]{#it{p}_{T,jet} > 5 GeV/c}"));
                                Tl.DrawLatex(0.575, 0.45, "#scale[1.0]{2.5 < #it{y}_{jet} < 4}");

                                TLegend *leg_pub = new TLegend(0.2, 0.5, 0.3, 0.85);
                                leg_pub->SetTextSize(0.05);
                                leg_pub->SetBorderSize(0);
                                leg_pub->AddEntry("B_mass", "Data", "PE");
                                leg_pub->AddEntry("model", "Fit model", "L");
                                leg_pub->AddEntry("sig", "Signal", "L");
                                leg_pub->AddEntry("bkg_cheb", "Comb. bkg.", "L");

                                linelow->SetY2(5000);
                                linehigh->SetY2(5000);
                                sideline1low->SetY2(5000);
                                sideline1high->SetY2(5000);
                                sideline2low->SetY2(5000);
                                sideline2high->SetY2(5000);
                        
                                linelow->Draw("SAME");
                                linehigh->Draw("SAME");
                                sideline1low->Draw("SAME");
                                sideline1high->Draw("SAME");
                                sideline2low->Draw("SAME");
                                sideline2high->Draw("SAME");
                                leg_pub->Draw("SAME");

                                c0->SaveAs("./plots/BMass.pdf");
                        }
                }
                
                if (!DoSystematic)
                        SigParams << mu_sig->getVal() << "," << sigma1->getVal() << "," << sigma2->getVal() << \
                         "," << alpha1_sig->getVal() << "," << p1_sig->getVal() << "," << alpha2_sig->getVal() << \
                         "," << p2_sig->getVal() << "," << nsig->getVal() << ", " << sig_frac << endl;
                else
                        SigParams << sig2_m->getVal() << "," << sig2_t->getVal() << "," << sig2_n->getVal() << \
                        "," << nsig->getVal() << ", " << sig_frac << endl;

                BkgParams << a1->getVal() << "," << nbkg->getVal() << endl;
                ResonantParams << mu->getVal() << "," << width->getVal() << "," << alpha1->getVal() << "," << p1->getVal() << "," << alpha2->getVal() << "," << p2->getVal() << "," << nsig->getVal() << endl;
        } // End of mass fit loop

        file_workspace->Close();
        file_workspace_misid->Close();

        latexFile.close();

        ican = 0;
        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        
        ican = 1;
        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        
        if (ican > -1) {
                cout << " You plotted " << ican + 1 << " canvasses......." << endl;
                ccan[ican]->Print(plotfileC.Data());
        }

        latexFile.open(latexFileName, std::ios::app);
        latexFile << setprecision(3);
        latexFile << "\\begin{table}[H]\n";
        latexFile << "    \\centering\n";
        latexFile << "    \\caption{Fraction of combinatorial and resonant background in each \\pthf bin}\n";
        latexFile << "    \\label{tab:fit_fractions}\n";
        latexFile << "    \\begin{tabular}{ccc}\n";
        latexFile << "        \\pthf & Comb bkg fraction & Res bkg fraction \\\\\n";
        
        for (int i = 0; i < vec_bkg_frac.size(); i++) {
                latexFile << "        " << ptHF_binedges[i] << " $< \\pt <$ " << ptHF_binedges[i + 1] << " & "
                          << vec_bkg_frac[i] * 100 << "\\%" << " & " << vec_res_frac[i] * 100 << "\\% \\\\\n";
        }

        latexFile << "    \\end{tabular}\n";
        latexFile << "\\end{table}\n";
        latexFile.close();
        
        float sigyield = std::accumulate(vec_sig_yield.begin(), vec_sig_yield.end(), 0);
        
        cout << "Total Signal Yield = " << sigyield << endl;

        // cout<<exp_c->getVal()<<endl;

        SigParams.close();
        BkgParams.close();
        ResonantParams.close();

        f.Write();
        f.Close();
}
