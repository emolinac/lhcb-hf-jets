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

#include "../Bjets/Settings.h"
#include "../Helpers_IC.h"

#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"

using namespace RooFit;

double ptMin_user = pTLow;
double ptMax_user = 250;
        
void MassFit()
{
        ptMin = ptMin_user;
        ptMax = ptMax_user;

        if(gSystem->AccessPathName((output_folder + "ntuple_bjets_misid.root").c_str())) {
                std::cout<<"MisID ntuple file not found! Check existence or location."<<std::endl;

                return;
        }
        
        TFile fread((output_folder + "ntuple_bjets_misid.root").c_str(), "READ");
        
        TTree* BTree = (TTree*) fread.Get("BTree");

        int NumEvts = BTree->GetEntries();
        
        float mup_px, mup_py, mup_pz, mup_e;
        float mum_px, mum_py, mum_pz, mum_e;
        float pi_px, pi_py, pi_pz, pi_e, pi_PIDK;
        float jet_pt, bmass_dtf, chi2ndf_dtf, tau_dtf;
        float HF_px, HF_py, HF_pz, HF_e;
        float tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e;
        float dphi;
        float jet_rap;
        bool isTrueBjet, TOS;
        int nSV;
        
        BTree->SetBranchAddress("HF_px", &HF_px);
        BTree->SetBranchAddress("HF_py", &HF_py);
        BTree->SetBranchAddress("HF_pz", &HF_pz);
        BTree->SetBranchAddress("HF_e", &HF_e);

        BTree->SetBranchAddress("tr_HF_px", &tr_HF_px);
        BTree->SetBranchAddress("tr_HF_py", &tr_HF_py);
        BTree->SetBranchAddress("tr_HF_pz", &tr_HF_pz);
        BTree->SetBranchAddress("tr_HF_e", &tr_HF_e);

        BTree->SetBranchAddress("mum_px", &mum_px);
        BTree->SetBranchAddress("mum_py", &mum_py);
        BTree->SetBranchAddress("mum_pz", &mum_pz);
        BTree->SetBranchAddress("mum_e", &mum_e);

        BTree->SetBranchAddress("mup_px", &mup_px);
        BTree->SetBranchAddress("mup_py", &mup_py);
        BTree->SetBranchAddress("mup_pz", &mup_pz);
        BTree->SetBranchAddress("mup_e", &mup_e);

        BTree->SetBranchAddress("pi_px", &pi_px);
        BTree->SetBranchAddress("pi_py", &pi_py);
        BTree->SetBranchAddress("pi_pz", &pi_pz);
        BTree->SetBranchAddress("pi_e", &pi_e);
        BTree->SetBranchAddress("pi_PIDK", &pi_PIDK);

        BTree->SetBranchAddress("jet_pt", &jet_pt);
        BTree->SetBranchAddress("jet_rap", &jet_rap);

        BTree->SetBranchAddress("dphi", &dphi);
        BTree->SetBranchAddress("nSV", &nSV);
        BTree->SetBranchAddress("bmass_dtf", &bmass_dtf);
        BTree->SetBranchAddress("chi2ndf_dtf", &chi2ndf_dtf);
        BTree->SetBranchAddress("tau_dtf", &tau_dtf);
        BTree->SetBranchAddress("isTrueBjet", &isTrueBjet);

        BTree->SetBranchAddress("TOS", &TOS);

        float mass_low = 5.15;
        float mass_high = 5.55; // note: originally 5.6, must not be that important
        double jpsiMass = 3.096916;

        vector<TH1D *> h1_mass_HFpt;
        vector<float> vec_bkg_frac, vec_bkg_yield, vec_sig_yield;

        int nBins = 200;
        float binsize;
        for (int i = 0; i < ptHFbinsize; i++) {
                if (i >= ptHFbinsize - 3)
                        nBins = 60;
                else
                        nBins = 80;
                
                TH1D *h1_temp = new TH1D(Form("h1_mass%d", i), "", nBins, mass_low, mass_high);
                
                h1_mass_HFpt.push_back(h1_temp);
                
                binsize = (mass_high - mass_low) / (float)nBins;
        }

        TH1D *h1_mass     = new TH1D("h1_mass"    , "", nBins, mass_low, mass_high);
        TH1D *h1_mass_all = new TH1D("h1_mass_all", "", nBins + 17, 4.8, mass_high);
        
        TLorentzVector HFmeson, HFjet, mum, mup, pion, jpsi, tr_HFmeson;

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

                bool PID_cond = (pi_PIDK > 0);
                bool DTF_cond = (chi2ndf_dtf < 9) && (tau_dtf > 0.3);
                
                tr_HFmeson.SetPxPyPzE(tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e);
                HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);

                // if (fabs(HFmeson.Pt() - tr_HFmeson.Pt()) / tr_HFmeson.Pt() > 0.03)
                //     continue;

                mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
                mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);

                jpsi = mup + mum;
                jpsi.SetE(sqrt(jpsi.E() * jpsi.E() - jpsi.M() * jpsi.M() + jpsiMass * jpsiMass));
                
                pion.SetPxPyPzE(pi_px, pi_py, pi_pz, pi_e);

                HFmeson = jpsi + pion;

                h1_mass->Fill(HFmeson.M());
                h1_mass_all->Fill(HFmeson.M());

                for (int i = 0; i < ptHFbinsize; i++) {
                        if (HFmeson.Pt() > ptHF_binedges[i] && HFmeson.Pt() < ptHF_binedges[i + 1]) {
                                h1_mass_HFpt[i]->Fill(HFmeson.M());
                                
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

        TString extension_plots = "./plots/";
        TString extension       = "massfits_missid";
        
        TString rootfile    = output_folder +  extension + TString(".root");
        TString plotfile    = extension_plots + extension + TString(".ps");
        TString plotfilePDF = extension_plots + extension + TString(".pdf");
        TString plotfileO   = plotfilePDF + TString("(");
        TString plotfileC   = plotfilePDF + TString("]");

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

        std::cout<<"S e t u p   c o m p o n e n t   p d f s"<<std::endl;

        // Declare observable x
        RooWorkspace *w_read;
        TFile *file_workspace;

        for (int i = 0; i < ptHFbinsize; i++) {
                if (i >= ptHFbinsize - 3)
                        nBins = 60;
                else
                        nBins = 80;
                
                binsize = (mass_high - mass_low) / (float)nBins;

                w_read = new RooWorkspace(Form("w%d", i));
                
                RooRealVar *mu_sig, *alpha1_sig, *alpha2_sig, *p1_sig, *p2_sig;
                
                RooRealVar HFMass("HFMass", "HFMass", mass_low, mass_high);
                RooDataHist B_mass("B_mass", "B_mass", RooArgList(HFMass), h1_mass_HFpt[i], 1.);

                HFMass.setBins(nBins);

                ////////////////////////////////////////////////////
                // build siganl
                ///////////////////////////////////////////////////

                RooRealVar* sigma_ratio = new RooRealVar("sigma_ratio", "sigma_ratio"      , 1.78528 / 3.15958);
                RooRealVar* mean        = new RooRealVar("mean"       ,"mean of gaussians" , 5.279, 5.27, 5.282);
                RooRealVar* width       = new RooRealVar("width"      ,"width of gaussians", 0.002, 0.001, 0.05);
                RooRealVar* sigma2      = new RooRealVar("sigma2"     ,"width of gaussians", 0.003, 0.001, 0.05); // CHANGE
                                
                // These parameters have been derived from simulation
                RooRealVar* mu     = new RooRealVar("mu", "mu", 5.32, 5.2, 5.4);
                RooRealVar* alpha1 = new RooRealVar("alpha1", "alpha1", 2., 0.01, 10.);
                RooRealVar* alpha2 = new RooRealVar("alpha2", "alpha2", 2., 0.01, 10.);
                RooRealVar* p1     = new RooRealVar("p1", "p1", 3., 1., 4.); // note: starting value 3, allowed to range between 1 and 4 
                RooRealVar* p2     = new RooRealVar("p2", "p2", 2., 1., 4.); // note: starting value 2, allowed to range between 1 and 4 
                // p1 = new RooRealVar("p1", "p1", 2.); // note: why commenting these? maybe to rigid ?
                // p2 = new RooRealVar("p2", "p2", 3.); // note: why commenting these? maybe to rigid ?
                // p1->setConstant(kTRUE);
                // p2->setConstant(kTRUE);
                
                // Define the CB functions
                RooCrystalBall sig("sig", "floatSidedCB_sig1", HFMass, *mu, *width, *alpha1, *p1, *alpha2, *p2);

                /////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Build linear background PDF, for combinatoric background
                ///////////////////////////////////////////////////

                RooRealVar *a1 = new RooRealVar("exp_c", "exp_c", -2., -4, -0.1);
                RooRealVar *a2 = new RooRealVar("c2", "coefficient #2", -0.1, -1., -0.001);
                
                RooExponential bkg("bkg", "bkg", HFMass, *a1);

                RooRealVar *a0_nosec = new RooRealVar("a0_nosec", "a0_nosec", 0.5, -1, 1.);
                RooRealVar *a1_nosec = new RooRealVar("a1_nosec", "a1_nosec", 0.2, 0., 1.);
                
                RooChebychev bkg_nosec("bkg_nosec", "Background_nosec", HFMass, RooArgSet(*a0_nosec, *a1_nosec));

                /////////////////////////////////////////////////////////////////////////////////////////////////////////
                // Fractions
                ///////////////////////////////////////////////////

                RooRealVar *nsig       = new RooRealVar("nsig", "fraction in signal", 10000, 0., 10000000.);
                RooRealVar *nsig_nosec = new RooRealVar("nsig_nosec", "fraction in signal", 500, 0., 1000000.);

                RooRealVar *nbkg       = new RooRealVar("nbkg", "fraction of background", 5000, 0., 1000000);
                RooRealVar *nbkg_nosec = new RooRealVar("nbkg_nosec", "fraction of background", 200, 0., 1000000);

                RooRealVar *ntanh      = new RooRealVar("ntanh", "fraction of background", 200, 0., 1000000);

                RooRealVar *nres       = new RooRealVar("nres", "fraction of background", 1000, 0., 0.1 * nsig->getVal());
                RooRealVar *nres_nosec = new RooRealVar("nres_nosec", "fraction of background", 200, 0., 1000000);
                
                // S a m p l e ,   f i t   a n d   p l o t   m o d e l
                // ---------------------------------------------------
                HFMass.setRange("comb1", 5.55, 5.78);
                HFMass.setRange("comb2", RooConst(mass_low), RooConst(mass_low + 0.1));
                HFMass.setRange("tanh", 5.05, 5.2);
                HFMass.setRange("noSec", 5.18, 5.6);
                HFMass.setRange("Full", mass_low, mass_high);

                // Fit model to data
                bkg_nosec.fitTo(B_mass, Range("comb1"));
                
                // Model in this case has no background
                RooAddPdf* model = new RooAddPdf("model", "g1+g2+a", RooArgList(sig), RooArgList(*nsig));

                model->fitTo(B_mass, PrintEvalErrors(-1), Save(true));

                // Begin plotting
                ican = 0;
                ccan[ican]->cd(i + 1);
                gPad->SetRightMargin(0.12);
                gPad->SetLeftMargin(0.12);
                ccan[ican]->SetFillStyle(4070);
                RooPlot *xframe = HFMass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
                xframe->SetTitleOffset(0.7, "Y");
                xframe->SetTitleOffset(0.9, "X");
                // model->paramOn(xframe, Layout(0.6, 0.9, 0.9));
                // xframe->getAttText()->SetTextSize(0.02);
                B_mass.plotOn(xframe, Name("B_mass"));
                // model->plotOn(xframe);

                // Overlay the background component of model with a dashed line
                model->plotOn(xframe, Name("sig"), Components(sig), LineStyle(5), LineColor(kGreen + 2), LineWidth(1.));
                // sig.plotOn(xframe, Name("sig"), LineStyle(kSolid), LineColor(kGreen), Range("Full"));
                // model->plotOn(xframe, Components(sig2), LineStyle(kSolid), LineColor(kMagenta));
                model->plotOn(xframe, Name("bkg"), Components(bkg), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
                // model->plotOn(xframe, Name("exp_bkg"), Components(exp_bkg), LineStyle(kDashed), LineColor(kBlue));
                // model->plotOn(xframe, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));
                model->plotOn(xframe, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1.));
                RooHist *hpull = xframe->pullHist();
                
                float chi2 = xframe->chiSquare(10);
                
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
                // sig.plotOn(xframe2, Name("sig"), LineStyle(kSolid), LineColor(kGreen), Range("Full"));
                // model->plotOn(xframe2, Components(sig2), LineStyle(kSolid), LineColor(kMagenta));
                model->plotOn(xframe2, Name("bkg"), Components(bkg), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
                // model->plotOn(xframe2, Name("exp_bkg"), Components(exp_bkg), LineStyle(kDashed), LineColor(kBlue));
                // model->plotOn(xframe2, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));
                model->plotOn(xframe2, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1.));

                // exp_bkg.plotOn(xframe2, LineColor(kMagenta));

                // Overlay the background+sig2 components of model with a dotted line
                // model->plotOn(xframe2, Components(RooArgSet(bkg, sig2)), LineStyle(kDotted));
                // h1_mass_all->Draw("SAME");
                xframe2->Draw("SAME");
                // gPad->SetLogy();
                xframe2->SetMinimum(1.);

                RooWorkspace *w = new RooWorkspace(Form("w%d", i), Form("workspace%d", i));
                // w->import(sign);
                // w->import(bkgn);
                w->import(*model);
                w->import(B_mass);
                w->Print();
                w->writeToFile(output_folder_massfits + Form("workspace%d_", i) + extension + ".root");

                cout << "Chi2/dof = " << chi2 << endl;
                // cout << sigma_arg->getVal() << endl;
                cout << "a1 = " << a1->getVal() << endl;
                cout << "a2 = " << a2->getVal() << endl;
                cout << "width = " << width->getVal() << endl;
                cout << "sigma2 = " << sigma2->getVal() << endl;

                float MassMu = mu->getVal();
                float MassSigma = width->getVal();
                float MassUp = MassMu + 2 * MassSigma;
                float MassDown = MassMu - 2 * MassSigma;
                cout << MassDown << ", " << MassUp << endl;
                // HFMass.setRange("signal", MassDown, MassUp);
                HFMass.setRange("signal", MassMin, MassMax);

                // intsigX = sig->createIntegral(*HFMass, NormSet(*HFMass), Range("signal"));
                // intbkgX = bkg->createIntegral(*HFMass, NormSet(*HFMass), Range("signal"));
                RooAbsReal *intsigX;
                RooAbsReal *intbkgX, *intresX;

                intsigX = sig.createIntegral(HFMass, NormSet(HFMass), Range("signal"));
                intbkgX = bkg.createIntegral(HFMass, NormSet(HFMass), Range("signal"));
                intresX = 0;
                
                float sig_yield = intsigX->getVal() * nsig->getVal();
                float bkg_yield = intbkgX->getVal() * nbkg->getVal();
                float res_yield = 0;
                float bkg_frac = (res_yield + bkg_yield) / (sig_yield + bkg_yield + res_yield);
                
                // cout << "Bkg Yield = " << intbkgX->getVal() << endl;
                // cout << "Sig Yield = " << intsigX->getVal() << endl;
                cout << "Signal Yield = " << sig_yield << endl;
                cout << "Bkg Fraction = " << bkg_frac << endl;
                vec_bkg_frac.push_back(bkg_frac);
                vec_bkg_yield.push_back((res_yield + bkg_yield));
                vec_sig_yield.push_back(sig_yield);
        }

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

        for (int i = 0; i < vec_bkg_frac.size(); i++) {
                cout << vec_bkg_frac[i] * 100 << " %"
                << ", " << vec_bkg_yield[i] << endl;
        }
        
        float sigyield = std::accumulate(vec_sig_yield.begin(), vec_sig_yield.end(), 0);
        cout << "Signal Yield = " << sigyield << endl;
}

