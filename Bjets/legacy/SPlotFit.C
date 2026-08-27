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
#include "Settings.h"
#include "../Helpers.h"
using namespace RooFit;

void SPlotFit(int NumEvts = -1,
              int dataset = 91599,
              bool isData = true,
              bool UseDTF = true,
              bool DoRecSelEff = 0,
              bool DoSystematic = 0,
              float ptmin_user = pTLow,
              float ptmax_user = 250,
              bool L0MuonDiMuon = false)
{
    bool MCflag = !isData;
    followHardest = false;
    truthLevel = false;
    chargedJetCut = false;
    ptMin = ptmin_user;
    ptMax = ptmax_user;
    if (truthLevel)
    {
        ghostCut = false;
    }
    int year = (dataset / 10000) % 10;
    int JetMeth = (dataset / 1000) % 10;
    int flavor = (dataset / 100) % 10;
    int ptRange = (dataset / 10) % 10;
    int Mag = (dataset / 1) % 10;
    int HF_pdgcode = -99;
    if (flavor == 5)
    {
        mass_num = 4.2;
        HF_pdgcode = 521;
    }
    else if (flavor == 4)
    {
        mass_num = 1.25;
        HF_pdgcode = 421;
    }
    else if (flavor == 1)
    {
        mass_num = 0.001;
        followHardest = true;
    }
    TString str_Mag = "";
    TString str_pT = "";
    TString str_level = "";
    TString str_followHard = "";
    TString str_flavor = "";
    TString str_ghost = "";
    TString str_charged = "";
    TString str_DTF = "";
    TString str_PID = "";
    TString str_year = "";
    if (year == 6)
        str_year = "2016";
    else if (year == 7)
        str_year = "2017";
    else if (year == 8)
        str_year = "2018";
    if (DTF_cut)
        str_DTF = "_DTF";
    if (PID_cut)
        str_PID = "_PID";

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
    if (isData)
        str_level = "data";
    else
    {
        if (truthLevel)
            str_level = "truth";
        else
            str_level = "reco";
    }
    if (flavor == 1)
        str_followHard = "_hard";
    else
    {
        if (followHardest)
            str_followHard = "_hard";
        else
            str_followHard = "_HF";
    }

    if (ghostCut)
        str_ghost = Form("_ghost_%.1f", ghostProb);

    if (chargedJetCut)
        str_charged = "_charge";


    TString str_L0 = "";
    if (L0MuonDiMuon)
        str_L0 = "_L0MuonDiMuon";

    TString extension, extension_reco, extension_misid;
    TString extension_read, extension_RootFiles, extension_RootFilesMC;
    
    extension_RootFiles = isData ? TString("../../root_files/Bjets/") : TString("../../root_files/BjetsMC/");
    extension_RootFilesMC = TString("../../root_files/BjetsMC/");
    
    //extension_reco = TString("splotfit_") + "reco" + Form("_ev_%d", -1) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_ghost + str_charged + str_Mag + str_flavor + str_DTF + str_PID + Form("_%d", dataset);
    extension_reco = TString("massfit_") + "reco" + Form("_ev_%d", -1) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_ghost + str_charged + str_Mag + str_flavor + str_DTF + str_PID + str_L0 + Form("_%d", dataset);
    
    //extension_misid = TString("splotfit_") + "reco" + Form("_ev_%d", -1) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_ghost + str_charged + str_flavor + str_DTF + str_PID + Form("_%d", 91599);
    extension_misid = TString("massfit_") + "reco" + Form("_ev_%d", -1) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_ghost + str_charged + str_flavor + str_DTF + str_PID + str_L0 + Form("_%d", 91599);

    extension = TString("splotfit_") + str_level + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_ghost + str_charged + str_Mag + str_flavor + str_DTF + str_PID + str_L0 + Form("_%d", dataset);

    if (DoRecSelEff)
    {
        extension = "recselsys_" + extension;
        extension_reco = "recselsys_" + extension_reco;
    }
    if (DoSystematic)
    {
        extension = "sys_" + extension;
        // extension_reco = "recselsys_" + extension_reco;
    }
    cout << extension << endl;
    cout << extension_reco << endl;
    // Setup Tree
    //extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + Form("_%d", dataset);

    // TFile fread(dir_deadcone + "hists/" + extension_read + ".root", "READ");
    //TFile fread(extension_RootFiles + extension_read + ".root", "READ");

    TChain *BTree = new TChain("BTree", "");
    vector<int> vec_datasets;
    if (Mag == 0)
        vec_datasets = {61590, 71590, 81590};
    else if (Mag == 1)
        vec_datasets = {61591, 71591, 81591};
    else
        vec_datasets = {61590, 61591, 71590, 71591, 81590, 81591};

    if (year == 9 && JetMeth != 9)
    {
        for (int i = 0; i < vec_datasets.size(); i++)
        {
            Mag = (vec_datasets[i] / 1) % 10;
            if (Mag == 0)
                str_Mag = "_MD";
            else if (Mag == 1)
                str_Mag = "_MU";
            extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + str_L0 + Form("_%d", vec_datasets[i]);
            cout << extension_read << endl;
            BTree->Add(extension_RootFiles + extension_read + ".root/BTree");
        }
    }
    else
    {
        extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + str_L0 + Form("_%d", dataset);
        BTree->Add(extension_RootFiles + extension_read + ".root/BTree");
    }
    if (NumEvts > BTree->GetEntries())
        NumEvts = BTree->GetEntries();
    if (NumEvts == -1)
        NumEvts = BTree->GetEntries();
    cout << BTree->GetEntries() << endl;
    float mup_px, mup_py, mup_pz, mup_e;
    float mum_px, mum_py, mum_pz, mum_e;
    float K_px, K_py, K_pz, K_e, K_PIDK;
    float jet_pt, bmass_dtf, chi2ndf_dtf, tau_dtf;
    float HF_px, HF_py, HF_pz, HF_e;
    float tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e;
    float Jpsi_CHI2NDOF, Jpsi_BPVDLS, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_IPCHI2;
    float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
    float dphi;
    float jet_rap;
    bool isTrueBjet, TOS;
    int nSV;
    float sweight;
    BTree->SetBranchAddress("HF_px", &HF_px);
    BTree->SetBranchAddress("HF_py", &HF_py);
    BTree->SetBranchAddress("HF_pz", &HF_pz);
    BTree->SetBranchAddress("HF_e", &HF_e);

    BTree->SetBranchAddress("tr_HF_px", &tr_HF_px);
    BTree->SetBranchAddress("tr_HF_py", &tr_HF_py);
    BTree->SetBranchAddress("tr_HF_pz", &tr_HF_pz);
    BTree->SetBranchAddress("tr_HF_e", &tr_HF_e);

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
    BTree->SetBranchAddress("jet_rap", &jet_rap);

    BTree->SetBranchAddress("dphi", &dphi);
    BTree->SetBranchAddress("nSV", &nSV);
    BTree->SetBranchAddress("bmass_dtf", &bmass_dtf);
    BTree->SetBranchAddress("chi2ndf_dtf", &chi2ndf_dtf);
    BTree->SetBranchAddress("tau_dtf", &tau_dtf);
    BTree->SetBranchAddress("isTrueBjet", &isTrueBjet);

    BTree->SetBranchAddress("TOS", &TOS);

    float mass_low = 5.15;
    float mass_high = 5.55;
    // if (UseDTF)


    vector<TH1D *> h1_mass_HFpt;
    vector<float> vec_bkg_frac, vec_res_frac, vec_bkg_yield, vec_sig_yield;

    TFile f(extension_RootFiles + extension + ".root", "RECREATE");

    int nBins = 80;
    float binsize;
    // int nBins[ptHFbinsize] = {100, 100, 100, 100, 100, 100};
    for (int i = 0; i < ptHFbinsize; i++)
    {
        if (i >= ptHFbinsize - 2)
            nBins = 60;
        else
            nBins = 80;
        TH1D *h1_temp = new TH1D(Form("h1_mass%d", i), "", nBins, mass_low, mass_high);
        h1_mass_HFpt.push_back(h1_temp);
        binsize = (mass_high - mass_low) / (float)nBins;
    }
    TH1D *h1_mass = new TH1D("h1_mass", "", nBins, mass_low, mass_high);
    TH1D *h1_mass_all = new TH1D("h1_mass_all", "", nBins + 17, 4.8, mass_high);
    TH1D *h1_MassMin = new TH1D("h1_MassMin", "", ptHFbinsize, ptHF_binedges);
    TH1D *h1_MassMax = new TH1D("h1_MassMax", "", ptHFbinsize, ptHF_binedges);
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
    TH1D *h1_BkgScale_forSysNear = new TH1D("h1_BkgScale_forSysNear", "", ptHFbinsize, ptHF_binedges);
    TH1D *h1_BkgScale_forSysFar = new TH1D("h1_BkgScale_forSysFar", "", ptHFbinsize, ptHF_binedges);

    // Define variables
    RooRealVar HFmass("HFmass", "HFmass", mass_low, mass_high);
    RooRealVar HFpt("HFpt", "Transverse Momentum", 0, 200);
    std::vector<double> sWeights(NumEvts, 0); // Store weights
    TTree *tree = new TTree("sWeightTree", "Tree with sWeights");
    TBranch *newBranch = tree->Branch("sweight", &sweight);

    // Create an unbinned dataset with mass and HFpt
    RooArgSet variables(HFmass, HFpt);
    RooDataSet dataSet("dataSet", "Unbinned Data", variables);

    // NumEvts = 100000;
    for (int ev = 0; ev < NumEvts; ev++)
    {
        BTree->GetEntry(ev);

        if (ev % 100000 == 0)
            cout << "Executing event " << ev << endl;
        if (jet_rap < etaMin || jet_rap > etaMax)
            continue;
        if (jet_pt > ptMax)
            continue;
        if (jet_pt < ptMin)
            continue;
        if (!TOS)
            continue;
        if (DoRecSelEff)
        {
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

        TLorentzVector HFmeson, HFjet, mum, mup, Kmeson, Jpsi, tr_HFmeson;
        mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
        mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);
        Kmeson.SetPxPyPzE(K_px, K_py, K_pz, K_e);
        // HFmeson = mup + mum + Kmeson;
        HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
        tr_HFmeson.SetPxPyPzE(tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e);

        // if(dtf_chi2ndf > 5) continue;
        // cout<<bmass_dtf<<endl;
        float bmass = HFmeson.M();
        if (UseDTF)
            bmass = bmass_dtf;
        bool PID_cond = (K_PIDK > 0);
        bool DTF_cond = (chi2ndf_dtf < 9) && (tau_dtf > 0.3);

        if (!DTF_cond && DTF_cut)
            continue;
        if (!PID_cond && PID_cut)
            continue;
        if (!isData)
        {
            if (!isTrueBjet)
                continue;
        }

        if (bmass < mass_low || bmass > mass_high)
            continue;

        HFmass.setVal(bmass);
        HFpt.setVal(HFmeson.Pt());
        dataSet.add(variables);

    }

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
    TString rootfile;
    TString plotfile;
    TString plotfilePDF;
    TString plotfileO;
    TString plotfileC;
    // TString OutputFileBase	= outbase+outinfo;
    TString plotextension = isData ? TString("../../plots/Bjets/") : TString("../../plots/BjetsMC/");   
    rootfile = extension_RootFiles + extension + TString(".root");
    plotfile = plotextension + extension + TString(".ps");
    plotfilePDF = plotextension + extension + TString(".pdf");
    plotfileO = plotfilePDF + TString("(");
    plotfileC = plotfilePDF + TString("]");

    // Plot data and PDF overlaid
    ican = 0;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 2100, 1800);
    ccan[ican]->SetFillColor(10);

    // gPad->SetBottomMargin(0.06);
    ccan[ican]->cd();
    ccan[ican]->Divide(3, 5, 0.0001, 0.0001);

    ican = 1;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 2100, 1800);
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
    {
        SigParams << "m,t,n,nsig,fraction" << endl;
    }
    else
        SigParams << "mu,sigma1,sigma2,alpha1,n1,n1_2,alpha2,n2,n2_2,nsig" << endl;

    BkgParams << "a0,a1,a2,nbkg,fraction" << endl;
    ResonantParams << "mu,sigma,alpha1,n1,alpha2,n2,nres,fraction" << endl;

    // S e t u p   c o m p o n e n t   p d f s
    // ---------------------------------------

    // Declare observable x
    for (int i = 0; i < ptHFbinsize; i++)
    {

        TFile *file_workspace = new TFile(Form( extension_RootFilesMC + "workspace%d_", i) + extension_reco + ".root", "READ");
        TFile *file_workspace_misid = new TFile(extension_RootFilesMC + TString("MisID/") + Form("workspace%d_", i) + extension_misid + TString(".root"), "READ");    

        std::cout <<  extension_RootFilesMC << Form("workspace%d_", i) << extension_reco << ".root" << std::endl;
        std::cout << extension_RootFilesMC << TString("MisID/") << extension_misid << TString(".root");

        if (file_workspace_misid == NULL)
        {
            cout << "no file at " << extension_misid << endl;
        }
        RooWorkspace *w_read = (RooWorkspace *)file_workspace->Get(Form("w%d", i));
        RooWorkspace *w_read_misid = (RooWorkspace *)file_workspace_misid->Get(Form("w%d", i));
        RooRealVar *sigma_ratio, *mean, *sigma1, *sigma2;
        RooRealVar *mu_sig, *alpha1_sig, *alpha2_sig, *p1_sig, *p2_sig;
        RooRealVar *mu_sig2, *alpha1_sig2, *alpha2_sig2, *p1_sig2, *p2_sig2;
        RooRealVar *mu, *width, *alpha1, *alpha2, *p1, *p2;

        // Select data in this pt range
        RooDataSet *data_pt_bin = (RooDataSet *)dataSet.reduce(Form("HFpt >= %f && HFpt < %f", ptHF_binedges[i], ptHF_binedges[i + 1]));
        RooDataHist B_mass("B_mass", "B_mass", RooArgSet(HFmass), *data_pt_bin);

        // Create two Gaussian PDFs g1(HFmass,mean1,sigma) anf g2(HFmass,mean2,sigma) and their parameters
        ////////////////////////////////////////////////////
        // build siganl
        ///////////////////////////////////////////////////

        if (UseDTF)
        {
            sigma_ratio = new RooRealVar("sigma_ratio", "sigma_ratio", 1.52571 / 0.856549);
            mean = new RooRealVar("mean", "mean of gaussians", 5.27966, 5.27, 5.282);
            sigma1 = new RooRealVar("sigma1", "width of gaussians", 0.011, 0.001, 0.03);
            sigma2 = new RooRealVar("sigma2", "width of gaussians", 0.007, 0.001, 0.03); // CHANGE
        }
        else
        {
            sigma_ratio = new RooRealVar("sigma_ratio", "sigma_ratio", 1.78528 / 3.15958);
            mean = new RooRealVar("mean", "mean of gaussians", 5.279, 5.27, 5.282);
            sigma1 = new RooRealVar("sigma1", "width of gaussians", 0.02, 0.001, 0.5);
            sigma2 = new RooRealVar("sigma2", "width of gaussians", 0.03, 0.001, 0.5); // CHANGE
        }

        RooGaussian gauss1("gauss_sig1", "Signal component 1", HFmass, *mean, *sigma1);
        RooGaussian gauss2("gauss_sig2", "Signal component 2", HFmass, *mean, *sigma2);

        mu_sig = new RooRealVar("mu_sig", "mu", 5.279, 5.27, 5.282);
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

        // else
        //     return;

        if (isData && w_read != NULL)
        {
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
        RooCrystalBall dcbPdf_sig1("dcbPdf_sig1", "floatSidedCB_sig1", HFmass, *mu_sig, *sigma1, *alpha1_sig, *p1_sig, *alpha2_sig, *p2_sig);
        // RooCrystalBall dcbPdf_sig2("dcbPdf_sig2", "floatSidedCB_sig2", HFmass, *mu_sig, *sigma2, *alpha1_sig, *p1_sig2, *alpha2_sig, *p2_sig2);
        RooCrystalBall dcbPdf_sig2("dcbPdf_sig2", "floatSidedCB_sig2", HFmass, *mu_sig, *sigma2, *alpha1_sig, *p1_sig, *alpha2_sig, *p2_sig);

        // Number of signal events
        RooRealVar *nsig1 = new RooRealVar("nsig1", "fraction of component 1 in signal", 0.3, 0., 1.);
        // RooRealVar *nsig1 = new RooRealVar("nsig1", "fraction of component 1 in signal", 0., 0., 0.);
        //if (isData && w_read != NULL)
        //{
        //    nsig1 = (RooRealVar *)w_read->obj("nsig1");
            //nsig1->setConstant(kTRUE);
        //}
        // Choose if you want two CB functions as signal, or two Gaussians.

        // RooAddPdf sig("sig", "Signal", RooArgList(gauss1, gauss2), RooArgList(*nsig1));
        RooAddPdf sig("sig", "Signal", RooArgList(dcbPdf_sig1, dcbPdf_sig2), RooArgList(*nsig1));

        RooRealVar *sig2_m = new RooRealVar("sig2_m", "sig2_m", 5.279, 5.27, 5.282);
        RooRealVar *sig2_n = new RooRealVar("sig2_n", "sig2_n", 5.);
        sig2_n->setConstant(kTRUE);
        RooRealVar *sig2_t = new RooRealVar("sig2_t", "sig2_t", 0.02, 0.001, 0.3);

        RooGenericPdf sig2("student-t", "", "tgamma((sig2_n+1)/2)/(tgamma(sig2_n/2)*sqrt(TMath::Pi()*sig2_n)*sig2_t) * pow(1+(1/sig2_n * pow((HFmass-sig2_m)/sig2_t, 2)), (-1-sig2_n)/2)", RooArgList(HFmass, *sig2_m, *sig2_n, *sig2_t));

        //
        // RooRealVar mean_nosec("mean_nosec", "mean of gaussians", 5.279, 5.25, 5.285);
        // RooRealVar sigma1_nosec("sigma1_nosec", "width of gaussians", 0.05, 0, 0.1);
        // RooRealVar sigma2_nosec("sigma2_nosec", "width of gaussians", 0.01, 0, 0.1);
        //
        // RooGaussian sig1_nosec("sig1", "Signal component 1", HFmass, mean_nosec, sigma1_nosec);
        // RooGaussian sig2_nosec("sig2", "Signal component 2", HFmass, mean_nosec, sigma2_nosec);
        //
        // RooRealVar nsig1_nosec("nsig1_nosec", "fraction of component 1 in signal", 0.3, 0., 1.);
        // RooAddPdf sig_nosec("sig_nosec", "Signal", RooArgList(sig1_nosec, sig2_nosec), RooArgList(nsig1_nosec));

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // build resonant structure, values from simulation with Jpsi fixed to its mass
        ///////////////////////////////////////////////////
        // EXT PARAMETER                                INTERNAL      INTERNAL        │
        // NO.   NAME      VALUE            ERROR       STEP SIZE       VALUE         │p.d.f.s
        //  1  a1           5.82768e+00   1.70615e-08   5.16495e-07  -1.08317e+00     │-------
        //  2  a2          -9.99999e-01   6.19314e-02   1.47440e-01  -1.57207e+00     │RooCrystalBall::dcbPdf_sig1[ x=HFmass x0=mu_sig sigmaL=sigma1
        //  3  alpha1       2.01441e+00   3.24539e-02   1.02160e-03  -6.39903e-01     │ sigmaR=sigma1 alphaL=alpha1_sig nL=p1_sig alphaR=alpha2_sig
        //  4  alpha2       5.70189e-01   5.50017e-03   5.10050e-04  -1.08856e+00     │nR=p2_sig ] = 0.00214893
        //  5  mu           5.32697e+00   1.04802e-04   8.24657e-05   2.56748e-01     │RooAddPdf::model[ nsig * dcbPdf_sig1 ] = 0.00214893
        //  6  nbkg         9.68910e+01   3.49620e+02   2.81877e-03  -1.55111e+00     │
        //  7  nsig         9.17769e+04   4.78325e+02   4.91049e-04  -9.55229e-01     │datasets
        //  8  p1           1.31970e+00   1.23509e-01   2.03443e-03  -8.27233e-01     │--------
        //  9  p2           2.27406e+00   2.53814e-02   1.32025e-03  -5.76613e-01     │RooDataHist::B_mass(HFmass)
        // 10  width        1.59333e-02   1.00715e-04   2.54854e-04  -1.31767e+00

        if (UseDTF && w_read_misid != NULL)
        {

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
        }
        else
        {
            // return;
            mu = new RooRealVar("mu", "mu", 5.33183e+00);
            width = new RooRealVar("width", "width", 2.62897e-02);
            alpha1 = new RooRealVar("alpha1", "alpha1", 1.72233e+00);
            alpha2 = new RooRealVar("alpha2", "alpha2", 7.67677e-01);
            p1 = new RooRealVar("p1", "p1", 1.73799e+00);
            p2 = new RooRealVar("p2", "p2", 2.09840e+00);
        }

        RooCrystalBall dcbPdf("dcbPdf", "floatSidedCB", HFmass, *mu, *width, *alpha1, *p1, *alpha2, *p2);

        //
        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Build argus background PDF, I can't get this to work, so DISREGARD FOR NOW!
        /////////////////////////////////////////////////////////////////////////////////

        RooRealVar *m_arg = new RooRealVar("m_arg", "m", 5.15, 5.1, 5.17); // The peak position of the Argus function
        RooRealVar *c_arg = new RooRealVar("c_arg", "c", -20, -30, -1);    // The curvature parameter of the Argus function
        RooRealVar *n_arg = new RooRealVar("n_arg", "n", 0.5);             // The power-law index of the Argus function

        RooRealVar *mean_arg = new RooRealVar("mean_arg", "mean", 0); // The mean of the Gaussian function
        // RooFormulaVar sigma_arg = new RooRealVar("sigma_arg", "sigma", "sqrt(sigma1*sigma1*nsig1 + (1-nsig1)*sigma2*sigma2)", RooArgList(sigma1, sigma2, nsig1));  // The standard deviation of the Gaussian function
        RooRealVar *sigma_arg = new RooRealVar("sigma_arg", "sigma", 0.01, 0.001, 0.05);

        // Create the Argus function and the Gaussian function
        RooArgusBG argus("argus", "argus", HFmass, *m_arg, *c_arg, *n_arg);
        RooGaussian gauss_arg("gauss_arg", "Signal component 2", HFmass, *mean_arg, *sigma_arg);

        RooFFTConvPdf argxgauss("argxgauss", "argus (X) gauss", HFmass, argus, gauss_arg);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Build tanh background PDF, left shoulder background
        ///////////////////////////////////////////////////

        RooRealVar *b1 = new RooRealVar("b1", "scale", 0.015, 0.001, 0.04);
        RooRealVar *b2 = new RooRealVar("b2", "scale", 0.015, 0.001, 0.07);
        RooRealVar *s1 = new RooRealVar("s1", "shift", 5.14, 5.12, 5.16);
        RooRealVar *s2 = new RooRealVar("s2", "shift", 4.95, 4.92, 4.96);
        RooGenericPdf tanhpdf("tanhpdf", "tanhpdf", "1 - tanh((HFmass-s1)/b1)", RooArgSet(HFmass, *b1, *s1));
        // RooGenericPdf tanhpdf("tanhpdf", "tanhpdf", "tanh((x-s2)/b2) - tanh((x-s1)/b1)", RooArgSet(x, b1, b2, s1, s2));

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Build linear background PDF, for combinatoric background
        ///////////////////////////////////////////////////

        RooRealVar *a2 = new RooRealVar("c2", "coefficient #2", -0.1, -1., -0.001);
        RooRealVar *a1 = new RooRealVar("exp_c", "exp_c", -2., -4, -0.1);
        RooExponential bkg("bkg", "bkg", HFmass, *a1);

        //  RooRealVar * a2("a2", "a2", -1, -2, -0.1);
        // RooRealVar * a1("a1", "a1", 2, 0., 9.);
        // RooPolynomial bkg("bkg", "Background", HFmass, RooArgList(a1, a2));

        RooRealVar *a0_cheb = new RooRealVar("a0_cheb", "a0_cheb", -0.5, -3, -0.0001);
        RooRealVar *a1_cheb = new RooRealVar("a1_cheb", "a1_cheb", -0.2, -3, -0.0001);
        RooRealVar *a2_cheb = new RooRealVar("a2_cheb", "a2_cheb", -0.2, -3, -0.0001);
        // RooChebychev bkg_cheb("bkg_cheb", "Background_cheb", HFmass, RooArgSet(*a0_cheb, *a1_cheb, *a2_cheb));
        RooChebychev bkg_cheb("bkg_cheb", "Background_cheb", HFmass, RooArgSet(*a0_cheb, *a1_cheb));

        /////////////////////////////////////////////////////////////////////////////////////////////////////////
        // Fractions
        ///////////////////////////////////////////////////

        // RooRealVar nsig2("nsig2", "fraction of component 2 in signal", 0.6, 0., 1.);
        RooRealVar *nsig = new RooRealVar("nsig", "fraction in signal", 10000, 0., 10000000.);

        // RooRealVar nsig2("nsig2", "fraction of component 2 in signal", 0.6, 0., 1.);
        RooRealVar *nsig_nosec = new RooRealVar("nsig_nosec", "fraction in signal", 500, 0., 1000000.);
        RooRealVar *nbkg = new RooRealVar("nbkg", "fraction of background", 1000, 0., 1000000);
        RooRealVar *nbkg_nosec = new RooRealVar("nbkg_nosec", "fraction of background", 200, 0., 1000000);
        RooRealVar *ntanh = new RooRealVar("ntanh", "fraction of background", 200, 0., 1000000);
        RooRealVar *nres = new RooRealVar("nres", "fraction of background", 20, 0., 0.07 * nsig->getVal());
        // RooFormulaVar nres("nres", "resonant bkg", "0.0384*nsig", RooArgList(nsig)); // CHANGE
        RooRealVar *nres_nosec = new RooRealVar("nres_nosec", "fraction of background", 200, 0., 1000000);
        // S a m p l e ,   f i t   a n d   p l o t   m o d e l
        // ---------------------------------------------------
        HFmass.setRange("comb1", 5.55, 5.78);
        HFmass.setRange("comb2", RooConst(mass_low), RooConst(mass_low + 0.1));
        HFmass.setRange("tanh", 5.05, 5.2);
        HFmass.setRange("noSec", 5.18, 5.6);
        HFmass.setRange("Full", mass_low, mass_high);

        // Fit model to data
        RooAddPdf *model, *model_nosec;
        // if(isData) model =new RooAddPdf("model", "g1+g2+a", RooArgList(bkg, sig, tanhpdf), RooArgList(bkgfrac,tanhfrac));
        if (isData)
        {
            if (UseDTF)
            {
                if (DoSystematic)
                {
                    model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg_cheb, sig2, dcbPdf), RooArgList(*nbkg, *nsig, *nres));
                }
                else
                {
                    model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg_cheb, sig, dcbPdf), RooArgList(*nbkg, *nsig, *nres));
                    // model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg_cheb, sig), RooArgList(nbkg, nsig));
                }
            }
            else
                model = new RooAddPdf("model", "g1+g2+a", RooArgList(bkg, sig, dcbPdf), RooArgList(*nbkg, *nsig, *nres));
        }
        else
        {
            model = new RooAddPdf("model", "g1+g2+a", RooArgList(sig), RooArgList(*nsig));
            // model = new RooAddPdf("model", "g1+g2+a", RooArgList(sig2), RooArgList(*nsig));
        }

        //model->fitTo(B_mass, PrintEvalErrors(-1), Save(true));
        RooFitResult *fitresult = model->fitTo(*data_pt_bin, RooFit::Extended(), RooFit::Save(), PrintEvalErrors(-1));
        RooStats::SPlot sPlot(Form("sPlot_pt_%d", i), "SPlot", *data_pt_bin, model, RooArgList(*nbkg, *nsig, *nres));
        int counter = 0;
        for (int j = 0; j < NumEvts; j++)
        {
            BTree->GetEntry(j);
            if (jet_rap < etaMin || jet_rap > etaMax)
                continue;
            if (jet_pt > ptMax)
                continue;
            if (jet_pt < ptMin)
                continue;
            if (!TOS)
                continue;
            if (DoRecSelEff)
            {
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

            TLorentzVector HFmeson, HFjet, mum, mup, Kmeson, Jpsi, tr_HFmeson;
            mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
            mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);
            Kmeson.SetPxPyPzE(K_px, K_py, K_pz, K_e);
            // HFmeson = mup + mum + Kmeson;
            HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
            tr_HFmeson.SetPxPyPzE(tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e);

            // if(dtf_chi2ndf > 5) continue;
            // cout<<bmass_dtf<<endl;
            float bmass = HFmeson.M();
            if (UseDTF)
                bmass = bmass_dtf;
            bool PID_cond = (K_PIDK > 0);
            bool DTF_cond = (chi2ndf_dtf < 9) && (tau_dtf > 0.3);

            if (!DTF_cond && DTF_cut)
                continue;
            if (!PID_cond && PID_cut)
                continue;
            if (!isData)
            {
                if (!isTrueBjet)
                    continue;
            }

            if (bmass < mass_low || bmass > mass_high)
                continue;

            if (HFmeson.Pt() > ptHF_binedges[i] && HFmeson.Pt() < ptHF_binedges[i + 1])
            {
                // cout << HFmeson.Pt() << ", " << "Filling" << ", " << data_pt_bin->numEntries() << ", " << counter << endl;
                sWeights[j] = sPlot.GetSWeight(counter, "nsig"); // Signal weight
                counter++;
            }
        }
        // RooStats::SPlot sData("sData", "SPlot with signal and background weights", *data, &model, RooArgList(nsig, nbkg));

        // RooMsgService::instance().setSilentMode(true);
        // RooStats::SPlot splot("splot", "splot", *data, &sum, RooArgSet(c1, c2));
        // Begin plotting
        //
        ican = 0;
        ccan[ican]->cd(i + 1);
        gPad->SetRightMargin(0.12);
        gPad->SetLeftMargin(0.15);
        ccan[ican]->SetFillStyle(4070);
        RooPlot *xframe = HFmass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
        xframe->SetTitleOffset(0.9, "Y");
        xframe->SetTitleOffset(0.9, "X");
        // model->paramOn(xframe, Layout(0.6, 0.9, 0.9));
        // xframe->getAttText()->SetTextSize(0.02);
        B_mass.plotOn(xframe, Name("B_mass"));
        // model->plotOn(xframe);

        // Overlay the background component of model with a dashed line
        model->plotOn(xframe, Name("sig"), Components(sig), LineStyle(5), LineColor(kGreen + 2), LineWidth(1.));
        // sig.plotOn(xframe, Name("sig"), LineStyle(kSolid), LineColor(kGreen), Range("Full"));
        if (DoSystematic)
            model->plotOn(xframe, Name("sig2"), Components(sig2), LineStyle(5), LineColor(kGreen + 2), LineWidth(1.));
        model->plotOn(xframe, Name("tanhpdf"), Components(tanhpdf), FillColor(TColor::GetColorTransparent(kCyan, 0.7)), FillStyle(4070), LineWidth(0), DrawOption("F"));
        model->plotOn(xframe, Name("argxgauss"), Components(argxgauss), FillColor(kCyan), DrawOption("F"));
        model->plotOn(xframe, Name("bkg"), Components(bkg), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
        model->plotOn(xframe, Name("bkg_cheb"), Components(bkg_cheb), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
        // model->plotOn(xframe, Name("exp_bkg"), Components(exp_bkg), LineStyle(kDashed), LineColor(kBlue));
        // model->plotOn(xframe, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));
        model->plotOn(xframe, Name("dcbPdf"), Components(dcbPdf), LineStyle(kDashed), LineColor(kMagenta), LineWidth(2.));

        model->plotOn(xframe, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1.));
        RooHist *hpull = xframe->pullHist();
        xframe->SetMaximum(30000);
        float chi2;
        if (isData)
            chi2 = xframe->chiSquare(10);
        else
            chi2 = xframe->chiSquare(10);

        // exp_bkg.plotOn(xframe, LineColor(kMagenta));

        // Overlay the background+sig2 components of model with a dotted line
        // model->plotOn(xframe, Components(RooArgSet(bkg, sig2)), LineStyle(kDotted));
        // h1_mass_all->Draw("SAME");
        xframe->Draw("SAME");
        gPad->SetLogy();
        xframe->SetMinimum(1.);

        TLegend *leg1 = new TLegend(0.5, 0.45, 0.6, 0.87);
        leg1->SetTextSize(0.07);
        // leg1->SetFillColor(kWhite);
        // leg1->SetLineColor(kWhite);
        leg1->SetBorderSize(0);
        leg1->AddEntry("B_mass", "Data", "PE");
        leg1->AddEntry("model", "Fit model", "L");
        leg1->AddEntry("sig", "Signal", "L");
        leg1->AddEntry("bkg", "Comb bkg", "L");
        if (isData)
        {
            // leg1->AddEntry("tanhpdf", "Part. Rec. bkg.", "F");
            leg1->AddEntry("dcbPdf", "misID bkg", "L");
        }
        // leg1->AddEntry((TObject*)0, "", "");
        leg1->AddEntry((TObject *)0, Form("#chi^{2}/ndf = %.2f", chi2), "");
        leg1->AddEntry((TObject *)0, Form("%.1f < p_{T}^{HF} < %.1f GeV", ptHF_binedges[i], ptHF_binedges[i + 1]), "");

        leg1->Draw("SAME");

        ican = 1;
        ccan[ican]->cd(i + 1);
        gPad->SetRightMargin(0.12);
        gPad->SetLeftMargin(0.12);
        ccan[ican]->SetFillStyle(4070);
        RooPlot *xframe2 = HFmass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
        xframe2->SetTitleOffset(0.7, "Y");
        xframe2->SetTitleOffset(0.9, "X");
        B_mass.plotOn(xframe2, Name("B_mass"));
        // model->plotOn(xframe2);

        // Overlay the background component of model with a dashed line
        model->plotOn(xframe2, Name("sig"), Components(sig), LineStyle(5), LineColor(kGreen + 2), LineWidth(1.));
        // sig.plotOn(xframe2, Name("sig"), LineStyle(kSolid), LineColor(kGreen), Range("Full"));
        // model->plotOn(xframe2, Components(sig2), LineStyle(kSolid), LineColor(kMagenta));
        model->plotOn(xframe2, Name("tanhpdf"), Components(tanhpdf), FillColor(TColor::GetColorTransparent(kCyan, 0.7)), FillStyle(4070), LineWidth(0), DrawOption("F"));
        model->plotOn(xframe2, Name("argxgauss"), Components(argxgauss), FillColor(kCyan), DrawOption("F"));
        model->plotOn(xframe2, Name("bkg"), Components(bkg), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
        model->plotOn(xframe2, Name("bkg_cheb"), Components(bkg_cheb), LineStyle(kDashed), LineColor(kBlue), LineWidth(1.));
        // model->plotOn(xframe2, Name("exp_bkg"), Components(exp_bkg), LineStyle(kDashed), LineColor(kBlue));
        // model->plotOn(xframe2, Name("dcbPdf"), Components(dcbPdf), FillColor(TColor::GetColorTransparent(kMagenta, 0.7)), FillStyle(3165), LineWidth(0), DrawOption("F"));

        model->plotOn(xframe2, Name("dcbPdf"), Components(dcbPdf), LineStyle(kDashed), LineColor(kMagenta), LineWidth(2.));

        model->plotOn(xframe2, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1.));

        // exp_bkg.plotOn(xframe2, LineColor(kMagenta));

        // Overlay the background+sig2 components of model with a dotted line
        // model->plotOn(xframe2, Components(RooArgSet(bkg, sig2)), LineStyle(kDotted));
        // h1_mass_all->Draw("SAME");
        xframe2->Draw("SAME");
        // gPad->SetLogy();
        xframe2->SetMinimum(0.);
        leg1->Draw("SAME");

        RooWorkspace *w = new RooWorkspace(Form("w%d", i), Form("workspace%d", i));
        // w->import(sign);
        // w->import(bkgn);
        w->import(*model);
        w->import(B_mass);
        w->import(sPlot);
        w->Print();
        w->writeToFile(extension_RootFiles + Form("workspace%d_", i) + extension + ".root");

        cout << "Chi2/dof = " << chi2 << endl;
        cout << sigma_arg->getVal() << endl;
        cout << "a1 = " << a1->getVal() << endl;
        cout << "a2 = " << a2->getVal() << endl;
        cout << "sigma1 = " << sigma1->getVal() << endl;
        cout << "sigma2 = " << sigma2->getVal() << endl;

        float MassMu;

        float MassSigma;
        // cout << "Nsig1 = " << nsig1->getVal() << endl;
        if (DoSystematic)
        {
            MassMu = sig2_m->getVal();
            MassSigma = sig2_t->getVal();
            cout << "MassMu = " << sig2_m->getVal() << endl;
            cout << "MassSigma = " << sig2_t->getVal() << endl;
        }
        else
        {
            MassMu = mu_sig->getVal();
            if (nsig1->getVal() > 0.5)
                MassSigma = sigma1->getVal();
            else
                MassSigma = sigma2->getVal();
            // if (sigma1->getVal() > sigma2->getVal())
            //     MassSigma = sigma1->getVal();
            // else
            //     MassSigma = sigma2->getVal();

            cout
                << "Sigma values = " << MassSigma << ", " << sigma1->getVal() << ", " << sigma2->getVal() << endl;
        }

        float MassUp, MassDown;
        if (i < 2)
        {
            MassDown = MassMu - 2 * MassSigma;
            MassUp = MassMu + 2 * MassSigma;
        }
        else if (i < 5)
        {
            MassDown = MassMu - 3 * MassSigma;
            MassUp = MassMu + 3 * MassSigma;
        }
        else
        {
            MassDown = MassMu - 4 * MassSigma;
            MassUp = MassMu + 4 * MassSigma;
        }
        
        Sideband1_Min = (MassMu - 14 * MassSigma) > 5.15 ? (MassMu - 14 * MassSigma) : 5.15;
        Sideband1_Max = (MassMu - 7 * MassSigma) > Sideband1_Min ? (MassMu - 7 * MassSigma) : 5.17;
        // Sideband2_Min = MassMu + 10 * MassSigma;
        // Sideband2_Max = MassMu + 14 * MassSigma;
        Sideband2_Min = MassMu + 10 * MassSigma;
        Sideband2_Max = (MassMu + 25 * MassSigma) < 5.55 ? (MassMu + 25 * MassSigma) : 5.55;
        

        float Sideband1_Min_forSysNear = Sideband1_Min + MassSigma;
        float Sideband1_Max_forSysNear = Sideband1_Max;
        float Sideband2_Min_forSysNear = Sideband2_Min;
        float Sideband2_Max_forSysNear = Sideband2_Max - MassSigma;

        float Sideband1_Min_forSysFar = Sideband1_Min;
        float Sideband1_Max_forSysFar = Sideband1_Max - MassSigma;
        float Sideband2_Min_forSysFar = Sideband2_Min + MassSigma;
        float Sideband2_Max_forSysFar = Sideband2_Max;

        cout << MassDown << ", " << MassUp << endl;
        cout << "Sideband1: " << Sideband1_Min << ", " << Sideband1_Max << endl;
        cout << "Sideband2: " << Sideband2_Min << ", " << Sideband2_Max << endl;
        HFmass.setRange("signal", MassDown, MassUp);
        HFmass.setRange("sideband1", Sideband1_Min, Sideband1_Max);
        HFmass.setRange("sideband2", Sideband2_Min, Sideband2_Max);

        HFmass.setRange("sideband1_forSysNear", Sideband1_Min_forSysNear, Sideband1_Max_forSysNear);
        HFmass.setRange("sideband2_forSysNear", Sideband2_Min_forSysNear, Sideband2_Max_forSysNear);
        HFmass.setRange("sideband1_forSysFar", Sideband1_Min_forSysFar, Sideband1_Max_forSysFar);
        HFmass.setRange("sideband2_forSysFar", Sideband2_Min_forSysFar, Sideband2_Max_forSysFar);
        // HFmass.setRange("signal", MassMin, MassMax);

        // intsigX = sig->createIntegral(*HFmass, NormSet(*HFmass), Range("signal"));
        // intbkgX = bkg->createIntegral(*HFmass, NormSet(*HFmass), Range("signal"));
        RooAbsReal *intsigX;
        RooAbsReal *intbkgX, *intresX, *intbkgSide, *intbkgSideLeft, *intbkgSideRight;
        intsigX = sig.createIntegral(HFmass, NormSet(HFmass), Range("signal"));
        intbkgX = bkg.createIntegral(HFmass, NormSet(HFmass), Range("signal"));

        intbkgSide = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband1", "sideband2"));
        intbkgSideLeft = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband1"));
        intbkgSideRight = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband2"));

        RooAbsReal *intbkgSide_forSysNear = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband1_forSysNear", "sideband2_forSysNear"));
        RooAbsReal *intbkgSideLeft_forSysNear = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband1_forSysNear"));
        RooAbsReal *intbkgSideRight_forSysNear = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband2_forSysNear"));
        RooAbsReal *intbkgSide_forSysFar = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband1_forSysFar", "sideband2_forSysFar"));
        RooAbsReal *intbkgSideLeft_forSysFar = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband1_forSysFar"));
        RooAbsReal *intbkgSideRight_forSysFar = bkg.createIntegral(HFmass, NormSet(HFmass), Range("sideband2_forSysFar"));

        intresX = dcbPdf.createIntegral(HFmass, NormSet(HFmass), Range("signal"));
        float sig_yield = intsigX->getVal() * nsig->getVal();
        float bkg_yield = intbkgX->getVal() * nbkg->getVal();
        float res_yield = intresX->getVal() * nres->getVal();
        float bkg_frac = (bkg_yield) / (sig_yield + bkg_yield + res_yield);
        float res_frac = (res_yield) / (sig_yield + bkg_yield + res_yield);
        float sig_frac = (sig_yield) / (sig_yield + bkg_yield + res_yield);

        float bkg_scale = intbkgX->getVal() / intbkgSide->getVal();
        float bkg_scale_forSysNear = intbkgX->getVal() / intbkgSide_forSysNear->getVal();
        float bkg_scale_forSysFar = intbkgX->getVal() / intbkgSide_forSysFar->getVal();
        // cout << "Bkg Yield = " << intbkgX->getVal() << endl;
        // cout << "Sig Yield = " << intsigX->getVal() << endl;
        cout << "Signal Yield = " << sig_yield << endl;
        cout << "Bkg Fraction = " << bkg_frac << endl;
        cout << "Bkg Fraction = " << res_frac << endl;
        cout << "Area_Sig / Area_Sideband = " << bkg_scale << endl;
        cout << "Area left side / Area right side = " << intbkgSideLeft->getVal() / intbkgSideRight->getVal() << endl;

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

        // Add a vertical dashed line at x = 2

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

        if (!DoSystematic)
        {
            SigParams
                << mu_sig->getVal() << "," << sigma1->getVal() << "," << sigma2->getVal() << "," << alpha1_sig->getVal() << "," << p1_sig->getVal() << "," << "," << alpha2_sig->getVal() << "," << p2_sig->getVal() << "," << "," << nsig->getVal() << ", " << sig_frac << endl;
        }
        else
            SigParams
                << sig2_m->getVal() << "," << sig2_t->getVal() << "," << sig2_n->getVal() << "," << nsig->getVal() << ", " << sig_frac << endl;

        BkgParams << a0_cheb->getVal() << "," << a1_cheb->getVal() << "," << a2_cheb->getVal() << "," << nbkg->getVal() << ", " << bkg_frac << endl;
        ResonantParams << mu->getVal() << "," << width->getVal() << "," << alpha1->getVal() << "," << p1->getVal() << "," << alpha2->getVal() << "," << p2->getVal() << "," << nres->getVal() << "," << res_frac << endl;
    }

    ican = 0;
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

    ican = 1;
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

    if (ican > -1)
    {
        cout << " You plotted " << ican + 1 << " canvasses......." << endl;
        ccan[ican]->Print(plotfileC.Data());
    }

    for (int i = 0; i < vec_bkg_frac.size(); i++)
    {
        cout << ptHF_binedges[i] << " < pT < " << ptHF_binedges[i + 1] << ": Signal = " << vec_sig_yield[i] << ", Background = " << vec_bkg_yield[i]
             << ", Bkg Frac = " << vec_bkg_frac[i] * 100 << " %"
             << ", Res Frac = " << vec_res_frac[i] * 100 << " %" << endl;
    }
    float sigyield = std::accumulate(vec_sig_yield.begin(), vec_sig_yield.end(), 0);
    cout << "Total Signal Yield = " << sigyield << endl;

    // cout<<exp_c->getVal()<<endl;

    SigParams.close();
    BkgParams.close();
    ResonantParams.close();

    f.cd();
    // Assign weights back to original TTree
    for (int j = 0; j < NumEvts; j++)
    {
        BTree->GetEntry(j);
        sweight = sWeights[j];
        tree->Fill();
    }
    tree->Write();
    // BTree->Write();

    f.Write();
    f.Close();
}
