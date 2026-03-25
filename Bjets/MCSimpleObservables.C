#include <iostream>
#include <TCanvas.h>
#include <vector>
#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"

//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/contrib/SoftDrop.hh"
//#include "../LundGen.hh"
//using namespace fastjet;
using namespace std;

void MCSimpleObservables(int NumEvts = -1, int dataset = 91599,
                  bool chargedJetCut_user = false,
                  bool SubtractGS = false,
                  bool onlysim9 = false)
{
        // This should be the file from MCMakeVarTree
        std::string extension = "simpleobservable_mc_";
        
        std::string extension_RootFilesMC = output_folder + "bjets-mc/";

        TFile fread((output_folder + "ntuple_bjets_mc.root").c_str(), "READ");
        TTree *BTree = (TTree *)fread.Get("BTree");
        
        if (NumEvts > BTree->GetEntries())
                NumEvts = BTree->GetEntries();
        
        if (NumEvts == -1)
                NumEvts = BTree->GetEntries();
        
        cout<<BTree->GetEntries()<<endl;

        TFile f((output_folder + "bjets_simpleobservable_mc.root").c_str(), "RECREATE");

        TH2D *h2_ptz  = new TH2D("ptz", ";z;p_{T,jet} [GeV/c]", zbinsize, z_binedges, ptbinsize, pt_binedges);
        TH2D *h2_ptjt = new TH2D("ptjt", ";j_{T} [GeV/c]; p_{T,jet} [GeV/c]", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
        TH2D *h2_ptr  = new TH2D("ptr", ";r; p_{T,jet} [GeV/c]", rbinsize, r_binedges, ptbinsize, pt_binedges);
        
        TH2D *h2_ptz_gluon        = new TH2D("h2_ptz_gluon", ";z;p_{T,jet} [GeV/c]", zbinsize, z_binedges, ptbinsize, pt_binedges);  
        TH2D *h2_ptz_gluon_ratio  = (TH2D*)h2_ptz_gluon->Clone("h2_ptz_gluon_ratio");
        
        TH2D *h2_ptjt_gluon       = new TH2D("h2_ptjt_gluon", ";j_{T} [GeV/c]; p_{T,jet} [GeV/c]", jtbinsize, jt_binedges, ptbinsize, pt_binedges);  
        TH2D *h2_ptjt_gluon_ratio = (TH2D*)h2_ptjt_gluon->Clone("h2_ptjt_gluon_ratio"); 
        
        TH2D *h2_ptr_gluon        = new TH2D("h2_ptr_gluon", ";r; p_{T,jet} [GeV/c]", rbinsize, r_binedges, ptbinsize, pt_binedges);  
        TH2D *h2_ptr_gluon_ratio  = (TH2D*)h2_ptr_gluon->Clone("h2_ptr_gluon_ratio");         
        
        TH3D *h3_ptzjt = new TH3D("ptzjt", ";z;j_{T} [GeV/c];p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_ptzr  = new TH3D("ptzr", ";z;r;p_{T,jet} [GeV/c]",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_ptjtr = new TH3D("ptjtr", ";j_{T} [GeV/c];r;p_{T,jet} [GeV/c]",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        
        TH3D *h3_ptzjt_gluon = new TH3D("ptzjt_gluon", ";z;j_{T} [GeV/c];p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_ptzr_gluon  = new TH3D("ptzr_gluon", ";z;r;p_{T,jet} [GeV/c]",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_ptjtr_gluon = new TH3D("ptjtr_gluon", ";j_{T} [GeV/c];r;p_{T,jet} [GeV/c]",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        
        //    TH1D *h1_frag_jt = new TH1D("frag_jt", "", 30, 0, 12);
        //    TH1D *h1_frag_z = new TH1D("frag_z", "", zbinsize, z_binedges);
        //    TH1D *h1_frag_r = new TH1D("frag_r", "", 30, 0, 1);

        TH1D *h1_frag_jt_gluon = new TH1D("frag_jt_gluon", "", 30, 0, 12);
        TH1D *h1_frag_z_gluon  = new TH1D("frag_z_gluon", "", zbinsize, z_binedges);
        TH1D *h1_frag_r_gluon  = new TH1D("frag_r_gluon", "", 30, 0, 1);

        TH2D *h2_frag_z_jetpt = new TH2D("frag_z_jetpt", "", zbinsize, z_binedges, ptbinsize, pt_binedges);

        TH1D *h1_z  = new TH1D("z", "", zbinsize, z_binedges);
        TH1D *h1_jt = new TH1D("jt", "", jtbinsize, jt_binedges);
        TH1D *h1_r  = new TH1D("r", "", rbinsize, r_binedges);
        
        TH1D *h1_meas_z  = new TH1D("meas_z", "", zbinsize, z_binedges);
        TH1D *h1_meas_jt = new TH1D("meas_jt", "", jtbinsize, jt_binedges);
        TH1D *h1_meas_r  = new TH1D("meas_r", "", rbinsize, r_binedges);
        
        
        TH2D *h2_zjt = new TH2D("zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        TH2D *h2_zr  = new TH2D("zr", "", zbinsize, z_binedges,  rbinsize, r_binedges);
        TH2D *h2_jtr = new TH2D("jtr", "", jtbinsize, jt_binedges,  rbinsize, r_binedges);
        
        TH2D *h2_zjt_gluon = new TH2D("zjt_gluon", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        TH2D *h2_zr_gluon  = new TH2D("zr_gluon", "", zbinsize, z_binedges,  rbinsize, r_binedges);
        TH2D *h2_jtr_gluon = new TH2D("jtr_gluon", "", jtbinsize, jt_binedges,  rbinsize, r_binedges);
        
        TH2D *h2_meas_zjt = new TH2D("meas_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        TH2D *h2_meas_zr  = new TH2D("meas_zr", "", zbinsize, z_binedges,  rbinsize, r_binedges);
        TH2D *h2_meas_jtr = new TH2D("meas_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
        

        // 2D Truth-Reco Correspondence (219 - 224)
        TH2D *h2_truthreco_z  = new TH2D("truthreco_z", ";Reco z; Truth z", zbinsize, z_binedges, zbinsize, z_binedges);
        TH2D *h2_truthreco_jt = new TH2D("truthreco_jt", ";Reco jT; Truth jT", jtbinsize, jt_binedges, jtbinsize, jt_binedges);
        TH2D *h2_truthreco_r  = new TH2D("truthreco_r", ";Reco r; Truth r", rbinsize, r_binedges, rbinsize, r_binedges);
        
        TH1D *h1_ratio_z0  = new TH1D("ratio_z0", ";truth z / reco z;", 30, 0, 2);
        TH1D *h1_ratio_jt0 = new TH1D("ratio_jt0", ";truth jt / reco jt;", 30, 0, 2);
        TH1D *h1_ratio_r0  = new TH1D("ratio_r0", ";truth r / reco r;", 30, 0, 2);
        TH1D *h1_ratio_z1  = new TH1D("ratio_z1", ";reco z / truth z;", 30, 0, 2);
        TH1D *h1_ratio_jt1 = new TH1D("ratio_jt1", ";reco jt / truth jt;", 30, 0, 2);
        TH1D *h1_ratio_r1  = new TH1D("ratio_r1", ";reco r / truth r;", 30, 0, 2);
        
        TH1D *h1_jet_flav = new TH1D("Jet_Flav", "", 7, -0.5, 6.5);
        TH1D *h1_jet_pt   = new TH1D("Jet_pT", "", ptbinsize, pt_binedges);
        TH1D *h1_jet_eta  = new TH1D("Jet_eta", "", 12, etaMin, etaMax);
        TH1D *h1_jet_rap  = new TH1D("Jet_rap", "", 12, etaMin, etaMax);
        TH1D *h1_jet_phi  = new TH1D("Jet_phi", "", 20, -3.14, 3.14);
        
        TH1D *h1_jet_ptbalance0 = new TH1D("jet_ptbalance0", "", 20, 0, 2);
        TH1D *h1_jet_ptbalance1 = new TH1D("jet_ptbalance1", "", 20, 0, 2);

        TH2D *h2_jetpteta             = new TH2D("h2_jetpteta", ";p_{T,jet} [GeV/c]; #eta", ptbinsize, pt_binedges, etabinsize, eta_binedges);
        TH2D *h2_jetpteta_gluon       = new TH2D("h2_jetpteta_gluon", ";p_{T,jet} [GeV/c]; #eta", ptbinsize, pt_binedges, etabinsize, eta_binedges);  
        TH2D *h2_jetpteta_gluon_ratio = (TH2D*)h2_jetpteta_gluon->Clone("h2_jetpteta_gluon_ratio");
        
        TH2D *h2_jetptp             = new TH2D("h2_jetptp", ";p_{T,jet} [GeV/c]; p_{jet} [GeV/c]", ptbinsize, pt_binedges, pbinsize, p_binedges);
        TH2D *h2_jetptp_gluon       = new TH2D("h2_jetptp_gluon", "p_{T,jet} [GeV/c]; p_{jet} [GeV/c]", ptbinsize, pt_binedges, pbinsize, p_binedges);  
        TH2D *h2_jetptp_gluon_ratio = (TH2D*)h2_jetptp_gluon->Clone("h2_jetptp_gluon_ratio");  
        
        TH1D *h1_meas_jet_pt  = new TH1D("meas_Jet_pT", "", ptbinsize, pt_binedges);
        TH1D *h1_meas_jet_eta = new TH1D("meas_Jet_eta", "", 12, etaMin, etaMax);
        TH1D *h1_meas_jet_rap = new TH1D("meas_Jet_rap", "", 12, etaMin, etaMax);
        TH1D *h1_meas_jet_phi = new TH1D("meas_Jet_phi", "", 20, -3.14, 3.14);

        TH1D *h1_jet_pt_noghost  = new TH1D("Jet_pT_noghost", "", 50, ptMin, ptMax);
        TH1D *h1_jet_eta_noghost = new TH1D("Jet_eta_noghost", "", 12, etaMin, etaMax);
        TH1D *h1_jet_phi_noghost = new TH1D("Jet_phi_noghost", "", 20, -3.14, 3.14);

        TH1D *h1_Jpsi_rap  = new TH1D("Jpsi_rap", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_Jpsi_pt   = new TH1D("Jpsi_pT", "", 50, 0, 100);
        TH1D *h1_Jpsi_phi  = new TH1D("Jpsi_phi", "", 20, -3.14, 3.14);
        TH1D *h1_Jpsi_mass = new TH1D("Jpsi_mass", "", 30, 3.1 - 0.1, 3.1 + 0.1);

        TH1D *h1_meas_Jpsi_rap  = new TH1D("meas_Jpsi_rap", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_meas_Jpsi_pt   = new TH1D("meas_Jpsi_pT", "", 50, 0, 100);
        TH1D *h1_meas_Jpsi_phi  = new TH1D("meas_Jpsi_phi", "", 20, -3.14, 3.14);
        TH1D *h1_meas_Jpsi_mass = new TH1D("meas_Jpsi_mass", "", 30, 3.1 - 0.1, 3.1 + 0.1);

        TH1D *h1_mup_eta = new TH1D("mup_eta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_mup_pt  = new TH1D("mup_pt", "", 40, 0, 20);
        TH1D *h1_mum_eta = new TH1D("mum_eta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_mum_pt  = new TH1D("mum_pt", "", 40, 0, 20);
        TH1D *h1_K_eta   = new TH1D("K_eta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_K_pt    = new TH1D("K_pt", "", 40, 0, 20);

        TH1D *h1_meas_mup_eta = new TH1D("meas_mup_eta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_meas_mup_pt  = new TH1D("meas_mup_pt", "", 40, 0, 20);
        TH1D *h1_meas_mum_eta = new TH1D("meas_mum_eta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_meas_mum_pt  = new TH1D("meas_mum_pt", "", 40, 0, 20);
        TH1D *h1_meas_K_eta   = new TH1D("meas_K_eta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_meas_K_pt    = new TH1D("meas_K_pt", "", 40, 0, 20);

        TH1D *h1_HF_rap          = new TH1D("HF_rap", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_HF_pt           = new TH1D("HF_pT", "", 50, 0, 100);
        TH1D *h1_HF_phi          = new TH1D("HF_phi", "", 20, -3.14, 3.14);
        TH1D *h1_HF_mass         = new TH1D("HF_mass", "", 30, 5.279 - 0.3, 5.279 + 0.3);
        TH1D *h1_HFjet_ptbalance = new TH1D("HFjet_ptbalance", "", 20, 0, 2);

        TH1D *h1_HFpt_GS    = new TH1D("HFpt_GS", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_HFpt_FC    = new TH1D("HFpt_FC", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_HFpt_Total = new TH1D("HFpt_Total", "", ptHFbinsize, ptHF_binedges);

        TH1D *h1_jetpt_GS    = new TH1D("jetpt_GS", "", ptbinsize, pt_binedges);
        TH1D *h1_jetpt_FC    = new TH1D("jetpt_FC", "", ptbinsize, pt_binedges);
        TH1D *h1_jetpt_Total = new TH1D("jetpt_Total", "", ptbinsize, pt_binedges);

        TH1D *h1_HFpt      = new TH1D("h1_HFpt", "", ptHFbinsize, ptHF_binedges);
        TH2D *h2_HFptjetpt = new TH2D("h2_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);

        TH1D *h1_HF_measeta = new TH1D("HF_measeta", "", 12, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_HF_measpt  = new TH1D("HF_measpT", "", 50, 0, 80);
        TH1D *h1_HF_measphi = new TH1D("HF_measphi", "", 20, -3.14, 3.14);

        TH1D *h1_Dtr_Eta = new TH1D("Dtr_Eta", "Dtr #eta", 50, etaMin - jetradius, etaMax + jetradius);
        TH1D *h1_dtr_phi = new TH1D("Dtr_phi", "Dtr #phi", 50, -3.14, 3.14);
        TH1D *h1_dtr_pT  = new TH1D("Dtr_pT", "Dtr pT", 50, 0, 10);
        TH1D *h1_dtr_p   = new TH1D("Dtr_p", "Dtr p", 50, 0, 150);
        TH1D *h1_keys    = new TH1D("Keys", "Keys", 100, 0, 100);
        
        TH1D *h1_hfdtrs_e_diff          = new TH1D("HFDtrsEDiff", "HFDtrsEDiff", 40, -100, 100);
        TH1D *h1_HFHadronMass           = new TH1D("HFHadronMass", "HFHadronMass", 100, 60, 120);
        TH1D *h1_HF_jet_dR              = new TH1D("HF_jet_dR", "", 24, 0, 4.);
        TH1D *h1_HF_jet_truedR          = new TH1D("HF_jet_truedR", "", 24, 0, 4.);
        TH1D *h1_HFdtr1_mass            = new TH1D("HFdtr1mass", "", 100, 1, 7);
        TH1D *h1_recoHFHadronMass       = new TH1D("recoHFHadronMass", "recoHFHadronMass", 40, 0, 15);
        TH1D *h1_nsplittings            = new TH1D("nSplittings", "nSplittings", 15, -0.5, 14.5);
        TH1D *h1_ptHF                   = new TH1D("ptHF", "pT_{HF}", 30, 0, 25);
        TH1D *h1_probNNghost            = new TH1D("probNNghost", "probNNghost", 20, 0, 1);
        TH1D *h1_nJetDtrs_noghost       = new TH1D("nJetDtrs_noghost", "nJetDtrs_noghost", 40, -0.5, 39.5);
        TH1D *h1_nJetDtrs_noghost_gluon = new TH1D("nJetDtrs_noghost_gluon", "nJetDtrs_noghost", 40, -0.5, 39.5);
        TH1F *h_nJetDaughters           = new TH1F("h_nJetDaughters", "Number of jet daughters", 40, -0.5, 39.5);

        TH2D *HFmesonyphi_h = new TH2D("HFmesonyphi_h", ";#phi [rad]; y", 100, -1 * TMath::Pi(), TMath::Pi(), 50, 1, 5);
        TH2D *HFjetetaphi_h = new TH2D("truejetetaphi_h", ";#phi [rad]; #eta", 100, -1 * TMath::Pi(), TMath::Pi(), 50, 1, 5);
        TH1D *HFjetpt_h     = new TH1D("truejetpt_h", ";p_{T}^{jet,MC} [GeV]", 50, 0, 200);
        TH1D *zmass         = new TH1D("zmass", ";M_{#mu#mu} [GeV]", 150, 0, 150);
        TH1D *truedphi_h    = new TH1D("truedphi_h", ";#Delta#phi [rad]", 50, 0, 3.14159);
        TH1D *truexj_h      = new TH1D("truexj_h", "x_{J}", 50, 0, 2);
        TH1D *Zpt_h         = new TH1D("Zpt_h", "p_{T}^{Z} [GeV]", 50, 0, 200);
        TH1D *truedy_h      = new TH1D("truedy_h", "#Deltay", 50, -3, 3);
        TH1D *truenconst_h = new TH1D("truenconst_h", ";N_{const}^{chg}", 40, -0.5, 39.5);

        //
        // Event loop
        //
        unsigned long long last_eventNum = 0;
        double last_jetpT = 0.;
        int event_counter = 0;

        double jet_pt, jet_eta, meas_jet_pt, meas_jet_eta;
        double jet_rap, meas_jet_rap;
        double jet_px, jet_py, jet_pz, jet_e;
        
        double mup_px, mup_py, mup_pz, mup_e;
        double mum_px, mum_py, mum_pz, mum_e;
        double K_px, K_py, K_pz, K_e;
        double HF_px, HF_py, HF_pz, HF_e, HF_pt;

        double meas_jet_px, meas_jet_py, meas_jet_pz, meas_jet_e;
        double meas_HF_px, meas_HF_py, meas_HF_pz, meas_HF_e;
        double meas_HF_pt;
        double meas_mup_px, meas_mup_py, meas_mup_pz, meas_mup_e;
        double meas_mum_px, meas_mum_py, meas_mum_pz, meas_mum_e;
        double meas_K_px, meas_K_py, meas_K_pz, meas_K_e;
        
        double truth_z, truth_jt, truth_r;
        double truth_z_b, truth_jt_b, truth_r_b;
        double truth_zg, truth_jtg, truth_rg;
        double meas_z, meas_jt, meas_r;
        double meas_inJet_z, meas_inJet_jt, meas_inJet_r;
        double jet_pt_recotruthratio, HF_pt_recotruthratio;

        int nsplits, ndtrs;
        int meas_nsplits, meas_ndtrs;
        int GluonTag, nTracks, NumHFHads;
        int eventNumber;
        bool Hasbbbar, hasRecoHF;
        // TLorentzVector
        
        BTree->SetBranchAddress("truth_z", &truth_z);
        BTree->SetBranchAddress("truth_jt", &truth_jt);
        BTree->SetBranchAddress("truth_r", &truth_r);
        BTree->SetBranchAddress("truth_z_b", &truth_z_b);
        BTree->SetBranchAddress("truth_jt_b", &truth_jt_b);
        BTree->SetBranchAddress("truth_r_b", &truth_r_b);
        BTree->SetBranchAddress("truth_zg", &truth_zg);
        BTree->SetBranchAddress("truth_jtg", &truth_jtg);
        BTree->SetBranchAddress("truth_rg", &truth_rg);
        
        BTree->SetBranchAddress("jet_pt", &jet_pt);
        BTree->SetBranchAddress("jet_eta", &jet_eta);
        BTree->SetBranchAddress("jet_rap", &jet_rap);
        
        BTree->SetBranchAddress("jet_px", &jet_px);
        BTree->SetBranchAddress("jet_py", &jet_py);
        BTree->SetBranchAddress("jet_pz", &jet_pz);
        BTree->SetBranchAddress("jet_e", &jet_e);
        
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

        BTree->SetBranchAddress("meas_z", &meas_z);
        BTree->SetBranchAddress("meas_jt", &meas_jt);
        BTree->SetBranchAddress("meas_r", &meas_r);
        BTree->SetBranchAddress("meas_inJet_z", &meas_inJet_z);
        BTree->SetBranchAddress("meas_inJet_jt", &meas_inJet_jt);
        BTree->SetBranchAddress("meas_inJet_r", &meas_inJet_r);
        
        BTree->SetBranchAddress("meas_jet_pt", &meas_jet_pt);
        BTree->SetBranchAddress("meas_jet_eta", &meas_jet_eta);
        BTree->SetBranchAddress("meas_jet_rap", &meas_jet_rap);
        
        BTree->SetBranchAddress("meas_jet_px", &meas_jet_px);
        BTree->SetBranchAddress("meas_jet_py", &meas_jet_py);
        BTree->SetBranchAddress("meas_jet_pz", &meas_jet_pz);
        BTree->SetBranchAddress("meas_jet_e", &meas_jet_e);
        
        BTree->SetBranchAddress("meas_HF_px", &meas_HF_px);
        BTree->SetBranchAddress("meas_HF_py", &meas_HF_py);
        BTree->SetBranchAddress("meas_HF_pz", &meas_HF_pz);
        BTree->SetBranchAddress("meas_HF_e", &meas_HF_e);
        BTree->SetBranchAddress("meas_HF_pt", &meas_HF_pt);
        
        BTree->SetBranchAddress("NumHFHads", &NumHFHads);
        BTree->SetBranchAddress("hasRecoHF", &hasRecoHF);

        BTree->SetBranchAddress("HF_px", &HF_px);
        BTree->SetBranchAddress("HF_py", &HF_py);
        BTree->SetBranchAddress("HF_pz", &HF_pz);
        BTree->SetBranchAddress("HF_e", &HF_e);
        BTree->SetBranchAddress("HF_pt", &HF_pt);
        
        BTree->SetBranchAddress("Hasbbbar", &Hasbbbar);
        BTree->SetBranchAddress("eventNumber", &eventNumber);

        BTree->Branch("meas_mum_px", &meas_mum_px);
        BTree->Branch("meas_mum_py", &meas_mum_py);
        BTree->Branch("meas_mum_pz", &meas_mum_pz);
        BTree->Branch("meas_mum_e", &meas_mum_e);
        BTree->Branch("meas_mup_px", &meas_mup_px);
        BTree->Branch("meas_mup_py", &meas_mup_py);
        BTree->Branch("meas_mup_pz", &meas_mup_pz);
        BTree->Branch("meas_mup_e", &meas_mup_e);
        BTree->Branch("meas_K_px", &meas_K_px);
        BTree->Branch("meas_K_py", &meas_K_py);
        BTree->Branch("meas_K_pz", &meas_K_pz);
        BTree->Branch("meas_K_e", &meas_K_e);
        
        BTree->SetBranchAddress("GluonTag", &GluonTag);
        BTree->SetBranchAddress("nTracks", &nTracks);
        //    BTree->SetBranchAddress("NumDtrRecoHF", &NumDtrRecoHF);
        BTree->SetBranchAddress("jet_pt_recotruthratio", &jet_pt_recotruthratio);
        BTree->SetBranchAddress("HF_pt_recotruthratio", &HF_pt_recotruthratio);

        int eventNum;
        int NumJets = 0;
        int NumJets_zdR = 0;
        int NumJets_ktdR = 0;
        int NumBJets = 0;

        TLorentzVector HFmeson, HFjet, mum, mup, K, Jpsi;
        TLorentzVector meas_HFmeson, meas_HFjet, meas_mum, meas_mup, meas_K, meas_Jpsi;
                
        cout << "Requested # of events" << NumEvts << endl;
        for (int ev = 0; ev < NumEvts; ev++) {
                BTree->GetEntry(ev);

                if (ev%10000 == 0) {
                        double percentage = 100.*ev/NumEvts;
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                HFjet.SetPxPyPzE(jet_px, jet_py, jet_pz, jet_e);
                mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
                mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);
                K.SetPxPyPzE(K_px, K_py, K_pz, K_e);
                HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
                Jpsi = mup + mum;

                meas_HFjet.SetPxPyPzE(meas_jet_px, meas_jet_py, meas_jet_pz, meas_jet_e);
                meas_mup.SetPxPyPzE(meas_mup_px, meas_mup_py, meas_mup_pz, meas_mup_e);
                meas_mum.SetPxPyPzE(meas_mum_px, meas_mum_py, meas_mum_pz, meas_mum_e);
                meas_K.SetPxPyPzE(meas_K_px, meas_K_py, meas_K_pz, meas_K_e);
                meas_HFmeson.SetPxPyPzE(meas_HF_px, meas_HF_py, meas_HF_pz, meas_HF_e);
                meas_Jpsi = meas_mup + meas_mum;
                
                if (jet_eta > etaMin && jet_eta < etaMax)
                        h2_HFptjetpt->Fill(HFmeson.Pt(), HFjet.Pt());

                bool pt_cond  = jet_pt > pTLow;
                bool rap_cond = jet_rap > etaMin && jet_rap < etaMax;
                
                if (!rap_cond)
                        continue;

                if (!pt_cond)
                        continue;

                if (Hasbbbar) {
                        h1_HFpt_GS->Fill(HFmeson.Pt());
                        h1_jetpt_GS->Fill(HFjet.Pt());
                } else {
                        h1_HFpt_FC->Fill(HFmeson.Pt());
                        h1_jetpt_FC->Fill(HFjet.Pt());
                }

                h1_HFpt_Total->Fill(HFmeson.Pt());
                h1_jetpt_Total->Fill(HFjet.Pt());

                if (SubtractGS && Hasbbbar) {
                        std::cout << " I am subtracting gluon splitting..." << std::endl;
                        continue;
                }

                event_counter++;
                // if(HFmeson.Pt() < 2) cout<<HFmeson.Pt()<<endl;

                // cout<<eventNumber<<": "<<jet_pt <<" -- "<< meas_jet_pt<<endl;
                // if(HFmeson.M() < 5.279) continue;
                // if(kts->size() < 1) cout<<"filled with 0 kts";
                NumBJets++;

                h1_jet_pt->Fill(jet_pt);
                h1_jet_eta->Fill(jet_eta);
                h1_jet_rap->Fill(jet_rap);

                h1_meas_jet_pt->Fill(meas_jet_pt);
                h1_meas_jet_eta->Fill(meas_jet_eta);
                h1_meas_jet_rap->Fill(meas_jet_rap);

                h1_HFpt->Fill(HFmeson.Pt());
                h1_HF_rap->Fill(HFmeson.Rapidity());
                h1_HF_pt->Fill(HFmeson.Pt());
                h1_HF_mass->Fill(HFmeson.M());

                h1_Jpsi_mass->Fill(Jpsi.M());
                h1_Jpsi_pt->Fill(Jpsi.Pt());
                h1_Jpsi_rap->Fill(Jpsi.Rapidity());

                h1_jet_ptbalance0->Fill(jet_pt / meas_jet_pt); //// Truth over Reco is the standard....
                h1_jet_ptbalance1->Fill(meas_jet_pt/ jet_pt );
                
                h1_HFjet_ptbalance->Fill(jet_pt / HFmeson.Pt());
                
                TVector3 HF_jet   = HFjet.Vect();
                TVector3 HF_meson = HFmeson.Vect();
                
                double frag_z = (HF_meson.Dot(HF_jet)) / (HF_jet.Mag2());  
                
                if (Hasbbbar) {
                        double frag_z_gluon = HFmeson.Vect().Dot(HFjet.Vect()) / (HFjet.Vect().Mag2());
                        double frag_jt_gluon = HFmeson.Vect().Cross(HFjet.Vect()).Mag() / HFjet.Vect().Mag();
                        double frag_r_gluon = HFmeson.DeltaR(HFjet, true);
                                
                        h1_frag_z_gluon->Fill(frag_z_gluon);
                        h1_frag_jt_gluon->Fill(frag_jt_gluon);
                        h1_frag_r_gluon->Fill(frag_r_gluon);
                        
                        h2_zjt_gluon->Fill(frag_z_gluon, frag_jt_gluon);
                        h2_zr_gluon->Fill(frag_z_gluon, frag_r_gluon);
                        h2_jtr_gluon->Fill(frag_jt_gluon, frag_r_gluon);
                                
                        h3_ptzjt_gluon->Fill(frag_z_gluon, frag_jt_gluon, jet_pt);
                        h3_ptzr_gluon->Fill(frag_z_gluon, frag_r_gluon, jet_pt);
                        h3_ptjtr_gluon->Fill(frag_jt_gluon, frag_r_gluon, jet_pt);
                }

                h1_z->Fill(truth_z);
                h1_jt->Fill(truth_jt);
                h1_r->Fill(truth_r);
                
                h1_meas_z->Fill(meas_z);
                h1_meas_jt->Fill(meas_jt);
                h1_meas_r->Fill(meas_r);
                
                h1_ratio_z0->Fill(truth_z / meas_z);
                h1_ratio_z1->Fill(meas_z/ truth_z);
                h1_ratio_jt0->Fill( truth_jt / meas_jt );
                h1_ratio_jt1->Fill(meas_jt/ truth_jt);
                h1_ratio_r0->Fill( truth_r / meas_r );
                h1_ratio_r1->Fill(meas_r/ truth_r);
                
                h2_frag_z_jetpt->Fill(frag_z, HFjet.Pt());
                
                if (Hasbbbar) {
                        //h2_frag_z_jetpt_gluon->Fill(HFmeson.Pt() / HFjet.Pt(), HFjet.Pt());
                        h2_jetpteta_gluon->Fill(HFjet.Pt(), HFjet.Eta());
                        h2_jetptp_gluon->Fill(HFjet.Pt(), HFjet.P());
                        h2_ptz_gluon->Fill(truth_z, jet_pt);
                        h2_ptjt_gluon->Fill(truth_jt, jet_pt);
                        h2_ptr_gluon->Fill(truth_r, jet_pt);                  
                }

                h2_jetpteta->Fill(HFjet.Pt(), HFjet.Eta());  
                h2_jetptp->Fill(HFjet.Pt(), HFjet.P());            

                h2_ptz->Fill(truth_z, jet_pt);
                h2_ptjt->Fill(truth_jt, jet_pt);
                h2_ptr->Fill(truth_r, jet_pt);
                
                h3_ptzjt->Fill(truth_z, truth_jt, jet_pt);
                h3_ptzr->Fill(truth_z, truth_r, jet_pt);
                h3_ptjtr->Fill(truth_jt, truth_r, jet_pt);    
                
                h2_zjt->Fill(truth_z, truth_jt);
                h2_zr->Fill(truth_z, truth_r);
                h2_jtr->Fill(truth_jt, truth_r);
                
                h2_meas_zjt->Fill(meas_z, meas_jt);
                h2_meas_zr->Fill(meas_z, meas_r);
                h2_meas_jtr->Fill(meas_jt, meas_r);

                bool hasEmission_ktdR = false;
                bool hasEmission_zdR = false;
        } // End of BTree entry loop
                
        h2_jetpteta_gluon_ratio->Divide(h2_jetpteta_gluon, h2_jetpteta);  
        h2_jetptp_gluon_ratio->Divide(h2_jetptp_gluon, h2_jetptp);
        h2_ptz_gluon_ratio->Divide(h2_ptz_gluon, h2_ptz);
        h2_ptjt_gluon_ratio->Divide(h2_ptjt_gluon, h2_ptjt);
        h2_ptr_gluon_ratio->Divide(h2_ptr_gluon, h2_ptr);   


        THStack *hs_ptz = new THStack("z_truth_all", ";z;#frac{1}{N_{jets}}#frac{dN}{dz}");
        THStack *hs_ptjt = new THStack("jt_truth_all", ";j_{T} [GeV/c];#frac{1}{N_{jets}}#frac{dN}{dj_{T}}");
        THStack *hs_ptr = new THStack("r_truth_all", ";r;#frac{1}{N_{jets}}#frac{dN}{dr}");         
        TH2D *h2_zjt_ptbinned[ptbinsize-1];
        TH2D *h2_zr_ptbinned[ptbinsize-1];
        TH2D *h2_jtr_ptbinned[ptbinsize-1];       
        
        for (int j = 1; j < ptbinsize; j++) {   
                TH1D *h1_ptz_temp = (TH1D *)h2_ptz->ProjectionX(Form("z_pt%d", j), j + 1, j + 1); 
                
                NormalizeHist(h1_ptz_temp);
                
                h1_ptz_temp->SetStats(0);
                h1_ptz_temp->SetMarkerStyle(j + 20);
                
                if (j!=5) {
                        h1_ptz_temp->SetMarkerColor(j);
                        h1_ptz_temp->SetLineColor(j);
                } else {
                        h1_ptz_temp->SetMarkerColor(j*j+3);
                        h1_ptz_temp->SetLineColor(j*j+3);
                }

                h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
                //h1_ptz_temp->SetOption("PE HIST");
                hs_ptz->Add(h1_ptz_temp);

                TH1D *h1_ptjt_temp = (TH1D *)h2_ptjt->ProjectionX(Form("jt_pt%d", j), j + 1, j + 1); 
                NormalizeHist(h1_ptjt_temp);        
                h1_ptjt_temp->SetStats(0);
                h1_ptjt_temp->SetMarkerStyle(j + 20);
                
                if (j!=5) {
                        h1_ptjt_temp->SetMarkerColor(j);
                        h1_ptjt_temp->SetLineColor(j);
                } else {
                        h1_ptjt_temp->SetMarkerColor(j*j+3);
                        h1_ptjt_temp->SetLineColor(j*j+3);
                }

                h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));        
                hs_ptjt->Add(h1_ptjt_temp);
                        
                TH1D *h1_ptr_temp = (TH1D *)h2_ptr->ProjectionX(Form("r_pt%d", j), j + 1, j + 1); 
                
                NormalizeHist(h1_ptr_temp);        
                
                h1_ptr_temp->SetStats(0);
                h1_ptr_temp->SetMarkerStyle(j + 20);
                
                if (j!=5) {
                        h1_ptr_temp->SetMarkerColor(j);
                        h1_ptr_temp->SetLineColor(j);
                } else {
                        h1_ptr_temp->SetMarkerColor(j*j+3);
                        h1_ptr_temp->SetLineColor(j*j+3);
                }

                h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));        
                hs_ptr->Add(h1_ptr_temp);   
                        
                h3_ptzjt->GetZaxis()->SetRange(j+1, j+1);      
                h2_zjt_ptbinned[j-1] = (TH2D *)h3_ptzjt->Project3D("yx");
                h2_zjt_ptbinned[j-1]->SetName(Form("zjt_truth_pt%d",j)); 
                NormalizeHist(h2_zjt_ptbinned[j-1]);	

                h3_ptzr->GetZaxis()->SetRange(j+1, j+1);        
                h2_zr_ptbinned[j-1] = (TH2D *)h3_ptzr->Project3D("yx");
                h2_zr_ptbinned[j-1]->SetName(Form("zr_truth_pt%d",j));
                NormalizeHist(h2_zr_ptbinned[j-1]);	
                        
                h3_ptjtr->GetZaxis()->SetRange(j+1, j+1);        
                h2_jtr_ptbinned[j-1] = (TH2D *)h3_ptjtr->Project3D("yx");
                h2_jtr_ptbinned[j-1]->SetName(Form("jtr_truth_pt%d",j));
                NormalizeHist(h2_jtr_ptbinned[j-1]);	      
        }  
        
        hs_ptz->Write();
        hs_ptjt->Write();
        hs_ptr->Write();
        
        h3_ptzjt->GetZaxis()->SetRange(1, ptbinsize+1);
        h3_ptzr->GetZaxis()->SetRange(1, ptbinsize+1);
        h3_ptjtr->GetZaxis()->SetRange(1, ptbinsize+1);    
                
        cout << event_counter << " events processed" << endl;

        TH1D *h1_GS_frac = (TH1D *)h1_HFpt_GS->Clone("h1_GS_frac");
        h1_GS_frac->Divide(h1_GS_frac, h1_HFpt_Total, 1, 1, "B");

        TH1D *h1_FC_frac = (TH1D *)h1_HFpt_FC->Clone("h1_FC_frac");
        h1_FC_frac->Divide(h1_FC_frac, h1_HFpt_Total, 1, 1, "B");

        SetRecoStyle(h1_GS_frac); //Blue
        SetTruthStyle(h1_FC_frac); //Green

        TH1D *h1_jetpt_GS_frac = (TH1D *)h1_jetpt_GS->Clone("h1_jetpt_GS_frac");
        h1_jetpt_GS_frac->Divide(h1_jetpt_GS_frac, h1_jetpt_Total, 1, 1, "B");

        TH1D *h1_jetpt_FC_frac = (TH1D *)h1_jetpt_FC->Clone("h1_jetpt_FC_frac");
        h1_jetpt_FC_frac->Divide(h1_jetpt_FC_frac, h1_jetpt_Total, 1, 1, "B");

        SetRecoStyle(h1_jetpt_GS_frac);
        SetTruthStyle(h1_jetpt_FC_frac);

        // h2_ktdR->Scale(1./ NumJets_ktdR);
        // h2_zdR->Scale(1./ NumJets_zdR);

        SetTruthStyle(h1_jet_eta);
        SetDataStyle(h1_jet_rap);
        SetTruthStyle(h1_jet_pt);

        SetRecoStyle(h1_meas_jet_eta);
        SetDataStyle(h1_meas_jet_rap);
        SetRecoStyle(h1_meas_jet_pt);

        SetTruthStyle(h1_HF_rap);
        SetTruthStyle(h1_HF_pt);
        SetTruthStyle(h1_HF_mass);

        SetTruthStyle(h1_z);
        SetTruthStyle(h1_jt);
        SetTruthStyle(h1_r);
        //  SetTruthStyle(h1_theta);

        SetRecoStyle(h1_meas_z);
        SetRecoStyle(h1_meas_jt);
        SetRecoStyle(h1_meas_r);
        //  SetRecoStyle(h1_meas_theta);

        SetTruthStyle(h1_Jpsi_rap);
        SetTruthStyle(h1_Jpsi_pt);
        SetTruthStyle(h1_Jpsi_mass);

        //---- paint setup...
        //
        int ican = -1, iframe = -1, itext = -1;
        TCanvas *ccan[1000];
        TH1F *frame[1000];
        TLatex *text[1000];
        for (int i = 0; i < 1000; i++) {
                text[i] = new TLatex();
                text[i]->SetNDC(kTRUE);
                text[i]->SetTextSize(0.06);
        }

        TLatex Tl;
        Tl.SetNDC(kTRUE);
        Tl.SetTextSize(0.04);
        //
        // gStyle->SetOptStat(0);
        // gStyle->SetPaperSize(TStyle::kUSLetter);
        // gStyle->SetPadBottomMargin(0.08);
        // gStyle->SetPadTopMargin(0.005);
        gStyle->SetPadLeftMargin(0.13);
        gStyle->SetPadRightMargin(0.13);
        gStyle->SetLabelSize(0.05, "X");
        gStyle->SetLabelSize(0.05, "Y");
        gStyle->SetTitleXSize(0.055);
        gStyle->SetTitleYSize(0.055);
        gStyle->SetTitleOffset(0.85, "X");
        gStyle->SetTitleOffset(1.2, "Y");
        gStyle->SetStatW(0.2);
        gStyle->SetPalette(kBird);
        gStyle->SetNumberContours(100);
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
        TString plotextension = TString("./plots/");
        rootfile = extension_RootFilesMC + extension + TString(".root");
        plotfile = plotextension + extension + TString(".ps");

        plotfilePDF = plotextension + extension + TString(".pdf");
        plotfileO = plotfilePDF + TString("(");
        plotfileC = plotfilePDF + TString("]");
        
        // c->SaveAs("plots/"+extension+".pdf");

        // c1->SaveAs(Form("plots_misc/Misc_%s"+extension+".pdf", file_name.c_str()));
        // cout<<"...outbase   = "<<outbase.Data()<<endl;
        // cout<<"...rootfile  = "<<rootfile.Data()<<endl;
        // cout<<"...plotfile  = "<<plotfile.Data()<<endl;
        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        
        h2_zjt->GetYaxis()->SetRangeUser(0, 12);
        h2_zjt->GetXaxis()->SetRangeUser(0, 1.05);
        h2_zjt->SetStats(0);
        h2_zjt->SetXTitle("z");
        h2_zjt->SetYTitle("j_{T}");
        h2_zjt->Draw("COLZ");

        // h2_lundplane->SetStats(0);
        // h2_lundplane->SetXTitle("ln(1/#theta)");
        // h2_lundplane->SetYTitle("ln(k_{T})");
        // h2_lundplane->Sethflocimum(.01);
        // h2_lundplane->SetMinimum(0.00001);
        //  gPad->SetLogz();
        Tl.DrawLatex(.58, 0.8, Form("%.1f<p_{T,jet}<%.1f GeV", ptMin, ptMax));
        // Tl.DrawLatex(0.7, 0.8, "No cuts");
        ccan[ican]->cd(2);
        // h2_thetaErad->Scale(1./NumJets);
        
        h2_zr->GetYaxis()->SetRangeUser(0, 1.0);
        h2_zr->GetXaxis()->SetRangeUser(0, 1.05);
        h2_zr->SetStats(0);
        h2_zr->SetXTitle("z");
        h2_zr->SetYTitle("r");
        // h2_thetaErad->Setmaximum(.01);
        // h2_thetaErad->SetMinimum(0.00001);
        //  gPad->SetLogz();
        h2_zr->Draw("COLZ");
        h2_zr->SetLineColor(1);
        h2_zr->SetLineStyle(9);

        Tl.DrawLatex(.5, 0.8, Form("%.1f<p_{T,jet}<%.1f GeV", ptMin, ptMax));
        Tl.DrawLatex(.6, 0.64, Form("%.1f<#eta<%.1f", etaMin, etaMax));
        ccan[ican]->cd(3);
        
        h2_jtr->SetStats(0);
        h2_jtr->SetXTitle("j_{T}");
        h2_jtr->SetYTitle("r");
        // h2_ktdR->SetMaximum(.01);
        // h2_ktdR->SetMinimum(0.00001);
        //  gPad->SetLogz();
        h2_jtr->Draw("COLZ");

        //  ccan[ican]->cd(4);
        //  h2_zdR->SetStats(0);
        //  h2_zdR->SetXTitle("ln(R/#DeltaR)");
        //  h2_zdR->SetYTitle("ln(1/z)");
        //  // h2_zdR->SetMaximum(.01);
        //  // h2_zdR->SetMinimum(0.00001);
        //  //  gPad->SetLogz();
        //  h2_zdR->Draw("COLZ");

        ccan[ican]->cd();
        ccan[ican]->Update();
        
        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        //

        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        h1_z->SetStats(0);
        h1_z->SetXTitle("z");

        h1_z->Draw("P E SAME");
        h1_z->SetMinimum(0.);
        // Tl.DrawLatex(.5, 0.75, Form("%.1f<p_{T,jet}<%.1f GeV", ptMin, ptMax));
        // Tl.DrawLatex(0.7, 0.8, "No cuts");
        h1_meas_z->SetStats(0);
        h1_meas_z->Draw("P E SAME");

        ccan[ican]->cd(2);
        // h2_thetaErad->Scale(1./NumJets);
        gPad->SetLogy();
        h1_jt->SetStats(0);
        h1_jt->SetXTitle("j_{T}");
        h1_jt->SetMinimum(0.1);
        h1_jt->SetMaximum(1e6);
        // h2_lundplane->Sethflocimum(.01);
        // h2_lundplane->SetMinimum(0.00001);
        h1_jt->Draw("P E SAME");
        h1_meas_jt->Draw("P E SAME");
        ccan[ican]->cd(3);
        
        gPad->SetLogy();
        h1_r->SetStats(0);
        h1_r->SetXTitle("r");
        h1_r->Draw("P E SAME");
        h1_r->SetMinimum(0.9);
        h1_r->SetMaximum(1e6);

        h1_meas_r->Draw("P E SAME");

        //  ccan[ican]->cd(4);
        //  h1_theta->SetStats(0);
        //  h1_theta->SetXTitle("ln(1/#theta)");
        //  h1_theta->SetMinimum(0.);
        //  // h2_lundplane->Sethflocimum(.01);
        //  // h2_lundplane->SetMinimum(0.00001);
        //  h1_theta->Draw("P E SAME");
        //  h1_meas_theta->Draw("P E SAME");

        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());

        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        h1_z->SetStats(0);
        h1_z->SetXTitle("z");

        THStack *stack_frag = new THStack("stack_frag", "");
        stack_frag->Add(h1_z);
        stack_frag->Add(h1_frag_z_gluon);

        SetTruthStyle(h1_z);
        SetDataStyle(h1_frag_z_gluon);
        stack_frag->Draw("P NOSTACK ");
        // h1_frag_z->Draw("P E SAME");

        // h1_frag_z_gluon->Draw("P E SAME");

        // SetTruthStyle(h1_tr_frag_z);
        // h1_tr_frag_z->Draw("P E SAME");

        ccan[ican]->cd(2);
        // h2_thetaErad->Scale(1./NumJets);
        gPad->SetLogy();
        h1_jt->SetStats(0);
        h1_jt->SetXTitle("jt");
        SetTruthStyle(h1_jt);
        h1_jt->Draw("P E SAME");

        // SetTruthStyle(h1_tr_frag_jt);
        // h1_tr_frag_jt->Draw("P E SAME");

        ccan[ican]->cd(3);
        gPad->SetLogy();
        SetTruthStyle(h1_r);
        h1_r->SetStats(0);
        h1_r->SetXTitle("r");  
        h1_r->Draw("P E SAME");

        // SetTruthStyle(h1_tr_frag_r);
        // h1_tr_frag_r->Draw("P E SAME");

        ccan[ican]->cd(4);
        NormalizeHist(h1_z);
        NormalizeHist(h1_frag_z_gluon);
        // stack_frag->Draw("P NOSTACK ");
        h1_z->Draw("PE SAME");
        h1_frag_z_gluon->Draw("PE SAME");

        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        
        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);

        for (int i = 1; i < ptbinsize; i++) {
                TH1D *h1_temp = (TH1D *)h2_frag_z_jetpt->ProjectionX(Form("htemp%d_fragzjetp", i), i + 1, i + 1);
                NormalizeHist(h1_temp);
                h1_temp ->SetStats(0);
                h1_temp->SetMarkerStyle(i + 20);
                h1_temp->SetMarkerColor(i + 1);
                h1_temp->SetLineColor(i + 1);
                h1_temp->Draw("P E SAME");
                h1_temp->Draw("HIST SAME");
                // h1_temp->SetMinimum(0.);
                // h1_temp->SetMaximum(1.05);
        }

        ccan[ican]->cd(2);
        /*
        auto legend_zg = new TLegend(0.55, 0.65, 0.9, 0.85);
        legend_zg ->SetTextSize(0.027);
        legend_zg ->SetBorderSize(0);
        legend_zg ->SetFillStyle(0);
        legend_zg ->SetFillColor(3);
        legend_zg->SetHeader("PYTHIA LHCb Unofficial","C");
        
        for (int i = 1; i < ptbinsize; i++)
        {
        TH1D *h1_temp = (TH1D *)h2_frag_z_jetpt_gluon->ProjectionX(Form("htemp%d_fragzjetpt_gluon", i), i + 1, i + 1);
        NormalizeHist(h1_temp);
        h1_temp ->SetStats(0);
        h1_temp->GetYaxis()->SetRangeUser(0, 2.1);
        h1_temp->SetMarkerStyle(i + 20);
        h1_temp->SetMarkerColor(i + 1);
        h1_temp->SetLineColor(i + 1);
        h1_temp->Draw("P E SAME");
        h1_temp->Draw("HIST SAME");
        h1_temp->SetMaximum(3.0);
        h1_temp->SetXTitle("z");
        h1_temp->SetYTitle(" 1/N dN/dz ");
        // h1_temp->SetMinimum(0.);
        // h1_temp->SetMaximum(1.05);
        legend_zg->AddEntry(h1_temp, Form(" %.1f < p_{T}^{jet} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));
        }
        */
        for (int i = 1; i < ptbinsize; i++) {
                TH1D *h1_temp = (TH1D *)h2_ptz->ProjectionX(Form("z_truth_pt%d", i), i + 1, i + 1);
                NormalizeHist(h1_temp);
                h1_temp ->SetStats(0);
                h1_temp->SetMarkerStyle(i + 20);
                h1_temp->SetMarkerColor(i + 1);
                h1_temp->SetLineColor(i + 1);
                h1_temp->Draw("P E SAME");
                h1_temp->Draw("HIST SAME");
                // h1_temp->SetMinimum(0.);
                // h1_temp->SetMaximum(1.05);
        }
        //  legend_zg ->Draw("SAME");

        // h1_frag_z_gluon->Draw("P E SAME");

        // SetTruthStyle(h1_tr_frag_z);
        // h1_tr_frag_z->Draw("P E SAME");
        /*
        ccan[ican]->cd(3);
        for (int i = 1; i < ptbinsize; i++)
        {
        TH1D *h1_temp = (TH1D *)h2_ptz->ProjectionX(Form("z_truth_pt%d", i), i + 1, i + 1);
        NormalizeHist(h1_temp);
        h1_temp ->SetStats(0);
        h1_temp->SetMarkerStyle(i + 20);
        h1_temp->SetMarkerColor(i + 1);
        h1_temp->SetLineColor(i + 1);
        h1_temp->Draw("P E SAME");
        h1_temp->Draw("HIST SAME");
        // h1_temp->SetMinimum(0.);
        // h1_temp->SetMaximum(1.05);
        }
        */
        ccan[ican]->cd();
        ccan[ican]->Update();
        
        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        //
        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        gPad->SetLogy();
        h1_jet_pt->SetStats(0);
        h1_jet_pt->SetXTitle("pT(jet)");
        h1_jet_pt->Draw("P E SAME");
        h1_meas_jet_pt->Draw("P E SAME");

        ccan[ican]->cd(2);
        // h2_thetaErad->Scale(1./NumJets);
        h1_jet_eta->SetStats(0);
        h1_jet_eta->SetXTitle("#eta(jet)");
        h1_jet_eta->Draw("P E SAME");
        h1_meas_jet_eta->Draw("P E SAME");

        ccan[ican]->cd(3);
        // h2_thetaErad->Scale(1./NumJets);
        h1_jet_rap->SetStats(0);
        h1_jet_rap->SetXTitle("#eta(jet)");
        h1_jet_rap->Draw("P E SAME");
        h1_meas_jet_rap->Draw("P E SAME");

        //  ccan[ican]->cd(4);
        //  // h1_nsplittings->Draw("P E SAME");
        //  h1_nJetDtrs_noghost->SetStats(0);
        //  h1_nJetDtrs_noghost_gluon->SetStats(0);
        //  SetTruthStyle(h1_nJetDtrs_noghost);
        //  SetDataStyle(h1_nJetDtrs_noghost_gluon);
        //  NormalizeHist(h1_nJetDtrs_noghost);
        //  NormalizeHist(h1_nJetDtrs_noghost_gluon);
        //  h1_nJetDtrs_noghost_gluon->Draw("P E SAME");
        //  h1_nJetDtrs_noghost->Draw("P E SAME");

        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        //

        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        gPad->SetLogy();
        h1_HF_pt->SetStats(0);
        h1_HF_pt->SetXTitle("pT(HF)");

        h1_HF_pt->Draw("P E SAME");
        // Tl.DrawLatex(.5, 0.75, Form("%.1f<p_{T,jet}<%.1f GeV", ptMin, ptMax));
        // Tl.DrawLatex(0.7, 0.8, "No cuts");
        // h1_tr_dR->SetStats(0);
        // h1_tr_dR->SetMarkerStyle(21);
        // h1_tr_dR->SetMarkerColor(kRed);
        // h1_tr_dR->SetLineColor(kRed);
        // h1_tr_dR->Draw("P E SAME");

        ccan[ican]->cd(2);
        // h2_thetaErad->Scale(1./NumJets);
        h1_HF_rap->SetStats(0);
        h1_HF_rap->SetXTitle("y(HF)");
        // h2_lundplane->Sethflocimum(.01);
        // h2_lundplane->SetMinimum(0.00001);
        h1_HF_rap->Draw("P E SAME");
        ccan[ican]->cd(3);
        h1_HF_mass->SetStats(0);
        // h1_dphi_HF_jet->SetXTitle("#Delta#Phi");
        // h2_lundplane->Sethflocimum(.01);
        // h2_lundplane->SetMinimum(0.00001);
        h1_HF_mass->Draw("P E SAME");

        // ccan[ican]->cd(4);
        // h1_FC_frac->Draw("PE1SAME");
        // h1_GS_frac->Draw("PE1SAME");

        // h1_FC_frac->SetMinimum(0.);
        // h1_FC_frac->SetMaximum(1.05);

        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());

        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        //  gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        //  gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        h1_FC_frac->SetStats(0);
        h1_GS_frac->SetStats(0);
        h1_FC_frac->Draw("PE1SAME");
        h1_GS_frac->Draw("PE1SAME");

        h1_FC_frac->SetMinimum(0.);
        h1_FC_frac->SetMaximum(1.05);

        ccan[ican]->cd(2);
        h1_jetpt_FC_frac->SetStats(0);
        h1_jetpt_GS_frac->SetStats(0);
        //  h1_jetpt_FC_frac->Draw("PE1SAME");
        h1_jetpt_GS_frac->Draw("PE1SAME");

        //  h1_jetpt_FC_frac->SetMinimum(0.);
        h1_jetpt_GS_frac->SetMaximum(1.05);
        h1_jetpt_GS_frac->SetXTitle("p_{T}^{Jet}");
        h1_jetpt_GS_frac->SetYTitle("Fraction of gluon splitting");
        h1_jetpt_GS_frac->GetXaxis()->SetLabelSize(0.03);
        h1_jetpt_GS_frac->GetYaxis()->SetLabelSize(0.03);
        
        auto legend_gs = new TLegend(0.45, 0.6, 0.75, 0.9);
        legend_gs ->SetTextSize(0.05);
        legend_gs ->SetBorderSize(0);
        legend_gs ->SetFillStyle(0);
        legend_gs ->SetFillColor(3);
        legend_gs->Draw("SAME");

        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        
        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        gPad->SetLogy();
        h1_Jpsi_pt->SetStats(0);
        h1_Jpsi_pt->SetXTitle("pT(J/#psi)");

        h1_Jpsi_pt->Draw("P E SAME");
        // Tl.DrawLatex(.5, 0.75, Form("%.1f<p_{T,jet}<%.1f GeV", ptMin, ptMax));
        // Tl.DrawLatex(0.7, 0.8, "No cuts");
        // h1_tr_dR->SetStats(0);
        // h1_tr_dR->SetMarkerStyle(21);
        // h1_tr_dR->SetMarkerColor(kRed);
        // h1_tr_dR->SetLineColor(kRed);
        // h1_tr_dR->Draw("P E SAME");

        ccan[ican]->cd(2);
        // h2_thetaErad->Scale(1./NumJets);
        h1_Jpsi_rap->SetStats(0);
        h1_Jpsi_rap->SetXTitle("y(J/#psi)");
        // h2_lundplane->Sethflocimum(.01);
        // h2_lundplane->SetMinimum(0.00001);
        h1_Jpsi_rap->Draw("P E SAME");
        ccan[ican]->cd(3);
        h1_Jpsi_mass->SetStats(0);
        // h1_dphi_Jpsi_jet->SetXTitle("#Delta#Phi");
        // h2_lundplane->Sethflocimum(.01);
        // h2_lundplane->SetMinimum(0.00001);
        h1_Jpsi_mass->Draw("P E SAME");

        ccan[ican]->cd();
        ccan[ican]->Update();
        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        
        ccan[ican]->cd(1);
        /// "crystalball" parameters: [alpha], [N], [Sigma], [mean]
        TF1 *fitFunc = new TF1("fitFunc", "crystalball", 0.0, 2.0);
        fitFunc->SetParameters(1.0, 1.0, 0.5, 1.5, 2.0);
        
        SetTruthStyle(h1_ratio_z0);
        SetRecoStyle(h1_ratio_z1);
        
        h1_ratio_z0->Fit(fitFunc, "RQ");
        
        h1_ratio_z0->SetStats(0);
        h1_ratio_z0->SetXTitle("z balance");
        h1_ratio_z0->Draw("P E SAME");
        fitFunc->Draw("P E SAME");
        h1_ratio_z1->SetStats(0);
        h1_ratio_z1->SetXTitle("z balance");
        h1_ratio_z1->Draw("P E SAME");
        
        auto legend_zbalance = new TLegend(0.2, 0.6, 0.4, 0.9);
        legend_zbalance ->SetTextSize(0.04);
        legend_zbalance ->SetBorderSize(0);
        legend_zbalance ->SetFillStyle(0);
        legend_zbalance ->SetFillColor(3);
        legend_zbalance->AddEntry(h1_ratio_z0, "z_{True}/z_{Reco}");
        legend_zbalance->AddEntry(h1_ratio_z1, "z_{Reco}/z_{Truth}");
        legend_zbalance->Draw("SAME");

        ccan[ican]->cd(2);
        SetTruthStyle(h1_ratio_jt0);
        SetRecoStyle(h1_ratio_jt1);
        h1_ratio_jt0->SetStats(0);
        h1_ratio_jt0->SetXTitle("jt balance");
        h1_ratio_jt0->Draw("P E SAME");
        h1_ratio_jt1->SetStats(0);
        h1_ratio_jt1->SetXTitle(" jt balance");
        h1_ratio_jt1->Draw("P E SAME");
        
        auto legend_jtbalance = new TLegend(0.6, 0.6, 0.9, 0.9);
        legend_jtbalance ->SetTextSize(0.04);
        legend_jtbalance ->SetBorderSize(0);
        legend_jtbalance ->SetFillStyle(0);
        legend_jtbalance ->SetFillColor(3);
        legend_jtbalance->AddEntry(h1_ratio_jt0, "j_{T; Truth}/j_{T; Reco}");
        legend_jtbalance->AddEntry(h1_ratio_jt1, "j_{T; Reco}/j_{T; Truth}");
        legend_jtbalance->Draw("SAME");


        ccan[ican]->cd(3);
        SetTruthStyle(h1_ratio_r0);
        SetRecoStyle(h1_ratio_r1);
        h1_ratio_r0->SetStats(0);
        h1_ratio_r0->SetXTitle(" r balance");
        h1_ratio_r0->Draw("P E SAME");
        h1_ratio_r1->SetStats(0);
        h1_ratio_r1->SetXTitle(" r balance");
        h1_ratio_r1->Draw("P E SAME");
        
        auto legend_rbalance = new TLegend(0.7, 0.6, 0.9, 0.9);
        legend_rbalance ->SetTextSize(0.04);
        legend_rbalance ->SetBorderSize(0);
        legend_rbalance ->SetFillStyle(0);
        legend_rbalance ->SetFillColor(3);
        legend_rbalance->AddEntry(h1_ratio_r0, "r_{Truth}/r_{Reco}");
        legend_rbalance->AddEntry(h1_ratio_r1, "r_{Reco}/r_{Truth}");
        legend_rbalance->Draw("SAME");
        
        ccan[ican]->cd(4);
        SetTruthStyle(h1_jet_ptbalance0);
        SetRecoStyle(h1_jet_ptbalance1);
        h1_jet_ptbalance0->SetStats(0);
        h1_jet_ptbalance0->SetXTitle(" pt balance");
        h1_jet_ptbalance0->Draw("P E SAME");
        h1_jet_ptbalance1->SetStats(0);
        h1_jet_ptbalance1->SetXTitle(" pt balance");
        h1_jet_ptbalance1->Draw("P E SAME");
        
        auto legend_ptbalance = new TLegend(0.6, 0.6, 0.9, 0.9);
        legend_ptbalance ->SetTextSize(0.04);
        legend_ptbalance ->SetBorderSize(0);
        legend_ptbalance ->SetFillStyle(0);
        legend_ptbalance ->SetFillColor(3);
        legend_ptbalance->AddEntry(h1_jet_ptbalance0, "p_{T; Truth}/p_{T; Reco}");
        legend_ptbalance->AddEntry(h1_jet_ptbalance1, "p_{T; Reco}/p_{T; Truth}");
        legend_ptbalance->Draw("SAME");


        ccan[ican]->cd();
        ccan[ican]->Update();

        if (ican == 0)
                ccan[ican]->Print(plotfileO.Data());
        else
                ccan[ican]->Print(plotfilePDF.Data());
        
        //
        ++ican;
        sprintf(buf, "ccan%d", ican);
        ccan[ican] = new TCanvas(buf, buf, 30 * ican, 30 * ican, 800, (8.5 / 11.) * 800);
        ccan[ican]->SetFillColor(10);
        gPad->SetLeftMargin(0.16);
        gPad->SetBottomMargin(0.06);
        gPad->SetRightMargin(0.15);
        ccan[ican]->cd();
        ccan[ican]->Divide(2, 2, 0.0001, 0.0001);
        ccan[ican]->cd(1);
        gPad->SetLogy();
        h1_mup_pt->SetStats(0);
        h1_mup_pt->SetXTitle("p_{T} [GeV]");

        SetRecoStyle(h1_mup_pt);
        SetTruthStyle(h1_mum_pt);
        SetDataStyle(h1_K_pt);
        h1_mup_pt->Draw("P E SAME");
        h1_mum_pt->Draw("P E SAME");
        h1_K_pt->Draw("P E SAME");

        ccan[ican]->cd(2);
        SetRecoStyle(h1_meas_mup_pt);
        SetTruthStyle(h1_meas_mum_pt);
        SetDataStyle(h1_meas_K_pt);
        h1_meas_mup_pt->Draw("P E SAME");
        h1_meas_mum_pt->Draw("P E SAME");
        h1_meas_K_pt->Draw("P E SAME");

        ccan[ican]->cd(3);
        // gPad->SetLogy();
        h1_mup_eta->SetStats(0);
        h1_mup_eta->SetXTitle("p_{T} [GeV]");

        SetRecoStyle(h1_mup_eta);
        SetTruthStyle(h1_mum_eta);
        SetDataStyle(h1_K_eta);
        h1_mup_eta->Draw("P E SAME");
        h1_mum_eta->Draw("P E SAME");
        h1_K_eta->Draw("P E SAME");

        ccan[ican]->cd(4);
        SetRecoStyle(h1_meas_mup_eta);
        SetTruthStyle(h1_meas_mum_eta);
        SetDataStyle(h1_meas_K_eta);
        h1_meas_mup_eta->Draw("P E SAME");
        h1_meas_mum_eta->Draw("P E SAME");
        h1_meas_K_eta->Draw("P E SAME");

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

        //
        f.Write();
        f.Close();

        cout << "Num of True B jets = " << NumBJets << endl;
}
