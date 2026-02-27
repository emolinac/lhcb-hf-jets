#include <TCanvas.h>
#include <vector>
#include <iostream>

#include "fastjet/ClusterSequence.hh"
// #include "fastjet/contrib/SoftDrop.hh"

#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/TBJetsMCReco.h"
#include "../include/TBJetsMCReco.C"

using namespace fastjet;
using namespace std;

void MakeVarTreeMCReco(int NumEvts_user = -1,
                       int year = 6, // 2016
                       int Mag = 1, // 1 up, -1 down
                       bool DoJESJER = false,
                       bool DoJetID = false,
                       bool L0MuonDiMuon = false)
{
        TBenchmark* benchmark = new TBenchmark();
        benchmark->Start("MakeVarTreeMCReco");

        int NumEvts      = NumEvts_user;
        int NumEvtsTruth = NumEvts_user;
        
        int HF_pdgcode = 521;

        TString str_followHard = "_HF";
        TString str_flavor     = "_b";
        TString str_level      = "reco";
        
        TString str_year = "2016";
        if (year == 6)
                str_year = "2016";
        else if (year == 7)
                str_year = "2017";
        else if (year == 8)
                str_year = "2018";

        TString str_Mag = "";
        if (Mag == -1)
                str_Mag = "_MD";
        else if (Mag == 1)
                str_Mag = "_MU";

        TString str_L0 = "";
        if (L0MuonDiMuon)
                str_L0 = "_L0MuonDiMuon";

        TString extension = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + 
                            str_followHard + str_Mag + str_flavor + str_L0;
        
        if (DoJESJER)
                extension = TString("JESJER_") + extension;
        if (DoJetID)
                extension = TString("jetid_") + extension;

        ////////////////////////////////////////////////////
        ///              Tracking and PID
        //////////////////////////////////////////////////
        TString eff_path = TString( "./efficiencies/");
        
        TString extension_pideffK       = "effhists-Turbo" + str_year + str_Mag + "-K-Brunel_DLLK>0-P.ETA.nTracks_Brunel";
        TString extension_pideffMu      = "effhists-Turbo" + str_year + str_Mag + "-Mu-IsMuon&Brunel_DLLmu>0-P.ETA.nTracks_Brunel";
        TString extension_trackeff      = "trackEff_" + str_year + "_Ratio_Full_Long_method";
        TString extension_trackeff_Muon = "trackEff_Muon_" + str_year + "_Ratio_Full_Long_method";
        TString extension_trackeff_Data = "trackEff_" + str_year + "_Data_Full_Long_method";
        TString extension_trackeff_MC   = "trackEff_" + str_year + "_MC_Full_Long_method";
        
        TFile file_trackeff(eff_path+ "TrackEff/" + extension_trackeff + ".root", "READ");
        TFile file_trackeff_Muon(eff_path+ "TrackEff/" + extension_trackeff_Muon + ".root", "READ");
        TFile file_trackeff_Data(eff_path+ "TrackEff/" + extension_trackeff_Data + ".root", "READ");
        TFile file_trackeff_MC(eff_path+ "TrackEff/" + extension_trackeff_MC + ".root", "READ");
        TFile file_pideffK(eff_path + "PIDEff/" + extension_pideffK + ".root", "READ");
        TFile file_pideffMu(eff_path + "PIDEff/" + extension_pideffMu + ".root", "READ");

        std::cout << extension_pideffK << std::endl;
        std::cout << extension_pideffMu << std::endl;
        std::cout << extension_trackeff << std::endl;

        TH2D* h2_ratio_trkeff_P_ETA_Muon       = (TH2D *)file_trackeff_Muon.Get("hP_ETA");
        TH2D* h2_ratio_trkeff_P_ETA_Muon_ERRHI = (TH2D *)file_trackeff_Muon.Get("hP_ETA_errhi");
        TH2D* h2_ratio_trkeff_P_ETA_Muon_ERRLO = (TH2D *)file_trackeff_Muon.Get("hP_ETA_errlo");

        TH2D* h2_ratio_trkeff_P_ETA       = (TH2D *)file_trackeff.Get("hP_ETA");
        TH2D* h2_ratio_trkeff_P_ETA_DATA  = (TH2D *)file_trackeff_Data.Get("hP_ETA");
        TH2D* h2_ratio_trkeff_P_ETA_MC    = (TH2D *)file_trackeff_MC.Get("hP_ETA");
        TH2D* h2_ratio_trkeff_P_ETA_ERRHI = (TH2D *)file_trackeff.Get("hP_ETA_errhi");
        TH2D* h2_ratio_trkeff_P_ETA_ERRLO = (TH2D *)file_trackeff.Get("hP_ETA_errlo");
        TH3D* h3_pideff_K_P_ETA_nTracks   = (TH3D *)file_pideffK.Get("eff_Brunel_DLLK>0");
        TH3D* h3_pideff_Mu_P_ETA_nTracks  = (TH3D *)file_pideffMu.Get("eff_IsMuon&Brunel_DLLmu>0");

        //////////////////////////////////////////////////
        ///              Trigger
        ////////////////////////////////////////////////
        TString extension_RootFilesTrig = TString("./efficiencies/TrigEff/");
        
        TString extension_trig_MC   = "PhotonHadronElectronTIS_jpsieff_reco_ev_-1_b_PID_91599.root";
        TString extension_trig_Data = "PhotonHadronElectronTIS_jpsieff_data_ev_-1_b_PID_91599.root";

        TFile file_trigeffMC(extension_RootFilesTrig + extension_trig_MC, "READ");
        TFile file_trigeffData(extension_RootFilesTrig + extension_trig_Data, "READ");

        TH2D* h2_trigeff_Data  = (TH2D *)file_trigeffData.Get("efficiency_Jpsiptrap");
        TH2D* h2_trigeff_MC    = (TH2D *)file_trigeffMC.Get("efficiency_Jpsiptrap");
        TH2D *h2_trigeff_ratio = (TH2D *)h2_trigeff_Data->Clone("h2_trigeff_ratio");
        
        h2_trigeff_ratio->Divide(h2_trigeff_MC);

        // WTA related stuff
        JetDefinition jet_def(cambridge_algorithm, JetDefinition::max_allowable_R);
        JetDefinition WTA(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme);

        PseudoJet dtr_pseudojet1, dtr_pseudojet2;

        vector<PseudoJet> jetdtrs, true_jetdtrs;

        double WTA_true_dist;
        double WTA_reco_dist;

        TLorentzVector WTA_true_axis;
        TLorentzVector WTA_reco_axis;
                
        TBJetsMCReco Tree;

        cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;
        
        if (NumEvts == -1)
                NumEvts = Tree.fChain->GetEntries();
        
        cout << "Executing CAJetAlgo" << endl;

        TFile f((output_folder + "ntuple_bjets_mcreco.root").c_str(), "RECREATE");
        
        TH1F *h1_TIS    = new TH1F("h1_TIS"   , "", ptJpsibinsize, ptJpsi_binedges);
        TH1F *h1_TISTOS = new TH1F("h1_TISTOS", "", ptJpsibinsize, ptJpsi_binedges);

        TH2F *h2_TIS_ptrap    = new TH2F("h2_TIS_ptrap"   , "", ptJpsibinsize, ptJpsi_binedges, HFetabinsize, HFeta_binedges);
        TH2F *h2_TISTOS_ptrap = new TH2F("h2_TISTOS_ptrap", "", ptJpsibinsize, ptJpsi_binedges, HFetabinsize, HFeta_binedges);

        TH2D *h2_HFpt_RM = new TH2D("h2_HFpt_RM", "", 30, 0, 100, 30, 0, 100);

        float jet_pt, jet_eta, tr_jet_pt, tr_jet_eta;
        float jet_rap, tr_jet_rap;
        float jet_px, jet_py, jet_pz, jet_e, jet_charge;
        float HF_px, HF_py, HF_pz, HF_e, HF_pt;
        float Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e;
        float mup_px, mup_py, mup_pz, mup_e;
        float mum_px, mum_py, mum_pz, mum_e;
        float K_px, K_py, K_pz, K_e, K_p, K_eta;
        float mup_CHI2NDOF, mup_GHOSTPROB, mup_IPCHI2;
        float mum_CHI2NDOF, mum_GHOSTPROB, mum_IPCHI2;
        float K_CHI2NDOF, K_GHOSTPROB, K_IPCHI2;
        float Jpsi_CHI2NDOF, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_BPVDLS, Jpsi_IPCHI2;
        float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
        int nTracks, nSPDHits;
        float jpsi_ipchi2, k_ipchi2;
        float K_PIDK;

        float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
        float tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e, tr_HF_pt;
        float tr_mup_px, tr_mup_py, tr_mup_pz, tr_mup_e;
        float tr_mum_px, tr_mum_py, tr_mum_pz, tr_mum_e;
        float tr_K_px, tr_K_py, tr_K_pz, tr_K_e;
        
        double z, jt, r;
        double tr_z, tr_jt, tr_r;
        double zg, jtg, rg;

        int nSV;
        bool isTrueBjet, Hasbbbar;
        float bmass_dtf, chi2ndf_dtf, tau_dtf;
        int NumBHads_tr, eventNumber;
        
        float pideff_K(1.0), pideff_mup(1.0), pideff_mum(1.0);
        float pideff_K_err(1.0), pideff_mup_err(1.0), pideff_mum_err(1.0);
        float trkeff_ratio_K(1.0), trkeff_ratio_mup(1.0), trkeff_ratio_mum(1.0);
        float trkeff_ratio_K_errhi(1.0), trkeff_ratio_mup_errhi(1.0), trkeff_ratio_mum_errhi(1.0);
        float trkeff_ratio_K_errlo(1.0), trkeff_ratio_mup_errlo(1.0), trkeff_ratio_mum_errlo(1.0);
        float trigeff_Data(1.0), trigeff_MC(1.0), trigeff_ratio(1.0);
        
        vector<float> dtr_pt, dtr_rap, dtr_id, dtr_3charge;
        
        float sv_mass, sv_chi2, sv_cosine, sv_ntrks;
        int SVTag;

        bool mup_L0, mum_L0;
        bool jpsi_L0, jpsi_L0Muon, jpsi_L0DiMuon, jpsi_Hlt1, jpsi_Hlt2, jpsi_Hlt2_Detached;
        bool Trig, TIS, TOS;
        
        TTree *BTree = new TTree("BTree", "B-jets Tree Variables");

        BTree->Branch("eventNumber", &eventNumber);

        BTree->Branch("dtr_pt", &dtr_pt);
        BTree->Branch("dtr_rap", &dtr_rap);
        BTree->Branch("dtr_id", &dtr_id);
        BTree->Branch("dtr_3charge", &dtr_3charge);
        
        BTree->Branch("jet_pt", &jet_pt);
        BTree->Branch("jet_eta", &jet_eta);
        BTree->Branch("jet_rap", &jet_rap);

        BTree->Branch("jet_px", &jet_px);
        BTree->Branch("jet_py", &jet_py);
        BTree->Branch("jet_pz", &jet_pz);
        BTree->Branch("jet_e", &jet_e);

        BTree->Branch("HF_px", &HF_px);
        BTree->Branch("HF_py", &HF_py);
        BTree->Branch("HF_pz", &HF_pz);
        BTree->Branch("HF_e", &HF_e);
        BTree->Branch("HF_pt", &HF_pt);

        BTree->Branch("Bu_IPCHI2", &Bu_IPCHI2);
        BTree->Branch("Bu_CHI2", &Bu_CHI2);
        BTree->Branch("Bu_CHI2NDOF", &Bu_CHI2NDOF);

        BTree->Branch("mum_px", &mum_px);
        BTree->Branch("mum_py", &mum_py);
        BTree->Branch("mum_pz", &mum_pz);
        BTree->Branch("mum_e", &mum_e);
        BTree->Branch("mum_IPCHI2", &mum_IPCHI2);
        BTree->Branch("mum_CHI2NDOF", &mum_CHI2NDOF);
        BTree->Branch("mum_GHOSTPROB", &mum_GHOSTPROB);

        BTree->Branch("mup_px", &mup_px);
        BTree->Branch("mup_py", &mup_py);
        BTree->Branch("mup_pz", &mup_pz);
        BTree->Branch("mup_e", &mup_e);
        BTree->Branch("mup_IPCHI2", &mup_IPCHI2);
        BTree->Branch("mup_CHI2NDOF", &mup_CHI2NDOF);
        BTree->Branch("mup_GHOSTPROB", &mup_GHOSTPROB);

        BTree->Branch("K_px", &K_px);
        BTree->Branch("K_py", &K_py);
        BTree->Branch("K_pz", &K_pz);
        BTree->Branch("K_e", &K_e);
        BTree->Branch("K_p", &K_p);
        BTree->Branch("K_eta", &K_eta);
        BTree->Branch("K_IPCHI2", &K_IPCHI2);
        BTree->Branch("K_CHI2NDOF", &K_CHI2NDOF);
        BTree->Branch("K_GHOSTPROB", &K_GHOSTPROB);
        BTree->Branch("k_ipchi2", &k_ipchi2);

        BTree->Branch("nTracks", &nTracks);
        BTree->Branch("nSPDHits", &nSPDHits);

        BTree->Branch("Jpsi_px", &Jpsi_px);
        BTree->Branch("Jpsi_py", &Jpsi_py);
        BTree->Branch("Jpsi_pz", &Jpsi_pz);
        BTree->Branch("Jpsi_e", &Jpsi_e);
        BTree->Branch("Jpsi_FDCHI2", &Jpsi_FDCHI2);
        BTree->Branch("Jpsi_BPVDLS", &Jpsi_BPVDLS);    
        BTree->Branch("Jpsi_CHI2", &Jpsi_CHI2);
        BTree->Branch("Jpsi_CHI2NDOF", &Jpsi_CHI2NDOF);
        BTree->Branch("jpsi_ipchi2", &jpsi_ipchi2);

        BTree->Branch("tr_jet_pt", &tr_jet_pt);
        BTree->Branch("tr_jet_eta", &tr_jet_eta);
        BTree->Branch("tr_jet_rap", &tr_jet_rap);

        BTree->Branch("tr_jet_px", &tr_jet_px);
        BTree->Branch("tr_jet_py", &tr_jet_py);
        BTree->Branch("tr_jet_pz", &tr_jet_pz);
        BTree->Branch("tr_jet_e", &tr_jet_e);

        BTree->Branch("tr_HF_px", &tr_HF_px);
        BTree->Branch("tr_HF_py", &tr_HF_py);
        BTree->Branch("tr_HF_pz", &tr_HF_pz);
        BTree->Branch("tr_HF_e", &tr_HF_e);
        BTree->Branch("tr_HF_pt", &tr_HF_pt);

        BTree->Branch("tr_mum_px", &tr_mum_px);
        BTree->Branch("tr_mum_py", &tr_mum_py);
        BTree->Branch("tr_mum_pz", &tr_mum_pz);
        BTree->Branch("tr_mum_e", &tr_mum_e);

        BTree->Branch("tr_mup_px", &tr_mup_px);
        BTree->Branch("tr_mup_py", &tr_mup_py);
        BTree->Branch("tr_mup_pz", &tr_mup_pz);
        BTree->Branch("tr_mup_e", &tr_mup_e);

        BTree->Branch("tr_K_px", &tr_K_px);
        BTree->Branch("tr_K_py", &tr_K_py);
        BTree->Branch("tr_K_pz", &tr_K_pz);
        BTree->Branch("tr_K_e", &tr_K_e);

        BTree->Branch("isTrueBjet", &isTrueBjet);

        BTree->Branch("nSV", &nSV);
        BTree->Branch("sv_mass", &sv_mass);
        BTree->Branch("sv_chi2", &sv_chi2);
        BTree->Branch("sv_ntrks", &sv_ntrks);
        BTree->Branch("sv_cosine", &sv_cosine);
        
        BTree->Branch("SVTag", &SVTag);

        BTree->Branch("bmass_dtf", &bmass_dtf);
        BTree->Branch("chi2ndf_dtf", &chi2ndf_dtf);
        BTree->Branch("tau_dtf", &tau_dtf);
        BTree->Branch("NumBHads_tr", &NumBHads_tr);
        BTree->Branch("Hasbbbar", &Hasbbbar);
        BTree->Branch("K_PIDK", &K_PIDK);
        BTree->Branch("pideff_K", &pideff_K);
        BTree->Branch("pideff_mum", &pideff_mum);
        BTree->Branch("pideff_mup", &pideff_mup);
        BTree->Branch("pideff_K_err", &pideff_K_err);
        BTree->Branch("pideff_mum_err", &pideff_mum_err);
        BTree->Branch("pideff_mup_err", &pideff_mup_err);

        BTree->Branch("trkeff_ratio_K", &trkeff_ratio_K);
        BTree->Branch("trkeff_ratio_mup", &trkeff_ratio_mup);
        BTree->Branch("trkeff_ratio_mum", &trkeff_ratio_mum);

        BTree->Branch("trkeff_ratio_K_errhi", &trkeff_ratio_K_errhi);
        BTree->Branch("trkeff_ratio_mup_errhi", &trkeff_ratio_mup_errhi);
        BTree->Branch("trkeff_ratio_mum_errhi", &trkeff_ratio_mum_errhi);
        BTree->Branch("trkeff_ratio_K_errlo", &trkeff_ratio_K_errlo);
        BTree->Branch("trkeff_ratio_mup_errlo", &trkeff_ratio_mup_errlo);
        BTree->Branch("trkeff_ratio_mum_errlo", &trkeff_ratio_mum_errlo);

        BTree->Branch("trigeff_Data", &trigeff_Data);
        BTree->Branch("trigeff_MC", &trigeff_MC);
        BTree->Branch("trigeff_ratio", &trigeff_ratio);

        BTree->Branch("mup_L0", &mup_L0);
        BTree->Branch("mum_L0", &mum_L0);
        
        BTree->Branch("jpsi_L0", &jpsi_L0);
        BTree->Branch("jpsi_L0Muon", &jpsi_L0Muon);
        BTree->Branch("jpsi_L0DiMuon", &jpsi_L0DiMuon);
        
        BTree->Branch("jpsi_Hlt1", &jpsi_Hlt1);
        BTree->Branch("jpsi_Hlt2", &jpsi_Hlt2);
        BTree->Branch("jpsi_Hlt2_Detached", &jpsi_Hlt2_Detached);
        
        BTree->Branch("Trig", &Trig);
        BTree->Branch("TIS", &TIS);
        BTree->Branch("TOS", &TOS);

        BTree->Branch("z", &z);
        BTree->Branch("jt", &jt);
        BTree->Branch("r", &r);
        BTree->Branch("zg", &zg);
        BTree->Branch("jtg", &jtg);
        BTree->Branch("rg", &rg);
        BTree->Branch("tr_z", &tr_z);
        BTree->Branch("tr_jt", &tr_jt);
        BTree->Branch("tr_r", &tr_r);
        
        BTree->Branch("WTA_true_dist", &WTA_true_dist);
        BTree->Branch("WTA_reco_dist", &WTA_reco_dist);

        // Event loop
        unsigned long long last_eventNum = 0;
        int events = 0;

        int ev_min = 0;
        TRandom3 *myRNG = new TRandom3();

        TLorentzVector HFjet, recojet, tr_truthjet, HFmeson, mup, mum, Jpsi, Kmeson;
        TLorentzVector tr_HFjet, tr_Kmeson, tr_mum, tr_mup, tr_HFmeson, tr_HFmeson_injet;
        TLorentzVector dtr, tr_dtr;
        TLorentzVector HFmatch;
        TLorentzVector dtr_matchtruthjet;
                
        ClusterSequence WTA_true_jets;
        ClusterSequence WTA_reco_jets;
        PseudoJet WTA_true_jet;
        PseudoJet WTA_reco_jet;
                
        for (int ev = ev_min; ev < NumEvts + ev_min; ev++) {
                dtr_pt.clear();
                dtr_id.clear();
                dtr_rap.clear();
                dtr_3charge.clear();
                jetdtrs.clear();
                true_jetdtrs.clear();

                Tree.GetEntry(ev);

                if (ev%10000 == 0) {
                        double percentage = 100.*ev/NumEvts;
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (ev != 0)
                if (Tree.eventNumber == last_eventNum)
                        continue;
                
                if (Tree.nPVs > 1)
                continue;
                
                bool WrongB = false;
                
                for (int dtrs0 = 0; dtrs0 < Tree.Jet_Dtr_nrecodtr; dtrs0++) {
                        if (abs(Tree.Jet_Dtr_ID[dtrs0]) == HF_pdgcode) {
                                if (fabs(Tree.Jet_Dtr_E[dtrs0] - Tree.Bu_PE) / (Tree.Bu_PE) > 0.001)
                                        WrongB = true;
                        
                                break;
                        }
                }
                
                if (WrongB)
                        continue;

                // Ibrahim trigger lines
                mup_L0 = -999; //Tree.mup_L0MuonDecision_TOS || Tree.mup_L0DiMuonDecision_TOS;
                mum_L0 = -999; //Tree.mum_L0MuonDecision_TOS || Tree.mum_L0DiMuonDecision_TOS;
        
                jpsi_L0 = Tree.Jpsi_L0MuonDecision_TOS || Tree.Jpsi_L0DiMuonDecision_TOS;
                jpsi_L0Muon = Tree.Jpsi_L0MuonDecision_TOS;
                jpsi_L0DiMuon = Tree.Jpsi_L0DiMuonDecision_TOS;

                jpsi_Hlt1 = Tree.Jpsi_Hlt1DiMuonHighMassDecision_TOS;
                jpsi_Hlt2 = Tree.Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS;
                jpsi_Hlt2_Detached = Tree.Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS;
        
                TIS = (Tree.Jpsi_L0Global_TIS && Tree.Jpsi_Hlt1Global_TIS && Tree.Jpsi_Hlt2Global_TIS);
                TOS = jpsi_L0 && jpsi_Hlt1 && jpsi_Hlt2;

                HFjet.SetPxPyPzE(Tree.Jet_PX / 1000.,
                                 Tree.Jet_PY / 1000.,
                                 Tree.Jet_PZ / 1000.,
                                 Tree.Jet_PE / 1000.);
                                
                mup.SetPxPyPzE(Tree.mup_PX / 1000., 
                               Tree.mup_PY / 1000., 
                               Tree.mup_PZ / 1000., 
                               Tree.mup_PE / 1000.);
                
                mum.SetPxPyPzE(Tree.mum_PX / 1000., 
                               Tree.mum_PY / 1000., 
                               Tree.mum_PZ / 1000., 
                               Tree.mum_PE / 1000.);
                
                Kmeson.SetPxPyPzE(Tree.K_PX / 1000., 
                                  Tree.K_PY / 1000., 
                                  Tree.K_PZ / 1000., 
                                  Tree.K_PE / 1000.);

                Jpsi.SetPxPyPzE(Tree.Jpsi_PX / 1000., 
                                Tree.Jpsi_PY / 1000., 
                                Tree.Jpsi_PZ / 1000., 
                                Tree.Jpsi_PE / 1000.);
                
                HFmeson = mup + mum + Kmeson;

                tr_HFjet.SetPxPyPzE(Tree.Jet_mcjet_PX / 1000.,
                                    Tree.Jet_mcjet_PY / 1000.,
                                    Tree.Jet_mcjet_PZ / 1000.,
                                    Tree.Jet_mcjet_PE / 1000.);

                tr_mup.SetPxPyPzE(Tree.mup_TRUEP_X / 1000., 
                                  Tree.mup_TRUEP_Y / 1000.,
                                  Tree.mup_TRUEP_Z / 1000., 
                                  Tree.mup_TRUEP_E / 1000.);

                tr_mum.SetPxPyPzE(Tree.mum_TRUEP_X / 1000., 
                                  Tree.mum_TRUEP_Y / 1000.,
                                  Tree.mum_TRUEP_Z / 1000., 
                                  Tree.mum_TRUEP_E / 1000.);

                tr_Kmeson.SetPxPyPzE(Tree.K_TRUEP_X / 1000., 
                                     Tree.K_TRUEP_Y / 1000.,
                                     Tree.K_TRUEP_Z / 1000., 
                                     Tree.K_TRUEP_E / 1000.);

                tr_HFmeson = tr_mup + tr_mum + tr_Kmeson;

                if (DoJESJER) {
                        const int n_iters = n_smearing_iter;
                        
                        for (int i_iter = 0; i_iter < n_iters; i_iter++) {
                                double rand = get_JES_JER(HFjet.Pt(), myRNG);  // Standard inclusive Z+jet values
                                // double rand = get_JES_JER(HFjet.Pt(), myRNG, DoJESJER);  // Low-multiplicity Z+jet values
                                
                                // Temp subtraction of HFmeson to perform JESJER.
                                HFjet -= HFmeson;
                                
                                //double newE2 = HFjet.E()*HFjet.E() + (rand*rand - 1) * HFjet.Pt()*HFjet.Pt();
                                double newE2 = HFjet.E()*HFjet.E() + (rand*rand - 1) * HFjet.P()*HFjet.P();
                                double newE = (newE2 < 0) ? 0 : std::sqrt(newE2);

                                HFjet.SetPxPyPzE(HFjet.Px()*rand, HFjet.Py()*rand, HFjet.Pz()*rand, newE);
                                
                                HFjet += HFmeson;
                        }
                }

                if (DoJetID) {
                        double mpt = 0;
                        double mtf = 0;
                        
                        int num_trk  = 0;
                        int num_neut = 0;
                        int num_part = 0;
                        
                        for (int dtrs0 = 0; dtrs0 < Tree.Jet_Dtr_nrecodtr; dtrs0++) {
                                if (fabs(Tree.Jet_Dtr_ThreeCharge[dtrs0]) == 0)
                                        num_neut++;
                                else
                                        num_trk++;
                        }

                        if (num_trk < 2)
                                continue;
                }


                bmass_dtf   = Tree.Bu_ConsBu_M[0] / 1000.;
                chi2ndf_dtf = Tree.Bu_ConsBu_chi2[0] / Tree.Bu_ConsBu_nDOF[0];
                tau_dtf     = Tree.Bu_ConsBu_ctau[0];

                nSV         = Tree.Jet_SVTag_Nvertices;
                jpsi_ipchi2 = log10(Tree.Jpsi_IPCHI2_OWNPV);
                k_ipchi2    = (Tree.K_IPCHI2_OWNPV);

                float leading_pT = 0;

                int n_maxpT_cand  = -999;
                int n_maxpT_entry = -999;
                int n_HFpt_entry  = -999;
                int HF_counter    = 0;

                bool hasHFhadron = false;

                int NumBHads = 0;

                for (int dtrs0 = 0; dtrs0 < Tree.Jet_Dtr_nrecodtr; dtrs0++) {
                        float trchi2ndf  = Tree.Jet_Dtr_TrackChi2[dtrs0] / Tree.Jet_Dtr_TrackNDF[dtrs0];
                        float dtr_charge = Tree.Jet_Dtr_ThreeCharge[dtrs0] / 3.;

                        dtr.SetPxPyPzE(Tree.Jet_Dtr_PX[dtrs0] / 1000.,
                                       Tree.Jet_Dtr_PY[dtrs0] / 1000.,
                                       Tree.Jet_Dtr_PZ[dtrs0] / 1000.,
                                       Tree.Jet_Dtr_E[dtrs0] / 1000.);
                        
                        // EFMC: I might want to check this for the PURITIES
                        tr_dtr.SetPxPyPzE(Tree.Jet_Dtr_TRUE_PX[dtrs0] / 1000.,
                                          Tree.Jet_Dtr_TRUE_PY[dtrs0] / 1000.,
                                          Tree.Jet_Dtr_TRUE_PZ[dtrs0] / 1000.,
                                          Tree.Jet_Dtr_TRUE_E[dtrs0] / 1000.);

                        if (abs(Tree.Jet_Dtr_ID[dtrs0]) != HF_pdgcode && 
                            !apply_chargedtrack_cuts(Tree.Jet_Dtr_ThreeCharge[dtrs0], 
                                                     dtr.P(), 
                                                     dtr.Pt(), 
                                                     trchi2ndf, 
                                                     Tree.Jet_Dtr_ProbNNghost[dtrs0], 
                                                     dtr.Rapidity()))
                                continue;

                        // EFMC: in contrast to MCMakeVar, here we know the B meson decays into a JPsi
                        // due to the trigger lines used.
                        jetdtrs.push_back(PseudoJet(Tree.Jet_Dtr_PX[dtrs0] / 1000.,
                                                    Tree.Jet_Dtr_PY[dtrs0] / 1000.,
                                                    Tree.Jet_Dtr_PZ[dtrs0] / 1000.,
                                                    Tree.Jet_Dtr_E[dtrs0] / 1000.));
                        
                        jetdtrs.back().set_user_info(new MyInfo(Tree.Jet_Dtr_ID[dtrs0]));
                        
                        if (abs(Tree.Jet_Dtr_ID[dtrs0]) == HF_pdgcode) {
                                HFmeson.SetPxPyPzE(dtr.Px(), dtr.Py(), dtr.Pz(), dtr.E());

                                // EFMC: I might want to check this for the PURITIES
                                HFmatch.SetPxPyPzE(Tree.Jet_Dtr_TRUE_PX[dtrs0] / 1000.,
                                                   Tree.Jet_Dtr_TRUE_PY[dtrs0] / 1000.,
                                                   Tree.Jet_Dtr_TRUE_PZ[dtrs0] / 1000.,
                                                   Tree.Jet_Dtr_TRUE_E[dtrs0] / 1000.);

                                hasHFhadron = true;
                                NumBHads++;
                        }

                        dtr_id.push_back(Tree.Jet_Dtr_ID[dtrs0]);
                        dtr_pt.push_back(dtr.Pt());
                        dtr_rap.push_back(dtr.Rapidity());
                        dtr_3charge.push_back(Tree.Jet_Dtr_ThreeCharge[dtrs0]);
                }

                if (!hasHFhadron)
                        continue;

                NumBHads_tr = 0;
                bool hasHFhadron_matched = false;
                
                for (int dtrs0 = 0; dtrs0 < Tree.Jet_mcjet_nmcdtrs; dtrs0++) {
                        float trchi2ndf = 0;

                        dtr_matchtruthjet.SetPxPyPzE(Tree.Jet_mcjet_dtrPX[dtrs0] / 1000.,
                                                     Tree.Jet_mcjet_dtrPY[dtrs0] / 1000.,
                                                     Tree.Jet_mcjet_dtrPZ[dtrs0] / 1000.,
                                                     Tree.Jet_mcjet_dtrE[dtrs0] / 1000.);

                        if (abs(Tree.Jet_mcjet_dtrID[dtrs0]) != HF_pdgcode && 
                                !apply_chargedtrack_momentum_cuts(Tree.Jet_mcjet_dtrThreeCharge[dtrs0], 
                                                                  dtr_matchtruthjet.P(), 
                                                                  dtr_matchtruthjet.Pt(),
                                                                  dtr_matchtruthjet.Rapidity()))
                                continue;

                        if (abs(Tree.Jet_mcjet_dtrID[dtrs0]) == HF_pdgcode) {
                                NumBHads_tr++;

                                if (fabs(dtr_matchtruthjet.Pt() - Tree.Bu_TRUEPT / 1000.) < 0.01) {
                                        tr_HFmeson_injet.SetPxPyPzE(dtr_matchtruthjet.Px(), 
                                                                    dtr_matchtruthjet.Py(), 
                                                                    dtr_matchtruthjet.Pz(), 
                                                                    dtr_matchtruthjet.E());

                                        h2_HFpt_RM->Fill(dtr.Pt(), HFmeson.Pt());

                                        hasHFhadron_matched = true;
                                } else {
                                        true_jetdtrs.push_back(PseudoJet(Tree.Jet_mcjet_dtrPX[dtrs0] / 1000.,
                                                                         Tree.Jet_mcjet_dtrPY[dtrs0] / 1000.,
                                                                         Tree.Jet_mcjet_dtrPZ[dtrs0] / 1000.,
                                                                         Tree.Jet_mcjet_dtrE[dtrs0] / 1000.));
                                        
                                        true_jetdtrs.back().set_user_info(new MyInfo(-999));
                                }
                        } else {
                                true_jetdtrs.push_back(PseudoJet(Tree.Jet_mcjet_dtrPX[dtrs0] / 1000.,
                                                                 Tree.Jet_mcjet_dtrPY[dtrs0] / 1000.,
                                                                 Tree.Jet_mcjet_dtrPZ[dtrs0] / 1000.,
                                                                 Tree.Jet_mcjet_dtrE[dtrs0] / 1000.));
                
                                true_jetdtrs.back().set_user_info(new MyInfo(Tree.Jet_mcjet_dtrID[dtrs0]));
                        }
                }

                if (hasHFhadron_matched) {
                        true_jetdtrs.push_back(PseudoJet(tr_HFmeson_injet.Px(),
                                                         tr_HFmeson_injet.Py(),
                                                         tr_HFmeson_injet.Pz(),
                                                         tr_HFmeson_injet.E()));

                        true_jetdtrs.back().set_user_info(new MyInfo(HF_pdgcode));
                }

                SVTag = 0;
                Hasbbbar = false;

                if (Tree.hasb && Tree.hasbbar) {
                        Hasbbbar = true;
                        
                        if (Tree.Jet_SVTag_Tag)
                                SVTag = 1;
                }

                if (Tree.Jet_SVTag_Tag && !(Tree.hasb && Tree.hasbbar)) // && !(Tree.tr_hasb && Tree.tr_hasbbar))
                        SVTag = 2;

                TVector3 HF_meson = HFmeson.Vect();
                TVector3 HF_jet   = HFjet.Vect();

                TVector3 tr_HF_meson = tr_HFmeson.Vect();
                TVector3 tr_HF_jet   = tr_HFjet.Vect();

                jt = (HF_jet.Cross(HF_meson).Mag()) / (HF_jet.Mag());
                z  = (HF_meson.Dot(HF_jet) ) / (HF_jet.Mag2() );
                r  = static_cast < TLorentzVector > (HFmeson).DeltaR(HFjet, kTRUE);
                
                if(hasHFhadron_matched ) {
                        tr_jt = (tr_HF_jet.Cross(tr_HF_meson).Mag()) / (tr_HF_jet.Mag());
                        tr_z  = (tr_HF_meson.Dot(tr_HF_jet) ) / (tr_HF_jet.Mag2() );
                        tr_r  =  static_cast < TLorentzVector > (tr_HFmeson).DeltaR(tr_HFjet, kTRUE);

                        if(Hasbbbar) {
                                jtg = (tr_HF_jet.Cross(tr_HF_meson).Mag()) / (tr_HF_jet.Mag());
                                zg  = (tr_HF_meson.Dot(tr_HF_jet) ) / (tr_HF_jet.Mag2() );
                                rg  = static_cast < TLorentzVector > (tr_HFmeson).DeltaR(tr_HFjet, kTRUE);
                        } else {
                                jtg = -999.;
                                zg  = -999.;
                                rg  = -999.;
                        }
                } else {
                        tr_jt = -999.;
                        tr_z  = -999.;
                        tr_r  = -999.;
                }

                // Check WTA Truth
                if (true_jetdtrs.size() > 0){
                        WTA_true_jets = ClusterSequence(true_jetdtrs, WTA);

                        WTA_true_jet = sorted_by_pt(WTA_true_jets.inclusive_jets())[0];
                        
                        WTA_true_axis.SetPxPyPzE(WTA_true_jet.px(), WTA_true_jet.py(), WTA_true_jet.pz(), WTA_true_jet.e());
                        WTA_true_dist = tr_HFmeson_injet.DeltaR(WTA_true_axis, true);
                } else {
                        WTA_true_dist = -999;
                }
                // Check WTA Reco
                WTA_reco_jets = ClusterSequence(jetdtrs, WTA);

                WTA_reco_jet = sorted_by_pt(WTA_reco_jets.inclusive_jets())[0];
                
                WTA_reco_axis.SetPxPyPzE(WTA_reco_jet.px(), WTA_reco_jet.py(), WTA_reco_jet.pz(), WTA_reco_jet.e());
                WTA_reco_dist = HFmeson.DeltaR(WTA_reco_axis, true);
                
                jet_pt  = HFjet.Pt();
                jet_eta = HFjet.Eta();
                jet_rap = HFjet.Rapidity();
                
                jet_px = HFjet.Px();
                jet_py = HFjet.Py();
                jet_pz = HFjet.Pz();
                jet_e  = HFjet.E();
                
                HF_px = HFmeson.Px();
                HF_py = HFmeson.Py();
                HF_pz = HFmeson.Pz();
                HF_e  = HFmeson.E();
                HF_pt = HFmeson.Pt();

                Bu_CHI2NDOF = Tree.Bu_ENDVERTEX_CHI2 / Tree.Bu_ENDVERTEX_NDOF;
                Bu_IPCHI2   = Tree.Bu_IPCHI2_OWNPV;
                Bu_CHI2     = Tree.Bu_ENDVERTEX_CHI2;
                // Bu_CHI2 = Tree.Bu_OWNPV_CHI2;
                
                mum_px = mum.Px();
                mum_py = mum.Py();
                mum_pz = mum.Pz();
                mum_e  = mum.E();
                
                mum_CHI2NDOF  = Tree.mum_TRACK_CHI2NDOF;
                mum_GHOSTPROB = Tree.mum_TRACK_GhostProb;
                mum_IPCHI2    = Tree.mum_IPCHI2_OWNPV;

                mup_px = mup.Px();
                mup_py = mup.Py();
                mup_pz = mup.Pz();
                mup_e  = mup.E();
                
                mup_CHI2NDOF  = Tree.mup_TRACK_CHI2NDOF;
                mup_GHOSTPROB = Tree.mup_TRACK_GhostProb;
                mup_IPCHI2    = Tree.mup_IPCHI2_OWNPV;

                K_px = Kmeson.Px();
                K_py = Kmeson.Py();
                K_pz = Kmeson.Pz();
                K_e  = Kmeson.E();
                K_p  = Kmeson.P() * 1000; // DUDE WHAT
                K_eta = Kmeson.Eta();
                
                K_CHI2NDOF  = Tree.K_TRACK_CHI2NDOF;
                K_GHOSTPROB = Tree.K_TRACK_GhostProb;
                K_IPCHI2    = Tree.K_IPCHI2_OWNPV;
                K_PIDK      = Tree.K_PIDK;

                Jpsi_px = Jpsi.Px();
                Jpsi_py = Jpsi.Py();
                Jpsi_pz = Jpsi.Pz();
                Jpsi_e  = Jpsi.E();

                Jpsi_CHI2NDOF = Tree.Jpsi_ENDVERTEX_CHI2 / Tree.Jpsi_ENDVERTEX_NDOF;
                Jpsi_IPCHI2   = Tree.Jpsi_IPCHI2_OWNPV;
                Jpsi_CHI2     = Tree.Jpsi_ENDVERTEX_CHI2;
                Jpsi_FDCHI2   = Tree.Jpsi_FDCHI2_OWNPV;
                Jpsi_BPVDLS   = Tree.Jpsi_BPVDLS;
                
                nTracks  = Tree.nTracks;
                nSPDHits = Tree.nSPDHits;

                tr_jet_eta = tr_HFjet.Eta();
                tr_jet_rap = tr_HFjet.Rapidity();
                tr_jet_pt  = tr_HFjet.Pt();
                tr_jet_px  = tr_HFjet.Px();
                tr_jet_py  = tr_HFjet.Py();
                tr_jet_pz  = tr_HFjet.Pz();
                tr_jet_e   = tr_HFjet.E();

                tr_HF_pt = tr_HFmeson.Pt();
                tr_HF_px = tr_HFmeson.Px();
                tr_HF_py = tr_HFmeson.Py();
                tr_HF_pz = tr_HFmeson.Pz();
                tr_HF_e  = tr_HFmeson.E();

                tr_mum_px = tr_mum.Px();
                tr_mum_py = tr_mum.Py();
                tr_mum_pz = tr_mum.Pz();
                tr_mum_e  = tr_mum.E();
                
                tr_mup_px = tr_mup.Px();
                tr_mup_py = tr_mup.Py();
                tr_mup_pz = tr_mup.Pz();
                tr_mup_e  = tr_mup.E();
                
                tr_K_px = tr_Kmeson.Px();
                tr_K_py = tr_Kmeson.Py();
                tr_K_pz = tr_Kmeson.Pz();
                tr_K_e  = tr_Kmeson.E();

                sv_mass   = Tree.Jet_SVTag_SigMaxMass / 1000.;
                sv_chi2   = log10(Tree.Jet_SVTag_SigMaxVtxChi2NDF);
                sv_ntrks  = Tree.Jet_SVTag_SigMaxNtracks;
                sv_cosine = Tree.Jet_SVTag_SigMaxDirAngleS2S;

                // MCReco case set to 1. What about when doing the closure test?
                pideff_K = pideff_mum = pideff_mup = trkeff_ratio_K = trkeff_ratio_mup = trkeff_ratio_mum = 1.0;
                
                if (hasHFhadron_matched)
                        isTrueBjet = true;
                else
                        isTrueBjet = false;

                last_eventNum = Tree.eventNumber;
                eventNumber   = Tree.eventNumber;

                events++;

                BTree->Fill();
        }

        cout << "Number of saved events = " << events << endl;
        f.Write();
        f.Close();

        benchmark->Show("MakeVarTreeMCReco");
}
