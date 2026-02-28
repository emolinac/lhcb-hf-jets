/*
   This macro creates a tree from a MC ROOT file that stores "TRUTH" Monte Carlo variables
*/

#include <iostream>
#include <TCanvas.h>
#include <vector>

#include "fastjet/ClusterSequence.hh"
// #include "fastjet/contrib/SoftDrop.hh"

#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/TBJetsMC.h"
#include "../include/TBJetsMC.C"

using namespace fastjet;
using namespace std;

// All events: -1
void MakeVarTreeMC(int NumEvts_user = -1) 
{    
        TBenchmark* benchmark = new TBenchmark();
        benchmark->Start("MakeVarTreeMC");

        int NumEvts      = NumEvts_user;
        int NumEvtsTruth = NumEvts_user;
        
        int flavor     = 5; // beauty
        int HF_pdgcode = -999;

        double mass_num;

        TString RecoHF = "OnlyOneRecoHF";

        if (flavor == 5) {
                mass_num   = 4.2;
                HF_pdgcode = 521;
        } else if (flavor == 4) {
                mass_num   = 1.25;
                HF_pdgcode = 421;
        } else if (flavor == 1) {
                mass_num      = 0.001;
                followHardest = true; 
        }

        JetDefinition jet_def(cambridge_algorithm, JetDefinition::max_allowable_R);
        JetDefinition WTA(cambridge_algorithm, JetDefinition::max_allowable_R, WTA_pt_scheme);

        vector<PseudoJet> jetdtrs, meas_jetdtrs;

        double WTA_true_dist;
        double WTA_reco_dist;

        TBJetsMC Tree;
    
        cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;
        
        if (NumEvts == -1)
                NumEvts = Tree.fChain->GetEntries();

        TFile f((output_folder + "ntuple_bjets_mc.root").c_str(), "RECREATE");
        
        // FF(z) Histograms 
        TH1D *h1_z_truth   = new TH1D("z_truth"  , "", 80, 0.0, 1.01);
        TH1D *h1_z_truth_b = new TH1D("z_truth_b", "", 80, 0.0, 1.01);
        TH1D *h1_z_truth_g = new TH1D("z_truth_g", "", 80, 0.0, 1.01);

        // Truth Variables (from Truth Tree)
        double jet_pt, jet_eta, meas_jet_pt, meas_jet_eta;
        double jet_px, jet_py, jet_pz, jet_e;
        double HF_px, HF_py, HF_pz, HF_e;
        double HF_pt;
        double jet_rap, meas_jet_rap;
        double K_px, K_py, K_pz, K_e;
        double mup_px, mup_py, mup_pz, mup_e;
        double mum_px, mum_py, mum_pz, mum_e;

        // Reco Variables (from Truth Tree)
        double meas_jet_px, meas_jet_py, meas_jet_pz, meas_jet_e;
        double meas_mup_px, meas_mup_py, meas_mup_pz, meas_mup_e;
        double meas_mum_px, meas_mum_py, meas_mum_pz, meas_mum_e;
        double meas_K_px, meas_K_py, meas_K_pz, meas_K_e;
        double meas_HF_px, meas_HF_py, meas_HF_pz, meas_HF_e;
        double meas_HF_pt;

        double truth_z, truth_jt, truth_r;
        double truth_z_b, truth_jt_b, truth_r_b;
        double truth_zg, truth_jtg, truth_rg;
        double meas_z, meas_jt, meas_r;
        double meas_inJet_z, meas_inJet_jt, meas_inJet_r;
        double jet_pt_recotruthratio, HF_pt_recotruthratio;

        int nsplits, ndtrs;
        int meas_nsplits, meas_ndtrs;
        int GluonTag, nTracks;
        int NumHFHads, eventNumber;
        int NumDtrRecoHF;

        double ndtrs_mc;
        double ndtrs_mcreco;

        bool hasRecoHF;
        bool hasbbbar;
        bool isSingle_tr_Bjet;
        bool isPhoton_tr_Bjet;

        // TLorentzVector
        TTree *BTree = new TTree("BTree", "B-jets Tree Variables");
        TNtuple* ntuple_true_jet  = new TNtuple("ntuple_true_jet" ,"","jet_pt:wta_distance");
        TNtuple* ntuple_true_pair = new TNtuple("ntuple_true_pair","","jet_pt:R_L:weight_pt:wta_distance");

        TNtuple* ntuple_reco_jet  = new TNtuple("ntuple_reco_jet" ,"","jet_pt:wta_distance");
        TNtuple* ntuple_reco_pair = new TNtuple("ntuple_reco_pair","","jet_pt:R_L:weight_pt:wta_distance");

        BTree->Branch("eventNumber", &eventNumber);

        BTree->Branch("jet_pt", &jet_pt);
        BTree->Branch("jet_eta", &jet_eta);
        BTree->Branch("jet_rap", &jet_rap);
        BTree->Branch("jet_px", &jet_px);
        BTree->Branch("jet_py", &jet_py);
        BTree->Branch("jet_pz", &jet_pz);
        BTree->Branch("jet_e" , &jet_e);

        BTree->Branch("HF_px", &HF_px);
        BTree->Branch("HF_py", &HF_py);
        BTree->Branch("HF_pz", &HF_pz);
        BTree->Branch("HF_e", &HF_e);
        BTree->Branch("HF_pt", &HF_pt);

        BTree->Branch("mum_px", &mum_px);
        BTree->Branch("mum_py", &mum_py);
        BTree->Branch("mum_pz", &mum_pz);
        BTree->Branch("mum_e", &mum_e);
        BTree->Branch("mup_px", &mup_px);
        BTree->Branch("mup_py", &mup_py);
        BTree->Branch("mup_pz", &mup_pz);
        BTree->Branch("mup_e", &mup_e);
        BTree->Branch("K_px", &K_px);
        BTree->Branch("K_py", &K_py);
        BTree->Branch("K_pz", &K_pz);
        BTree->Branch("K_e", &K_e);

        BTree->Branch("meas_jet_pt", &meas_jet_pt);
        BTree->Branch("meas_jet_eta", &meas_jet_eta);
        BTree->Branch("meas_jet_rap", &meas_jet_rap);
        BTree->Branch("meas_jet_px", &meas_jet_px);
        BTree->Branch("meas_jet_py", &meas_jet_py);
        BTree->Branch("meas_jet_pz", &meas_jet_pz);
        BTree->Branch("meas_jet_e", &meas_jet_e);

        BTree->Branch("meas_HF_px", &meas_HF_px);
        BTree->Branch("meas_HF_py", &meas_HF_py);
        BTree->Branch("meas_HF_pz", &meas_HF_pz);
        BTree->Branch("meas_HF_e", &meas_HF_e);
        BTree->Branch("meas_HF_pt", &meas_HF_pt);

        // BTree->Branch("meas_mum_px", &meas_mum_px);
        // BTree->Branch("meas_mum_py", &meas_mum_py);
        // BTree->Branch("meas_mum_pz", &meas_mum_pz);
        // BTree->Branch("meas_mum_e", &meas_mum_e);
        // BTree->Branch("meas_mup_px", &meas_mup_px);
        // BTree->Branch("meas_mup_py", &meas_mup_py);
        // BTree->Branch("meas_mup_pz", &meas_mup_pz);
        // BTree->Branch("meas_mup_e", &meas_mup_e);
        // BTree->Branch("meas_K_px", &meas_K_px);
        // BTree->Branch("meas_K_py", &meas_K_py);
        // BTree->Branch("meas_K_pz", &meas_K_pz);
        // BTree->Branch("meas_K_e", &meas_K_e);

        BTree->Branch("meas_z", &meas_z);
        BTree->Branch("meas_jt", &meas_jt);
        BTree->Branch("meas_r", &meas_r);

        BTree->Branch("meas_inJet_z", &meas_inJet_z);
        BTree->Branch("meas_inJet_jt", &meas_inJet_jt);
        BTree->Branch("meas_inJet_r", &meas_inJet_r);
        
        BTree->Branch("jet_pt_recotruthratio", &jet_pt_recotruthratio);
        BTree->Branch("HF_pt_recotruthratio", &HF_pt_recotruthratio);
        BTree->Branch("GluonTag", &GluonTag);
        BTree->Branch("nTracks", &nTracks);
        BTree->Branch("hasbbbar", &hasbbbar);
        BTree->Branch("NumHFHads", &NumHFHads);
        BTree->Branch("NumDtrRecoHF", &NumDtrRecoHF);
        BTree->Branch("hasRecoHF", &hasRecoHF);
        BTree->Branch("isSingle_tr_Bjet", &isSingle_tr_Bjet);
        BTree->Branch("isPhoton_tr_Bjet", &isPhoton_tr_Bjet);
        BTree->Branch("ndtrs_mc",     &ndtrs_mc);
        BTree->Branch("ndtrs_mcreco", &ndtrs_mcreco);
        
        BTree->Branch("truth_z", &truth_z);
        BTree->Branch("truth_jt", &truth_jt);
        BTree->Branch("truth_r", &truth_r);
        BTree->Branch("truth_z_b", &truth_z_b);
        BTree->Branch("truth_jt_b", &truth_jt_b);
        BTree->Branch("truth_r_b", &truth_r_b);
        BTree->Branch("truth_zg", &truth_zg);
        BTree->Branch("truth_jtg", &truth_jtg);
        BTree->Branch("truth_rg", &truth_rg);

        BTree->Branch("WTA_true_dist", &WTA_true_dist);
        BTree->Branch("WTA_reco_dist", &WTA_reco_dist);

        // Event loop
        int eventNum;
        unsigned long long last_eventNum = 0;
        int events = 0;

        bool maxjetpT_found = false;
        int NumRecoHF = 0;

        int cut_HFdR = 0;
        
        int Num_Single_tr_Bjet = 0;
        int Num_Photon_tr_Bjet = 0;

        TLorentzVector HFjet, recojet, meas_HFjet, HFmeson, mup, mum, Kmeson, Jpsi;
        TLorentzVector meas_Zdtr3, meas_mup, meas_mum, meas_HFmeson, meas_Kmeson;
        TLorentzVector dtr0, dtr1;
        TLorentzVector WTA_true_axis;
        TLorentzVector WTA_reco_axis;

        ClusterSequence WTA_true_jets;
        ClusterSequence WTA_reco_jets;
        PseudoJet WTA_true_jet;
        PseudoJet WTA_reco_jet;
                
        for (int ev = 0; ev < NumEvts; ev++) {
                jetdtrs.clear();
                meas_jetdtrs.clear();
                
                Tree.GetEntry(ev);

                if (ev%10000 == 0) {
                        double percentage = 100.*ev/NumEvts;
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (ev != 0)
                        if (last_eventNum == Tree.eventNumber)
                                continue;

                events++;

                if (Tree.nPVs > 1) 
                        continue;

                mup.SetPxPyPzE(Tree.MCJet_truth_mup_PX / 1000.,
                               Tree.MCJet_truth_mup_PY / 1000.,
                               Tree.MCJet_truth_mup_PZ / 1000.,
                               Tree.MCJet_truth_mup_PE / 1000.);

                if (!apply_muon_cuts(mup.Pt()))
                        continue;
                
                mum.SetPxPyPzE(Tree.MCJet_truth_mum_PX / 1000.,
                               Tree.MCJet_truth_mum_PY / 1000.,
                               Tree.MCJet_truth_mum_PZ / 1000.,
                               Tree.MCJet_truth_mum_PE / 1000.);

                if (!apply_muon_cuts(mum.Pt()))
                        continue;
                
                Kmeson.SetPxPyPzE(Tree.MCJet_truth_K_PX / 1000.,
                                  Tree.MCJet_truth_K_PY / 1000.,
                                  Tree.MCJet_truth_K_PZ / 1000.,
                                  Tree.MCJet_truth_K_PE / 1000.);

                HFjet.SetPxPyPzE(Tree.MCJet_PX / 1000.,
                                 Tree.MCJet_PY / 1000.,
                                 Tree.MCJet_PZ / 1000.,
                                 Tree.MCJet_PE / 1000.);

                if (!apply_jet_cuts(HFjet.Rapidity(), HFjet.Pt()))
                        continue;

                //HFmeson = mup + mum + Kmeson;
                HFmeson.SetPxPyPzE(Tree.MCJet_truth_B_PX / 1000.,
                                   Tree.MCJet_truth_B_PY / 1000.,
                                   Tree.MCJet_truth_B_PZ / 1000.,
                                   Tree.MCJet_truth_B_PE / 1000.);        
                
                Jpsi = mup + mum;

                double HF_jet_truedR = static_cast < TLorentzVector > (HFmeson).DeltaR(HFjet, true);

                if (HF_jet_truedR > jetradius) 
                        continue;
        
                meas_HFjet.SetPxPyPzE(Tree.MCJet_recojet_PX / 1000.,
                                      Tree.MCJet_recojet_PY / 1000.,
                                      Tree.MCJet_recojet_PZ / 1000.,
                                      Tree.MCJet_recojet_PE / 1000.);

                if (!apply_jet_cuts(meas_HFjet.Rapidity(), meas_HFjet.Pt()))
                        continue;

                double matched_mup_e = std::sqrt(pow(Tree.MCJet_truth_match_mup_PX, 2) +
                                                 pow(Tree.MCJet_truth_match_mup_PY, 2) +
                                                 pow(Tree.MCJet_truth_match_mup_PZ, 2));
                
                meas_mup.SetPxPyPzE(Tree.MCJet_truth_match_mup_PX / 1000.,
                                    Tree.MCJet_truth_match_mup_PY / 1000.,
                                    Tree.MCJet_truth_match_mup_PZ / 1000.,
                                    matched_mup_e / 1000.);
                
                double matched_mum_e = std::sqrt(pow(Tree.MCJet_truth_match_mum_PX, 2) +
                                                 pow(Tree.MCJet_truth_match_mum_PY, 2) +
                                                 pow(Tree.MCJet_truth_match_mum_PZ, 2));
                
                meas_mum.SetPxPyPzE(Tree.MCJet_truth_match_mum_PX / 1000.,
                                    Tree.MCJet_truth_match_mum_PY / 1000.,
                                    Tree.MCJet_truth_match_mum_PZ / 1000.,
                                    matched_mum_e / 1000.);
                
                double matched_K_e = std::sqrt(pow(Tree.MCJet_truth_match_K_PX, 2) +
                                               pow(Tree.MCJet_truth_match_K_PY, 2) +
                                               pow(Tree.MCJet_truth_match_K_PZ, 2));
                
                meas_Kmeson.SetPxPyPzE(Tree.MCJet_truth_match_K_PX / 1000.,
                                       Tree.MCJet_truth_match_K_PY / 1000.,
                                       Tree.MCJet_truth_match_K_PZ / 1000.,
                                       matched_K_e / 1000.);

                meas_HFmeson = meas_mup + meas_mum + meas_Kmeson;

                double HF_jet_measdR = static_cast < TLorentzVector > (meas_HFmeson).DeltaR(meas_HFjet, true); 
                
                if (HF_jet_measdR > jetradius) 
                        continue;
                
                hasbbbar = false;
                // GluonTag = false;
                
                if (Tree.hasb && Tree.hasbbar) {
                        hasbbbar = true;
                        // GluonTag = true;
                }
                
                ///// Experimenting with Stripping Line cuts //// //
                // if (Kmeson.Pt() < 0.25 || mup.Pt() < 0.25 || mum.Pt() < 0.25)
                //     continue;

                // if (HF_jet_measdR > jetradius) {
                //      cut_HFdR++;
                //      continue;
                // }

                // jet_Nmcdtrs = recojetNdtrs;
                
                bool hasHFhadron = false;
                bool hasJpsi     = false;
                bool isGluonJet  = false;
                bool hasbquark   = false;
                bool hasbbar     = false;
                
                GluonTag     = false;
                NumHFHads    = 0;
                NumDtrRecoHF = 0;
                
                int quark_pdg = -999;

                // Loop to check if HFMeson is in the jet
                for (int dtrs0 = 0; dtrs0 < Tree.MCJet_Dtr_nmcdtrs; dtrs0++) {
                        dtr0.SetPxPyPzE(Tree.MCJet_Dtr_PX[dtrs0] / 1000.,
                                        Tree.MCJet_Dtr_PY[dtrs0] / 1000.,
                                        Tree.MCJet_Dtr_PZ[dtrs0] / 1000.,
                                        Tree.MCJet_Dtr_E[dtrs0]  / 1000.);

                        if (std::abs(Tree.MCJet_Dtr_ID[dtrs0]) != HF_pdgcode && 
                            !apply_particle_momentum_cuts(dtr0.P(),
                                                          dtr0.Pt(),
                                                          dtr0.Rapidity()))
                                continue;

                        // Check for HF branch and jet dtr
                        if (std::abs(Tree.MCJet_Dtr_ID[dtrs0]) == HF_pdgcode) {
                                if (std::fabs(dtr0.Px() - HFmeson.Px()) < 1e-2 && std::fabs(dtr0.Py() - HFmeson.Py()) < 1e-2) {
                                        NumHFHads++;
                                        hasHFhadron = true;

                                        // This is the signal, follow it during reclustering
                                        jetdtrs.push_back(PseudoJet(Tree.MCJet_Dtr_PX[dtrs0] / 1000.,
                                                                    Tree.MCJet_Dtr_PY[dtrs0] / 1000.,
                                                                    Tree.MCJet_Dtr_PZ[dtrs0] / 1000.,
                                                                    Tree.MCJet_Dtr_E[dtrs0]  / 1000.));

                                        jetdtrs.back().set_user_info(new MyInfo(Tree.MCJet_Dtr_ID[dtrs0]));
                                } else {
                                        NumHFHads++;
                                        
                                        // If this is another 521 (Bmeson) that didn't decay to Jpsi K
                                        // do not follow it during reclustering, set ID to -999
                                        jetdtrs.push_back(PseudoJet(Tree.MCJet_Dtr_PX[dtrs0] / 1000.,
                                                                    Tree.MCJet_Dtr_PY[dtrs0] / 1000.,
                                                                    Tree.MCJet_Dtr_PZ[dtrs0] / 1000.,
                                                                    Tree.MCJet_Dtr_E[dtrs0]  / 1000.));
                                        
                                        jetdtrs.back().set_user_info(new MyInfo(-999));
                                }
                        } else {
                                jetdtrs.push_back(PseudoJet(Tree.MCJet_Dtr_PX[dtrs0] / 1000.,
                                                            Tree.MCJet_Dtr_PY[dtrs0] / 1000.,
                                                            Tree.MCJet_Dtr_PZ[dtrs0] / 1000.,
                                                            Tree.MCJet_Dtr_E[dtrs0]  / 1000.));
                                
                                jetdtrs.back().set_user_info(new MyInfo(Tree.MCJet_Dtr_ID[dtrs0]));
                        }

                        if (Tree.MCJet_Dtr_TopMotherID[dtrs0] == 5)
                                hasbquark = true;

                        if (Tree.MCJet_Dtr_TopMotherID[dtrs0] == -5)
                                hasbbar = true;
                }
                
                // Veto events that don't have a B meson
                if (!hasHFhadron)
                        continue;

                // Check WTA Truth
                // ClusterSequence WTA_true_jets(jetdtrs, WTA);

                // PseudoJet WTA_true_jet = sorted_by_pt(WTA_true_jets.inclusive_jets())[0];

                WTA_true_jets = ClusterSequence(jetdtrs, WTA);

                WTA_true_jet = sorted_by_pt(WTA_true_jets.inclusive_jets())[0];
                
                WTA_true_axis.SetPxPyPzE(WTA_true_jet.px(), WTA_true_jet.py(), WTA_true_jet.pz(), WTA_true_jet.e());
                WTA_true_dist = HFmeson.DeltaR(WTA_true_axis, true);

                for (int dtrs0 = 0; dtrs0 < Tree.MCJet_Dtr_nmcdtrs; dtrs0++) {
                        if (std::abs(Tree.MCJet_Dtr_ID[dtrs0]) < 100)
                                continue;

                        dtr0.SetPxPyPzE(Tree.MCJet_Dtr_PX[dtrs0] / 1000.,
                                        Tree.MCJet_Dtr_PY[dtrs0] / 1000.,
                                        Tree.MCJet_Dtr_PZ[dtrs0] / 1000.,
                                        Tree.MCJet_Dtr_E[dtrs0]  / 1000.);

                        if (std::abs(Tree.MCJet_Dtr_ID[dtrs0]) != HF_pdgcode && 
                            !apply_particle_momentum_cuts(dtr0.P(),
                                                          dtr0.Pt(),
                                                          dtr0.Rapidity()))
                                continue;
                        
                        for (int dtrs1 = dtrs0 + 1; dtrs1 < Tree.MCJet_Dtr_nmcdtrs; dtrs1++) {
                                if (std::abs(Tree.MCJet_Dtr_ID[dtrs1]) < 100)
                                        continue;

                                dtr1.SetPxPyPzE(Tree.MCJet_Dtr_PX[dtrs1] / 1000.,
                                                Tree.MCJet_Dtr_PY[dtrs1] / 1000.,
                                                Tree.MCJet_Dtr_PZ[dtrs1] / 1000.,
                                                Tree.MCJet_Dtr_E[dtrs1]  / 1000.);
                        
                                if (std::abs(Tree.MCJet_Dtr_ID[dtrs1]) != HF_pdgcode && 
                                    !apply_particle_momentum_cuts(dtr1.P(),
                                                                  dtr1.Pt(),
                                                                  dtr1.Rapidity()))
                                        continue;

                                ntuple_true_pair->Fill(HFjet.Pt(), dtr0.DeltaR(dtr1, true), (dtr0.Pt() * dtr1.Pt()) / std::pow(HFjet.Pt(), 2),WTA_true_dist);
                        }
                }

                // if (NumHFHads > 1) 
                //      continue;

                if (hasbquark && hasbbar)
                        GluonTag = true;
        
                /// Experimenting Gluon Jet Daughters
                //  for (int i = 0; i <  Tree.MCJet_Dtr_nmcdtrs; i++) {
                //      if (GluonTag) {
                //          cout << "Gluon jet Dtrs: " <<  endl;
                //          cout << Tree.MCJet_Dtr_ID[i] <<endl;
                //      }
                //  }
                
                hasRecoHF = false;
                
                int NumHFinRecoJet = 0;
                
                // NOTE : Check if the reco HFJet has the HF Meson inside
                for (int dtrs0 = 0; dtrs0 < Tree.MCJet_recojet_nrecodtrs; dtrs0++) {
                        double trchi2ndf = Tree.MCJet_recojet_Dtr_TrackChi2[dtrs0] / Tree.MCJet_recojet_Dtr_TrackNDF[dtrs0];

                        dtr0.SetPxPyPzE(Tree.MCJet_recojet_Dtr_PX[dtrs0] / 1000.,
                                        Tree.MCJet_recojet_Dtr_PY[dtrs0] / 1000.,
                                        Tree.MCJet_recojet_Dtr_PZ[dtrs0] / 1000.,
                                        Tree.MCJet_recojet_Dtr_E[dtrs0]  / 1000.);

                        if (std::abs(Tree.MCJet_recojet_Dtr_ID[dtrs0]) != HF_pdgcode && 
                            !apply_particle_cuts(dtr0.P(),
                                                 dtr0.Pt(),
                                                 trchi2ndf,
                                                 Tree.MCJet_recojet_Dtr_ProbNNghost[dtrs0],
                                                 dtr0.Rapidity()))
                                continue;

                        meas_jetdtrs.push_back(PseudoJet(Tree.MCJet_recojet_Dtr_PX[dtrs0] / 1000.,
                                                         Tree.MCJet_recojet_Dtr_PY[dtrs0] / 1000.,
                                                         Tree.MCJet_recojet_Dtr_PZ[dtrs0] / 1000.,
                                                         Tree.MCJet_recojet_Dtr_E[dtrs0]  / 1000.));
                        
                        meas_jetdtrs.back().set_user_info(new MyInfo(Tree.MCJet_recojet_Dtr_ID[dtrs0]));

                        if (abs(Tree.MCJet_recojet_Dtr_ID[dtrs0]) == HF_pdgcode) {
                                meas_HFmeson.SetPxPyPzE(dtr0.Px(), dtr0.Py(), dtr0.Pz(), dtr0.E());
                                
                                NumDtrRecoHF++;
                                NumHFinRecoJet++;
                                
                                hasRecoHF = true;
                                
                                if (NumHFinRecoJet == 2)
                                        cout<< "Found two +-521 in reco jet: " << Tree.MCJet_recojet_Dtr_ID[dtrs0] << endl;
                        }
                }

                // If reconstructed HF meson (from reco Daughters) is not inside the reco jet, then the Reco HF jet end up without the HF meson inside //// 
                if (static_cast < TLorentzVector > (meas_HFmeson).DeltaR(meas_HFjet, true) > jetradius) 
                        hasRecoHF = false;

                // if (meas_HFmeson.DeltaR(meas_HFjet) > jetradius) 
                //      hasRecoHF = false;

                if (hasRecoHF)
                        NumRecoHF++;

                // Check WTA reco
                // ClusterSequence WTA_reco_jets(meas_jetdtrs, WTA);

                // PseudoJet WTA_reco_jet = sorted_by_pt(WTA_reco_jets.inclusive_jets())[0];

                WTA_reco_jets = ClusterSequence(meas_jetdtrs, WTA);

                WTA_reco_jet = sorted_by_pt(WTA_reco_jets.inclusive_jets())[0];
                
                WTA_reco_axis.SetPxPyPzE(WTA_reco_jet.px(), WTA_reco_jet.py(), WTA_reco_jet.pz(), WTA_reco_jet.e());
                WTA_reco_dist = meas_HFmeson.DeltaR(WTA_reco_axis, true);

                for (int dtrs0 = 0; dtrs0 < Tree.MCJet_recojet_nrecodtrs; dtrs0++) {
                        if (std::abs(Tree.MCJet_recojet_Dtr_ID[dtrs0]) < 100)
                                continue;

                        dtr0.SetPxPyPzE(Tree.MCJet_recojet_Dtr_PX[dtrs0] / 1000.,
                                        Tree.MCJet_recojet_Dtr_PY[dtrs0] / 1000.,
                                        Tree.MCJet_recojet_Dtr_PZ[dtrs0] / 1000.,
                                        Tree.MCJet_recojet_Dtr_E[dtrs0]  / 1000.);

                        double trchi2ndf = Tree.MCJet_recojet_Dtr_TrackChi2[dtrs0] / Tree.MCJet_recojet_Dtr_TrackNDF[dtrs0];

                        if (std::abs(Tree.MCJet_recojet_Dtr_ID[dtrs0]) != HF_pdgcode && 
                            !apply_particle_cuts(dtr0.P(),
                                                 dtr0.Pt(),
                                                 trchi2ndf,
                                                 Tree.MCJet_recojet_Dtr_ProbNNghost[dtrs0],
                                                 dtr0.Rapidity()))
                                continue;

                        for (int dtrs1 = dtrs0 + 1; dtrs1 < Tree.MCJet_recojet_nrecodtrs; dtrs1++) {
                                if (std::abs(Tree.MCJet_recojet_Dtr_ID[dtrs1]) < 100)
                                        continue;

                                dtr1.SetPxPyPzE(Tree.MCJet_recojet_Dtr_PX[dtrs1] / 1000.,
                                                Tree.MCJet_recojet_Dtr_PY[dtrs1] / 1000.,
                                                Tree.MCJet_recojet_Dtr_PZ[dtrs1] / 1000.,
                                                Tree.MCJet_recojet_Dtr_E[dtrs1]  / 1000.);
                        
                                double trchi2ndf_1 = Tree.MCJet_recojet_Dtr_TrackChi2[dtrs1] / Tree.MCJet_recojet_Dtr_TrackNDF[dtrs1];
                        
                                if (std::abs(Tree.MCJet_recojet_Dtr_ID[dtrs1]) != HF_pdgcode && 
                                    !apply_particle_cuts(dtr1.P(),
                                                         dtr1.Pt(),
                                                         trchi2ndf_1,
                                                         Tree.MCJet_recojet_Dtr_ProbNNghost[dtrs1],
                                                         dtr1.Rapidity()))
                                        continue;

                                ntuple_reco_pair->Fill(meas_HFjet.Pt(), dtr0.DeltaR(dtr1, true), (dtr0.Pt() * dtr1.Pt()) / std::pow(meas_HFjet.Pt(), 2),WTA_reco_dist);
                        }
                }

                // Calculate jet substructure kinematical variables
                TVector3 HF_meson = HFmeson.Vect();
                TVector3 HF_jet   = HFjet.Vect();
                TVector3 meas_HF_meson = meas_HFmeson.Vect();
                TVector3 meas_HF_jet   = meas_HFjet.Vect();
                
                if (!hasbbbar) {
                        truth_z_b  = (HF_meson.Dot(HF_jet) ) / (HF_jet.Mag2() );
                        truth_jt_b = (HF_jet.Cross(HF_meson).Mag()) / (HF_jet.Mag());
                        truth_r_b  = static_cast < TLorentzVector > (HFmeson).DeltaR(HFjet, true);
                        
                        h1_z_truth_b->Fill(truth_z_b);
                }

                meas_z  = (meas_HF_meson.Dot(meas_HF_jet) ) / (meas_HF_jet.Mag2() );
                meas_jt = (meas_HF_jet.Cross(meas_HF_meson).Mag()) / (meas_HF_jet.Mag());
                meas_r  = static_cast < TLorentzVector > (meas_HFmeson).DeltaR(meas_HFjet, true);

                // Below: Not Needed --- Just for curiosity //
                if (hasRecoHF) {
                        meas_inJet_z  = (meas_HF_meson.Dot(meas_HF_jet) ) / (meas_HF_jet.Mag2() );
                        meas_inJet_jt = (meas_HF_jet.Cross(meas_HF_meson).Mag()) / (meas_HF_jet.Mag());
                        meas_inJet_r  = static_cast < TLorentzVector > (meas_HFmeson).DeltaR(meas_HFjet, true);
                }
                
                if (hasbbbar) {
                        truth_zg  = (HF_meson.Dot(HF_jet) ) / (HF_jet.Mag2() );
                        truth_jtg = (HF_jet.Cross(HF_meson).Mag()) / (HF_jet.Mag());
                        truth_rg  = static_cast < TLorentzVector > (HFmeson).DeltaR(HFjet, true);
                        
                        h1_z_truth_g->Fill(truth_zg);
                }

                jet_pt_recotruthratio = meas_HFjet.Pt()/HFjet.Pt();
                HF_pt_recotruthratio  = meas_HFmeson.Pt()/HFmeson.Pt();

                jet_pt  = HFjet.Pt();
                jet_eta = HFjet.Eta();
                jet_rap = HFjet.Rapidity();
                
                jet_px = HFjet.Px();
                jet_py = HFjet.Py();
                jet_pz = HFjet.Pz();
                jet_e  = HFjet.E();

                // mum_px = mum.Px();
                // mum_py = mum.Py();
                // mum_pz = mum.Pz();
                // mum_e  = mum.E();
                
                // mup_px = mup.Px();
                // mup_py = mup.Py();
                // mup_pz = mup.Pz();
                // mup_e  = mup.E();
                
                // K_px = Kmeson.Px();
                // K_py = Kmeson.Py();
                // K_pz = Kmeson.Pz();
                // K_e  = Kmeson.E();

                HF_px = HFmeson.Px();
                HF_py = HFmeson.Py();
                HF_pz = HFmeson.Pz();
                HF_e  = HFmeson.E();
                HF_pt = HFmeson.Pt();
                
                truth_z  = HF_meson.Dot(HF_jet) / (HF_jet.Mag2());
                truth_jt = (HF_jet.Cross(HF_meson).Mag()) / (HF_jet.Mag());
                truth_r  = static_cast < TLorentzVector > (HFmeson).DeltaR(HFjet, true);                
                
                meas_jet_pt  = meas_HFjet.Pt();
                meas_jet_eta = meas_HFjet.Eta();
                meas_jet_rap = meas_HFjet.Rapidity();

                meas_jet_px = meas_HFjet.Px();
                meas_jet_py = meas_HFjet.Py();
                meas_jet_pz = meas_HFjet.Pz();
                meas_jet_e  = meas_HFjet.E();

                meas_HF_pt = meas_HFmeson.Pt();
                meas_HF_px = meas_HFmeson.Px();
                meas_HF_py = meas_HFmeson.Py();
                meas_HF_pz = meas_HFmeson.Pz();
                meas_HF_e  = meas_HFmeson.E();

                // meas_mum_px = meas_mum.Px();
                // meas_mum_py = meas_mum.Py();
                // meas_mum_pz = meas_mum.Pz();
                // meas_mum_e  = meas_mum.E();
                
                // meas_mup_px = meas_mup.Px();
                // meas_mup_py = meas_mup.Py();
                // meas_mup_pz = meas_mup.Pz();
                // meas_mup_e  = meas_mup.E();
                
                // meas_K_px = meas_Kmeson.Px();
                // meas_K_py = meas_Kmeson.Py();
                // meas_K_pz = meas_Kmeson.Pz();
                // meas_K_e  = meas_Kmeson.E();

                ndtrs_mc = Tree.MCJet_Dtr_nmcdtrs;
                ndtrs_mcreco = Tree.MCJet_recojet_nrecodtrs;
                
                nTracks = Tree.nTracks; 

                last_eventNum = Tree.eventNumber;
                eventNumber   = Tree.eventNumber;

                ntuple_reco_jet->Fill(meas_HFjet.Pt(),WTA_reco_dist);
                ntuple_true_jet->Fill(HFjet.Pt(),WTA_true_dist);

                BTree->Fill();
        }

        cout << "Total number of events processed = " << events << endl;
        cout << "NumRecoHF = " << NumRecoHF << endl;
        cout << "NumHFHads = " << NumHFHads << endl;
        cout << "Num Single Bjets    = " << Num_Single_tr_Bjet << endl;
        cout << "Num Photon + B jets = " << Num_Photon_tr_Bjet << endl;

        cout << "Events blocked: " << endl;
        cout << "Dinjet = " << cut_HFdR << endl;
        
        f.Write();
        f.Close();

        benchmark->Show("MakeVarTreeMC");
}
