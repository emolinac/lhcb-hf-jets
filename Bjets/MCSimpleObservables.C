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
#include "../include/utils.cpp"
#include "../include/utils.h"
#include "../include/utils-visual.cpp"
#include "../include/utils-visual.h"

//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/contrib/SoftDrop.hh"
//#include "../LundGen.hh"
//using namespace fastjet;
using namespace std;

void MCSimpleObservables(int NumEvts = -1)
{
        // This should be the file from MCMakeVarTree
        std::string extension = "simpleobservable_mc_";
        
        std::string extension_RootFilesMC = output_folder + "bjets-mc/";
        std::cout<<"Starting"<<std::endl;


        TFile fread((output_folder + "ntuple_bjets_mc.root").c_str(), "READ");
        TTree *BTree = (TTree *)fread.Get("BTree");
        
        std::cout<<"Files opened"<<std::endl;

        if (NumEvts > BTree->GetEntries())
                NumEvts = BTree->GetEntries();
        
        if (NumEvts == -1)
                NumEvts = BTree->GetEntries();
        
        cout<<BTree->GetEntries()<<endl;

        TFile f((output_folder + "bjets_simpleobservable_mc.root").c_str(), "RECREATE");

        // 2D Truth-Reco Correspondence (219 - 224)
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

        // EEC RELATED PLOTS
        TH3D* h_npair_mc          = new TH3D("h_npair_mc"         , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h_npair_eqcharge_mc = new TH3D("h_npair_eqcharge_mc", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h_npair_opcharge_mc = new TH3D("h_npair_opcharge_mc", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);

        TH2D* h_eec_mc_rl_jetpt          = new TH2D("h_eec_mc_rl_jetpt"         , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h_eec_eqcharge_mc_rl_jetpt = new TH2D("h_eec_eqcharge_mc_rl_jetpt", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h_eec_opcharge_mc_rl_jetpt = new TH2D("h_eec_opcharge_mc_rl_jetpt", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);

        TH1F* hmc_eec[ptbinsize]; 
        TH1F* hmc_eqcheec[ptbinsize]; 
        TH1F* hmc_neqcheec[ptbinsize]; 
        TH1F* hmc_npair[ptbinsize]; 

        std::cout<<"Histos defined"<<std::endl;

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

        vector<float> *pair_rl = 0, *pair_weight = 0, *pair_chargeprod = 0;

        double WTA_true_dist;
        double WTA_reco_dist;

        std::cout<<"Setting branches"<<std::endl;
        
        BTree->SetBranchAddress("pair_rl"        , &pair_rl);
        BTree->SetBranchAddress("pair_weight"    , &pair_weight);
        BTree->SetBranchAddress("pair_chargeprod", &pair_chargeprod);
        
        BTree->SetBranchAddress("jet_pt", &jet_pt);
        BTree->SetBranchAddress("jet_eta", &jet_eta);
        BTree->SetBranchAddress("jet_rap", &jet_rap);
        BTree->SetBranchAddress("jet_px", &jet_px);
        BTree->SetBranchAddress("jet_py", &jet_py);
        BTree->SetBranchAddress("jet_pz", &jet_pz);
        BTree->SetBranchAddress("jet_e" , &jet_e);

        BTree->SetBranchAddress("HF_px", &HF_px);
        BTree->SetBranchAddress("HF_py", &HF_py);
        BTree->SetBranchAddress("HF_pz", &HF_pz);
        BTree->SetBranchAddress("HF_e", &HF_e);
        BTree->SetBranchAddress("HF_pt", &HF_pt);

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

        BTree->SetBranchAddress("meas_mum_px", &meas_mum_px);
        BTree->SetBranchAddress("meas_mum_py", &meas_mum_py);
        BTree->SetBranchAddress("meas_mum_pz", &meas_mum_pz);
        BTree->SetBranchAddress("meas_mum_e", &meas_mum_e);
        BTree->SetBranchAddress("meas_mup_px", &meas_mup_px);
        BTree->SetBranchAddress("meas_mup_py", &meas_mup_py);
        BTree->SetBranchAddress("meas_mup_pz", &meas_mup_pz);
        BTree->SetBranchAddress("meas_mup_e", &meas_mup_e);
        BTree->SetBranchAddress("meas_K_px", &meas_K_px);
        BTree->SetBranchAddress("meas_K_py", &meas_K_py);
        BTree->SetBranchAddress("meas_K_pz", &meas_K_pz);
        BTree->SetBranchAddress("meas_K_e", &meas_K_e);

        BTree->SetBranchAddress("jet_pt_recotruthratio", &jet_pt_recotruthratio);
        BTree->SetBranchAddress("HF_pt_recotruthratio", &HF_pt_recotruthratio);
        BTree->SetBranchAddress("GluonTag", &GluonTag);
        BTree->SetBranchAddress("nTracks", &nTracks);
        BTree->SetBranchAddress("hasbbbar", &hasbbbar);
        BTree->SetBranchAddress("NumHFHads", &NumHFHads);
        BTree->SetBranchAddress("NumDtrRecoHF", &NumDtrRecoHF);
        BTree->SetBranchAddress("hasRecoHF", &hasRecoHF);
        BTree->SetBranchAddress("isSingle_tr_Bjet", &isSingle_tr_Bjet);
        BTree->SetBranchAddress("isPhoton_tr_Bjet", &isPhoton_tr_Bjet);
        BTree->SetBranchAddress("ndtrs_mc",     &ndtrs_mc);
        BTree->SetBranchAddress("ndtrs_mcreco", &ndtrs_mcreco);
        
        BTree->SetBranchAddress("WTA_true_dist", &WTA_true_dist);
        BTree->SetBranchAddress("WTA_reco_dist", &WTA_reco_dist);

        int eventNum;
        int NumJets = 0;
        int NumJets_zdR = 0;
        int NumJets_ktdR = 0;
        int NumBJets = 0;

        TLorentzVector HFmeson, HFjet, mum, mup, K, Jpsi;
        TLorentzVector meas_HFmeson, meas_HFjet, meas_mum, meas_mup, meas_K, meas_Jpsi;
                
        cout << "Requested # of events " << NumEvts << endl;
        
        // Event loop
        unsigned long long last_eventNum = 0;
        double last_jetpT = 0.;

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
                
                h2_jetpteta->Fill(HFjet.Pt(), HFjet.Eta());  
                h2_jetptp->Fill(HFjet.Pt(), HFjet.P());            

                if (!pair_rl->empty()) {
                        ULong_t vector_size = pair_rl->size();

                        float *rl_info         = pair_rl->data();
                        float *weight_info     = pair_weight->data();
                        float *chargeprod_info = pair_chargeprod->data();
                        
                        for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                h_npair_mc->Fill(rl_info[vector_index], HFjet.Pt(), weight_info[vector_index]);

                                if (chargeprod_info[vector_index] > 0)
                                        h_npair_eqcharge_mc->Fill(rl_info[vector_index], HFjet.Pt(), weight_info[vector_index]);
                                else if (chargeprod_info[vector_index] < 0)
                                        h_npair_opcharge_mc->Fill(rl_info[vector_index], HFjet.Pt(), weight_info[vector_index]);
                        }
                }

                
                
                bool hasEmission_ktdR = false;
                bool hasEmission_zdR = false;
        } // End of BTree entry loop

        // Operate the EEC related plots
        
        apply_unfolded_weights(h_npair_mc, h_eec_mc_rl_jetpt);
        apply_unfolded_weights(h_npair_eqcharge_mc, h_eec_eqcharge_mc_rl_jetpt);
        apply_unfolded_weights(h_npair_opcharge_mc, h_eec_opcharge_mc_rl_jetpt);

        TH2D* h_npair_mc_rl_jetpt = (TH2D*) h_npair_mc->Project3D("yx");

        for (int bin = 0 ; bin < ptbinsize ; bin++) {
                int nominal_jet_pt_bin = bin + 3;

                // Pseudodata operations
                hmc_eec[bin]      = new TH1F(Form("hmc_eec%i",bin)     , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                hmc_eqcheec[bin]  = new TH1F(Form("hmc_eqcheec%i",bin) , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hmc_neqcheec[bin] = new TH1F(Form("hmc_neqcheec%i",bin), "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
                hmc_npair[bin]    = new TH1F(Form("hmc_npair%i",bin)   , "", nbin_rl_nominal_unfolding,unfolding_rl_nominal_binning);
                
                set_histogram_style(hmc_eec[bin]     , corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hmc_eqcheec[bin] , corr_marker_color_jet_pt[bin], std_line_width-1, std_marker_style_jet_pt[bin] , std_marker_size+1);
                set_histogram_style(hmc_neqcheec[bin], corr_marker_color_jet_pt[bin], std_line_width-1, corr_marker_style_jet_pt[bin], std_marker_size+1);
                set_histogram_style(hmc_npair[bin]   , corr_marker_color_jet_pt[bin], std_line_width, corr_marker_style_jet_pt[bin], std_marker_size+1);

                project_nominal_phase_space(h_eec_mc_rl_jetpt          , hmc_eec[bin]     , bin + 1);
                project_nominal_phase_space(h_eec_eqcharge_mc_rl_jetpt , hmc_eqcheec[bin] , bin + 1);
                project_nominal_phase_space(h_eec_opcharge_mc_rl_jetpt , hmc_neqcheec[bin], bin + 1);
                project_nominal_phase_space(h_npair_mc_rl_jetpt        , hmc_npair[bin]   , bin + 1);

                hmc_eqcheec[bin]->Divide(hmc_eec[bin]);
                hmc_neqcheec[bin]->Divide(hmc_eec[bin]);

                hmc_eec[bin]->Scale(1./h1_jet_pt->GetBinContent(bin + 1),"width");
                hmc_npair[bin]->Scale(1./h1_jet_pt->GetBinContent(bin + 1),"width");
        }
        
        //
        f.Write();
        f.Close();

        cout << "Num of True B jets = " << NumBJets << endl;
}
