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

using namespace std;

// Description: Processes output of MCMakeVarTree.C to obtain the denominator for the efficiency calculation (efficiency = # of matched B-jets / # of truth B-jets)

void MCSimpleEff(int NumEvts = -1, int dataset = 91599, int flavor = 5,
                 bool chargedJetCut_user = false,
                 bool SubtractGS = false,
                 bool onlysim9 = false) 
{    
        // This should be the file from MCMakeVarTree
        std::string extension = "efficiency_";
        
        std::string extension_RootFilesMC = output_folder + "bjets-mc/";

        TFile fread((output_folder + "ntuple_bjets_mc.root").c_str(), "READ");
        TTree *BTree = (TTree *)fread.Get("BTree");
        
        if (NumEvts > BTree->GetEntries())
                NumEvts = BTree->GetEntries();
        
        if (NumEvts == -1)
                NumEvts = BTree->GetEntries();
        
        cout<<BTree->GetEntries()<<endl;

        TFile f((output_folder + "bjets_efficiencies.root").c_str(), "RECREATE");

        TH1D *h1_denom_efficiency_jetpt = new TH1D("denom_efficiency_jetpt", "", ptbinsize, pt_binedges); //stay 
        TH1D *h1_denom_efficiency_HFpt  = new TH1D("denom_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
        
        TH2D *h2_denom_efficiency_jetpteta    = new TH2D("denom_efficiency_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
        TH2D *h2_denom_efficiency_HFpteta     = new TH2D("denom_efficiency_HFpteta"     , "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges); 
        TH2D *h2_denom_efficiency_HFptjetpt   = new TH2D("denom_efficiency_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
        TH2D *h2_denom_efficiency_HFptnTracks = new TH2D("denom_efficiency_HFptnTracks","",ptHFbinsize, ptHF_binedges, nTracksbinsize, nTrack_binedges);
        
        TH3D *h3_denom_efficiency_HFptetajetpt = new TH3D("denom_efficiency_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, customptbinsize, custompt_binedges);

        // EEC-related 
        TH3D *h3_denom_efficiency_rl_jetpt_weight       = new TH3D("denom_efficiency_rl_jetpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_denom_efficiency_rl_jetpt_weight_eqch  = new TH3D("denom_efficiency_rl_jetpt_weight_eqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_denom_efficiency_rl_jetpt_weight_neqch = new TH3D("denom_efficiency_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH2D *h2_SVTag_eff_denom   = new TH2D("h2_SVTag_eff_denom","", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
        TH2D *h2_SVTag_eff_denom_z = new TH2D("h2_SVTag_eff_denom_z","", zbinsize, z_binedges, customptbinsize, custompt_binedges); 
        
        // Event loop
        unsigned long long last_eventNum = 0;
        double last_jetpT = 0.;
        int event_counter = 0;

        // Truth Variables (from Truth Tree)
        double jet_pt, jet_rap, meas_jet_pt, meas_jet_rap;
        double jet_px, jet_py, jet_pz, jet_e;
        double HF_px, HF_py, HF_pz, HF_e;
        double HF_pt;
        double mup_px, mup_py, mup_pz, mup_e;
        double mum_px, mum_py, mum_pz, mum_e;
        
        double jet_eta, meas_jet_eta;
        double K_px, K_py, K_pz, K_e;

        // Reco Variables (from Truth Tree)
        double meas_jet_px, meas_jet_py, meas_jet_pz, meas_jet_e;
        double meas_HF_px, meas_HF_py, meas_HF_pz, meas_HF_e;
        
        double meas_HF_pt;

        double jet_pt_recotruthratio, HF_pt_recotruthratio;

        int nsplits, ndtrs, NumHFHads;
        int meas_nsplits, meas_ndtrs;
        int eventNumber;
        int GluonTag, nTracks;
        //    bool GluonTag, Hasbbbar;
        
        int NumDtrRecoHF;
        bool hasRecoHF, Hasbbbar;

        double WTA_true_dist;

        vector<float> *pair_rl = 0, *pair_weight = 0, *pair_chargeprod = 0;

        BTree->SetBranchAddress("WTA_true_dist", &WTA_true_dist);

        BTree->SetBranchAddress("pair_rl"        , &pair_rl);
        BTree->SetBranchAddress("pair_weight"    , &pair_weight);
        BTree->SetBranchAddress("pair_chargeprod", &pair_chargeprod);
        
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

        BTree->SetBranchAddress("GluonTag", &GluonTag);
        BTree->SetBranchAddress("nTracks", &nTracks);
        
        BTree->SetBranchAddress("NumDtrRecoHF", &NumDtrRecoHF);
        BTree->SetBranchAddress("jet_pt_recotruthratio", &jet_pt_recotruthratio);
        BTree->SetBranchAddress("HF_pt_recotruthratio", &HF_pt_recotruthratio);


        int eventNum;
        int NumBjets = 0;
        int NumRecoBjets = 0;

        TLorentzVector HFmeson, meas_HFmeson, mup, mum;

        cout << "Requested # of events " << NumEvts << endl;
        for (int ev = 0; ev < NumEvts; ev++) 
        {
                BTree->GetEntry(ev);

                if (ev%10000 == 0) {
                        double percentage = 100.*ev/NumEvts;
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }

                if (WTA_true_dist > 0.005)
                        continue;
                
                HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
                meas_HFmeson.SetPxPyPzE(meas_HF_px, meas_HF_py, meas_HF_pz, meas_HF_e);
                mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
                mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);

                // new
                //    if(NumDtrRecoHF > 1) continue;
                
                bool rap_cond      = (jet_rap > etaMin && jet_rap < etaMax);
                bool pt_cond       = (jet_pt > pTLow);
                bool meas_rap_cond = (meas_jet_rap > etaMin && meas_jet_rap < etaMax);
                bool meas_pt_cond  = (meas_jet_pt > pTLow);

                //        if (NumHFHads > 1) {
                //            continue;
                //        }
                //
                //        if(Hasbbbar) continue;
                
                if (SubtractGS && Hasbbbar)
                        continue;
                
                NumBjets++;
                
                if (pt_cond) 
                        h2_denom_efficiency_jetpteta->Fill(jet_pt, jet_eta);

                if (pt_cond && rap_cond) {
                        h1_denom_efficiency_HFpt->Fill(HF_pt);
                        
                        h2_denom_efficiency_HFpteta->Fill(HFmeson.Pt(), HFmeson.Rapidity());
                        h3_denom_efficiency_HFptetajetpt->Fill(HFmeson.Pt(), HFmeson.Rapidity(), jet_pt);
                        h1_denom_efficiency_jetpt->Fill(jet_pt);
                        h2_denom_efficiency_HFptjetpt->Fill(HFmeson.Pt(), jet_pt);
                        
                        h2_denom_efficiency_HFptnTracks->Fill(HFmeson.Pt(), nTracks);
                        
                        if (!pair_rl->empty()) {
                                ULong_t vector_size = pair_rl->size();

                                float *rl_info         = pair_rl->data();
                                float *weight_info     = pair_weight->data();
                                float *chargeprod_info = pair_chargeprod->data();
                                
                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        h3_denom_efficiency_rl_jetpt_weight->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);

                                        if (chargeprod_info[vector_index] > 0)
                                                h3_denom_efficiency_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                        else if (chargeprod_info[vector_index] < 0)
                                                h3_denom_efficiency_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                }
                        }
                }

                event_counter++;
        }
        
        cout << event_counter << " events processed" << endl;

        f.Write();
        f.Close();
        
        cout << "Num of True B jets = " << NumBjets << endl;
        cout << "Num of Reco B jets = " << NumRecoBjets << endl;
}
