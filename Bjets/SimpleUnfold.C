#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <TH3.h>
#include <TF1.h>
#include <TLatex.h>
#include <THStack.h>
#include <TChain.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TEfficiency.h>
#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"

using namespace std;

void SimpleUnfold(int NumEvts = -1,
                  bool SubtractGS = false,                  
                  bool DoJESJER = false,
                  bool DoRecSelEff = false,
                  bool DoSignalSys = false,
                  bool DoUnfoldPrior = false,
                  bool sPlotFit = false)
{
        TFile *file_eff = new TFile((output_folder + "bjets_efficiencies.root").c_str(), "READ");
        
        TChain *BTree = new TChain("BTree", "");
        BTree->Add((output_folder + "ntuple_bjets_mcreco.root/BTree").c_str());
        
        if (NumEvts > BTree->GetEntries())
                NumEvts = BTree->GetEntries();
        if (NumEvts == -1)
                NumEvts = BTree->GetEntries();
        
        TH1D *h1_denom_efficiency_HFpt   = (TH1D*) file_eff->Get("denom_efficiency_HFpt");
        TH1D *h1_denom_efficiency_jetpt  = (TH1D*) file_eff->Get("denom_efficiency_jetpt");
        
        TH2D *h2_denom_efficiency_HFptjetpt = (TH2D*) file_eff->Get("denom_efficiency_HFptjetpt");
        TH2D *h2_denom_efficiency_HFpteta   = (TH2D*) file_eff->Get("denom_efficiency_HFpteta");
        TH2D *h2_denom_efficiency_jetpteta  = (TH2D*) file_eff->Get("denom_efficiency_jetpteta");
        
        TH3D *h3_denom_efficiency_HFptetajetpt = (TH3D*) file_eff->Get("denom_efficiency_HFptetajetpt");
        
        TH3D *h3_denom_efficiency_rl_jetpt_weight       = (TH3D*) file_eff->Get("denom_efficiency_rl_jetpt_weight");
        TH3D *h3_denom_efficiency_rl_jetpt_weight_eqch  = (TH3D*) file_eff->Get("denom_efficiency_rl_jetpt_weight_eqch");
        TH3D *h3_denom_efficiency_rl_jetpt_weight_neqch = (TH3D*) file_eff->Get("denom_efficiency_rl_jetpt_weight_neqch");  
        
        //    /////////////////// Mass Fit Parameters /////////////////////////////////
        // TString extension_mass("");
        // if (sPlotFit)
        //         extension_mass = TString("splotfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.4_b" + str_PID + str_L0 + "_91599.root";
        // else
        //         extension_mass = TString("massfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.4_b" + str_PID + str_L0 + "_91599.root";
        
        // if (DoRecSelEff)
        //         extension_mass = "recselsys_" + extension_mass;

        // if (DoSignalSys)
        //         extension_mass = "sys_" + extension_mass;

        // extension_mass = extension_RootFilesData + extension_mass; // EFMC: currently, not in use
        
        TFile f_massfit((output_folder + "mass-fits/results_mass_fit_mcreco.root").c_str(), "READ");
        
        TH1D *h1_MassMin      = (TH1D*) f_massfit.Get("h1_MassMin");
        TH1D *h1_MassMax      = (TH1D*) f_massfit.Get("h1_MassMax");
        TH1D *h1_BkgScale     = (TH1D*) f_massfit.Get("h1_BkgScale");
        TH1D *h1_BkgScaleNear = (TH1D*) f_massfit.Get("h1_BkgScale_forSysNear");
        TH1D *h1_BkgScaleFar  = (TH1D*) f_massfit.Get("h1_BkgScale_forSysFar");
        
        
        if (h1_BkgScaleNear == NULL)
                cout << "NULL NEAR!" << endl;
        if (h1_BkgScaleFar == NULL)
                cout << "NULL FAR!" << endl;

        TFile* file_reco_weights;
        TH3D*  h3_ptzjt_ratio;

        if (DoUnfoldPrior) {
                file_reco_weights = new TFile((output_folder + "MC_DATA_WEIGHTS.root").c_str(), "READ"); 

                h3_ptzjt_ratio = (TH3D *) file_reco_weights->Get("ptzjt_ratio");
        }
        
        TFile *f = TFile::Open((output_folder + "bjets_corrections.root").c_str(), "RECREATE");
        
        //Reco 1D Observables
        TH1D *h1_z  = new TH1D("z", "", zbinsize, z_binedges);
        TH1D *h1_jt = new TH1D("jt", "", jtbinsize, jt_binedges);
        TH1D *h1_r  = new TH1D("r", "", rbinsize, r_binedges);
        
        // 2D Truth-Reco Correspondence (219 - 224)
        TH2D *h2_truthreco_z  = new TH2D("truthreco_z", ";Reco z; Truth z", zbinsize, z_binedges, zbinsize, z_binedges);
        TH2D *h2_truthreco_jt = new TH2D("truthreco_jt", ";Reco jT; Truth jT", jtbinsize, jt_binedges, jtbinsize, jt_binedges);
        TH2D *h2_truthreco_r  = new TH2D("truthreco_r", ";Reco r; Truth r", rbinsize, r_binedges, rbinsize, r_binedges);

        TH2D *h2_truthreco_jetpt = new TH2D("truthreco_jetpt", "", longptbinsize, longpt_binedges, longptbinsize, longpt_binedges);
        
        // 1D Denom Efficiencies and Purities of Observables (237 - 246)
        TH1D *h1_num_purity_jetpt   = new TH1D("num_purity_jetpt", "", ptbinsize, pt_binedges);
        TH1D *h1_denom_purity_jetpt = new TH1D("denom_purity_jetpt", "", ptbinsize, pt_binedges);
        
        TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", ptbinsize, pt_binedges);
        
        TH1D *h1_jetpt       = new TH1D("jetpt"      , "", ptbinsize, pt_binedges);
        TH1D *h1_jetpt_truth = new TH1D("jetpt_truth", "", ptbinsize, pt_binedges);
       
        RooUnfoldResponse *response_jetpt = new RooUnfoldResponse(h1_jetpt, h1_jetpt_truth, "response_jetpt");
        
        // 1D RMs (for visualization purposes)
        TH1D *h1_form_rl     = new TH1D("h1_form_rl"    , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
        TH1D *h1_form_weight = new TH1D("h1_form_weight", "", nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_rl     = new RooUnfoldResponse(h1_form_rl    , h1_form_rl    , "response_rl");
        RooUnfoldResponse *response_weight = new RooUnfoldResponse(h1_form_weight, h1_form_weight, "response_weight");
        
        // 3D RM
        TH3D *h3_meas_rl_jetpt_weight = new TH3D("h3_meas_rl_jetpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight = new TH3D("h3_true_rl_jetpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair = new RooUnfoldResponse(h3_meas_rl_jetpt_weight, h3_true_rl_jetpt_weight, "response_npair");
        
        TH3D *h3_meas_rl_jetpt_weight_eqch = new TH3D("h3_meas_rl_jetpt_weight_eqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight_eqch = new TH3D("h3_true_rl_jetpt_weight_eqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair_eqch = new RooUnfoldResponse(h3_meas_rl_jetpt_weight_eqch, h3_true_rl_jetpt_weight_eqch, "response_npair_eqch");
        
        TH3D *h3_meas_rl_jetpt_weight_neqch = new TH3D("h3_meas_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight_neqch = new TH3D("h3_true_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair_neqch = new RooUnfoldResponse(h3_meas_rl_jetpt_weight_neqch, h3_true_rl_jetpt_weight_neqch, "response_npair_neqch");
        
        // 3D Efficiencies
        TH3D* h3_num_efficiency_rl_jetpt_weight = new TH3D("h3_num_efficiency_rl_jetpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight     = new TH3D("h3_efficiency_rl_jetpt_weight"    , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_efficiency_rl_jetpt_weight_eqch = new TH3D("h3_num_efficiency_rl_jetpt_weight_eqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight_eqch     = new TH3D("h3_efficiency_rl_jetpt_weight_eqch"    , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_efficiency_rl_jetpt_weight_neqch = new TH3D("h3_num_efficiency_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight_neqch     = new TH3D("h3_efficiency_rl_jetpt_weight_neqch"    , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        // 3D Purities
        TH3D* h3_num_purity_rl_jetpt_weight     = new TH3D("h3_num_purity_rl_jetpt_weight"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight   = new TH3D("h3_denom_purity_rl_jetpt_weight" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight         = new TH3D("h3_purity_rl_jetpt_weight"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_purity_rl_jetpt_weight_eqch   = new TH3D("h3_num_purity_rl_jetpt_weight_eqch"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight_eqch = new TH3D("h3_denom_purity_rl_jetpt_weight_eqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight_eqch       = new TH3D("h3_purity_rl_jetpt_weight_eqch"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_purity_rl_jetpt_weight_neqch   = new TH3D("h3_num_purity_rl_jetpt_weight_neqch"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight_neqch = new TH3D("h3_denom_purity_rl_jetpt_weight_neqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight_neqch       = new TH3D("h3_purity_rl_jetpt_weight_neqch"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        // Prior-weighted observables
        TH3D *h3_rl_jetpt_weight_weighted = new TH3D("h3_rl_jetpt_weight_weighted", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        // Event loop
        unsigned long long last_eventNum = 0;
        float last_jetpT = 0.;
        int event_counter = 0;
        
        cout << BTree->GetEntries() << endl;
        
        vector<float> *zs(0), *tr_zs(0);
        vector<float> *sv_ipchi2(0);
        
        float jet_pt, HF_pt, jet_rap, tr_jet_pt, tr_jet_rap, tr_HF_pt, bmass_dtf;
        float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
        float jet_px, jet_py, jet_pz, jet_e;
        float HF_px, HF_py, HF_pz, HF_e;
        bool isTrueBjet;
        int NumBHads_tr;
        int eventNumber, nSV;
        bool Hasbbbar, TOS;
        float Jpsi_CHI2NDOF, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_IPCHI2, Jpsi_BPVDLS;;
        float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
        float sv_mass, sv_chi2, sv_ntrks, sv_cosine;
        int SVTag;
        float K_PIDK, chi2ndf_dtf, tau_dtf;
        
        vector<float> *pair_rl = 0, *pair_weight = 0, *pair_chargeprod = 0;
        vector<float> *truthmatched_pair_rl = 0, *truthmatched_pair_weight = 0;
        
        BTree->SetBranchAddress("pair_rl"        , &pair_rl);
        BTree->SetBranchAddress("pair_weight"    , &pair_weight);
        BTree->SetBranchAddress("pair_chargeprod", &pair_chargeprod);
        
        BTree->SetBranchAddress("truthmatched_pair_rl"    , &truthmatched_pair_rl);
        BTree->SetBranchAddress("truthmatched_pair_weight", &truthmatched_pair_weight);
        
        BTree->SetBranchAddress("HF_px", &HF_px);
        BTree->SetBranchAddress("HF_py", &HF_py);
        BTree->SetBranchAddress("HF_pz", &HF_pz);
        BTree->SetBranchAddress("HF_e", &HF_e);
        BTree->SetBranchAddress("HF_pt", &HF_pt);
        
        BTree->SetBranchAddress("Bu_IPCHI2", &Bu_IPCHI2);
        BTree->SetBranchAddress("Bu_CHI2", &Bu_CHI2);
        BTree->SetBranchAddress("Bu_CHI2NDOF", &Bu_CHI2NDOF);
        
        BTree->SetBranchAddress("Jpsi_FDCHI2", &Jpsi_FDCHI2);
        BTree->SetBranchAddress("Jpsi_CHI2", &Jpsi_CHI2);
        BTree->SetBranchAddress("Jpsi_CHI2NDOF", &Jpsi_CHI2NDOF);
        BTree->SetBranchAddress("Jpsi_BPVDLS", &Jpsi_BPVDLS);    
        
        BTree->SetBranchAddress("jet_pt" , &jet_pt);
        BTree->SetBranchAddress("jet_rap", &jet_rap);
        BTree->SetBranchAddress("jet_px" , &jet_px);
        BTree->SetBranchAddress("jet_py" , &jet_py);
        BTree->SetBranchAddress("jet_pz" , &jet_pz);
        BTree->SetBranchAddress("jet_e"  , &jet_e);
        
        BTree->SetBranchAddress("tr_jet_px", &tr_jet_px);
        BTree->SetBranchAddress("tr_jet_py", &tr_jet_py);
        BTree->SetBranchAddress("tr_jet_pz", &tr_jet_pz);
        BTree->SetBranchAddress("tr_jet_e", &tr_jet_e);
        
        BTree->SetBranchAddress("tr_jet_pt", &tr_jet_pt);
        BTree->SetBranchAddress("tr_HF_pt", &tr_HF_pt);
        BTree->SetBranchAddress("tr_jet_rap", &tr_jet_rap);
        BTree->SetBranchAddress("isTrueBjet", &isTrueBjet);
        BTree->SetBranchAddress("NumBHads_tr", &NumBHads_tr);
        BTree->SetBranchAddress("bmass_dtf", &bmass_dtf);
        BTree->SetBranchAddress("eventNumber", &eventNumber);
        BTree->SetBranchAddress("chi2ndf_dtf", &chi2ndf_dtf);
        BTree->SetBranchAddress("tau_dtf", &tau_dtf);
        
        BTree->SetBranchAddress("Hasbbbar", &Hasbbbar);
        BTree->SetBranchAddress("nSV", &nSV);
        BTree->SetBranchAddress("sv_mass", &sv_mass);
        BTree->SetBranchAddress("sv_chi2", &sv_chi2);
        BTree->SetBranchAddress("sv_cosine", &sv_cosine);
        BTree->SetBranchAddress("sv_ipchi2", &sv_ipchi2);
        BTree->SetBranchAddress("sv_ntrks", &sv_ntrks);
        BTree->SetBranchAddress("SVTag", &SVTag);
        BTree->SetBranchAddress("K_PIDK", &K_PIDK);
        
        BTree->SetBranchAddress("TOS", &TOS);
        
        cout << "Requested # of events" << NumEvts << endl;
        
        
        int eventNum;
        int NumBjets = 0;
        int NumTrueBjets = 0;

        TLorentzVector HFmeson, HFjet, tr_HFjet;
        
        for (int ev = 0; ev < NumEvts; ev++) {
                BTree->GetEntry(ev);
                
                if (ev%10000 == 0) {
                        double percentage = 100.*ev/NumEvts;
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }
                
                HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
                HFjet.SetPxPyPzE(jet_px, jet_py, jet_pz, jet_e);
                tr_HFjet.SetPxPyPzE(tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e);
                
                float eff_weight = 1.0;
                float prior_weight = 1.0;
                
                // if (DoUnfoldPrior)
                //         prior_weight = h3_ptzjt_ratio->GetBinContent(h3_ptzjt_ratio->GetXaxis()->FindBin(z),
                //                                                      h3_ptzjt_ratio->GetYaxis()->FindBin(jt),
                //                                                      h3_ptzjt_ratio->GetZaxis()->FindBin(jet_pt));
                
                // Reweight by inverse of (number of smearing trials)
                if (DoJESJER)
                        eff_weight = 1. / n_smearing_iter;
                
                float bkg_weight = h1_BkgScale != NULL ? h1_BkgScale->GetBinContent(h1_BkgScale->FindBin(HFmeson.Pt())) : 1.0;
                float MassHigh   = h1_MassMax != NULL ? h1_MassMax->GetBinContent(h1_MassMax->FindBin(HFmeson.Pt())) : 5.31;
                float MassLow    = h1_MassMin != NULL ? h1_MassMin->GetBinContent(h1_MassMin->FindBin(HFmeson.Pt())) : 5.24;
                
                
                bool mass_cond   = (bmass_dtf > MassLow && bmass_dtf < MassHigh);
                bool PID_cond    = (K_PIDK > 0);
                bool rap_cond    = (jet_rap > etaMin && jet_rap < etaMax);
                bool pt_cond     = (jet_pt > pTLow);
                bool tr_rap_cond = (tr_jet_rap > etaMin && tr_jet_rap < etaMax);
                bool tr_pt_cond  = (tr_jet_pt > pTLow);
                bool SV_cond     = (nSV > 0) && mass_cond && sv_mass > 0.4;
                bool gluon_cond  = mass_cond && Hasbbbar;
                
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
                
                if (!TOS)
                        continue;

                if (!mass_cond)
                        continue;
                
                if (SubtractGS && Hasbbbar)
                        continue;
                
                float HF_rap = HFmeson.Rapidity();
        
                bool num_cond_eff = tr_pt_cond && tr_rap_cond && isTrueBjet && pt_cond && rap_cond;
                bool num_cond_pur = tr_pt_cond && tr_rap_cond && isTrueBjet && pt_cond && rap_cond;
                bool denom_cond   = pt_cond && rap_cond;
                
                if (num_cond_eff) {
                        h1_num_efficiency_jetpt->Fill(tr_jet_pt, eff_weight);

                        if (!truthmatched_pair_rl->empty()) {
                                ULong_t vector_size = truthmatched_pair_rl->size();

                                // Note: this vectors will have the same size always, dont worry
                                float *rl_info         = truthmatched_pair_rl->data();
                                float *weight_info     = truthmatched_pair_weight->data();
                                float *chargeprod_info = pair_chargeprod->data();
                                
                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        if (rl_info[vector_index] == -999)
                                                continue;

                                        h3_num_efficiency_rl_jetpt_weight->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        if (chargeprod_info[vector_index] > 0)
                                                h3_num_efficiency_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                        else if (chargeprod_info[vector_index] < 0)
                                                h3_num_efficiency_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                }
                        }
                }
                
                if (denom_cond) {
                        h1_denom_purity_jetpt->Fill(jet_pt, eff_weight);
                        
                        if (!truthmatched_pair_rl->empty()) {
                                ULong_t vector_size = truthmatched_pair_rl->size();

                                // Note: this vectors will have the same size always, dont worry
                                float *truthmatched_rl_info = truthmatched_pair_rl->data();

                                float *rl_info         = pair_rl->data();
                                float *weight_info     = pair_weight->data();
                                float *chargeprod_info = pair_chargeprod->data();
                                
                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        if (truthmatched_rl_info[vector_index] == -999)
                                                continue;
                                        
                                        h3_denom_purity_rl_jetpt_weight->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        if (chargeprod_info[vector_index] > 0)
                                                h3_denom_purity_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                        else if (chargeprod_info[vector_index] < 0)
                                                h3_denom_purity_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                }
                        }
                }    

                if (num_cond_pur) {
                        // Note: Two birds in one shot baby
                        h1_num_purity_jetpt->Fill(jet_pt, eff_weight);

                        response_jetpt->Fill(jet_pt, tr_jet_pt, prior_weight);
                        
                        if (!truthmatched_pair_rl->empty()) {
                                ULong_t vector_size = pair_rl->size();

                                // Note: this vectors will have the same size always, dont worry
                                float *rl_info         = pair_rl->data();
                                float *weight_info     = pair_weight->data();
                                float *chargeprod_info = pair_chargeprod->data();

                                float *truthmatched_rl_info     = truthmatched_pair_rl->data();
                                float *truthmatched_weight_info = truthmatched_pair_weight->data();
                                
                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        h3_num_purity_rl_jetpt_weight->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        if (truthmatched_rl_info[vector_index] != -999) {
                                                response_npair->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index]);
                                        
                                                response_rl->Fill(rl_info[vector_index], truthmatched_rl_info[vector_index]);

                                                response_weight->Fill(weight_info[vector_index], truthmatched_weight_info[vector_index]);
                                        }

                                        if (chargeprod_info[vector_index] > 0) {
                                                h3_num_purity_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                                if (truthmatched_rl_info[vector_index] != -999)
                                                        response_npair_eqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index]);

                                        } else if (chargeprod_info[vector_index] < 0) {
                                                h3_num_purity_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                                if (truthmatched_rl_info[vector_index] != -999)
                                                        response_npair_neqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index]);
                                        }
                                }
                        }

                        NumTrueBjets++;              
                }
                
                // if (!pt_cond || !tr_pt_cond || !isTrueBjet || !rap_cond || !tr_rap_cond )
                //         continue;

                // response_ptz->Fill( z,  jet_pt, tr_z, tr_jet_pt, prior_weight);        
                // response_ptjt->Fill( jt, jet_pt, tr_jt, tr_jet_pt, prior_weight);
                // response_ptr->Fill( r,  jet_pt, tr_r, tr_jet_pt, prior_weight);
                
                // response_ptzjt->Fill( z, jt, jet_pt, tr_z, tr_jt, tr_jet_pt, prior_weight);
                // response_ptzr->Fill( z, r, jet_pt, tr_z, tr_r, tr_jet_pt, prior_weight);
                // response_ptjtr->Fill( jt,  r, jet_pt, tr_jt, tr_r, tr_jet_pt, prior_weight);  
                
                // h2_truthreco_jetpt->Fill(jet_pt, tr_jet_pt, prior_weight);
                // h2_truthreco_z->Fill(z, tr_z, prior_weight);
                // h2_truthreco_jt->Fill(jt, tr_jt, prior_weight);
                // h2_truthreco_r->Fill(r, tr_r, prior_weight);
                
                // // For a cross check on the unfold prior systematic
                // h2_ptz_weighted->Fill(z, jet_pt, prior_weight);
                // h2_ptjt_weighted->Fill(jt, jet_pt, prior_weight);
                // h2_ptr_weighted->Fill(r, jet_pt, prior_weight);  
                // h1_z_weighted->Fill(z, prior_weight);
                // h1_jt_weighted->Fill(jt, prior_weight);
                // h1_r_weighted->Fill(r, prior_weight);       
                
                // h2_ptz->Fill(z, jet_pt);
                // h2_ptjt->Fill(jt, jet_pt);
                // h2_ptr->Fill(r, jet_pt);
                // h1_z->Fill(z);
                // h1_jt->Fill(jt);
                // h1_r->Fill(r);
                
                event_counter++;
        }
                
        TH2 *h2_response_jetpt  = response_jetpt->Hresponse();
        TH2 *h2_response_rl     = response_rl->Hresponse();
        TH2 *h2_response_weight = response_weight->Hresponse();
        
        TH2 *h3_response_npair       = response_npair->Hresponse();        
        TH2 *h3_response_npair_eqch  = response_npair_eqch->Hresponse();
        TH2 *h3_response_npair_neqch = response_npair_neqch->Hresponse();
                
        h2_response_jetpt->GetXaxis()->SetTitle("reco p_{T, jet} [GeV/c]");
        h2_response_jetpt->GetYaxis()->SetTitle("truth p_{T, jet} [GeV/c]");
        h2_response_rl->GetXaxis()->SetTitle("reco R_{L}");
        h2_response_rl->GetYaxis()->SetTitle("truth R_{L}");
        h2_response_weight->GetXaxis()->SetTitle("reco w");
        h2_response_weight->GetYaxis()->SetTitle("truth w");
        
        // h3_response_npair->GetXaxis()->SetTitle("reco bin # (flattened)");
        // h3_response_npair->GetYaxis()->SetTitle("truth bin # (flattened)");
        // h3_response_npair_eqch->GetXaxis()->SetTitle("reco bin # (flattened)");
        // h3_response_npair_eqch->GetYaxis()->SetTitle("truth bin # (flattened)");
        // h3_response_npair_neqch->GetXaxis()->SetTitle("reco bin # (flattened)");
        // h3_response_npair_neqch->GetYaxis()->SetTitle("truth bin # (flattened)");
        
        h2_response_jetpt->Write("response_jetpt");
        h2_response_rl->Write("response_rl");
        h2_response_weight->Write("response_weight");
        
        h3_response_npair->Write("response_npair");
        h3_response_npair_eqch->Write("response_npair_eqch");
        h3_response_npair_neqch->Write("response_npair_neqch");    

        response_jetpt->Write("Roo_response_jetpt");
        response_rl->Write("Roo_response_rl");
        response_weight->Write("Roo_response_weigh");
        
        response_npair->Write("Roo_response_npair" );        
        response_npair_eqch->Write("Roo_response_npair_eqch" );
        response_npair_neqch->Write( "Roo_response_npair_neqch");
                
        TH1D *h1_purity_jetpt = (TH1D *)h1_num_purity_jetpt->Clone("h1_purity_jetpt");
        h1_purity_jetpt->Divide(h1_num_purity_jetpt, h1_denom_purity_jetpt, 1, 1, "B");
                
        TH1D *h1_efficiency_jetpt = (TH1D *)h1_num_efficiency_jetpt->Clone("h1_efficiency_jetpt");
        h1_efficiency_jetpt->Divide(h1_num_efficiency_jetpt, h1_denom_efficiency_jetpt, 1, 1, "B");
                                                
        h3_efficiency_rl_jetpt_weight->Divide(h3_num_efficiency_rl_jetpt_weight, h3_denom_efficiency_rl_jetpt_weight, 1, 1, "B");
        h3_efficiency_rl_jetpt_weight_eqch->Divide(h3_num_efficiency_rl_jetpt_weight_eqch, h3_denom_efficiency_rl_jetpt_weight_eqch, 1, 1, "B");
        h3_efficiency_rl_jetpt_weight_neqch->Divide(h3_num_efficiency_rl_jetpt_weight_neqch, h3_denom_efficiency_rl_jetpt_weight_neqch, 1, 1, "B");    
        
        h3_purity_rl_jetpt_weight->Divide(h3_num_purity_rl_jetpt_weight, h3_denom_purity_rl_jetpt_weight, 1, 1);
        h3_purity_rl_jetpt_weight_eqch->Divide(h3_num_purity_rl_jetpt_weight_eqch, h3_denom_purity_rl_jetpt_weight_eqch, 1, 1);
        h3_purity_rl_jetpt_weight_neqch->Divide(h3_num_purity_rl_jetpt_weight_neqch, h3_denom_purity_rl_jetpt_weight_neqch, 1, 1);
        
        h1_efficiency_jetpt->Write("efficiency_jetpt");
        h1_purity_jetpt->Write("purity_jetpt");
        
        h3_efficiency_rl_jetpt_weight->Write("efficiency_rl_jetpt_weight");
        h3_efficiency_rl_jetpt_weight_eqch->Write("efficiency_rl_jetpt_weight_eqch");
        h3_efficiency_rl_jetpt_weight_neqch->Write("efficiency_rl_jetpt_weight_neqch");    
        
        h3_purity_rl_jetpt_weight->Write("purity_rl_jetpt_weight");
        h3_purity_rl_jetpt_weight_eqch->Write("purity_rl_jetpt_weight_eqch");
        h3_purity_rl_jetpt_weight_neqch->Write("purity_rl_jetpt_weight_neqch");    
        
        // TH2D *h2_efficiency_zjt_ptbinned[ptbinsize-1];
        // TH2D *h2_efficiency_zr_ptbinned[ptbinsize-1];
        // TH2D *h2_efficiency_jtr_ptbinned[ptbinsize-1];  
        // TH2D *h2_purity_zjt_ptbinned[ptbinsize-1];
        // TH2D *h2_purity_zr_ptbinned[ptbinsize-1];
        // TH2D *h2_purity_jtr_ptbinned[ptbinsize-1];      
        
        // THStack *hs_efficiency_ptz  = new THStack("efficiency_z_all", ";z;efficiency");
        // THStack *hs_efficiency_ptjt = new THStack("efficiency_jt_all", ";j_{T} [GeV/c];efficiency");
        // THStack *hs_efficiency_ptr  = new THStack("efficiency_r_all", ";r;efficiency");
        // THStack *hs_purity_ptz  = new THStack("purity_z_all", ";z;purity");
        // THStack *hs_purity_ptjt = new THStack("purity_jt_all", ";j_{T} [GeV/c];purity");
        // THStack *hs_purity_ptr  = new THStack("purity_r_all", ";r;purity");
                                
        // for (int i = 1; i < ptbinsize; i++) {   
        //         TH1D *h1_ptz_temp = (TH1D *)h2_efficiency_ptz->ProjectionX(Form("efficiency_z_pt%d", i), i + 1, i + 1); 
        //         h1_ptz_temp->SetStats(0);
        //         h1_ptz_temp->SetMarkerStyle(i + 20);
                
        //         if (i!=5) {
        //                 h1_ptz_temp->SetMarkerColor(i);
        //                 h1_ptz_temp->SetLineColor(i);
        //         } else {
        //                 h1_ptz_temp->SetMarkerColor(i*i+3);
        //                 h1_ptz_temp->SetLineColor(i*i+3);
        //         }

        //         h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));
        //         hs_efficiency_ptz->Add(h1_ptz_temp);
                
        //         TH1D *h1_ptjt_temp = (TH1D *)h2_efficiency_ptjt->ProjectionX(Form("efficiency_jt_pt%d", i), i + 1, i + 1); 
        //         h1_ptjt_temp->SetStats(0);
        //         h1_ptjt_temp->SetMarkerStyle(i + 20);
                
        //         if (i!=5) {
        //                 h1_ptjt_temp->SetMarkerColor(i);
        //                 h1_ptjt_temp->SetLineColor(i);
        //         } else {
        //                 h1_ptjt_temp->SetMarkerColor(i*i+3);
        //                 h1_ptjt_temp->SetLineColor(i*i+3);
        //         }

        //         h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
        //         hs_efficiency_ptjt->Add(h1_ptjt_temp);
                
        //         TH1D *h1_ptr_temp = (TH1D *)h2_efficiency_ptr->ProjectionX(Form("efficiency_r_pt%d", i), i + 1, i + 1); 
        //         h1_ptr_temp->SetStats(0);
        //         h1_ptr_temp->SetMarkerStyle(i + 20);

        //         if (i!=5) {
        //                 h1_ptr_temp->SetMarkerColor(i);
        //                 h1_ptr_temp->SetLineColor(i);
        //         } else {
        //                 h1_ptr_temp->SetMarkerColor(i*i+3);
        //                 h1_ptr_temp->SetLineColor(i*i+3);
        //         }

        //         h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
        //         hs_efficiency_ptr->Add(h1_ptr_temp);
                
                
        //         h1_ptz_temp = (TH1D *)h2_purity_ptz->ProjectionX(Form("purity_z_pt%d", i), i + 1, i + 1); 
        //         h1_ptz_temp->SetStats(0);
        //         h1_ptz_temp->SetMarkerStyle(i + 20);
                
        //         if (i!=5) {
        //                 h1_ptz_temp->SetMarkerColor(i);
        //                 h1_ptz_temp->SetLineColor(i);
        //         } else {
        //                 h1_ptz_temp->SetMarkerColor(i*i+3);
        //                 h1_ptz_temp->SetLineColor(i*i+3);
        //         }

        //         h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
        //         hs_purity_ptz->Add(h1_ptz_temp);
                
        //         h1_ptjt_temp = (TH1D *)h2_purity_ptjt->ProjectionX(Form("purity_jt_pt%d", i), i + 1, i + 1); 
        //         h1_ptjt_temp->SetStats(0);
        //         h1_ptjt_temp->SetMarkerStyle(i + 20);

        //         if (i!=5) {
        //                 h1_ptjt_temp->SetMarkerColor(i);
        //                 h1_ptjt_temp->SetLineColor(i);
        //         } else {
        //                 h1_ptjt_temp->SetMarkerColor(i*i+3);
        //                 h1_ptjt_temp->SetLineColor(i*i+3);
        //         }

        //         h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
        //         hs_purity_ptjt->Add(h1_ptjt_temp);
                
        //         h1_ptr_temp = (TH1D *)h2_purity_ptr->ProjectionX(Form("purity_r_pt%d", i), i + 1, i + 1); 
        //         h1_ptr_temp->SetStats(0);
        //         h1_ptr_temp->SetMarkerStyle(i + 20);

        //         if (i!=5) {
        //                 h1_ptr_temp->SetMarkerColor(i);
        //                 h1_ptr_temp->SetLineColor(i);
        //         } else {
        //                 h1_ptr_temp->SetMarkerColor(i*i+3);
        //                 h1_ptr_temp->SetLineColor(i*i+3);
        //         }

        //         h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
        //         hs_purity_ptr->Add(h1_ptr_temp); 
                                        
        //         h3_efficiency_ptzjt->GetZaxis()->SetRange(i+1, i+1);      
        //         h2_efficiency_zjt_ptbinned[i-1] = (TH2D *)h3_efficiency_ptzjt->Project3D("yx");
        //         h2_efficiency_zjt_ptbinned[i-1]->SetStats(0);
        //         h2_efficiency_zjt_ptbinned[i-1]->SetName(Form("efficiency_zjt_pt%d", i));
        //         h2_efficiency_zjt_ptbinned[i-1]->Write(); 

        //         h3_efficiency_ptzr->GetZaxis()->SetRange(i+1, i+1); 
        //         h2_efficiency_zr_ptbinned[i-1] = (TH2D *)h3_efficiency_ptzr->Project3D("yx");
        //         h2_efficiency_zr_ptbinned[i-1]->SetStats(0);       
        //         h2_efficiency_zr_ptbinned[i-1]->SetName(Form("efficiency_zr_pt%d",i)); 
        //         h2_efficiency_zr_ptbinned[i-1]->Write();
                
        //         h3_efficiency_ptjtr->GetZaxis()->SetRange(i+1, i+1);        
        //         h2_efficiency_jtr_ptbinned[i-1] = (TH2D *)h3_efficiency_ptjtr->Project3D("yx");
        //         h2_efficiency_jtr_ptbinned[i-1]->SetStats(0);        
        //         h2_efficiency_jtr_ptbinned[i-1]->SetName(Form("efficiency_jtr_pt%d",i));        
        //         h2_efficiency_jtr_ptbinned[i-1]->Write();
                
        //         h3_purity_ptzjt->GetZaxis()->SetRange(i+1, i+1);      
        //         h2_purity_zjt_ptbinned[i-1] = (TH2D *)h3_purity_ptzjt->Project3D("yx");
        //         h2_purity_zjt_ptbinned[i-1]->SetStats(0);                       
        //         h2_purity_zjt_ptbinned[i-1]->SetName(Form("purity_zjt_pt%d",i)); 
        //         h2_purity_zjt_ptbinned[i-1]->Write();         

        //         h3_purity_ptzr->GetZaxis()->SetRange(i+1, i+1);        
        //         h2_purity_zr_ptbinned[i-1] = (TH2D *)h3_purity_ptzr->Project3D("yx");
        //         h2_purity_zr_ptbinned[i-1]->SetStats(0);
        //         h2_purity_zr_ptbinned[i-1]->SetName(Form("purity_zr_pt%d",i));
        //         h2_purity_zr_ptbinned[i-1]->Write();        
                
        //         h3_purity_ptjtr->GetZaxis()->SetRange(i+1, i+1);        
        //         h2_purity_jtr_ptbinned[i-1] = (TH2D *)h3_purity_ptjtr->Project3D("yx");
        //         h2_purity_jtr_ptbinned[i-1]->SetStats(0);
        //         h2_purity_jtr_ptbinned[i-1]->SetName(Form("purity_jtr_pt%d",i));
        //         h2_purity_jtr_ptbinned[i-1]->Write();                        
        // }   
        
        // hs_efficiency_ptz->Write();
        // hs_efficiency_ptjt->Write();
        // hs_efficiency_ptr->Write();
        // hs_purity_ptz->Write();
        // hs_purity_ptjt->Write();
        // hs_purity_ptr->Write();       
                        
        // NormalizeHist(h2_ptz_weighted);
        // NormalizeHist(h2_ptjt_weighted);
        // NormalizeHist(h2_ptr_weighted);   
        
        // NormalizeHist(h1_z_weighted);
        // NormalizeHist(h1_jt_weighted);
        // NormalizeHist(h1_r_weighted);    
                
        // NormalizeHist(h2_ptz);
        // NormalizeHist(h2_ptjt);
        // NormalizeHist(h2_ptr);   
        
        // NormalizeHist(h1_z);
        // NormalizeHist(h1_jt);
        // NormalizeHist(h1_r); 
                        
        // h2_ptz_weighted->Write();
        // h2_ptjt_weighted->Write();
        // h2_ptr_weighted->Write();
        // h1_z_weighted->Write();
        // h1_jt_weighted->Write();
        // h1_r_weighted->Write();   
        
        // h2_ptz->Write();
        // h2_ptjt->Write();
        // h2_ptr->Write();
        // h1_z->Write();
        // h1_jt->Write();
        // h1_r->Write();       

        cout << event_counter << " events processed" << endl;
                
        /////////////////////////////////////////////////////
        TH1D *h1_jetpt_true = (TH1D *)h1_denom_efficiency_jetpt->Clone("h1_jetpt_true");
        TH1D *h1_jetpt_reco = (TH1D *)h1_denom_purity_jetpt->Clone("h1_jetpt_reco");
        
        h1_jetpt_reco->Multiply(h1_purity_jetpt);
        
        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_reco, 4);
        
        h1_jetpt_reco = (TH1D *)unfold_jetpt.Hreco();
        
        h1_jetpt_reco->Divide(h1_efficiency_jetpt);
        
        ////////////////////////////////////////////////////

        f->Close();     
        
}
