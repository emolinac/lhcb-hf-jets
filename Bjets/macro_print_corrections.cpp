#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"

void macro_print_corrections()
{
        TFile* f_mc     = new TFile((output_folder + "ntuple_bjets_mc.root").c_str());
        TFile* f_mcreco = new TFile((output_folder + "ntuple_bjets_mcreco.root").c_str());

        TTree* tree_mc = (TTree*) f_mc->Get("BTree");
        TTree* tree_mcreco = (TTree*) f_mcreco->Get("BTree");

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

        double truth_z, truth_jt, truth_r;
        double truth_z_b, truth_jt_b, truth_r_b;
        double truth_zg, truth_jtg, truth_rg;
        double meas_z, meas_jt, meas_r;
        double meas_inJet_z, meas_inJet_jt, meas_inJet_r;
        double jet_pt_recotruthratio, HF_pt_recotruthratio;

        int nsplits, ndtrs, NumHFHads;
        int meas_nsplits, meas_ndtrs;
        int eventNumber;
        int GluonTag, nTracks;
        //    bool GluonTag, Hasbbbar;
        
        int NumDtrRecoHF;
        bool hasRecoHF, Hasbbbar;

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

        BTree->SetBranchAddress("GluonTag", &GluonTag);
        BTree->SetBranchAddress("nTracks", &nTracks);
        
        BTree->SetBranchAddress("NumDtrRecoHF", &NumDtrRecoHF);
        BTree->SetBranchAddress("jet_pt_recotruthratio", &jet_pt_recotruthratio);
        BTree->SetBranchAddress("HF_pt_recotruthratio", &HF_pt_recotruthratio);

        
}