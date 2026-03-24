#include "./BMCDecayTree.C"
#include "./BDecayTree.C"
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"


using namespace std;

void HFRecoEff(int NumEvts_user = 10000, int dataset = 1510,
               bool chargedJetCut_user = false,
               bool DoRecSelEff = false,
               bool DoSignalSys = false,
               bool SubtractGS = false)
{

  int NumEvts = NumEvts_user;
  int NumEvtsMC = NumEvts_user;
  followHardest = false;
  truthLevel = true;
  chargedJetCut = chargedJetCut_user;
  if (truthLevel)
  {
    ghostCut = false;
  }

  int year = (dataset / 10000) % 10;
  int JetMeth = (dataset / 1000) % 10;
  int flavor = (dataset / 100) % 10;
  int ptRange = (dataset / 10) % 10;
  int Mag = (dataset / 1) % 10;

  int HF_pdgcode = -999;

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

  //
  //  Naming Convention:
  //  {$1}{$2}{$3}{$4}{$5}
  //  $1: Year: 2016 = 6, 2017 = 7, 2018 = 8, All = 9
  //  $2: jets or dijets: TaggedDijets = 1, else = 2
  //  $3: flavor: b = 5, c = 4, udsg = 1
  //  $4: pT: pt15pt20 = 0, pt20pt50 = 1, pt50 = 2, else = 3
  //  $5: Magnet: MD = 0, MU = 1, Both = 9

  TString str_Mag = "";
  TString str_pT = "";
  TString str_level = "";
  TString str_followHard = "";
  TString str_flavor = "";
  TString str_ghost = "";
  TString str_PID = "";
  TString str_charged = "";
  TString str_GS = "";

  TString str_year = "2016";

  if (SubtractGS)
  {
    str_GS = "_GSsub";
  }
  if (year == 6)
    str_year = "2016";
  else if (year == 7)
    str_year = "2017";
  else if (year == 8)
    str_year = "2018";

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

  if (PID_cut)
    str_PID = "_PID";

  if (truthLevel)
    str_level = "truth";
  else
    str_level = "reco";

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
  TString str_tree;
    

  TString extension_read, extension_RootFilesMC, extension, extension_RootFilesData;
  extension_RootFilesMC = TString("../../root_files/BjetsMC/");
  extension_RootFilesData = TString("../../root_files/Bjets/");  

  TString  extension_eff,  extension_prefix;
    
  extension = TString("HFeff_") + str_level + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + str_PID + str_GS + Form("_%d", dataset);
  extension_eff = TString("efficiency_truth") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_charged + str_Mag + str_flavor + str_GS + Form("_%d", dataset);

  cout << extension_eff << endl;

  TFile *file_eff = new TFile((output_folder + "bjets_efficiencies.root").c_str(), "READ");

  TH1D *h1_denom_efficiency_HFpt = (TH1D *)file_eff->Get("h1_denom_efficiency_HFpt");
  TH1D *h1_denom_efficiency_jetpt = (TH1D *)file_eff->Get("h1_denom_efficiency_jetpt");
//  TH2D *h2_denom_efficiency_HFpteta = (TH2D *)file_eff->Get("h2_denom_efficiency_HFpteta");
    TH2D *h2_denom_efficiency_HFpteta = (TH2D *)file_eff->Get("denom_efficiency_HFpteta");
//  TH3D *h3_denom_efficiency_HFptetajetpt = (TH3D *)file_eff->Get("h3_denom_efficiency_HFptetajetpt");
    TH3D *h3_denom_efficiency_HFptetajetpt = (TH3D *)file_eff->Get("denom_efficiency_HFptetajetpt");
  TH2D *h2_denom_efficiency_HFptjetpt = (TH2D *)file_eff->Get("h2_denom_efficiency_HFptjetpt");
//  TH2D *h2_denom_efficiency_HFptnTracks = (TH2D *)file_eff->Get("h2_denom_efficiency_HFptnTracks");
    TH2D *h2_denom_efficiency_HFptnTracks = (TH2D *)file_eff->Get("denom_efficiency_HFptnTracks");

  /////////////////// Mass Fit Parameters /////////////////////////////////
  //TString extension_mass(extension_RootFilesData + "massfit_data_ev_-1_ptj_12250_eta_2.54.0_ghost_0.5_b_PID_91599.root");
  TString extension_mass( extension_RootFilesData + TString("massfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.5_b" + str_PID + "_91599.root");  
  if (DoRecSelEff)
    extension_mass = "recselsys_" + extension_mass;
  if (DoSignalSys)
    extension_mass = "sys_" + extension_mass;

  TFile f_massfit((output_folder + "mass-fits/results_mass_fit_data.root").c_str(), "READ");
  TH1D *h1_MassMin = (TH1D *)f_massfit.Get("h1_MassMin");
  TH1D *h1_MassMax = (TH1D *)f_massfit.Get("h1_MassMax");
  TH1D *h1_BkgScale = (TH1D *)f_massfit.Get("h1_BkgScale");
  TH1D *h1_BkgScaleNear = (TH1D *)f_massfit.Get("h1_BkgScale_forSysNear");
  TH1D *h1_BkgScaleFar = (TH1D *)f_massfit.Get("h1_BkgScale_forSysFar");

  // BMCDecayTree MCTree(0, dataset, false);
   BDecayTree Tree(0, dataset, false);

  // cout << "Total number of MC events = " << MCTree.fChain->GetEntries() << endl;
  // cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;

  TChain *HFRecoTree = new TChain("HFRecoTree", "");
  vector<int> vec_datasets;
  if (Mag == 0)
    vec_datasets = {61590, 71590, 81590};
  else if (Mag == 1)
    vec_datasets = {61591, 71591, 81591};
  else
    vec_datasets = {61590, 61591, 71590, 71591, 81590, 81591};

  cout << "Loading dataset(s): " << endl;
  if (year == 9 && JetMeth != 9)
  {
    for (int i = 0; i < vec_datasets.size(); i++)
    {
      Mag = (vec_datasets[i] / 1) % 10;
      if (Mag == 0)
        str_Mag = "_MD";
      else if (Mag == 1)
        str_Mag = "_MU";
      extension_read = TString("tree_HFeff_") + str_level + Form("_ev_%d", NumEvts) + str_charged + str_Mag + str_flavor + Form("_%d", vec_datasets[i]);
      cout << extension_read << endl;
      HFRecoTree->Add(extension_RootFilesMC + extension_read + ".root/HFRecoTree");
    }
  }
  else
  {
    extension_read = TString("tree_HFeff_") + str_level + Form("_ev_%d", NumEvts) + str_charged + str_Mag + str_flavor + Form("_%d", dataset);
    // extension_read = extension_prefix + extension_read;
    cout << extension_read << endl;
    HFRecoTree->Add(extension_RootFilesMC  + extension_read + ".root/HFRecoTree");
  }

  if (NumEvts > HFRecoTree->GetEntries())
    NumEvts = HFRecoTree->GetEntries();
  if (NumEvts == -1)
    NumEvts = HFRecoTree->GetEntries();

  // int NumEvts = 0;
  // cout <<"Choose number of events (-1: All Events, or enter integer): ";
  // cin>> NumEvts;

//  cout << "Executing CAJetAlgo" << endl;

  TFile f((output_folder + "HF_efficiencies.root").c_str(), "RECREATE");

  TH1D *h1_HF_pt_reco = new TH1D("h1_HF_pt_reco", "", ptHFbinsize, ptHF_binedges);
  TH1D *h1_HF_eta_reco = new TH1D("h1_HF_eta_reco", "", 20, 2, 5);
  TH1D *h1_HF_rap_reco = new TH1D("h1_HF_rap_reco", "", 20, 2, 5);
  TH1D *h1_HF_mass_reco = new TH1D("h1_HF_mass_reco", "", 60, 5.18, 5.38);
  TH1D *h1_HF_dtfmass_reco = new TH1D("h1_HF_dtfmass_reco", "", 60, 5.18, 5.38);

  TH1D *h1_HF_pt = new TH1D("h1_HF_pt", "", ptHFbinsize, ptHF_binedges);
  TH1D *h1_HF_eta = new TH1D("h1_HF_eta", "", HFetabinsize, HFeta_binedges);

  TH2D *h2_HF_pteta = new TH2D("h2_HF_pteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
  TH2D *h2_HF_pteta_reco = new TH2D("h2_HF_pteta_reco", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
  TH2D *h2_jetpt_jetrap = new TH2D("h2_jetpt_jetrap", "", ptbinsize, pt_binedges, HFetabinsize, HFeta_binedges);
  TH2D *h2_jetpt_HFrap = new TH2D("h2_jetpt_HFrap", "", ptbinsize, pt_binedges, HFetabinsize, HFeta_binedges);

  TH1D *h1_num_efficiency_HFpt = new TH1D("num_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
  TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", ptbinsize, pt_binedges);
//  TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", longptbinsize, longpt_binedges);
  // TH1D *h1_denom_efficiency_HFpt = new TH1D("denom_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
  TH1D *h1_num_efficiency_HFeta = new TH1D("num_efficiency_HFeta", "", HFetabinsize, HFeta_binedges);
  TH1D *h1_denom_efficiency_HFeta = new TH1D("denom_efficiency_HFeta", "", HFetabinsize, HFeta_binedges);
    
  TH2D *h2_num_efficiency_HFpteta = new TH2D("num_efficiency_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
//    TH2D *h2_denom_efficiency_HFpteta = new TH2D("denom_efficiency_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
//    
    
  TH3D *h3_num_efficiency_HFptetajetpt = new TH3D("num_efficiency_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, customptbinsize, custompt_binedges);
  TH2D *h2_num_efficiency_HFptjetpt = new TH2D("num_efficiency_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
  TH2D *h2_num_efficiency_HFptnTracks = new TH2D("num_efficiency_HFptnTracks", "", ptHFbinsize, ptHF_binedges, nTracksbinsize, nTrack_binedges);
  // TH2D *h2_denom_efficiency_HFpteta = new TH2D("denom_efficiency_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);

  TH1D *h1_HF_IPCHI2 = new TH1D("h1_HF_IPCHI2", "", 60, -1, 1000);
  TH1D *h1_HF_LOGIPCHI2 = new TH1D("h1_HF_LOGIPCHI2", "", 60, -5, 10);

//  JetDefinition jet_def(cambridge_algorithm, jetradius);

  // TClonesArray *arr = new TClonesArray("TLorentzVector");
    float HF_px, HF_py, HF_pz, HF_e;
    float Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e;
    float K_px, K_py, K_pz, K_e;
    float mum_px, mum_py, mum_pz, mum_e;
    float mup_px, mup_py, mup_pz, mup_e;
    float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e; 
    //float tr_jet_rap, tr_jet_pt;
    float dtf_mass, dtf_chi2ndf, dtf_ctau, K_PIDK;
    int eventNumber, nTracks, NumHFHads;
    bool Hasbbbar;

      bool mup_L0, mum_L0;
      bool jpsi_L0, jpsi_Hlt1, jpsi_Hlt2;
      bool Trig, TIS, TOS;
      
      // Not needed -- Obsolete/Irrelevant Variables
      float mup_TRCHI2DOF, mup_IPCHI2, mup_TRGHOSTPROB;
      float mum_TRCHI2DOF, mum_IPCHI2, mum_TRGHOSTPROB;
      float K_TRCHI2DOF, K_IPCHI2, K_TRGHOSTPROB;
      bool K_ISLONG, mup_ISMUON, mum_ISMUON;
      float Jpsi_TRCHI2DOF, Jpsi_IPCHI2, Jpsi_TRGHOSTPROB, Jpsi_AMAXDOCA, Jpsi_BPVDIRA, Jpsi_BPVVDCHI2, Jpsi_VTXCHI2;
      float Bu_TRCHI2DOF, Bu_IPCHI2, Bu_TRGHOSTPROB, Bu_AMAXDOCA, Bu_BPVDIRA, Bu_BPVIPCHI2, Bu_VTXCHI2;

    
  HFRecoTree->SetBranchAddress("eventNumber", &eventNumber);

  HFRecoTree->SetBranchAddress("mup_px", &mup_px);
  HFRecoTree->SetBranchAddress("mup_py", &mup_py);
  HFRecoTree->SetBranchAddress("mup_pz", &mup_pz);
  HFRecoTree->SetBranchAddress("mup_pe", &mup_e);
  HFRecoTree->SetBranchAddress("mup_ISMUON", &mup_ISMUON);

  HFRecoTree->SetBranchAddress("mum_px", &mum_px);
  HFRecoTree->SetBranchAddress("mum_py", &mum_py);
  HFRecoTree->SetBranchAddress("mum_pz", &mum_pz);
  HFRecoTree->SetBranchAddress("mum_pe", &mum_e);
  HFRecoTree->SetBranchAddress("mum_ISMUON", &mum_ISMUON);

  HFRecoTree->SetBranchAddress("Jpsi_px", &Jpsi_px);
  HFRecoTree->SetBranchAddress("Jpsi_py", &Jpsi_py);
  HFRecoTree->SetBranchAddress("Jpsi_pz", &Jpsi_pz);
  HFRecoTree->SetBranchAddress("Jpsi_pe", &Jpsi_e);

  HFRecoTree->SetBranchAddress("K_px", &K_px);
  HFRecoTree->SetBranchAddress("K_py", &K_py);
  HFRecoTree->SetBranchAddress("K_pz", &K_pz);
  HFRecoTree->SetBranchAddress("K_pe", &K_e);
  HFRecoTree->SetBranchAddress("K_ISLONG", &K_ISLONG);
    
  HFRecoTree->SetBranchAddress("HF_px", &HF_px);
  HFRecoTree->SetBranchAddress("HF_py", &HF_py);
  HFRecoTree->SetBranchAddress("HF_pz", &HF_pz);
  HFRecoTree->SetBranchAddress("HF_pe", &HF_e);

  HFRecoTree->SetBranchAddress("tr_jet_px", &tr_jet_px);
  HFRecoTree->SetBranchAddress("tr_jet_py", &tr_jet_py);
  HFRecoTree->SetBranchAddress("tr_jet_pz", &tr_jet_pz);
  HFRecoTree->SetBranchAddress("tr_jet_pe", &tr_jet_e);
  //HFRecoTree->SetBranchAddress("tr_jet_pt", &tr_jet_pt);
  //HFRecoTree->SetBranchAddress("tr_jet_rap", &tr_jet_rap);
  

  HFRecoTree->SetBranchAddress("dtf_mass", &dtf_mass);
  HFRecoTree->SetBranchAddress("dtf_chi2ndf", &dtf_chi2ndf);
  HFRecoTree->SetBranchAddress("dtf_ctau", &dtf_ctau);
  HFRecoTree->SetBranchAddress("K_PIDK", &K_PIDK);
  HFRecoTree->SetBranchAddress("nTracks", &nTracks);
  HFRecoTree->SetBranchAddress("Hasbbbar", &Hasbbbar);
  HFRecoTree->SetBranchAddress("NumHFHads", &NumHFHads);
    
  HFRecoTree->SetBranchAddress("mup_L0", &mup_L0);
  HFRecoTree->SetBranchAddress("mum_L0", &mum_L0);
  HFRecoTree->SetBranchAddress("jpsi_L0", &jpsi_L0);
  HFRecoTree->SetBranchAddress("jpsi_Hlt1", &jpsi_Hlt1);
  HFRecoTree->SetBranchAddress("jpsi_Hlt2", &jpsi_Hlt2);
  HFRecoTree->SetBranchAddress("Trig", &Trig);
  HFRecoTree->SetBranchAddress("TIS", &TIS);
  HFRecoTree->SetBranchAddress("TOS", &TOS);

//  HFRecoTree->SetBranchAddress("mup_TRCHI2DOF", &mup_TRCHI2DOF);
//  HFRecoTree->SetBranchAddress("mup_IPCHI2", &mup_IPCHI2);
//  HFRecoTree->SetBranchAddress("mup_TRGHOSTPROB", &mup_TRGHOSTPROB);
//  HFRecoTree->SetBranchAddress("mum_TRCHI2DOF", &mum_TRCHI2DOF);
//  HFRecoTree->SetBranchAddress("mum_IPCHI2", &mum_IPCHI2);
//  HFRecoTree->SetBranchAddress("mum_TRGHOSTPROB", &mum_TRGHOSTPROB);
//  HFRecoTree->SetBranchAddress("K_TRCHI2DOF", &K_TRCHI2DOF);
//  HFRecoTree->SetBranchAddress("K_IPCHI2", &K_IPCHI2);
//  HFRecoTree->SetBranchAddress("K_TRGHOSTPROB", &K_TRGHOSTPROB);
//  HFRecoTree->SetBranchAddress("Jpsi_TRCHI2DOF", &Jpsi_TRCHI2DOF);
//  HFRecoTree->SetBranchAddress("Jpsi_IPCHI2", &Jpsi_IPCHI2);
//  HFRecoTree->SetBranchAddress("Jpsi_TRGHOSTPROB", &Jpsi_TRGHOSTPROB);
//  HFRecoTree->SetBranchAddress("Jpsi_AMAXDOCA", &Jpsi_AMAXDOCA);
//  HFRecoTree->SetBranchAddress("Jpsi_BPVDIRA", &Jpsi_BPVDIRA);
//  HFRecoTree->SetBranchAddress("Jpsi_BPVVDCHI2", &Jpsi_BPVVDCHI2);
//  HFRecoTree->SetBranchAddress("Jpsi_VTXCHI2", &Jpsi_VTXCHI2);
//  HFRecoTree->SetBranchAddress("Bu_TRCHI2DOF", &Bu_TRCHI2DOF);
//  HFRecoTree->SetBranchAddress("Bu_IPCHI2", &Bu_IPCHI2);
//  HFRecoTree->SetBranchAddress("Bu_TRGHOSTPROB", &Bu_TRGHOSTPROB);
//  HFRecoTree->SetBranchAddress("Bu_AMAXDOCA", &Bu_AMAXDOCA);
//  HFRecoTree->SetBranchAddress("Bu_BPVDIRA", &Bu_BPVDIRA);
//  HFRecoTree->SetBranchAddress("Bu_BPVIPCHI2", &Bu_BPVIPCHI2);
//  HFRecoTree->SetBranchAddress("Bu_VTXCHI2", &Bu_VTXCHI2);
 
  //
  // Event loop
  //

  int eventNum;
  unsigned long long last_eventNum = 0;
  int events = 0;

  bool maxjetpT_found = false;
  int NumRecoHF = 0;

  int cut_npvs = 0;
  int cut_eta = 0;
  int cut_pt = 0;
  int cut_HFdR = 0;

  for (int ev = 0; ev < NumEvts; ev++)
  {
    HFRecoTree->GetEntry(ev);
    // EventTree->GetEntry(ev);
    // cout<<recojet_e<<", ";
    if (ev % 10000 == 0)
      cout << "Executing event " << ev << endl;

    // if (last_eventNum == Tree.eventNumber)
    //   continue;

    TLorentzVector HFjet, HFmeson, mup, mum, Kmeson, Jpsi;

    HFjet.SetPxPyPzE(tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e);
    HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
    mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
    mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);
    Kmeson.SetPxPyPzE(K_px, K_py, K_pz, K_e);
    Jpsi.SetPxPyPzE(Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e);
      
    float JpsiMass = 3.0969; // GeV
    bool PID_cond = (K_PIDK > 0);

    float MassHigh = h1_MassMax != NULL ? h1_MassMax->GetBinContent(h1_MassMax->FindBin(HFmeson.Pt())) : 5.31;
    float MassLow = h1_MassMin != NULL ? h1_MassMin->GetBinContent(h1_MassMin->FindBin(HFmeson.Pt())) : 5.24;
      
//      if (h1_MassMax == NULL){cout << "No Mass file "<< endl;}
//      float MassHigh =  5.31;
//      float MassLow =  5.24;
    bool mass_cond = (dtf_mass > MassLow && dtf_mass < MassHigh);

    bool K_cuts = (K_ISLONG) && (K_TRCHI2DOF < 4) && (K_IPCHI2 > 25.) && (K_TRGHOSTPROB < 0.4) && (Kmeson.Pt() > 0.25);
    bool mum_cuts = (mum_ISMUON) && (mum_TRCHI2DOF < 4) && (mum_IPCHI2 > 25.) && (mum_TRGHOSTPROB < 0.4);
    bool mup_cuts = (mup_ISMUON) && (mup_TRCHI2DOF < 4) && (mup_IPCHI2 > 25.) && (mup_TRGHOSTPROB < 0.4);
//
    bool Jpsi_cuts = (fabs(Jpsi.M() - JpsiMass) < 0.1 && (Jpsi_AMAXDOCA < 0.3));
    bool Jpsi_mom_cuts = (Jpsi_VTXCHI2 < 9) && (Jpsi_BPVDIRA > 0) && (Jpsi_BPVVDCHI2 > 225);
    bool Bu_mom_cuts = (Bu_BPVIPCHI2 < 25) && (Bu_VTXCHI2 < 45);
//     bool Bu_mom_cuts = (Bu_VTXCHI2 < 45);

    // cout << Jpsi_AMAXDOCA << ", " << Jpsi_BPVDIRA << ", " << Bu_BPVIPCHI2 << ", " << Bu_VTXCHI2 << endl;
//     if (!K_cuts)
//       continue;
//     if (!mum_cuts || !mup_cuts)
//       continue;
//     if (!Jpsi_cuts)
//       continue;
//     if (!Jpsi_mom_cuts)
//       continue;
//     if (!Bu_mom_cuts)
//       continue;
      
    if (!TOS)
        continue;
    if (!mass_cond)
      continue;

    // Not implement in MC since PIDCalib has no MC/Data tables
    // if (PID_cut && !PID_cond)
    //   continue;

    if (SubtractGS && Hasbbbar)
      continue;
 
    
    //bool rap_cond = tr_jet_rap > etaMin && tr_jet_rap < etaMax;
    //bool pt_cond = tr_jet_pt > pTLow;
    //bool num_cond = pt_cond && rap_cond;

    // cout << HF_px << ", " << HF_py << ", " << HF_pz << ", " << HF_e << endl;

    Jpsi = mup + mum;
    h1_HF_IPCHI2->Fill(Bu_BPVIPCHI2);
    h1_HF_LOGIPCHI2->Fill(log(Bu_BPVIPCHI2));

    h1_HF_pt_reco->Fill(HFmeson.Pt());
    h1_HF_eta_reco->Fill(HFmeson.Eta());
    h1_HF_mass_reco->Fill(HFmeson.M());
    h1_HF_dtfmass_reco->Fill(dtf_mass);
    
    //if (num_cond)
    //{
      h1_num_efficiency_HFpt->Fill(HFmeson.Pt());
      // h2_num_efficiency_HFpteta->Fill(HFmeson.Pt(), HFmeson.Eta());
      h2_num_efficiency_HFpteta->Fill(HFmeson.Pt(), HFmeson.Rapidity());
      h3_num_efficiency_HFptetajetpt->Fill(HFmeson.Pt(), HFmeson.Rapidity(), HFjet.Pt());
      // cout << HFmeson.Eta() << ", ";
      h2_num_efficiency_HFptjetpt->Fill(HFmeson.Pt(), HFjet.Pt());
      h2_num_efficiency_HFptnTracks->Fill(HFmeson.Pt(), nTracks);
      h1_num_efficiency_jetpt->Fill(HFjet.Pt());
    //}
    h1_HF_rap_reco->Fill(HFmeson.Rapidity());
    // h1_HF_eta_reco->Fill(HFmeson.Eta());
    h2_jetpt_jetrap->Fill(HFjet.Pt(), HFjet.Rapidity());
    h2_jetpt_HFrap->Fill(HFjet.Pt(), HFmeson.Rapidity());

    events++;
    //  }
  }

  cout << h1_denom_efficiency_HFpt->Integral() << endl;
  cout << h1_num_efficiency_HFpt->Integral() << endl;

  TH1D *h1_efficiency_HFpt = (TH1D *)h1_num_efficiency_HFpt->Clone("h1_efficiency_HFpt");
  h1_efficiency_HFpt->Divide(h1_num_efficiency_HFpt, h1_denom_efficiency_HFpt, 1, 1, "B");

  TH1D *h1_efficiency_jetpt = (TH1D *)h1_num_efficiency_jetpt->Clone("h1_efficiency_jetpt");
  h1_efficiency_jetpt->Divide(h1_num_efficiency_jetpt, h1_denom_efficiency_jetpt, 1, 1, "B");
  // cout << h1_num_efficiency_HFpt->Integral() << ", " << h1_denom_efficiency_HFpt->Integral() << endl;

  TH2D *h2_efficiency_HFpteta = (TH2D *)h2_num_efficiency_HFpteta->Clone("h2_efficiency_HFpteta");
  h2_efficiency_HFpteta->Divide(h2_num_efficiency_HFpteta, h2_denom_efficiency_HFpteta, 1, 1, "B");

  TH3D *h3_efficiency_HFptetajetpt = (TH3D *)h3_num_efficiency_HFptetajetpt->Clone("h3_efficiency_HFptetajetpt");
  h3_efficiency_HFptetajetpt->Divide(h3_num_efficiency_HFptetajetpt, h3_denom_efficiency_HFptetajetpt, 1, 1, "B");

  TH2D *h2_efficiency_HFptjetpt = (TH2D *)h2_num_efficiency_HFptjetpt->Clone("h2_efficiency_HFptjetpt");
  h2_efficiency_HFptjetpt->Divide(h2_num_efficiency_HFptjetpt, h2_denom_efficiency_HFptjetpt, 1, 1, "B");

  TH2D *h2_efficiency_HFptnTracks = (TH2D *)h2_num_efficiency_HFptnTracks->Clone("h2_efficiency_HFptnTracks");
  h2_efficiency_HFptnTracks->Divide(h2_num_efficiency_HFptnTracks, h2_denom_efficiency_HFptnTracks, 1, 1, "B");

  h1_efficiency_HFpt->Write("efficiency_HFpt");
  h1_efficiency_jetpt->Write("efficiency_jetpt");
  h1_num_efficiency_HFpt->Write("num_efficiency_HFpt");
  h1_denom_efficiency_HFpt->Write("denom_efficiency_HFpt");

  h2_efficiency_HFpteta->Write("efficiency_HFpteta");
  h2_num_efficiency_HFpteta->Write("num_efficiency_HFpteta");
  h2_denom_efficiency_HFpteta->Write("denom_efficiency_HFpteta");

  h3_efficiency_HFptetajetpt->Write("efficiency_HFptetajetpt");
  h3_num_efficiency_HFptetajetpt->Write("num_efficiency_HFptetajetpt");
  h3_denom_efficiency_HFptetajetpt->Write("denom_efficiency_HFptetajetpt");

  h2_efficiency_HFptjetpt->Write("efficiency_HFptjetpt");
  h2_num_efficiency_HFptjetpt->Write("num_efficiency_HFptjetpt");
  h2_denom_efficiency_HFptjetpt->Write("denom_efficiency_HFptjetpt");

  h2_efficiency_HFptnTracks->Write("efficiency_HFptnTracks");
  h2_num_efficiency_HFptnTracks->Write("num_efficiency_HFptnTracks");
  h2_denom_efficiency_HFptnTracks->Write("denom_efficiency_HFptnTracks");

//  ////////////////////////////////////////////
   TF1 *f_eff1 = new TF1("f_eff1", "[0]*tanh([1]*(x-[3]))+[2]", 0, 150);
    TF1 *f_eff2 = new TF1("f_eff2", "[0]*tanh([1]*([3] - x))+[2]", 25, 100);
   f_eff1->SetParameters(0.2, 0.1, 0.2, 6.4);
   f_eff2->SetParameters(0.2, 0.1, 0.2, 100);
   f_eff2->SetParLimits(3, 80, 200);
//   h1_efficiency_HFpt->Fit("f_eff1", "E", "SAME", 2, 25);

//  h1_efficiency_HFpt->Fit("pol4", "E", "SAME", 0, 150);

  // f_eff1->Write("f_eff1");

  //////////////////////////////////////////

  cout << "Total number of events processed = " << events << endl;
  cout << "NumRecoHF = " << NumRecoHF << endl;

  cout << "Events blocked: " << endl;
  cout << "npvs = " << cut_npvs << endl;
  cout << "jetpt = " << cut_pt << endl;
  cout << "jeteta = " << cut_eta << endl;
  cout << "Dinjet = " << cut_HFdR << endl;
  // cout<<"hasHF = "<< cut_hasHF<<endl;
//
  //
  int ican = -1, iframe = -1, itext = -1;
  TCanvas *ccan[1000];
  TH1F *frame[1000];
  TLatex *text[1000];
  for (int i = 0; i < 1000; i++)
  {
    text[i] = new TLatex();
    text[i]->SetNDC(kTRUE);
    text[i]->SetTextSize(0.06);
  }
  TLatex Tl;
  Tl.SetNDC(kTRUE);
  Tl.SetTextSize(0.04);
  //
  Int_t lhcbFont = 132; // Old LHCb style: 62;
  // Line thickness
  float_t lhcbWidth = 2.00; // Old LHCb style: 3.00;
  // Text size
  float_t lhcbTSize = 0.06;

  // gROOT->SetStyle("Plain");
  gStyle->SetCanvasDefH(1414);
  gStyle->SetCanvasDefW(1000);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.12);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.2);
  gStyle->SetStatStyle(0);
  // gStyle->SetLineWidth(3);
  gStyle->SetOptStat(0);
  //---- paint...
  char buf[100];
  char bufb[100];
  TString rootfile;
  TString plotfile;
  TString plotfilePDF;
  TString plotfileO;
  TString plotfileC;
  // TString OutputFileBase	= outbase+outinfo;
  TString plotextension = TString("../../plots/BjetsMC/");
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

  // //
  //
  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 1);
  ccan[ican]->SetFillColor(10);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.06);
  gPad->SetRightMargin(0.15);
  ccan[ican]->cd();
  ccan[ican]->Divide(2, 3);
  ccan[ican]->cd(1);
  SetRecoStyle(h1_efficiency_HFpt);
  h1_efficiency_HFpt->Draw("PE1SAME");
  h1_efficiency_HFpt->SetTitle("; p_{T}^{HF} [GeV]; Efficiency");
  h1_efficiency_HFpt->SetTitleOffset(1.2);
  h1_efficiency_HFpt->GetXaxis()->CenterTitle(false);
  h1_efficiency_HFpt->SetMaximum(1.05);
  h1_efficiency_HFpt->SetMinimum(0.0);

  // Tl.drawLatex(0.7, 0.8, "No cuts");

  ccan[ican]->cd(2);
  // gPad->SetLogy();
  // SetTruthStyle(h1_HF_pt);
  // SetRecoStyle(h1_HF_pt_reco);
  // h1_HF_pt_reco->Draw("PE1SAME");
  // h1_HF_pt->Draw("PE1SAME");
  SetRecoStyle(h1_efficiency_jetpt);
  h1_efficiency_jetpt->Draw("PE SAME");
  h1_efficiency_jetpt->SetTitle("; p_{T}^{jet} [GeV]; Efficiency");
  h1_efficiency_jetpt->SetTitleOffset(1.2);
  h1_efficiency_jetpt->SetMaximum(1.05);
  h1_efficiency_jetpt->SetMinimum(0.00);
  h1_efficiency_jetpt->GetXaxis()->CenterTitle(false);

  ccan[ican]->cd(3);
  SetTruthStyle(h1_HF_eta);
  SetRecoStyle(h1_HF_eta_reco);
  h1_HF_eta_reco->Draw("PE1SAME");
  h1_HF_eta_reco->SetTitle("; #eta; Efficiency");
  h1_HF_eta_reco->SetTitleOffset(1.2);
  h1_HF_eta_reco->GetXaxis()->CenterTitle(false);
  h1_HF_eta->Draw("PE1SAME");
  // h1_HF_rap_reco->Draw("PE1SAME");

  ccan[ican]->cd(4);
    auto legendmass = new TLegend(0.6, 0.6, 0.85, 0.9);
    legendmass->SetTextSize(0.04);
    legendmass->SetBorderSize(0);
    legendmass->SetFillStyle(0);
    legendmass->SetFillColor(3);
    
    
  SetDataStyle(h1_HF_mass_reco);
  SetRecoStyle(h1_HF_dtfmass_reco);
  h1_HF_dtfmass_reco->Draw("PE1SAME");
  h1_HF_mass_reco->Draw("PE1SAME");
    
    legendmass->AddEntry(h1_HF_dtfmass_reco, "DTF");
    legendmass->AddEntry(h1_HF_mass_reco, "M_{#mu #mu K}");
    
    legendmass->Draw("SAME");
    
//
  ccan[ican]->cd(5);

  auto legendpteta = new TLegend(0.45, 0.55, 0.75, 0.85);
  legendpteta->SetTextSize(0.03);
  legendpteta->SetBorderSize(0);
  legendpteta->SetFillStyle(0);
  legendpteta->SetFillColor(3);
  for (int i = 0; i < HFetabinsize; i++)
  {
    TH1D *h1_temp = (TH1D *)h2_efficiency_HFpteta->ProjectionX(Form("htemp%d_hfpteta", i), i + 1, i + 1);
    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(i + 1);
    h1_temp->SetLineColor(i + 1);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMaximum(1.05);
    h1_temp->SetMinimum(0.);
    h1_temp->SetTitle("; p_{T}^{HF} [GeV]; Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);

    legendpteta->AddEntry(h1_temp, Form("%.1f < #eta < %.1f GeV", HFeta_binedges[i], HFeta_binedges[i + 1]));
  }
  legendpteta->Draw("SAME");
  // h2_HF_pteta_ratio->Draw("COLZ TEXT");

  ccan[ican]->cd(6);
//
  auto legendptjetpt = new TLegend(0.45, 0.55, 0.75, 0.85);
  legendptjetpt->SetTextSize(0.03);
  legendptjetpt->SetBorderSize(0);
  legendptjetpt->SetFillStyle(0);
  legendptjetpt->SetFillColor(3);
  for (int i = 0; i < customptbinsize; i++)
  {
    TH1D *h1_temp = (TH1D *)h2_efficiency_HFptjetpt->ProjectionX(Form("htemp%d_hfptjetpt", i), i + 1, i + 1);
    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerSize(0.3);
    h1_temp->SetMarkerColor(i + 1);
    h1_temp->SetLineColor(i + 1);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMaximum(1.05);
    h1_temp->SetMinimum(0.);
    h1_temp->SetTitle("; p_{T}^{HF} [GeV]; Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);

    legendptjetpt->AddEntry(h1_temp, Form("%.1f < p_{T}^{jet} < %.1f GeV", custompt_binedges[i], custompt_binedges[i + 1]));
  }
  legendptjetpt->Draw("SAME");
  // h2_HF_pteta_ratio->Draw("COLZ TEXT");

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
  //

  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 1);
  ccan[ican]->SetFillColor(10);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.12);
  gPad->SetRightMargin(0.15);
  ccan[ican]->cd();
  ccan[ican]->Divide(2, 3);
  ccan[ican]->cd(1);

  NormalizeHist(h1_HF_rap_reco);
  NormalizeHist(h1_HF_eta_reco);

  SetRecoStyle(h1_HF_rap_reco);
  SetTruthStyle(h1_HF_eta_reco);
  h1_HF_rap_reco->Draw("PE1SAME");
  h1_HF_eta_reco->Draw("PE1SAME");

  ccan[ican]->cd(2);
  SetRecoStyle(h1_HF_IPCHI2);
  h1_HF_IPCHI2->Draw("PE SAME");

  ccan[ican]->cd(3);
  // gPad->SetLogz();
  // h2_jetpt_HFrap->Draw("COLZ");
  SetRecoStyle(h1_HF_LOGIPCHI2);
  h1_HF_LOGIPCHI2->Draw("PE SAME");
  ccan[ican]->cd(4);
  h3_efficiency_HFptetajetpt->GetYaxis()->SetRange(2, 2);
  // h3_purity_HFptetajetpt->GetYaxis()->SetRange(2, 2);
  TH2D *h2_temp2 = (TH2D *)(h3_efficiency_HFptetajetpt->Project3D("zx_eff2"));
  // TH2D *h2_temp_pur = (TH2D *)(h3_purity_HFptetajetpt->Project3D("zx_pur"));
  for (int i = 0; i < customptbinsize; i++)
  {

    TH1D *h1_temp = (TH1D *)(h2_temp2->ProjectionX(Form("htemp%d_HFptetajetpt2", i), i + 1, i + 1));
    // TH1D *h1_temp_pur = (TH1D *)(h2_temp_pur->ProjectionX(Form("htemp_pur%d_HFptetajetpt", i), i + 1, i + 1))

    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(i + 1);
    h1_temp->SetLineColor(i + 1);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMinimum(0.);
    h1_temp->SetMaximum(1.05);
    h1_temp->SetTitle("; p_{T}^{HF} [GeV]; Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);
      
      
  }

  ccan[ican]->cd(5);
  h3_efficiency_HFptetajetpt->GetYaxis()->SetRange(3, 3);
  // h3_purity_HFptetajetpt->GetYaxis()->SetRange(2, 2);
  TH2D *h2_temp3 = (TH2D *)(h3_efficiency_HFptetajetpt->Project3D("zx_eff3"));
  // TH2D *h2_temp_pur = (TH2D *)(h3_purity_HFptetajetpt->Project3D("zx_pur"));
//    auto legendpteta = new TLegend(0.45, 0.55, 0.75, 0.85);
//    legendpteta->SetTextSize(0.03);
//    legendpteta->SetBorderSize(0);
//    legendpteta->SetFillStyle(0);
//    legendpteta->SetFillColor(3);
  for (int i = 0; i < customptbinsize; i++)
  {

    TH1D *h1_temp = (TH1D *)(h2_temp3->ProjectionX(Form("htemp%d_HFptetajetpt3", i), i + 1, i + 1));
    // TH1D *h1_temp_pur = (TH1D *)(h2_temp_pur->ProjectionX(Form("htemp_pur%d_HFptetajetpt", i), i + 1, i + 1))

    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(i + 1);
    h1_temp->SetLineColor(i + 1);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMinimum(0.);
    h1_temp->SetMaximum(1.05);
    h1_temp->SetTitle("; p_{T}^{HF} [GeV]; Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);
//    legendptnTracks->AddEntry(h1_temp, Form("%.1f < p_{T}^{Jet} < %.1f GeV", nTrack_binedges[i], nTrack_binedges[i + 1]));
      
      
  }

  ccan[ican]->cd(6);

  auto legendptnTracks = new TLegend(0.45, 0.55, 0.75, 0.85);
  legendptnTracks->SetTextSize(0.03);
  legendptnTracks->SetBorderSize(0);
  legendptnTracks->SetFillStyle(0);
  legendptnTracks->SetFillColor(3);
  for (int i = 0; i < nTracksbinsize; i++)
  {
    TH1D *h1_temp = (TH1D *)h2_efficiency_HFptnTracks->ProjectionX(Form("htemp%d_hfptnTracks", i), i + 1, i + 1);
    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(i + 1);
    h1_temp->SetLineColor(i + 1);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMaximum(1.05);
    h1_temp->SetMinimum(0.);
    h1_temp->SetTitle("; p_{T}^{HF} [GeV]; Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);

    legendptnTracks->AddEntry(h1_temp, Form("%.1f < nTracks < %.1f GeV", nTrack_binedges[i], nTrack_binedges[i + 1]));
  }
  legendptnTracks->Draw("SAME");
  // h2_HF_pteta_ratio->Draw("COLZ TEXT");

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

  f.Write();
  f.Close();
}
//
