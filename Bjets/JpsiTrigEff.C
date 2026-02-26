

#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"
#include "../Helpers.h"

using namespace std;
using namespace RooFit;

void JpsiTrigEff(int NumEvts_user = -1, int dataset = 91599,
                 bool isData = true,
                 bool L0MuonDiMuon = true)
{

  int NumEvts = NumEvts_user;
  int NumEvtsMC = NumEvts_user;
  followHardest = false;
  truthLevel = false;
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
  TString str_DTF = "";
  TString str_PID = "";
  TString str_charged = "";
  TString str_GS = "";
  TString str_L0 = "";

  TString str_year = "2016";

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

  if (DTF_cut)
    str_DTF = "_DTF";
  if (PID_cut)
    str_PID = "_PID";

  if (truthLevel)
    str_level = "truth";
  else
  {
    if (isData)
      str_level = "data";
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
  if (L0MuonDiMuon)
    str_L0 = "_L0MuonDiMuon";

  if (ghostCut)
    str_ghost = Form("_ghost_%.1f", ghostProb);

  TString str_tree;
  TString extension_reco;
  extension_reco = TString("jpsieff_") + "reco" + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + str_DTF + str_PID + str_GS + str_L0 + Form("_%d", dataset);

  TString extension, extension_read, extension_prefix;
  extension = TString("jpsieff_") + str_level + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + str_DTF + str_PID + str_GS + str_L0 + Form("_%d", dataset);

  cout << extension_reco << endl;

  // cout << "Total number of MC events = " << MCTree.fChain->GetEntries() << endl;
  // cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;

  TChain *JpsiRecoTree = new TChain("JpsiRecoTree", "");
  vector<int> vec_datasets;
  if (Mag == 0)
    vec_datasets = {61590, 71590, 81590};
  else if (Mag == 1)
    vec_datasets = {61591, 71591, 81591};
  else
    vec_datasets = {61590, 61591, 71590, 71591, 81590, 81591};


  TString extension_rootFiles = isData ? "./efficiencies/TrigEff/data/" : "./efficiencies/TrigEff/MC/";
  TString extension_tree = isData ? "../../root_files/Bjets/" : "../../root_files/BjetsMC/";
  TString extension_plots = isData ? "../../plots/Bjets/" : "../../plots/BjetsMC/";

/*
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
      extension_read = TString("tree_jpsieff_") + str_level + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + Form("_%d", vec_datasets[i]);
      cout << extension_read << endl;
      JpsiRecoTree->Add(extension_rootFiles + extension_read + ".root/JpsiRecoTree");
    }
  }
*/
 // else
 // {
    extension_read = TString("tree_jpsieff_") + str_level + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + Form("_%d", dataset);
    // extension_read = extension_prefix + extension_read;
    cout << extension_read << endl;
    JpsiRecoTree->Add(extension_tree + extension_read + ".root/JpsiRecoTree");
 // }
  std::cout << "FILE NAME" << std::endl;
  std::cout << extension_rootFiles + extension_read + ".root" << std::endl;
  if (NumEvts > JpsiRecoTree->GetEntries())
    NumEvts = JpsiRecoTree->GetEntries();
  if (NumEvts == -1)
    NumEvts = JpsiRecoTree->GetEntries();

  // int NumEvts = 0;
  // cout <<"Choose number of events (-1: All Events, or enter integer): ";
  // cin>> NumEvts;


  TFile f(extension_rootFiles + extension + ".root", "RECREATE");

  TH1D *h1_HF_pt_reco = new TH1D("h1_HF_pt_reco", "", ptJpsibinsize, ptJpsi_binedges);
  TH1D *h1_HF_eta_reco = new TH1D("h1_HF_eta_reco", "", 20, 2, 5);
  TH1D *h1_HF_rap_reco = new TH1D("h1_HF_rap_reco", "", 20, 2, 5);
  TH1D *h1_HF_mass_reco = new TH1D("h1_HF_mass_reco", "", 60, 5.18, 5.38);
  TH1D *h1_HF_dtfmass_reco = new TH1D("h1_HF_dtfmass_reco", "", 60, 5.18, 5.38);

  TH1D *h1_HF_pt = new TH1D("h1_HF_pt", "", ptJpsibinsize, ptJpsi_binedges);
  TH1D *h1_HF_eta = new TH1D("h1_HF_eta", "", Jpsietabinsize, Jpsieta_binedges);

  TH2F *h2_HF_pteta = new TH2F("h2_HF_pteta", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_HF_pteta_reco = new TH2F("h2_HF_pteta_reco", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_jetpt_jetrap = new TH2F("h2_jetpt_jetrap", "", ptbinsize, pt_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_jetpt_HFrap = new TH2F("h2_jetpt_HFrap", "", ptbinsize, pt_binedges, Jpsietabinsize, Jpsieta_binedges);

  TH1D *h1_num_efficiency_HFpt = new TH1D("num_efficiency_HFpt", "", ptJpsibinsize, ptJpsi_binedges);
  TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", ptbinsize, pt_binedges);
  TH1D *h1_efficiency_HFpt = new TH1D("efficiency_HFpt", "", ptJpsibinsize, ptJpsi_binedges);
  TH2F *h2_efficiency_Jpsiptrap = new TH2F("efficiency_Jpsiptrap", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);

  TH2F *h2_efficiency_ratio = (TH2F *)h2_efficiency_Jpsiptrap->Clone("h2_eff_ratio");

  // TH1D *h1_denom_efficiency_HFpt = new TH1D("denom_efficiency_HFpt", "", ptJpsibinsize, ptJpsi_binedges);
  TH1D *h1_num_efficiency_Jpsieta = new TH1D("num_efficiency_Jpsieta", "", Jpsietabinsize, Jpsieta_binedges);
  TH1D *h1_denom_efficiency_Jpsieta = new TH1D("denom_efficiency_Jpsieta", "", Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_num_efficiency_HFpteta = new TH2F("num_efficiency_HFpteta", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH3D *h3_num_efficiency_HFptetajetpt = new TH3D("num_efficiency_HFptetajetpt", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges, customptbinsize, custompt_binedges);
  TH2F *h2_num_efficiency_HFptjetpt = new TH2F("num_efficiency_HFptjetpt", "", ptJpsibinsize, ptJpsi_binedges, customptbinsize, custompt_binedges);
  TH2F *h2_num_efficiency_HFptnTracks = new TH2F("num_efficiency_HFptnTracks", "", ptJpsibinsize, ptJpsi_binedges, nTracksbinsize, nTrack_binedges);
  TH2F *h2_denom_efficiency_HFpteta = new TH2F("denom_efficiency_HFpteta", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);

  TH1D *h1_HF_IPCHI2 = new TH1D("h1_HF_IPCHI2", "", 60, -1, 1000);
  TH1D *h1_HF_LOGIPCHI2 = new TH1D("h1_HF_LOGIPCHI2", "", 60, -5, 10);

  TH2F *h2_TRIG = new TH2F("h2_TRIG", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_TIS = new TH2F("h2_TIS", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_TOS = new TH2F("h2_TOS", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);
  TH2F *h2_TISTOS = new TH2F("h2_TISTOS", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges);

  TH2F *h2_num_efficiency_jpsiptmass = new TH2F("num_efficiency_jpsiptmass", "", ptJpsibinsize, ptJpsi_binedges, 60, 3.0, 3.19);
  TH2F *h2_denom_efficiency_jpsiptmass = new TH2F("denom_efficiency_jpsiptmass", "", ptJpsibinsize, ptJpsi_binedges, 60, 3.0, 3.19);

  const int mJpsibinsize = 35;
  double mJpsi_binedges[mJpsibinsize + 1];
  for (int i = 0; i < mJpsibinsize + 1; i++)
  {
    double top = 3.19;
    double bottom = 3.0;
    double step = (top - bottom) / mJpsibinsize;
    mJpsi_binedges[i] = step * i + bottom;
    // cout << mJpsi_binedges[i] << endl;
  }

  TH3F *h3_num_efficiency_jpsiptrapmass = new TH3F("num_efficiency_jpsiptrapmass", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges, mJpsibinsize, mJpsi_binedges);
  TH3F *h3_denom_efficiency_jpsiptrapmass = new TH3F("denom_efficiency_jpsiptrapmass", "", ptJpsibinsize, ptJpsi_binedges, Jpsietabinsize, Jpsieta_binedges, mJpsibinsize, mJpsi_binedges);

  // TClonesArray *arr = new TClonesArray("TLorentzVector");
  double HF_px, HF_py, HF_pz, HF_e;
  double Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e;
  double K_px, K_py, K_pz, K_e;
  double mum_px, mum_py, mum_pz, mum_e;
  double mup_px, mup_py, mup_pz, mup_e;
  double tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
  double dtf_mass, dtf_chi2ndf, dtf_ctau, K_PIDK;
  int eventNumber, nTracks, NumHFHads;
  bool Hasbbbar;
  double mup_TRCHI2DOF, mup_IPCHI2, mup_TRGHOSTPROB;
  double mum_TRCHI2DOF, mum_IPCHI2, mum_TRGHOSTPROB;
  double K_TRCHI2DOF, K_IPCHI2, K_TRGHOSTPROB;
  bool K_ISLONG, mup_ISMUON, mum_ISMUON;

  double Jpsi_TRCHI2DOF, Jpsi_IPCHI2, Jpsi_TRGHOSTPROB, Jpsi_AMAXDOCA, Jpsi_BPVDIRA, Jpsi_BPVVDCHI2, Jpsi_VTXCHI2;
  double Bu_TRCHI2DOF, Bu_IPCHI2, Bu_TRGHOSTPROB, Bu_AMAXDOCA, Bu_BPVDIRA, Bu_BPVIPCHI2, Bu_VTXCHI2;

  bool mup_L0Muon_TIS, mup_L0DiMuon_TIS;
  bool mum_L0Muon_TIS, mum_L0DiMuon_TIS;
  bool Jpsi_L0Muon_TIS, Jpsi_L0DiMuon_TIS;
  bool L0Hlt1_TIS;
  bool Jpsi_Hlt1Phys_TIS;
  bool Jpsi_Hlt2Phys_TIS;
  bool mup_L0Muon_TOS, mup_L0DiMuon_TOS;
  bool mum_L0Muon_TOS, mum_L0DiMuon_TOS;
  bool Jpsi_L0Global_TIS, Jpsi_L0Global_TOS;
  bool Jpsi_L0Muon_TOS, Jpsi_L0DiMuon_TOS, Jpsi_Hlt1_TOS, Jpsi_Hlt2_TOS, Jpsi_Hlt2DiMuonJPsi_TOS;
  bool Jpsi_L0Photon_TIS, Jpsi_L0Hadron_TIS, Jpsi_L0Electron_TIS;
  bool Jpsi_Hlt2DiMuonDetachedJPsi_TOS, Jpsi_Hlt2DiMuonDetachedHeavy_TOS;
  bool Jpsi_Hlt2DiMuonJpsiHighPT_TOS;

  JpsiRecoTree->SetBranchAddress("eventNumber", &eventNumber);

  JpsiRecoTree->SetBranchAddress("mup_px", &mup_px);
  JpsiRecoTree->SetBranchAddress("mup_py", &mup_py);
  JpsiRecoTree->SetBranchAddress("mup_pz", &mup_pz);
  JpsiRecoTree->SetBranchAddress("mup_pe", &mup_e);
  JpsiRecoTree->SetBranchAddress("mup_ISMUON", &mup_ISMUON);

  JpsiRecoTree->SetBranchAddress("mum_px", &mum_px);
  JpsiRecoTree->SetBranchAddress("mum_py", &mum_py);
  JpsiRecoTree->SetBranchAddress("mum_pz", &mum_pz);
  JpsiRecoTree->SetBranchAddress("mum_pe", &mum_e);
  JpsiRecoTree->SetBranchAddress("mum_ISMUON", &mum_ISMUON);

  JpsiRecoTree->SetBranchAddress("Jpsi_px", &Jpsi_px);
  JpsiRecoTree->SetBranchAddress("Jpsi_py", &Jpsi_py);
  JpsiRecoTree->SetBranchAddress("Jpsi_pz", &Jpsi_pz);
  JpsiRecoTree->SetBranchAddress("Jpsi_pe", &Jpsi_e);

  JpsiRecoTree->SetBranchAddress("mup_TRCHI2DOF", &mup_TRCHI2DOF);
  JpsiRecoTree->SetBranchAddress("mup_IPCHI2", &mup_IPCHI2);
  JpsiRecoTree->SetBranchAddress("mup_TRGHOSTPROB", &mup_TRGHOSTPROB);

  JpsiRecoTree->SetBranchAddress("mum_TRCHI2DOF", &mum_TRCHI2DOF);
  JpsiRecoTree->SetBranchAddress("mum_IPCHI2", &mum_IPCHI2);
  JpsiRecoTree->SetBranchAddress("mum_TRGHOSTPROB", &mum_TRGHOSTPROB);

  JpsiRecoTree->SetBranchAddress("Jpsi_TRCHI2DOF", &Jpsi_TRCHI2DOF);
  JpsiRecoTree->SetBranchAddress("Jpsi_IPCHI2", &Jpsi_IPCHI2);
  JpsiRecoTree->SetBranchAddress("Jpsi_TRGHOSTPROB", &Jpsi_TRGHOSTPROB);
  JpsiRecoTree->SetBranchAddress("Jpsi_AMAXDOCA", &Jpsi_AMAXDOCA);
  JpsiRecoTree->SetBranchAddress("Jpsi_BPVDIRA", &Jpsi_BPVDIRA);
  JpsiRecoTree->SetBranchAddress("Jpsi_BPVVDCHI2", &Jpsi_BPVVDCHI2);
  JpsiRecoTree->SetBranchAddress("Jpsi_VTXCHI2", &Jpsi_VTXCHI2);

  JpsiRecoTree->SetBranchAddress("mup_L0Muon_TOS", &mup_L0Muon_TOS);
  JpsiRecoTree->SetBranchAddress("mup_L0DiMuon_TOS", &mup_L0DiMuon_TOS);

  JpsiRecoTree->SetBranchAddress("mum_L0Muon_TOS", &mum_L0Muon_TOS);
  JpsiRecoTree->SetBranchAddress("mum_L0DiMuon_TOS", &mum_L0DiMuon_TOS);

  JpsiRecoTree->SetBranchAddress("Jpsi_L0Muon_TOS", &Jpsi_L0Muon_TOS);
  JpsiRecoTree->SetBranchAddress("Jpsi_L0DiMuon_TOS", &Jpsi_L0DiMuon_TOS);

  JpsiRecoTree->SetBranchAddress("mup_L0Muon_TIS", &mup_L0Muon_TIS);
  JpsiRecoTree->SetBranchAddress("mup_L0DiMuon_TIS", &mup_L0DiMuon_TIS);

  JpsiRecoTree->SetBranchAddress("mum_L0Muon_TIS", &mum_L0Muon_TIS);
  JpsiRecoTree->SetBranchAddress("mum_L0DiMuon_TIS", &mum_L0DiMuon_TIS);

  JpsiRecoTree->SetBranchAddress("Jpsi_L0Muon_TIS", &Jpsi_L0Muon_TIS);
  JpsiRecoTree->SetBranchAddress("Jpsi_L0DiMuon_TIS", &Jpsi_L0DiMuon_TIS);

  JpsiRecoTree->SetBranchAddress("Jpsi_L0Electron_TIS", &Jpsi_L0Electron_TIS);
  JpsiRecoTree->SetBranchAddress("Jpsi_L0Photon_TIS", &Jpsi_L0Photon_TIS);
  JpsiRecoTree->SetBranchAddress("Jpsi_L0Hadron_TIS", &Jpsi_L0Hadron_TIS);

  JpsiRecoTree->SetBranchAddress("Jpsi_L0Global_TIS", &Jpsi_L0Global_TIS);
  JpsiRecoTree->SetBranchAddress("Jpsi_L0Global_TOS", &Jpsi_L0Global_TOS);

  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt1_TOS", &Jpsi_Hlt1_TOS);
  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt1Phys_TIS", &Jpsi_Hlt1Phys_TIS);
  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt2Phys_TIS", &Jpsi_Hlt2Phys_TIS);

  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt2DiMuonJpsiHighPT_TOS", &Jpsi_Hlt2DiMuonJpsiHighPT_TOS);
  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt2DiMuonJPsi_TOS", &Jpsi_Hlt2DiMuonJPsi_TOS);
  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsi_TOS", &Jpsi_Hlt2DiMuonDetachedJPsi_TOS);
  JpsiRecoTree->SetBranchAddress("Jpsi_Hlt2_TOS", &Jpsi_Hlt2_TOS);

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
  cout << "NumEvts = " << NumEvts << endl;
  for (int ev = 0; ev < NumEvts; ev++)
  {
    JpsiRecoTree->GetEntry(ev);
    // EventTree->GetEntry(ev);
    // cout<<recojet_e<<", ";
    if (ev % 100000 == 0)
      cout << "Executing event " << ev << endl;

    // cout << Jpsi_Hlt1_TOS << ", " << Jpsi_Hlt2_TOS << endl;

    // if (last_eventNum == Tree.eventNumber)
    //   continue;

    TLorentzVector HFjet, HFmeson, mup, mum, Kmeson, Jpsi;

    // mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
    // mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);
    Jpsi.SetPxPyPzE(Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e);
    // cout << Jpsi.Pt() << ", ";

    // bool TOS = (Jpsi_L0DiMuon_TOS || Jpsi_L0Muon_TOS) && Jpsi_Hlt1_TOS;
    bool TOS = (Jpsi_L0DiMuon_TOS || Jpsi_L0Muon_TOS) && Jpsi_Hlt1_TOS && (Jpsi_Hlt2DiMuonJPsi_TOS || Jpsi_Hlt2DiMuonJpsiHighPT_TOS || Jpsi_Hlt2DiMuonDetachedJPsi_TOS);
    bool TIS = (Jpsi_L0Electron_TIS || Jpsi_L0Hadron_TIS || Jpsi_L0Photon_TIS) && Jpsi_Hlt1Phys_TIS && Jpsi_Hlt2Phys_TIS;

    if (L0MuonDiMuon)
    {
      TOS = (Jpsi_L0DiMuon_TOS || Jpsi_L0Muon_TOS) && Jpsi_Hlt1_TOS && (Jpsi_Hlt2DiMuonJPsi_TOS || Jpsi_Hlt2DiMuonJpsiHighPT_TOS || Jpsi_Hlt2DiMuonDetachedJPsi_TOS);
    }
    else
    {
      TOS = (Jpsi_L0DiMuon_TOS) && Jpsi_Hlt1_TOS && (Jpsi_Hlt2DiMuonJPsi_TOS || Jpsi_Hlt2DiMuonJpsiHighPT_TOS || Jpsi_Hlt2DiMuonDetachedJPsi_TOS);
    }


    // if (isData)
    // {
    if (TIS)
    {
      if (TOS)
      {
        h3_num_efficiency_jpsiptrapmass->Fill(Jpsi.Pt(), Jpsi.Rapidity(), Jpsi.M());
        // h2_num_efficiency_HFpteta->Fill(Jpsi.Pt(), Jpsi.Rapidity());
        h2_num_efficiency_jpsiptmass->Fill(Jpsi.Pt(), Jpsi.M());
      }

      h3_denom_efficiency_jpsiptrapmass->Fill(Jpsi.Pt(), Jpsi.Rapidity(), Jpsi.M());
      // h2_denom_efficiency_HFpteta->Fill(Jpsi.Pt(), Jpsi.Rapidity());
      h2_denom_efficiency_jpsiptmass->Fill(Jpsi.Pt(), Jpsi.M());
    }
    // }

    if (TOS)
      h2_num_efficiency_HFpteta->Fill(Jpsi.Pt(), Jpsi.Rapidity());
    h2_denom_efficiency_HFpteta->Fill(Jpsi.Pt(), Jpsi.Rapidity());

    Jpsi = mup + mum;
    events++;
    //  }
  }

  vector<RooAddPdf *> vec_models;

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
  // gStyle->SetCanvasDefH(2000);
  // gStyle->SetCanvasDefW(3508);
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
  gStyle->SetTitleOffset(1., "Y");
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
  rootfile = extension_rootFiles + extension + extension + TString(".root");
  plotfile = extension_plots + extension + TString(".ps");
  plotfilePDF = extension_plots + extension + TString(".pdf");
  plotfileO = plotfilePDF + TString("(");
  plotfileC = plotfilePDF + TString("]");

  // Plot data and PDF overlaid
  ++ican;
  cout << "ican = " << ican << endl;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 1);
  ccan[ican]->SetFillColor(10);

  // gPad->SetBottomMargin(0.16);
  ccan[ican]->cd();
  // ccan[ican]->Divide(ptJpsibinsize, Jpsietabinsize, 0.0001, 0.0001);
  ccan[ican]->Divide(9, 9, 0.0001, 0.0001);

  ++ican;
  cout << "ican = " << ican << endl;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 1);
  ccan[ican]->SetFillColor(10);

  // gPad->SetBottomMargin(0.16);
  ccan[ican]->cd();
  // ccan[ican]->Divide(ptJpsibinsize, Jpsietabinsize, 0.0001, 0.0001);
  ccan[ican]->Divide(9, 9, 0.0001, 0.0001);
  // c->SaveAs("plots/"+extension+".pdf");
  // h1_mass->Scale(1./h1_mass->GetEntries());

  // c1->SaveAs(Form("plots_misc/Misc_%s"+extension+".pdf", file_name.c_str()));
  // cout<<"...outbase   = "<<outbase.Data()<<endl;
  // cout<<"...rootfile  = "<<rootfile.Data()<<endl;
  // cout<<"...plotfile  = "<<plotfile.Data()<<endl;
  RooRealVar JpsiMass("JpsiMass", "JpsiMass", 2.7, 3.3);

  // if (isData)
  // {

  // RooWorkspace *w = new RooWorkspace(Form("w%d", i), Form("workspace%d", i));

  // TFile *file_workspace = new TFile("hists/workspace_" + extension_reco + ".root", "READ");
  // cout << extension_reco << endl;
  // RooWorkspace *w = new RooWorkspace("workspace", "workspace");

  for (int j1 = 0; j1 < ptJpsibinsize; j1++)
  {
    for (int j = 0; j < Jpsietabinsize; j++)
    {
      int i = j1 * Jpsietabinsize + j;
      TFile *file_workspacenum = new TFile(Form("./efficiencies/TrigEff/MC/workspacenum%d_", i) + extension_reco + ".root", "READ");
      TFile *file_workspacedenom = new TFile(Form("./efficiencies/TrigEff/MC/workspacedenom%d_", i) + extension_reco + ".root", "READ");

      RooWorkspace *w_readnum = (RooWorkspace *)file_workspacenum->Get(Form("workspacenum%d", i));
      RooWorkspace *w_readdenom = (RooWorkspace *)file_workspacedenom->Get(Form("workspacedenom%d", i));

      if (w_readnum == NULL && isData)
        return;
      cout << " i = " << i << endl;
      // TH1D *h1_tempnum = (TH1D *)h2_num_efficiency_jpsiptmass->ProjectionY(Form("htempnum1%d_hfptmass", i), i + 1, i + 1);
      // TH1D *h1_tempdenom = (TH1D *)h2_denom_efficiency_jpsiptmass->ProjectionY(Form("htempdenom1%d_hfptmass", i), i + 1, i + 1);
      TH1D *h1_tempnum = (TH1D *)h3_num_efficiency_jpsiptrapmass->ProjectionZ(Form("htempnum1%d_hfptmass", i), j1 + 1, j1 + 1, j + 1, j + 1);
      TH1D *h1_tempdenom = (TH1D *)h3_denom_efficiency_jpsiptrapmass->ProjectionZ(Form("htempdenom1%d_hfptmass", i), j1 + 1, j1 + 1, j + 1, j + 1);

      RooDataHist Jpsimass_num("Jpsimass", "Jpsimass", RooArgList(JpsiMass), h1_tempnum, 1.);
      RooDataHist Jpsimass_denom("Jpsimass", "Jpsimass", RooArgList(JpsiMass), h1_tempdenom, 1.);

      JpsiMass.setBins(mJpsibinsize);

      RooRealVar *mu_num, *sigma_num, *alpha1_num, *alpha2_num, *p1_num, *p2_num;
      RooRealVar *mu_denom, *sigma_denom, *alpha1_denom, *alpha2_denom, *p1_denom, *p2_denom;

      mu_num = new RooRealVar("mu_num", "mu", 3.097, 3.09, 3.1);
      sigma_num = new RooRealVar("sigma_num", "sigma", 0.01, 0.0005, 0.1);
      alpha1_num = new RooRealVar("alpha1_num", "alpha1", 2., 1., 8.); // CHANGE
      alpha2_num = new RooRealVar("alpha2_num", "alpha2", 2., 1., 8.); // CHANGE
      p1_num = new RooRealVar("p1_num", "p1", 1., 1., 8.);             // CHANGE
      p2_num = new RooRealVar("p2_num", "p2", 1., 1., 8.);

      mu_denom = new RooRealVar("mu_denom", "mu", 3.097, 3.09, 3.12);
      sigma_denom = new RooRealVar("sigma_denom", "sigma", 0.01, 0.0005, 0.1);
      alpha1_denom = new RooRealVar("alpha1_denom", "alpha1", 2., 1, 8.); // CHANGE
      alpha2_denom = new RooRealVar("alpha2_denom", "alpha2", 2., 1, 8.); // CHANGE
      p1_denom = new RooRealVar("p1_denom", "p1", 1., 1, 8.);             // CHANGE
      p2_denom = new RooRealVar("p2_denom", "p2", 1., 1, 8.);
      if (isData && w_readnum != NULL)
      {
        cout << "w_readnum not null" << endl;
        mu_num = (RooRealVar *)w_readnum->obj("mu_num");
        sigma_num = (RooRealVar *)w_readnum->obj("sigma_num");
        alpha1_num = (RooRealVar *)w_readnum->obj("alpha1_num");
        alpha2_num = (RooRealVar *)w_readnum->obj("alpha2_num");
        p1_num = (RooRealVar *)w_readnum->obj("p1_num");
        p2_num = (RooRealVar *)w_readnum->obj("p2_num");
        alpha1_num->setConstant(kTRUE);
        alpha2_num->setConstant(kTRUE);
        p1_num->setConstant(kTRUE);
        p2_num->setConstant(kTRUE);

        cout << "alpha1 = " << alpha1_num->getVal() << endl;

        mu_denom = (RooRealVar *)w_readdenom->obj("mu_denom");
        sigma_denom = (RooRealVar *)w_readdenom->obj("sigma_denom");

        alpha1_denom = (RooRealVar *)w_readdenom->obj("alpha1_denom");
        alpha2_denom = (RooRealVar *)w_readdenom->obj("alpha2_denom");
        p1_denom = (RooRealVar *)w_readdenom->obj("p1_denom");
        p2_denom = (RooRealVar *)w_readdenom->obj("p2_denom");
        alpha1_denom->setConstant(kTRUE);
        alpha2_denom->setConstant(kTRUE);
        p1_denom->setConstant(kTRUE);
        p2_denom->setConstant(kTRUE);
      }

      RooCrystalBall dcbPdf_num("dcbPdf_num", "DoubleSidedCB_num", JpsiMass, *mu_num, *sigma_num, *alpha1_num, *p1_num, *alpha2_num, *p2_num);
      RooCrystalBall dcbPdf_denom("dcbPdf_denom", "DoubleSidedCB_denom", JpsiMass, *mu_denom, *sigma_denom, *alpha1_denom, *p1_denom, *alpha2_denom, *p2_denom);

      // RooCBShape dcbPdf_num("dcbPdf_num", "CB_num", JpsiMass, mu_num, sigma_num, alpha1_num, p1_num);
      // RooCBShape dcbPdf_denom("dcbPdf_denom", "CB_denom", JpsiMass, mu_denom, sigma_denom, alpha1_denom, p1_denom);

      RooRealVar *a1_num = new RooRealVar("exp_c_num", "exp_c_num", -2., -10, -0.001);
      RooRealVar *a1_denom = new RooRealVar("exp_c_denom", "exp_c_denom", -2., -10, -0.001);

      RooExponential bkg_num("bkg_num", "bkg_num", JpsiMass, *a1_num);
      RooExponential bkg_denom("bkg_denom", "bkg_denom", JpsiMass, *a1_denom);

      RooRealVar *nbkg_num = new RooRealVar("nbkg_num", "fraction of background", 1000, 0., 1000000);
      RooRealVar *nbkg_denom = new RooRealVar("nbkg_denom", "fraction of background", 1000, 0., 1000000);

      RooRealVar *nsig_num = new RooRealVar("nsig_num", "fraction of background", 1000, 0., 1000000);
      RooRealVar *nsig_denom = new RooRealVar("nsig_denom", "fraction of background", 1000, 0., 1000000);

      RooAddPdf *model_num;
      RooAddPdf *model_denom;

      if (isData)
      {
        model_num = new RooAddPdf("model_num", "g1+g2+a", RooArgList(dcbPdf_num, bkg_num), RooArgList(*nsig_num, *nbkg_num));
        model_denom = new RooAddPdf("model_denom", "g1+g2+a", RooArgList(dcbPdf_denom, bkg_denom), RooArgList(*nsig_denom, *nbkg_denom));
      }
      else
      {
        model_num = new RooAddPdf("model_num", "g1+g2+a", RooArgList(dcbPdf_num), RooArgList(*nsig_num));
        model_denom = new RooAddPdf("model_denom", "g1+g2+a", RooArgList(dcbPdf_denom), RooArgList(*nsig_denom));
      }
      // model_num = new RooAddPdf("model_num", "g1+g2+a", RooArgList(dcbPdf_num, bkg_num), RooArgList(nsig_num, nbkg_num));
      // model_denom = new RooAddPdf("model_denom", "g1+g2+a", RooArgList(dcbPdf_denom, bkg_denom), RooArgList(nsig_denom, nbkg_denom));
      // JpsiMass.setRange("comblow", 3., 3.03);
      // JpsiMass.setRange("combhigh", 3.17, 3.19);
      // Jpsimass_denom.setRange("denom_comblow", 3., 3.04);
      // Jpsimass_denom.setRange("denom_combhigh", 3.16, 3.19);

      // bkg_num.fitTo(Jpsimass_num, Range("comblow", "combhigh"));
      // a1_num.setConstant(kTRUE);

      // bkg_denom.fitTo(Jpsimass_denom, Range("comblow", "combhigh"));
      // a1_denom.setConstant(kTRUE);

      model_num->fitTo(Jpsimass_num, PrintEvalErrors(-1), Save(true));
      model_denom->fitTo(Jpsimass_denom, PrintEvalErrors(-1), Save(true));

      vec_models.push_back(model_num);

      double N_tisnottos = (nsig_denom->getVal() - nsig_num->getVal());
      double N_tis = nsig_denom->getVal();
      double N_tistos = nsig_num->getVal();
      double err_tistos = nsig_num->getError();
      double err_tis = nsig_denom->getError();
      double eff = nsig_num->getVal() / nsig_denom->getVal();
      double err_tisnottos = N_tisnottos * sqrt(pow(err_tistos / N_tistos, 2) + pow(err_tis / N_tis, 2));
      double eff_err = sqrt((((1 - eff) * (1 - eff) * err_tistos * err_tistos) + (eff * eff * err_tisnottos * err_tisnottos)) / (N_tis * N_tis));
      // double eff_err = eff * sqrt(pow(nsig_num->getError() / nsig_num->getVal(), 2) + pow(nsig_denom->getError() / nsig_denom->getVal(), 2));
      // h1_efficiency_HFpt->SetBinContent(i + 1, eff);
      // h1_efficiency_HFpt->SetBinError(i + 1, eff_err);

      h2_efficiency_Jpsiptrap->SetBinContent(j1 + 1, j + 1, eff);
      h2_efficiency_Jpsiptrap->SetBinError(j1 + 1, j + 1, eff_err);

      ican = 0;
      ccan[ican]->cd(i + 1);
      double binsize = (3.2 - 3.07) / 60;
      gPad->SetRightMargin(0.12);
      gPad->SetLeftMargin(0.12);
      ccan[ican]->SetFillStyle(4070);
      RooPlot *xframe = JpsiMass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
      xframe->SetTitleOffset(0.7, "Y");
      xframe->SetTitleOffset(0.9, "X");

      Jpsimass_num.plotOn(xframe, Name("Jpsi_mass"), RooFit::MarkerSize(0.5));

      vec_models[i]->plotOn(xframe, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1));
      vec_models[i]->plotOn(xframe, Name("signal_num"), Components(dcbPdf_num), LineStyle(kSolid), LineColor(kGreen), LineWidth(1));
      vec_models[i]->plotOn(xframe, Name("bkg_num"), Components(bkg_num), LineStyle(kDashed), LineColor(kBlue), LineWidth(1));
      xframe->Draw("SAME");

      ican = 1;
      ccan[ican]->cd(i + 1);
      gPad->SetRightMargin(0.12);
      gPad->SetLeftMargin(0.12);
      ccan[ican]->SetFillStyle(4070);
      RooPlot *xframe2 = JpsiMass.frame(Title(Form(";M_{#mu#muK} [GeV];Events/(%.1f MeV) ", binsize * 1000.)));
      xframe2->SetTitleOffset(0.7, "Y");
      xframe2->SetTitleOffset(0.9, "X");

      Jpsimass_denom.plotOn(xframe2, Name("Jpsi_mass"), MarkerSize(0.5));

      model_denom->plotOn(xframe2, Name("model"), LineStyle(kSolid), LineColor(kRed), LineWidth(1));
      model_denom->plotOn(xframe2, Name("signal_denom"), Components(dcbPdf_denom), LineStyle(kSolid), LineColor(kGreen), LineWidth(1));
      if (isData)
        model_denom->plotOn(xframe2, Name("bkg_denom"), Components(bkg_denom), LineStyle(kDashed), LineColor(kBlue), LineWidth(1));
      xframe2->Draw("SAME");

      RooWorkspace *wnum = new RooWorkspace(Form("workspacenum%d", i), Form("workspacenum%d", i));
      RooWorkspace *wdenom = new RooWorkspace(Form("workspacedenom%d", i), Form("workspacedenom%d", i));
      // w->import(sign);
      // w->import(bkgn);
      wnum->import(*model_num);
      // wnum->import(JpsiMass);
      // wnum->Print();
      wnum->writeToFile(Form(extension_rootFiles + "workspacenum%d_", i) + extension + ".root");

      wdenom->import(*model_denom);
      // wdenom->import(JpsiMass);
      // wdenom->Print();
      wdenom->writeToFile(Form(extension_rootFiles + "workspacedenom%d_", i) + extension + ".root");
    }
  }
  // }
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
  // cout << f1_num[0]->GetMaximum() << endl;
  TH2F *h2_efficiency_HFpteta = (TH2F *)h2_num_efficiency_HFpteta->Clone("h2_efficiency_HFpteta");
  h2_efficiency_HFpteta->Divide(h2_num_efficiency_HFpteta, h2_denom_efficiency_HFpteta, 1, 1, "B");

  // if (!isData)
  // {
  //   h2_efficiency_Jpsiptrap = (TH2F *)h2_num_efficiency_HFpteta->Clone("h2_efficiency_Jpsiptrap");
  //   h2_efficiency_Jpsiptrap->Divide(h2_num_efficiency_HFpteta, h2_denom_efficiency_HFpteta, 1, 1, "B");
  // }

  h2_efficiency_ratio->Divide(h2_efficiency_Jpsiptrap, h2_efficiency_HFpteta, 1, 1);

  // TH2F *h2_efficiency_HFpteta = (TH2F *)h2_TRIG->Clone("h2_efficiency_HFpteta");
  // h2_TIS->Multiply(h2_TOS);
  // h2_TIS->Divide(h2_TISTOS);
  // h2_efficiency_HFpteta->Divide(h2_TIS);
  // h3_efficiency_HFptetajetpt->Write("efficiency_HFptetajetpt");
  // h3_num_efficiency_HFptetajetpt->Write("num_efficiency_HFptetajetpt");
  // h3_denom_efficiency_HFptetajetpt->Write("denom_efficiency_HFptetajetpt");

  //////////////////////////////////////////

  cout
      << "Total number of events processed = " << events << endl;
  cout << "NumRecoHF = " << NumRecoHF << endl;

  cout << "Events blocked: " << endl;
  cout << "npvs = " << cut_npvs << endl;
  cout << "jetpt = " << cut_pt << endl;
  cout << "jeteta = " << cut_eta << endl;
  cout << "Dinjet = " << cut_HFdR << endl;
  // cout<<"hasHF = "<< cut_hasHF<<endl;

  // //
  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 1);
  ccan[ican]->SetFillColor(10);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  ccan[ican]->cd();
  // ccan[ican]->Divide(2, 2);
  ccan[ican]->cd(1);

  auto legendrap = new TLegend(0.45, 0.2, 0.75, 0.4);
  legendrap->SetTextSize(0.03);
  legendrap->SetBorderSize(0);
  legendrap->SetFillStyle(0);
  legendrap->SetFillColor(3);

  // h2_efficiency_Jpsiptrap->Draw("COLZ");
  for (int i = 0; i < Jpsietabinsize; i++)
  {
    // ccan[ican]->cd(i + 1);

    TH1D *h1_temp = (TH1D *)h2_efficiency_Jpsiptrap->ProjectionX(Form("htempeff%d_jpsiptrap", i), i + 1, i + 1);
    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(colors[i]);
    h1_temp->SetMarkerSize(1.);
    h1_temp->SetLineColor(colors[i]);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMaximum(1.05);
    h1_temp->SetMinimum(0.);
    if (isData)
      h1_temp->SetTitle("Trigger Eff. in Data; p_{T}^{J/#psi} [GeV]; Efficiency");
    else
      h1_temp->SetTitle("Trigger Eff. in MC; p_{T}^{J/#psi} [GeV]; Efficiency");

    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);

    legendrap->AddEntry(h1_temp, Form("%.1f < y^{J/#psi} < %.1f", Jpsieta_binedges[i], Jpsieta_binedges[i + 1]));
  }
  legendrap->Draw("SAME");

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
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  ccan[ican]->cd();
  // ccan[ican]->Divide(2, 2);
  ccan[ican]->cd(1);

  // h2_efficiency_Jpsiptrap->Draw("COLZ");
  for (int i = 0; i < Jpsietabinsize; i++)
  {
    // ccan[ican]->cd(i + 1);

    TH1D *h1_temp = (TH1D *)h2_efficiency_HFpteta->ProjectionX(Form("htempeff%d_HFpteta", i), i + 1, i + 1);
    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(colors[i]);
    h1_temp->SetMarkerSize(1.);
    h1_temp->SetLineColor(colors[i]);
    h1_temp->Draw("PE1SAME");
    h1_temp->SetMaximum(1.05);
    h1_temp->SetMinimum(0.);
    h1_temp->SetTitle("; p_{T}^{J/#psi} [GeV]; Trigger Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);

    // legendrap->AddEntry(h1_temp, Form("%.1f < y^{J/#psi} < %.1f", Jpsieta_binedges[i], Jpsieta_binedges[i + 1]));
  }
  legendrap->Draw("SAME");

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

  TString extension_trig_MC, extension_trig_Data;
  TH2D *h2_trigeff_Data;
  TH2D *h2_trigeff_MC;

  extension_trig_MC = "jpsieff_reco_ev_-1_b_PID" + str_L0 + Form("_%d", dataset);
  extension_trig_Data = "jpsieff_data_ev_-1_b_PID" + str_L0 + Form("_%d ", dataset);

  cout << extension_trig_MC << endl;

  TFile *file_trigeffMC;
  if (isData)
  {
    file_trigeffMC = new TFile("./efficiencies/TrigEff/MC/" + extension_trig_MC + ".root", "READ");
    // TFile file_trigeffData("hists/" + extension_trig_Data + ".root", "READ");

    // h2_trigeff_Data = (TH2D *)file_trigeffData.Get("efficiency_Jpsiptrap");
    h2_trigeff_MC = (TH2D *)file_trigeffMC->Get("efficiency_Jpsiptrap");
    TH2D *h2_trigeff_ratio = (TH2D *)h2_efficiency_Jpsiptrap->Clone("h2_trigeff_ratio");
    h2_trigeff_ratio->Divide(h2_trigeff_MC);
    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 1);
    ccan[ican]->SetFillColor(10);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.15);
    ccan[ican]->cd();
    // ccan[ican]->Divide(2, 2);
    ccan[ican]->cd(1);

    // h2_efficiency_Jpsiptrap->Draw("COLZ");
    for (int i = 0; i < Jpsietabinsize; i++)
    {
      // ccan[ican]->cd(i + 1);

      TH1D *h1_temp = (TH1D *)h2_trigeff_ratio->ProjectionX(Form("htempeff%d_dataMC", i), i + 1, i + 1);
      h1_temp->SetMarkerStyle(i + 20);
      h1_temp->SetMarkerColor(colors[i]);
      h1_temp->SetMarkerSize(1.);
      h1_temp->SetLineColor(colors[i]);
      h1_temp->Draw("PE1SAME");
      h1_temp->SetMaximum(1.3);
      h1_temp->SetMinimum(0.8);
      h1_temp->SetTitle("; p_{T}^{J/#psi} [GeV]; Data/MC");
      h1_temp->SetTitleOffset(1.2);
      h1_temp->GetXaxis()->CenterTitle(false);

      // legendrap->AddEntry(h1_temp, Form("%.1f < y^{J/#psi} < %.1f", Jpsieta_binedges[i], Jpsieta_binedges[i + 1]));
    }
    legendrap->Draw("SAME");

    // Tl.drawLatex(0.7, 0.8, "No cuts");

    // gPad->SetLogy();

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
  }

  //
  if (!isData)
  {

    ++ican;
    sprintf(buf, "ccan%d", ican);
    ccan[ican] = new TCanvas(buf, buf, 1);
    ccan[ican]->SetFillColor(10);
    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);
    gPad->SetRightMargin(0.15);
    ccan[ican]->cd();
    // ccan[ican]->Divide(2, 2);
    ccan[ican]->cd(1);

    // h2_efficiency_Jpsiptrap->Draw("COLZ");
    for (int i = 0; i < Jpsietabinsize; i++)
    {
      // ccan[ican]->cd(i + 1);

      TH1D *h1_temp = (TH1D *)h2_efficiency_ratio->ProjectionX(Form("htempeff%d_ratio", i), i + 1, i + 1);
      h1_temp->SetMarkerStyle(i + 20);
      h1_temp->SetMarkerColor(colors[i]);
      h1_temp->SetMarkerSize(1.);
      h1_temp->SetLineColor(colors[i]);
      h1_temp->Draw("PE1SAME");
      h1_temp->SetMaximum(1.1);
      h1_temp->SetMinimum(0.8);
      h1_temp->SetTitle("; p_{T}^{J/#psi} [GeV]; TISTOS/True MC");
      h1_temp->SetTitleOffset(1.2);
      h1_temp->GetXaxis()->CenterTitle(false);

      // legendrap->AddEntry(h1_temp, Form("%.1f < y^{J/#psi} < %.1f", Jpsieta_binedges[i], Jpsieta_binedges[i + 1]));
    }
    legendrap->Draw("SAME");

    // Tl.drawLatex(0.7, 0.8, "No cuts");

    // gPad->SetLogy();

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
  }
  //

  //
  ++ican;
  sprintf(buf, "ccan%d", ican);
  ccan[ican] = new TCanvas(buf, buf, 1);
  ccan[ican]->SetFillColor(10);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.16);
  gPad->SetRightMargin(0.15);
  ccan[ican]->cd();
  ccan[ican]->Divide(3, 3);

  // gPad->SetLogy();

  for (int i = 0; i < ptJpsibinsize; i++)
  {
    ccan[ican]->cd(i + 1);

    TH1D *h1_temp = (TH1D *)h2_num_efficiency_jpsiptmass->ProjectionY(Form("htempnum%d_hfptmass", i), i + 1, i + 1);
    h1_temp->SetMarkerStyle(i + 20);
    h1_temp->SetMarkerColor(colors[i]);
    h1_temp->SetMarkerSize(1.);
    h1_temp->SetLineColor(colors[i]);
    h1_temp->Draw("PE1SAME");
    // h1_temp->SetMaximum(1.05);
    // h1_temp->SetMinimum(0.);
    h1_temp->SetTitle("; p_{T}^{J/#psi} [GeV]; Efficiency");
    h1_temp->SetTitleOffset(1.2);
    h1_temp->GetXaxis()->CenterTitle(false);

    // legendrap->AddEntry(h1_temp, Form("%.1f < y^{J/#psi} < %.1f", Jpsieta_binedges[i], Jpsieta_binedges[i + 1]));
  }
  // legendrap->Draw("SAME");

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

  if (ican > -1)
  {
    cout << " You plotted " << ican + 1 << " canvasses......." << endl;
    ccan[ican]->Print(plotfileC.Data());
  }
  std::cout << extension_rootFiles + extension_read + ".root" << std::endl;


  f.Write();
  f.Close();
}
//

