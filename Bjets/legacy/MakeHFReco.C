#include "./BMCDecayTree.C"
#include "./BDecayTree.C"
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"
#include "../Helpers.h"

using namespace std;

void MakeHFReco(int NumEvts_user = -1, int dataset = 1510,
                bool chargedJetCut_user = false)
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

  TString str_charged = "";

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

  // if(ptRange == 0) str_pT = "pt15pt20";
  // else if(ptRange == 1) str_pT = "pt20pt50";
  // else if(ptRange == 2) str_pT = "pt50";

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
  
  TString extension_read, extension_RootFilesMC, extension;
  extension_RootFilesMC = TString("../../root_files/BjetsMC/");

  TString  extension_eff;
  extension = TString("tree_HFeff_") + str_level + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + Form("_%d", dataset);
    
  extension_eff = TString("efficiency_truth") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(ptMin), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_charged + str_Mag + str_flavor + Form("_%d", dataset);

    
//  cout << extension_eff << endl;
//  TFile *file_eff = new TFile(extension_eff + TString(".root"), "READ");

//  BMCDecayTree MCTree(0, dataset, false);
  BDecayTree Tree(0, dataset, false);

//  cout << "Total number of MC events = " << MCTree.fChain->GetEntries() << endl;
  cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;
//  if (NumEvtsMC == -1)
//  {
//    NumEvtsMC = MCTree.fChain->GetEntries();
//  }
  if (NumEvts == -1)
  {
    NumEvts = Tree.fChain->GetEntries();
  }

  // int NumEvts = 0;
  // cout <<"Choose number of events (-1: All Events, or enter integer): ";
  // cin>> NumEvts;
//
//  cout << "Executing CAJetAlgo" << endl;

  TFile f(extension_RootFilesMC + extension + ".root", "RECREATE");


  float HF_px, HF_py, HF_pz, HF_e;
  float Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e;
  float K_px, K_py, K_pz, K_e;
  float mum_px, mum_py, mum_pz, mum_e;
  float mup_px, mup_py, mup_pz, mup_e;
  float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
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

  TTree *HFRecoTree = new TTree("HFRecoTree", "HFReco Tree Variables");

  HFRecoTree->Branch("eventNumber", &eventNumber);

  HFRecoTree->Branch("mup_px", &mup_px);
  HFRecoTree->Branch("mup_py", &mup_py);
  HFRecoTree->Branch("mup_pz", &mup_pz);
  HFRecoTree->Branch("mup_pe", &mup_e);
  HFRecoTree->Branch("mup_ISMUON", &mup_ISMUON);

  HFRecoTree->Branch("mum_px", &mum_px);
  HFRecoTree->Branch("mum_py", &mum_py);
  HFRecoTree->Branch("mum_pz", &mum_pz);
  HFRecoTree->Branch("mum_pe", &mum_e);
  HFRecoTree->Branch("mum_ISMUON", &mum_ISMUON);

  HFRecoTree->Branch("Jpsi_px", &Jpsi_px);
  HFRecoTree->Branch("Jpsi_py", &Jpsi_py);
  HFRecoTree->Branch("Jpsi_pz", &Jpsi_pz);
  HFRecoTree->Branch("Jpsi_pe", &Jpsi_e);

  HFRecoTree->Branch("K_px", &K_px);
  HFRecoTree->Branch("K_py", &K_py);
  HFRecoTree->Branch("K_pz", &K_pz);
  HFRecoTree->Branch("K_pe", &K_e);
  HFRecoTree->Branch("K_ISLONG", &K_ISLONG);

  HFRecoTree->Branch("HF_px", &HF_px);
  HFRecoTree->Branch("HF_py", &HF_py);
  HFRecoTree->Branch("HF_pz", &HF_pz);
  HFRecoTree->Branch("HF_pe", &HF_e);

  HFRecoTree->Branch("tr_jet_px", &tr_jet_px);
  HFRecoTree->Branch("tr_jet_py", &tr_jet_py);
  HFRecoTree->Branch("tr_jet_pz", &tr_jet_pz);
  HFRecoTree->Branch("tr_jet_pe", &tr_jet_e);

  HFRecoTree->Branch("dtf_mass", &dtf_mass);
  HFRecoTree->Branch("dtf_chi2ndf", &dtf_chi2ndf);
  HFRecoTree->Branch("dtf_ctau", &dtf_ctau);
  HFRecoTree->Branch("K_PIDK", &K_PIDK);
  HFRecoTree->Branch("nTracks", &nTracks);
  HFRecoTree->Branch("Hasbbbar", &Hasbbbar);
  HFRecoTree->Branch("NumHFHads", &NumHFHads);
    
    HFRecoTree->Branch("mup_L0", &mup_L0);
    HFRecoTree->Branch("mum_L0", &mum_L0);
    HFRecoTree->Branch("jpsi_L0", &jpsi_L0);
    HFRecoTree->Branch("jpsi_Hlt1", &jpsi_Hlt1);
    HFRecoTree->Branch("jpsi_Hlt2", &jpsi_Hlt2);
    HFRecoTree->Branch("Trig", &Trig);
    HFRecoTree->Branch("TIS", &TIS);
    HFRecoTree->Branch("TOS", &TOS);

//  HFRecoTree->Branch("mup_TRCHI2DOF", &mup_TRCHI2DOF);
//  HFRecoTree->Branch("mup_IPCHI2", &mup_IPCHI2);
//  HFRecoTree->Branch("mup_TRGHOSTPROB", &mup_TRGHOSTPROB);
//  HFRecoTree->Branch("mum_TRCHI2DOF", &mum_TRCHI2DOF);
//  HFRecoTree->Branch("mum_IPCHI2", &mum_IPCHI2);
//  HFRecoTree->Branch("mum_TRGHOSTPROB", &mum_TRGHOSTPROB);
//  HFRecoTree->Branch("K_TRCHI2DOF", &K_TRCHI2DOF);
//  HFRecoTree->Branch("K_IPCHI2", &K_IPCHI2);
//  HFRecoTree->Branch("K_TRGHOSTPROB", &K_TRGHOSTPROB);
//  HFRecoTree->Branch("Jpsi_TRCHI2DOF", &Jpsi_TRCHI2DOF);
//  HFRecoTree->Branch("Jpsi_IPCHI2", &Jpsi_IPCHI2);
//  HFRecoTree->Branch("Jpsi_TRGHOSTPROB", &Jpsi_TRGHOSTPROB);
//  HFRecoTree->Branch("Jpsi_AMAXDOCA", &Jpsi_AMAXDOCA);
//  HFRecoTree->Branch("Jpsi_BPVDIRA", &Jpsi_BPVDIRA);
//  HFRecoTree->Branch("Jpsi_BPVVDCHI2", &Jpsi_BPVVDCHI2);
//  HFRecoTree->Branch("Jpsi_VTXCHI2", &Jpsi_VTXCHI2);
//  HFRecoTree->Branch("Bu_TRCHI2DOF", &Bu_TRCHI2DOF);
//  HFRecoTree->Branch("Bu_IPCHI2", &Bu_IPCHI2);
//  HFRecoTree->Branch("Bu_TRGHOSTPROB", &Bu_TRGHOSTPROB);
//  HFRecoTree->Branch("Bu_AMAXDOCA", &Bu_AMAXDOCA);
//  HFRecoTree->Branch("Bu_BPVDIRA", &Bu_BPVDIRA);
//  HFRecoTree->Branch("Bu_BPVIPCHI2", &Bu_BPVIPCHI2);
//  HFRecoTree->Branch("Bu_VTXCHI2", &Bu_VTXCHI2);
    

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
    Tree.GetEntry(ev);
    // EventTree->GetEntry(ev);
    // cout<<recojet_e<<", ";
    if (ev % 10000 == 0)
      cout << "Executing event " << ev << endl;

    if (ev != 0)
    {
      if (last_eventNum == Tree.eventNumber)
        continue;
    }
    TLorentzVector HFmeson;
    float tr_jet_rap, tr_jet_pt;
    bool FromJet;
    float leading_pT = 0;
    int n_maxpT_cand = -999;
    int n_maxpT_entry = -999;
    int n_HFpt_entry = -999;
    int HF_counter = 0;
    int ncand = Tree.nCandidate;

    if (Tree.totCandidates > 1)
    {
      for (int icand = 0; icand < Tree.totCandidates; icand++)
      {
        Tree.GetEntry(ev + icand - ncand);
        tr_jet_pt = Tree.Bu_MCJet_PT / 1000.;
        tr_jet_rap = Tree.Bu_MCJet_Rapidity;
        Hasbbbar = false;
        FromJet = tr_jet_pt > pTLow && tr_jet_rap > etaMin && tr_jet_rap < etaMax;

        if (Tree.Bu_hasb && Tree.Bu_hasbbar)
          Hasbbbar = true;

        if (!FromJet)
          continue;

        HFmeson.SetPxPyPzE(Tree.Bu_PX / 1000.,
                           Tree.Bu_PY / 1000.,
                           Tree.Bu_PZ / 1000.,
                           Tree.Bu_PE / 1000.);
        // if (Hasbbbar)
        //   cout << Tree.eventNumber << ", " << Tree.Bu_MCJet_PT / 1000. << ", " << HFmeson.Pt() << endl;
        if (HFmeson.Pt() > leading_pT)
        {
          // cout<<leading_pT<<", "<<cand_jetpt<<endl;
          leading_pT = HFmeson.Pt();
          n_maxpT_cand = Tree.nCandidate;
          n_maxpT_entry = ev + icand - ncand;
          // cout<<n_maxpT_entry<<endl;
        }
      }
    }
    Tree.GetEntry(n_maxpT_entry);

    TLorentzVector HFmeson2, mup, mum, Kmeson, Jpsi;

    HFmeson.SetPxPyPzE(Tree.Bu_PX / 1000.,
                       Tree.Bu_PY / 1000.,
                       Tree.Bu_PZ / 1000.,
                       Tree.Bu_PE / 1000.);
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
    Jpsi = mup + mum;
    // HFmeson = mup + mum + Kmeson;

    jpsi_L0 = Tree.Jpsi_L0MuonDecision_TOS || Tree.Jpsi_L0DiMuonDecision_TOS;
    mup_L0 = Tree.mup_L0MuonDecision_TOS || Tree.mup_L0DiMuonDecision_TOS;
    mum_L0 = Tree.mum_L0MuonDecision_TOS || Tree.mum_L0DiMuonDecision_TOS;
      
    dtf_mass = Tree.Bu_ConsBu_M[0] / 1000.;
    bool PID_cond = (Tree.K_PIDK > 0);
    bool mass_cond = (dtf_mass > 5.24 && dtf_mass < 5.31);
    bool DTF_cond = (Tree.Bu_ConsBu_chi2[0] / Tree.Bu_ConsBu_nDOF[0] < 9) && (Tree.Bu_ConsBu_ctau[0] > 0.3);
    bool eta_cond = HFmeson.Eta() > 2.0 && HFmeson.Eta() < 5.;
    tr_jet_pt = Tree.Bu_MCJet_PT / 1000.;
    tr_jet_rap = Tree.Bu_MCJet_Rapidity;
    Hasbbbar = false;
    FromJet = tr_jet_pt > 12.5 && tr_jet_rap > etaMin && tr_jet_rap < etaMax;


    // if (HFmeson.Pt() < 2)
    //   continue;
    // if (!mass_cond)
    //   continue;
//     if (!eta_cond)
//       continue;

    if (!FromJet)
      continue;

    // if (Tree.Bu_TRUEPT <= 0)
    //   continue;

    // if (Tree.Bu_MCJet_NumHFHads > 1)
    //   continue;
    // if (Tree.Bu_TRUEP_E <= 0)
    //   continue;

    // cout << Tree.Bu_MCHF_PT / 1000. << ", " << HFmeson.Pt() << endl;
    // cout << "new: " << Tree.eventNumber << endl;
    // cout << Tree.Bu_TRUEP_E / 1000. << ", " << HFmeson.E() << endl;

    // cout << Tree.mup_TRUEP_E / 1000. << ", " << mup.E() << endl;
    // cout << Tree.mup_TRUEP_X / 1000. << ", " << mup.Px() << endl;
    // cout << Tree.mup_TRUEP_Y / 1000. << ", " << mup.Py() << endl;
    // cout << Tree.mup_TRUEP_Z / 1000. << ", " << mup.Pz() << endl;

    // cout << Tree.mum_TRUEP_E / 1000. << ", " << mum.E() << endl;
    // cout << Tree.mum_TRUEP_X / 1000. << ", " << mum.Px() << endl;
    // cout << Tree.mum_TRUEP_Y / 1000. << ", " << mum.Py() << endl;
    // cout << Tree.mum_TRUEP_Z / 1000. << ", " << mum.Pz() << endl;

    // cout << Tree.K_TRUEP_E / 1000. << ", " << Kmeson.E() << endl;
    // cout << Tree.K_TRUEP_X / 1000. << ", " << Kmeson.Px() << endl;
    // cout << Tree.K_TRUEP_Y / 1000. << ", " << Kmeson.Py() << endl;
    // cout << Tree.K_TRUEP_Z / 1000. << ", " << Kmeson.Pz() << endl;
    // }

    // cout << Tree.Bu_MCJet_NumHFHads << endl;

    // if (tr_jet_pt > 40)
    //   cout << Tree.eventNumber << ",";
    HF_px = HFmeson.Px();
    HF_py = HFmeson.Py();
    HF_pz = HFmeson.Pz();
    HF_e = HFmeson.E();

    K_px = Kmeson.Px();
    K_py = Kmeson.Py();
    K_pz = Kmeson.Pz();
    K_e = Kmeson.E();

    mup_px = mup.Px();
    mup_py = mup.Py();
    mup_pz = mup.Pz();
    mup_e = mup.E();

    mum_px = mum.Px();
    mum_py = mum.Py();
    mum_pz = mum.Pz();
    mum_e = mum.E();

    Jpsi_px = Jpsi.Px();
    Jpsi_py = Jpsi.Py();
    Jpsi_pz = Jpsi.Pz();
    Jpsi_e = Jpsi.E();

    tr_jet_px = Tree.Bu_MCJet_Px / 1000.;
    tr_jet_py = Tree.Bu_MCJet_Py / 1000.;
    tr_jet_pz = Tree.Bu_MCJet_Pz / 1000.;
    tr_jet_e = Tree.Bu_MCJet_E / 1000.;

    K_PIDK = Tree.K_PIDK;
    dtf_chi2ndf = Tree.Bu_ConsBu_chi2[0] / Tree.Bu_ConsBu_nDOF[0];
    dtf_ctau = Tree.Bu_ConsBu_ctau[0];

      jpsi_Hlt1 = Tree.Jpsi_Hlt1DiMuonHighMassDecision_TOS;
      jpsi_Hlt2 = Tree.Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS;
      TIS = (Tree.Jpsi_L0Global_TIS && Tree.Jpsi_Hlt1Global_TIS && Tree.Jpsi_Hlt2Global_TIS);
      TOS = jpsi_L0 && jpsi_Hlt1 && jpsi_Hlt2;
      
      nTracks = Tree.nTracks;

//    mup_TRCHI2DOF = Tree.mup_TRCHI2DOF;
//      mup_IPCHI2 = Tree.mup_IPCHI2;
//    mup_TRGHOSTPROB = Tree.mup_TRGHOSTPROB;
//    mup_ISMUON = Tree.mup_isMuon;
//
//    mum_TRCHI2DOF = Tree.mum_TRCHI2DOF;
//    mum_IPCHI2 = Tree.mum_IPCHI2;
//    mum_TRGHOSTPROB = Tree.mum_TRGHOSTPROB;
//    mum_ISMUON = Tree.mum_isMuon;
//
//    K_TRCHI2DOF = Tree.K_TRCHI2DOF;
//    K_IPCHI2 = Tree.K_IPCHI2;
//    K_TRGHOSTPROB = Tree.K_TRGHOSTPROB;
//    K_ISLONG = (Tree.K_TRACK_Type == 3);
//
//    Jpsi_TRCHI2DOF = Tree.Jpsi_TRCHI2DOF;
//    Jpsi_IPCHI2 = Tree.Jpsi_IPCHI2;
//    Jpsi_TRGHOSTPROB = Tree.Jpsi_TRGHOSTPROB;
//    Jpsi_AMAXDOCA = Tree.Jpsi_AMAXDOCA;
//    Jpsi_BPVDIRA = Tree.Jpsi_BPVDIRA;
//    Jpsi_BPVVDCHI2 = Tree.Jpsi_BPVVDCHI2;
//    Jpsi_VTXCHI2 = Tree.Jpsi_VTXCHI2;
//
//    Bu_TRCHI2DOF = Tree.Bu_TRCHI2DOF;
//    Bu_IPCHI2 = Tree.Bu_IPCHI2;
//    Bu_TRGHOSTPROB = Tree.Bu_TRGHOSTPROB;
//    Bu_AMAXDOCA = Tree.Bu_AMAXDOCA;
//    Bu_BPVDIRA = Tree.Bu_BPVDIRA;
//    Bu_BPVIPCHI2 = Tree.Bu_BPVIPCHI2;
//    Bu_VTXCHI2 = Tree.Bu_VTXCHI2;
      
//      mup_TRCHI2DOF = Tree.mup_TRACK_CHI2NDOF;
//      mup_IPCHI2 = Tree.mup_TRACK_PCHI2;
//      mup_TRGHOSTPROB = Tree.mup_TRACK_GhostProb;
//      mup_ISMUON = Tree.mup_isMuon;
//
//      mum_TRCHI2DOF = Tree.mum_TRACK_CHI2NDOF;
//      mum_IPCHI2 = Tree.mum_TRACK_PCHI2;
//      mum_TRGHOSTPROB = Tree.mum_TRACK_GhostProb;
//      mum_ISMUON = Tree.mum_isMuon;
//
//      K_TRCHI2DOF = Tree.K_TRACK_CHI2NDOF;
//      K_IPCHI2 = Tree.K_TRACK_PCHI2;
//      K_TRGHOSTPROB = Tree.K_TRACK_GhostProb;
//      K_ISLONG = (Tree.K_TRACK_Type == 3);
////      Jpsi_TRCHI2DOF = Tree.Jpsi_TRACK_CHI2NDOF;
////      Jpsi_IPCHI2 = Tree.Jpsi_PCHI2;
////      Jpsi_TRGHOSTPROB = Tree.Jpsi_TRACK_GhostProb;
//      Jpsi_BPVDIRA = Tree.Jpsi_DIRA_OWNPV;
//      Jpsi_BPVVDCHI2 = Tree.Jpsi_IPCHI2_OWNPV;
////      Jpsi_AMAXDOCA = Tree.Jpsi_AMAXDOCA;
////      Jpsi_VTXCHI2 = Tree.Jpsi_VTXCHI2;
////      Bu_TRCHI2DOF = Tree.Bu_TRACK_CHI2NDOF;
////      Bu_IPCHI2 = Tree.Bu_PCHI2;
////      Bu_TRGHOSTPROB = Tree.Bu_TRACK_GhostProb;
////      Bu_AMAXDOCA = Tree.Bu_AMAXDOCA;
//      Bu_BPVDIRA = Tree.Bu_DIRA_OWNPV;
//      Bu_BPVIPCHI2 = Tree.Bu_IPCHI2_OWNPV;
////      Bu_VTXCHI2 = Tree.Bu_VTXCHI2;

    // cout << endl;
    // cout << Jpsi_BPVDIRA << ", " << Tree.Jpsi_DIRA_OWNPV << endl;
    // cout << Jpsi_BPVVDCHI2 << ", " << Tree.Jpsi_FDCHI2_OWNPV << endl;
    // cout << Jpsi_VTXCHI2 << ", " << Tree.Jpsi_FDCHI2_ORIVX << endl;
    // cout << Tree.K_TRACK_Type << endl;

    NumHFHads = Tree.Bu_MCJet_NumHFHads;
    last_eventNum = Tree.eventNumber;

    if (Tree.Bu_hasb && Tree.Bu_hasbbar)
      Hasbbbar = true;

    // if (Hasbbbar)
    //   cout << Tree.eventNumber << ", " << HFmeson.Pt() << ", " << Tree.Bu_MCHF_PT / 1000. << ", " << Tree.Bu_TRUEPT / 1000. << endl;

    // cout << Tree.eventNumber << ", " << Tree.totCandidates << ", " << Tree.Bu_PT << ", " << Tree.Bu_TRUEPT << ", " << Tree.Bu_MCHF_PT << endl;
    // if (Tree.Bu_TRUEP_E == 0)
    // {
    //   cout << "new: " << endl;
    //   cout << tr_jet_rap << ", " << HFmeson.Rapidity() << endl;
    //   cout << tr_jet_pt << ", " << HFmeson.Pt() << endl;
    //   continue;
    // }
    events++;
    HFRecoTree->Fill();
    //  }
  }

  f.Write();
  f.Close();
}
//
