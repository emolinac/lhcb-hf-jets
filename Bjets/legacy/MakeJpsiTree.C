

// #include "AnalyzeDijetData_withsubsets.h"
// #include "LundPlaneData.h"

// #include "BMCDecayTree.C"
#include "JpsiTree.C"
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"
#include "../Helpers.h"

using namespace std;

void MakeJpsiTree(int NumEvts_user = -1, int dataset = 91599,
                  bool isData = true,
                  bool chargedJetCut_user = false)
{

  int NumEvts = NumEvts_user;
  int NumEvtsMC = NumEvts_user;
  followHardest = false;
  truthLevel = false;
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
  if (isData)
    str_level = "data";
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

  TString extension, extension_eff;
  extension = TString("tree_jpsieff_") + str_level + Form("_ev_%d", NumEvts) + str_Mag + str_flavor + Form("_%d", dataset);

  JpsiTree Tree(0, dataset, isData);

  cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;

  if (NumEvts == -1)
  {
    NumEvts = Tree.fChain->GetEntries();
  }

  // int NumEvts = 0;
  // cout <<"Choose number of events (-1: All Events, or enter integer): ";
  // cin>> NumEvts;

  cout << "Executing CAJetAlgo" << endl;


  TString extension_rootFiles = isData ? "../../root_files/Bjets/" : "../../root_files/BjetsMC/";
  TFile f(extension_rootFiles + extension + ".root", "RECREATE");

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
  bool Jpsi_L0Photon_TIS, Jpsi_L0Electron_TIS, Jpsi_L0Hadron_TIS;
  bool L0Hlt1_TIS;
  bool Jpsi_Hlt1Phys_TIS;
  bool Jpsi_Hlt2Phys_TIS;
  bool mup_L0Muon_TOS, mup_L0DiMuon_TOS;
  bool mum_L0Muon_TOS, mum_L0DiMuon_TOS;
  bool Jpsi_L0Global_TIS, Jpsi_L0Global_TOS;
  bool Jpsi_L0Muon_TOS, Jpsi_L0DiMuon_TOS, Jpsi_Hlt1_TOS, Jpsi_Hlt2_TOS, Jpsi_Hlt2DiMuonJPsi_TOS, Jpsi_Hlt2DiMuonJpsiHighPT_TOS;
  bool Jpsi_Hlt2DiMuonDetachedJPsi_TOS, Jpsi_Hlt2DiMuonDetachedHeavy_TOS;
  TTree *JpsiRecoTree = new TTree("JpsiRecoTree", "JpsiReco Tree Variables");

  JpsiRecoTree->Branch("eventNumber", &eventNumber);

  JpsiRecoTree->Branch("mup_px", &mup_px);
  JpsiRecoTree->Branch("mup_py", &mup_py);
  JpsiRecoTree->Branch("mup_pz", &mup_pz);
  JpsiRecoTree->Branch("mup_pe", &mup_e);
  JpsiRecoTree->Branch("mup_ISMUON", &mup_ISMUON);

  JpsiRecoTree->Branch("mum_px", &mum_px);
  JpsiRecoTree->Branch("mum_py", &mum_py);
  JpsiRecoTree->Branch("mum_pz", &mum_pz);
  JpsiRecoTree->Branch("mum_pe", &mum_e);
  JpsiRecoTree->Branch("mum_ISMUON", &mum_ISMUON);

  JpsiRecoTree->Branch("Jpsi_px", &Jpsi_px);
  JpsiRecoTree->Branch("Jpsi_py", &Jpsi_py);
  JpsiRecoTree->Branch("Jpsi_pz", &Jpsi_pz);
  JpsiRecoTree->Branch("Jpsi_pe", &Jpsi_e);

  JpsiRecoTree->Branch("mup_TRCHI2DOF", &mup_TRCHI2DOF);
  JpsiRecoTree->Branch("mup_IPCHI2", &mup_IPCHI2);
  JpsiRecoTree->Branch("mup_TRGHOSTPROB", &mup_TRGHOSTPROB);

  JpsiRecoTree->Branch("mum_TRCHI2DOF", &mum_TRCHI2DOF);
  JpsiRecoTree->Branch("mum_IPCHI2", &mum_IPCHI2);
  JpsiRecoTree->Branch("mum_TRGHOSTPROB", &mum_TRGHOSTPROB);

  JpsiRecoTree->Branch("Jpsi_TRCHI2DOF", &Jpsi_TRCHI2DOF);
  JpsiRecoTree->Branch("Jpsi_IPCHI2", &Jpsi_IPCHI2);
  JpsiRecoTree->Branch("Jpsi_TRGHOSTPROB", &Jpsi_TRGHOSTPROB);
  JpsiRecoTree->Branch("Jpsi_AMAXDOCA", &Jpsi_AMAXDOCA);
  JpsiRecoTree->Branch("Jpsi_BPVDIRA", &Jpsi_BPVDIRA);
  JpsiRecoTree->Branch("Jpsi_BPVVDCHI2", &Jpsi_BPVVDCHI2);
  JpsiRecoTree->Branch("Jpsi_VTXCHI2", &Jpsi_VTXCHI2);

  JpsiRecoTree->Branch("mup_L0Muon_TOS", &mup_L0Muon_TOS);
  JpsiRecoTree->Branch("mup_L0DiMuon_TOS", &mup_L0DiMuon_TOS);

  JpsiRecoTree->Branch("mum_L0Muon_TOS", &mum_L0Muon_TOS);
  JpsiRecoTree->Branch("mum_L0DiMuon_TOS", &mum_L0DiMuon_TOS);

  JpsiRecoTree->Branch("Jpsi_L0Muon_TOS", &Jpsi_L0Muon_TOS);
  JpsiRecoTree->Branch("Jpsi_L0DiMuon_TOS", &Jpsi_L0DiMuon_TOS);

  JpsiRecoTree->Branch("mup_L0Muon_TIS", &mup_L0Muon_TIS);
  JpsiRecoTree->Branch("mup_L0DiMuon_TIS", &mup_L0DiMuon_TIS);

  JpsiRecoTree->Branch("mum_L0Muon_TIS", &mum_L0Muon_TIS);
  JpsiRecoTree->Branch("mum_L0DiMuon_TIS", &mum_L0DiMuon_TIS);

  JpsiRecoTree->Branch("Jpsi_L0Muon_TIS", &Jpsi_L0Muon_TIS);
  JpsiRecoTree->Branch("Jpsi_L0Photon_TIS", &Jpsi_L0Photon_TIS);
  JpsiRecoTree->Branch("Jpsi_L0Hadron_TIS", &Jpsi_L0Hadron_TIS);
  JpsiRecoTree->Branch("Jpsi_L0Electron_TIS", &Jpsi_L0Electron_TIS);

  JpsiRecoTree->Branch("Jpsi_L0DiMuon_TIS", &Jpsi_L0DiMuon_TIS);

  JpsiRecoTree->Branch("Jpsi_L0Global_TIS", &Jpsi_L0Global_TIS);
  JpsiRecoTree->Branch("Jpsi_L0Global_TOS", &Jpsi_L0Global_TOS);

  JpsiRecoTree->Branch("Jpsi_Hlt1_TOS", &Jpsi_Hlt1_TOS);
  JpsiRecoTree->Branch("Jpsi_Hlt1Phys_TIS", &Jpsi_Hlt1Phys_TIS);
  JpsiRecoTree->Branch("Jpsi_Hlt2Phys_TIS", &Jpsi_Hlt2Phys_TIS);

  JpsiRecoTree->Branch("Jpsi_Hlt2DiMuonJpsiHighPT_TOS", &Jpsi_Hlt2DiMuonJpsiHighPT_TOS);
  JpsiRecoTree->Branch("Jpsi_Hlt2DiMuonJPsi_TOS", &Jpsi_Hlt2DiMuonJPsi_TOS);
  JpsiRecoTree->Branch("Jpsi_Hlt2DiMuonDetachedJPsi_TOS", &Jpsi_Hlt2DiMuonDetachedJPsi_TOS);
  JpsiRecoTree->Branch("Jpsi_Hlt2_TOS", &Jpsi_Hlt2_TOS);

  // TClonesArray *arr = new TClonesArray("TLorentzVector");

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
    if (ev % 1000000 == 0)
      cout << "Executing event " << ev << endl;
    TLorentzVector HFmeson2, mup, mum, Jpsi;

    Jpsi.SetPxPyPzE(Tree.J_psi_1S0_PX / 1000.,
                    Tree.J_psi_1S0_PY / 1000.,
                    Tree.J_psi_1S0_PZ / 1000.,
                    Tree.J_psi_1S0_PE / 1000.);

    // cout << Jpsi.E() << endl;
    if (ev != 0)
    {

      if (last_eventNum == Tree.eventNumber)
      {
        // cout << "Skipping! " << Jpsi.E() << endl;
        continue;
      }
    }
    // mup_L0Muon_TIS = Tree.mup_L0MuonDecision_TIS;
    // mup_L0DiMuon_TIS = Tree.mup_L0DiMuonDecision_TIS;

    // mum_L0Muon_TIS = Tree.mum_L0MuonDecision_TIS;
    // mum_L0DiMuon_TIS = Tree.mum_L0DiMuonDecision_TIS;

    Jpsi_L0Muon_TIS = Tree.J_psi_1S_L0MuonDecision_TIS;
    Jpsi_L0Photon_TIS = Tree.J_psi_1S_L0PhotonDecision_TIS;
    Jpsi_L0Electron_TIS = Tree.J_psi_1S_L0ElectronDecision_TIS;
    Jpsi_L0Hadron_TIS = Tree.J_psi_1S_L0HadronDecision_TIS;
    Jpsi_L0DiMuon_TIS = Tree.J_psi_1S_L0DiMuonDecision_TIS;

    Jpsi_Hlt1Phys_TIS = Tree.J_psi_1S_Hlt1Phys_TIS;
    Jpsi_Hlt2Phys_TIS = Tree.J_psi_1S_Hlt2Phys_TIS;

    L0Hlt1_TIS = (Jpsi_L0Muon_TIS || Jpsi_L0DiMuon_TIS) && Jpsi_Hlt1Phys_TIS;

    // mup_L0Muon_TOS = Tree.mup_L0MuonDecision_TOS;
    // mup_L0DiMuon_TOS = Tree.mup_L0DiMuonDecision_TOS;

    // mum_L0Muon_TOS = Tree.mum_L0MuonDecision_TOS;
    // mum_L0DiMuon_TOS = Tree.mum_L0DiMuonDecision_TOS;

    Jpsi_L0Muon_TOS = Tree.J_psi_1S_L0MuonDecision_TOS;
    Jpsi_L0DiMuon_TOS = Tree.J_psi_1S_L0DiMuonDecision_TOS;

    Jpsi_L0Global_TIS = Tree.J_psi_1S_L0Global_TIS;
    Jpsi_L0Global_TOS = Tree.J_psi_1S_L0Global_TOS;

    Jpsi_Hlt1_TOS = Tree.J_psi_1S_Hlt1DiMuonHighMassDecision_TOS;
    Jpsi_Hlt2DiMuonJpsiHighPT_TOS = Tree.J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS;
    Jpsi_Hlt2DiMuonJPsi_TOS = Tree.J_psi_1S_Hlt2DiMuonJPsiDecision_TOS;
    Jpsi_Hlt2DiMuonDetachedJPsi_TOS = Tree.J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TOS;
    // cout << Jpsi_Hlt2DiMuonJpsiHighPT_TOS << ", " << Jpsi_Hlt2DiMuonJPsi_TOS << ", " << Jpsi_Hlt2DiMuonDetachedJPsi_TOS << endl;
    // Jpsi_Hlt2DiMuonDetachedHeavy_TOS = Tree.J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TOS;
    // bool TIS = (Jpsi_L0DiMuon_TIS || Jpsi_L0Muon_TIS) && (Jpsi_Hlt1Phys_TIS) && Jpsi_Hlt2Phys_TIS;
    bool TIS = (Jpsi_L0Global_TIS) && (Jpsi_Hlt1Phys_TIS) && Jpsi_Hlt2Phys_TIS;
    // cout << (Jpsi_L0DiMuon_TIS || Jpsi_L0Muon_TIS) << ", " << Jpsi_Hlt1Phys_TIS << ", " << Jpsi_Hlt2Phys_TIS << endl;
    // cout << TIS << endl;

    // if (!(Jpsi_L0Photon_TIS || Jpsi_L0Hadron_TIS || Jpsi_L0Electron_TIS || Jpsi_L0Muon_TIS || Jpsi_L0DiMuon_TIS) && Jpsi_L0Global_TIS)
    //   cout << Jpsi_L0Global_TIS << ", " << Jpsi_L0Photon_TIS << ", " << Jpsi_L0Hadron_TIS << ", " << Jpsi_L0Electron_TIS << ", " << Jpsi_L0Muon_TIS << ", " << Jpsi_L0DiMuon_TIS << endl;

    // if (!Jpsi_L0Global_TIS && ((Jpsi_L0Photon_TIS || Jpsi_L0Hadron_TIS || Jpsi_L0Electron_TIS || Jpsi_L0Muon_TIS || Jpsi_L0DiMuon_TIS)))
    //   cout << Jpsi_L0Global_TIS << ", " << Jpsi_L0Photon_TIS << ", " << Jpsi_L0Hadron_TIS << ", " << Jpsi_L0Electron_TIS << ", " << Jpsi_L0Muon_TIS << ", " << Jpsi_L0DiMuon_TIS << endl;

    if (!TIS && isData)
      continue;
    TLorentzVector HFmeson;
    double tr_jet_rap, tr_jet_pt;
    bool FromJet;
    // cout << "passed TIS";
    // cout << Tree.J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS << ", ";

    // Jpsi = mup + mum;
    Jpsi.SetPxPyPzE(Tree.J_psi_1S0_PX / 1000.,
                    Tree.J_psi_1S0_PY / 1000.,
                    Tree.J_psi_1S0_PZ / 1000.,
                    Tree.J_psi_1S0_PE / 1000.);
    // HFmeson = mup + mum + Kmeson;

    // if (!Tree.J_psi_1S0_Hlt2DiMuonJPsiDecision_TOS)
    //   cout << Jpsi.M() << endl;

    Jpsi_px = Jpsi.Px();
    Jpsi_py = Jpsi.Py();
    Jpsi_pz = Jpsi.Pz();
    Jpsi_e = Jpsi.E();

    // Jpsi_TRCHI2DOF = Tree.J_psi_1S_TRCHI2DOF;
    // Jpsi_IPCHI2 = Tree.J_psi_1S_IPCHI2;
    // Jpsi_TRGHOSTPROB = Tree.J_psi_1S_TRGHOSTPROB;
    // Jpsi_AMAXDOCA = Tree.J_psi_1S_AMAXDOCA;
    // Jpsi_BPVDIRA = Tree.J_psi_1S_BPVDIRA;
    // Jpsi_BPVVDCHI2 = Tree.J_psi_1S_BPVVDCHI2;
    // Jpsi_VTXCHI2 = Tree.J_psi_1S_VTXCHI2;

    nTracks = Tree.nTracks;

    last_eventNum = Tree.eventNumber;

    events++;
    JpsiRecoTree->Fill();
    //  }
  }
  cout << "Saved " << events << " events" << endl;

  f.Write();
  f.Close();
}
//
