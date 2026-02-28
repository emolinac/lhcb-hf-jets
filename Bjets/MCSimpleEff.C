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

void MCSimpleEff(int NumEvts = -1, int dataset = 91599,
                 bool chargedJetCut_user = false,
                 bool SubtractGS = false,
                 bool onlysim9 = false) {
    
    bool MCflag = true;
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

    if (flavor == 5)
    {
      mass_num = 4.2;
    }
    else if (flavor == 4)
    {
      mass_num = 1.25;
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
    TString str_year = "";
    TString str_GS = "";
    TString str_charged = "";

    if (SubtractGS)
      str_GS = "_GSsub";
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
//    else if (Mag == 9)
//      str_Mag = "_MU&MD";

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

    TString str_simversion = "";     
    if (onlysim9)
    {
      str_simversion = "_sim9";
    }    

    if (chargedJetCut)
      str_charged = "_charge";
    // TString str_trees[5];
    // str_trees[0] = "TaggedDijets/DecayTree";
    // str_trees[1] = "D0KPiJet/DecayTree";
    // str_trees[2] = "B0KPiJet/DecayTree";
    // str_trees[3] = "Jets/DecayTree";
    TString str_tree;
    
    TString extension_read, extension_RootFilesMC, extension;
    
    TString RecoHF = "OnlyOneRecoHF";

    // This should be the file from MCMakeVarTree
    
    extension = TString("efficiency_") + str_level + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + str_GS + str_simversion + Form("_%d", dataset);
    
//    extension_read = TString("MCtree_") + str_level + str_mag + Form("_ev_%d", NumEvtsTrue3) + str_flavor;
    extension_RootFilesMC = TString("../../root_files/BjetsMC/");

    extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + str_simversion + Form("_%d", dataset);

   TFile fread((output_folder + "ntuple_bjets_mc.root").c_str(), "READ");
   TTree *BTree = (TTree *)fread.Get("BTree");
    // TChain *BTree = new TChain("BTree", "");
    // vector<int> vec_datasets;
    // if (Mag ==0)
    //     /// MC MD : 2016, 2017, 2018
    //   vec_datasets = {61590, 71590, 81590};
    // else if (Mag == 1)
    //     /// MC MU : 2016, 2017, 2018
    //   vec_datasets = {61591, 71591, 81591};
    // else
    //     ///
    //   vec_datasets = {61590, 61591, 71590, 71591, 81590, 81591};

    // cout << "Loading ntuples: " << endl;
    // if (year == 9 && JetMeth != 9)
    // {
    //   for (int i = 0; i < vec_datasets.size(); i++)
    //   {
    //     Mag = (vec_datasets[i] / 1) % 10;
    //     if (Mag == 0)
    //       str_Mag = "_MD";
    //     else if (Mag == 1)
    //       str_Mag = "_MU";
    //     extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + str_simversion + Form("_%d", vec_datasets[i]);
    //     // extension_read = extension_prefix + extension_read;
    //     cout << extension_read << endl;
    //     BTree->Add(extension_RootFilesMC  + extension_read + ".root/BTree");
    //   }
    // }
    // else
    // {
    //   extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_charged + str_Mag + str_flavor + str_simversion + Form("_%d", dataset);
    //   // extension_read = extension_prefix + extension_read;
    //   cout << extension_read << endl;
    //   BTree->Add(extension_RootFilesMC  + extension_read + ".root/BTree");
    // }
    if (NumEvts > BTree->GetEntries())
      NumEvts = BTree->GetEntries();
    if (NumEvts == -1)
      NumEvts = BTree->GetEntries();
    cout << BTree->GetEntries() << endl;

    TFile f((output_folder + "bjets_efficiencies.root").c_str(), "RECREATE");


//    if (NumEvts > BTree->GetEntries()) {
//        NumEvts = BTree->GetEntries();
//    }
//
//    if (NumEvts == -1) {
//        NumEvts = BTree->GetEntries();
//    }
//
//    cout << BTree->GetEntries() << endl;
//    extension = TString("efficiency_") + str_level + Form("_ev_%d", NumEvts)  + str_flavor;
//    extension = TString("efficiency_") + str_level + str_mag + Form("_ev_%d", NumEvts)  + str_flavor;
//    extension = TString("efficiency_") + str_level + Form("_ev_%d", NumEvts)  + str_flavor + RecoHF;
//    TFile f(extension_RootFilesMC + extension + ".root", "RECREATE");

    //  Naming Convention:
    //  {$1}{$2}{$3}{$4}{$5}
    //  $1: Year: 2016 = 6, 2017 = 7, 2018 = 8, All = 9
    //  $2: jets or dijets: TaggedDijets = 1, else = 2
    //  $3: flavor: b = 5, c = 4, udsg = 1
    //  $4: pT: pt15pt20 = 0, pt20pt50 = 1, pt50 = 2, else = 3
    //  $5: Magnet: MD = 0, MU = 1, Both = 9

    // if(ptRange == 0) str_pT = "pt15pt20";
    // else if(ptRange == 1) str_pT = "pt20pt50";
    // else if(ptRange == 2) str_pT = "pt50";

//    cout << BTree->GetEntries() << endl;
    

    /// Integrated Observables
    
    //Reco 1D Observables
    TH1D *h1_z = new TH1D("z", "", zbinsize, z_binedges);
    TH1D *h1_jt = new TH1D("jt", "", jtbinsize, jt_binedges);
    TH1D *h1_r = new TH1D("r", "", rbinsize, r_binedges);
    
    //Truth 1D Observables
    TH1D *h1_z_truth = new TH1D("z_truth", "", zbinsize, z_binedges);
    TH1D *h1_jt_truth = new TH1D("jt_truth", "", jtbinsize, jt_binedges);
    TH1D *h1_r_truth = new TH1D("r_truth", "", rbinsize, r_binedges);
    
    //Reco 2D Observables (208 - 217)
    TH2D *h2_zjt = new TH2D("zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
    TH2D *h2_zr = new TH2D("zr", "", zbinsize, z_binedges,  rbinsize, r_binedges);
    TH2D *h2_jtr = new TH2D("jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
    
    //Truth 2D Observables
    TH2D *h2_zjt_truth = new TH2D("zjt_truth", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
    TH2D *h2_zr_truth = new TH2D("zr_truth", "", zbinsize, z_binedges,  rbinsize, r_binedges);
    TH2D *h2_jtr_truth = new TH2D("jtr_truth", "", jtbinsize, jt_binedges,  rbinsize, r_binedges);
    
    // 2D Truth-Reco Correspondence (219 - 224)
    TH2D *h2_truthreco_z = new TH2D("truthreco_z", ";Reco z; Truth z", zbinsize, z_binedges, zbinsize, z_binedges);
    TH2D *h2_truthreco_jt = new TH2D("truthreco_jt", ";Reco jT; Truth jT", jtbinsize, jt_binedges, jtbinsize, jt_binedges);
    TH2D *h2_truthreco_r = new TH2D("truthreco_r", ";Reco r; Truth r", rbinsize, r_binedges, rbinsize, r_binedges);

    // 1D jet_pt-Binned observables
    
    TH1D *h1_z_truth_b = new TH1D("z_truth_b", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_g = new TH1D("z_truth_g", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_12515 = new TH1D("z_truth_12515", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_1520 = new TH1D("z_truth_1520", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_2030 = new TH1D("z_truth_2030", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_3050 = new TH1D("z_truth_3050", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_50100 = new TH1D("z_truth_50100", "", zbinsize, z_binedges);
    TH1D *h1_z_truth_100150 = new TH1D("z_truth_100150", "", zbinsize, z_binedges);
    TH1D *h1_jt_truth_12515 = new TH1D("jt_truth_12515", "", jtbinsize, jt_binedges);
    TH1D *h1_jt_truth_1520 = new TH1D("jt_truth_1520", "", jtbinsize, jt_binedges);
    TH1D *h1_jt_truth_2030 = new TH1D("jt_truth_2030", "", jtbinsize, jt_binedges);
    TH1D *h1_jt_truth_3050 = new TH1D("jt_truth_3050", "", jtbinsize, jt_binedges);
    TH1D *h1_jt_truth_50100 = new TH1D("jt_truth_50100", "", jtbinsize, jt_binedges);
    TH1D *h1_jt_truth_100150 = new TH1D("jt_truth_100150", "", jtbinsize, jt_binedges);
    TH1D *h1_r_truth_12515 = new TH1D("r_truth_12515", "", rbinsize, r_binedges);
    TH1D *h1_r_truth_1520 = new TH1D("r_truth_1520", "", rbinsize, r_binedges);
    TH1D *h1_r_truth_2030 = new TH1D("r_truth_2030", "", rbinsize, r_binedges);
    TH1D *h1_r_truth_3050 = new TH1D("r_truth_3050", "", rbinsize, r_binedges);
    TH1D *h1_r_truth_50100 = new TH1D("r_truth_50100", "", rbinsize, r_binedges);
    TH1D *h1_r_truth_100150 = new TH1D("r_truth_100150", "", rbinsize, r_binedges);
    
    // 1D Num Efficiencies and Purities of Observables (226 - 235)
    TH1D *h1_num_purity_z = new TH1D("num_purity_z", "", zbinsize, z_binedges);
    TH1D *h1_num_purity_jt = new TH1D("num_purity_jt", "", jtbinsize, jt_binedges );
    TH1D *h1_num_purity_r = new TH1D("num_purity_r", "", rbinsize, r_binedges);
    TH1D *h1_num_efficiency_z = new TH1D("num_efficiency_z", "", zbinsize, z_binedges);
    TH1D *h1_num_efficiency_jt = new TH1D("num_efficiency_jt", "", jtbinsize, jt_binedges );
    TH1D *h1_num_efficiency_r = new TH1D("num_efficiency_r", "", rbinsize, r_binedges);
    
    // 1D Denom Efficiencies and Purities of Observables (237 - 246)
    TH1D *h1_denom_purity_z = new TH1D("denom_purity_z", "", zbinsize, z_binedges);
    TH1D *h1_denom_purity_jt = new TH1D("denom_purity_jt", "", jtbinsize, jt_binedges );
    TH1D *h1_denom_purity_r = new TH1D("denom_purity_r", "", rbinsize, r_binedges);
    TH1D *h1_denom_efficiency_z = new TH1D("denom_efficiency_z", "", zbinsize, z_binedges);
    TH1D *h1_denom_efficiency_jt = new TH1D("denom_efficiency_jt", "", jtbinsize, jt_binedges );
    TH1D *h1_denom_efficiency_r = new TH1D("denom_efficiency_r", "", rbinsize, r_binedges);

    // 1D Kinematic Efficiencies and Purities (248-251)
    TH1D *h1_num_purity_jetpt = new TH1D("num_purity_jetpt", "", longptbinsize, longpt_binedges);
    TH1D *h1_denom_purity_jetpt = new TH1D("denom_purity_jetpt", "", longptbinsize, longpt_binedges);
    TH1D *h1_num_purity_rap = new TH1D("num_purity_rap", "", etabinsize, eta_binedges);       /// Jet Rap or HF rap?
    TH1D *h1_denom_purity_rap = new TH1D("denom_purity_rap", "", etabinsize, eta_binedges);  /// Jet Rap or HF rap?
    
    // 2D Efficiencies and Purities Observables (257 - 273)
    TH2D *h2_num_purity_zjt = new TH2D("num_purity_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
    TH2D *h2_denom_purity_zjt = new TH2D("denom_purity_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
    
    TH2D *h2_num_purity_zr = new TH2D("num_purity_zr", "", zbinsize, z_binedges, rbinsize, r_binedges);
    TH2D *h2_denom_purity_zr = new TH2D("denom_purity_zr", "", zbinsize, z_binedges, rbinsize, r_binedges);
    
    TH2D *h2_num_purity_jtr = new TH2D("num_purity_jtr", "", jtbinsize, jt_binedges, rbinsize, r_binedges);
    TH2D *h2_denom_purity_jtr = new TH2D("denom_purity_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
    
    TH2D *h2_num_efficiency_zjt = new TH2D("num_efficiency_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
    TH2D *h2_denom_efficiency_zjt = new TH2D("denom_efficiency_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
    
    TH2D *h2_num_efficiency_zr = new TH2D("num_efficiency_zr", "",  zbinsize, z_binedges, rbinsize, r_binedges);
    TH2D *h2_denom_efficiency_zr = new TH2D("denom_efficiency_zr", "",  zbinsize, z_binedges, rbinsize, r_binedges);
    
    TH2D *h2_num_efficiency_jtr = new TH2D("num_efficiency_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
    TH2D *h2_denom_efficiency_jtr = new TH2D("denom_efficiency_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
    
    
    
    // 2D and 4D Observables response Matrices (286 - 292)
    RooUnfoldResponse *response_zjt = new RooUnfoldResponse( h2_zjt, h2_zjt_truth, "response_zjt" );
    RooUnfoldResponse *response_zr = new RooUnfoldResponse( h2_zr, h2_zr_truth, "respose_zr");
    RooUnfoldResponse *response_jtr = new RooUnfoldResponse( h2_jtr , h2_jtr_truth, "respose_jtr");
    RooUnfoldResponse *response_z = new RooUnfoldResponse(h1_denom_purity_z, h1_denom_efficiency_z, "response_z");
    RooUnfoldResponse *response_jt = new RooUnfoldResponse(h1_denom_purity_jt, h1_denom_efficiency_jt, "response_jt");
    RooUnfoldResponse *response_r = new RooUnfoldResponse(h1_denom_purity_r, h1_denom_efficiency_r, "response_r");
    
    // 1D Kinematic Efficiencies and Purities (299 - 306)
//    TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", longptbinsize, longpt_binedges);
    TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", ptbinsize, pt_binedges);
    
    TH1D *h1_denom_efficiency_jetpt = new TH1D("denom_efficiency_jetpt", "", ptbinsize, pt_binedges);
//    TH1D *h1_denom_efficiency_jetpt = new TH1D("denom_efficiency_jetpt", "", longptbinsize, longpt_binedges);
    
    TH1D *h1_num_efficiency_HFpt = new TH1D("num_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
    TH1D *h1_denom_efficiency_HFpt = new TH1D("denom_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
    TH1D *h1_num_efficiency_HFeta = new TH1D("num_efficiency_HFeta", "", etabinsize, eta_binedges);
    TH1D *h1_denom_efficiency_HFeta = new TH1D("denom_efficiency_HFeta", "", etabinsize, eta_binedges);
    TH1D *h1_num_efficiency_jetrap = new TH1D("num_efficiency_rap", "", etabinsize, eta_binedges);
    TH1D *h1_denom_efficiency_jetrap = new TH1D("denom_efficiency_rap", "", etabinsize, eta_binedges);
    
    // 2D jet_pt Efficiencies (308 - 309)
    TH2D *h2_num_efficiency_jetpteta = new TH2D("num_efficiency_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
    TH2D *h2_denom_efficiency_jetpteta = new TH2D("denom_efficiency_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
    
    // 2D HFpt-eta Efficiencies; 3D HFpt-eta-jet_pt Efficiencies (311 - 318)
    TH2D *h2_num_efficiency_HFpteta = new TH2D("num_efficiency_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
    TH2D *h2_denom_efficiency_HFpteta = new TH2D("denom_efficiency_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
    TH3D *h3_denom_efficiency_HFptetajetpt = new TH3D("denom_efficiency_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, customptbinsize, custompt_binedges);

    TH2D *h2_denom_efficiency_HFptjetpt = new TH2D("denom_efficiency_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
    TH2D *h2_denom_efficiency_HFptnTracks = new TH2D("denom_efficiency_HFptnTracks","",ptHFbinsize, ptHF_binedges, nTracksbinsize, nTrack_binedges);
    
    // Reco/Truth Jet pt 1D Hists; Jet Pt Response Matrix (320 - 324)
    TH1D *h1_jetpt = new TH1D("jetpt", "", ptbinsize, pt_binedges);
    TH1D *h1_jetpt_truth = new TH1D("jetpt_truth", "", ptbinsize, pt_binedges);
    TH1D *h1_jetpt_test = new TH1D("jetpt_test", "", ptbinsize, pt_binedges);
    
    RooUnfoldResponse *response_jetpt = new RooUnfoldResponse(h1_jetpt, h1_jetpt_truth, "response_jetpt", ";Reco Jet pT; Truth Jet pT");
    
    TH2D *h2_truthreco_jetpt = new TH2D("truthreco_jetpt", ";Reco Jet pT; Truth Jet pT", longptbinsize, longpt_binedges, longptbinsize, longpt_binedges);
    
    // 3D Efficiencies/Purities of 2D Observables vs jet_pt; 6D Response Matrices (324 - 357)
    /// Q: Is Ibrahim's "matched obsv." is my "truth obsv." ?

    TH3D *h3_ptzjt = new TH3D("ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_truthptzjt = new TH3D("truthptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_num_efficiency_ptzjt = new TH3D("num_efficiency_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_denom_efficiency_ptzjt = new TH3D("denom_efficiency_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_efficiency_ptzjt = new TH3D("efficiency_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_num_purity_ptzjt = new TH3D("num_purity_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_denom_purity_ptzjt = new TH3D("denom_purity_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_purity_ptzjt = new TH3D("purity_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
    
    RooUnfoldResponse *response_ptzjt = new RooUnfoldResponse(h3_ptzjt, h3_truthptzjt, "response_ptzjt");
    
    TH3D *h3_ptzr = new TH3D("ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_truthptzr = new TH3D("truthptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_num_efficiency_ptzr = new TH3D("num_efficiency_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_denom_efficiency_ptzr = new TH3D("denom_efficiency_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_efficiency_ptzr = new TH3D("efficiency_ptzr", "", zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_num_purity_ptzr = new TH3D("num_purity_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_denom_purity_ptzr = new TH3D("denom_purity_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_purity_ptzr = new TH3D("purity_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );

    RooUnfoldResponse *response_ptzr = new RooUnfoldResponse(h3_ptzr, h3_truthptzr, "response_ptzr");
    
    TH3D *h3_ptjtr = new TH3D("ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_truthptjtr = new TH3D("truthptjtr", "", jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_num_efficiency_ptjtr = new TH3D("num_efficiency_ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_denom_efficiency_ptjtr = new TH3D("denom_efficiency_ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize,  pt_binedges );
    TH3D *h3_efficiency_ptjtr = new TH3D("efficiency_ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_num_purity_ptjtr = new TH3D("num_purity_ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_denom_purity_ptjtr = new TH3D("denom_purity_ptjtr", "", jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    TH3D *h3_purity_ptjtr = new TH3D("purity_ptjtr", "", jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D,  ptbinsize, pt_binedges );

    RooUnfoldResponse *response_ptjtr = new RooUnfoldResponse(h3_ptjtr, h3_truthptjtr, "response_ptjtr");
    
    // 2D Jet_pt vs Observable histograms, Purities/Efficiencies, Responses (359 - 368)
    TH2D *h2_ptz = new TH2D("ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_truthptz = new TH2D("truthptz", "", zbinsize, z_binedges, ptbinsize, pt_binedges );
    TH2D *h2_num_efficiency_ptz = new TH2D("num_efficiency_ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_denom_efficiency_ptz = new TH2D("denom_efficiency_ptz", "", zbinsize, z_binedges, ptbinsize, pt_binedges );
    TH2D *h2_efficiency_ptz = new TH2D("efficiency_ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_num_purity_ptz = new TH2D("num_purity_ptz", "", zbinsize, z_binedges, ptbinsize, pt_binedges );
    TH2D *h2_denom_purity_ptz = new TH2D("denom_purity_ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_purity_ptz = new TH2D("purity_ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
    
    RooUnfoldResponse *response_ptz = new RooUnfoldResponse(h2_ptz, h2_truthptz, "response_ptz");
    
    TH2D *h2_ptjt = new TH2D("ptjt", "", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_truthptjt = new TH2D("truthptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges );
    TH2D *h2_num_efficiency_ptjt = new TH2D("num_efficiency_ptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
    TH2D *h2_denom_efficiency_ptjt = new TH2D("denom_efficiency_ptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
    TH2D *h2_efficiency_ptjt = new TH2D("efficiency_ptjt", "", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_num_purity_ptjt = new TH2D("num_purity_ptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges );
    TH2D *h2_denom_purity_ptjt = new TH2D("denom_purity_ptjt", "", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_purity_ptjt = new TH2D("purity_ptjt", "", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
    
    RooUnfoldResponse *response_ptjt = new RooUnfoldResponse(h2_ptjt, h2_truthptjt, "response_ptjt");
   
    TH2D *h2_ptr = new TH2D("ptr", "", rbinsize, r_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_truthptr = new TH2D("truthptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges );
    TH2D *h2_num_efficiency_ptr = new TH2D("num_efficiency_ptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges);
    TH2D *h2_denom_efficiency_ptr = new TH2D("denom_efficiency_ptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges);
    TH2D *h2_efficiency_ptr = new TH2D("efficiency_ptr", "", rbinsize, r_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_num_purity_ptr = new TH2D("num_purity_ptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges );
    TH2D *h2_denom_purity_ptr = new TH2D("denom_purity_ptr", "", rbinsize, r_binedges,  ptbinsize, pt_binedges );
    TH2D *h2_purity_ptr = new TH2D("purity_ptr", "", rbinsize, r_binedges,  ptbinsize, pt_binedges );
    
    RooUnfoldResponse *response_ptr = new RooUnfoldResponse(h2_ptr, h2_truthptr, "response_ptr");
    
    TH2D *h2_SVTag_eff_denom = new TH2D("h2_SVTag_eff_denom","", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
    TH2D *h2_SVTag_eff_denom_z = new TH2D("h2_SVTag_eff_denom_z","", zbinsize, z_binedges, customptbinsize, custompt_binedges);
    
    //// Other Stuff
    TH1D *h1_jetpt_reco_ext = new TH1D("jetpt_reco_ext", "", 20, 0.0, 100.);
    TH1D *h1_jetpt_truth_ext = new TH1D("jetpt_truth_ext", "", 20, 0.0, 100.);
   
    TH1D *h1_HFpt_reco_ext = new TH1D("HFpt_reco_ext", "", 20, 0.0, 100.);
    TH1D *h1_HFpt_truth_ext = new TH1D("HFpt_truth_ext", "", 20, 0.0, 100.);

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


    int eventNum;
    int NumBjets = 0;
    int NumRecoBjets = 0;
//  double etaMin = 2.5;
//  double etaMax =  4.0;
//  double ptMin = 10.0;
//  double ptMax = 120.0;

    cout << "Requested # of events " << NumEvts << endl;

    for (int ev = 0; ev < NumEvts; ev++) 
    {
        BTree->GetEntry(ev);

        if (ev % 10000 == 0)
            cout << "Executing event: " << ev << endl;

        TLorentzVector HFmeson(HF_px, HF_py, HF_pz, HF_e);
        TLorentzVector meas_HFmeson(meas_HF_px, meas_HF_py, meas_HF_pz, meas_HF_e);
        TLorentzVector mup(mup_px, mup_py, mup_pz, mup_e);
        TLorentzVector mum(mum_px, mum_py, mum_pz, mum_e);

        // new
//    if(NumDtrRecoHF > 1) continue;
        
        bool rap_cond = (jet_rap > etaMin && jet_rap < etaMax);
        bool pt_cond = (jet_pt > pTLow);
        bool meas_rap_cond = (meas_jet_rap > etaMin && meas_jet_rap < etaMax);
        bool meas_pt_cond = (meas_jet_pt > pTLow);

//        if (NumHFHads > 1) {
//            continue;
//        }
//
//        if(Hasbbbar) continue;
        
        if (SubtractGS && Hasbbbar)
        {
          continue;
        }
        
        h1_jetpt_reco_ext->Fill(meas_jet_pt);
        h1_jetpt_truth_ext->Fill(jet_pt);

        h1_HFpt_reco_ext->Fill(meas_HFmeson.Pt());
        h1_HFpt_truth_ext->Fill(HFmeson.Pt());
    

        if (hasRecoHF && meas_pt_cond && meas_rap_cond)
        {
            response_jetpt->Fill(meas_jet_pt, jet_pt);
            
            response_z->Fill(meas_z, truth_z);
            response_jt->Fill(meas_jt, truth_jt);
            response_r->Fill(meas_r, truth_r);
            
            response_zjt->Fill(meas_z, meas_jt, truth_z, truth_jt);
            response_zr->Fill(meas_z, meas_r, truth_z, truth_r);
            response_jtr->Fill(meas_jt, meas_r, truth_jt, truth_r);
            
            response_ptjt->Fill( meas_jt,  meas_jet_pt, truth_jt, jet_pt);
            response_ptr->Fill( meas_r,  meas_jet_pt, truth_r, jet_pt);
            response_ptz->Fill( meas_z,  meas_jet_pt, truth_z, jet_pt);
            
            response_ptzjt->Fill( meas_z, meas_jt, meas_jet_pt, truth_z, truth_jt, jet_pt);
            response_ptzr->Fill( meas_z, meas_r, meas_jet_pt, truth_z, truth_r, jet_pt);
            response_ptjtr->Fill( meas_jt,  meas_r, meas_jet_pt, truth_jt, truth_r, jet_pt);
            
        }
        
//        if (jet_pt > ptMin && jet_pt < ptMax)
//        {
//         
//            
////            RooUnfoldResponse * = new RooUnfoldResponse( h2_rjt , h2_rjt_truth, "respose_rjt");
//            
//
//        }

        NumBjets++;
        
        if (pt_cond) {
            h2_denom_efficiency_jetpteta->Fill(jet_pt, jet_eta);
        }

        if (pt_cond && rap_cond) 
        {
            
            h1_denom_efficiency_HFpt->Fill(HF_pt);
            h1_denom_efficiency_HFeta->Fill(HFmeson.Eta());
            h1_denom_efficiency_jetrap->Fill(jet_rap);
            
            h2_denom_efficiency_HFpteta->Fill(HFmeson.Pt(), HFmeson.Rapidity());
            h3_denom_efficiency_HFptetajetpt->Fill(HFmeson.Pt(), HFmeson.Rapidity(), jet_pt);
            h1_denom_efficiency_jetpt->Fill(jet_pt);
            h2_denom_efficiency_HFptjetpt->Fill(HFmeson.Pt(), jet_pt);
            
//            h2_denom_efficiency_fragzjetpt->Fill(HFmeson.Pt() / jet_pt, jet_pt);
            h2_denom_efficiency_HFptnTracks->Fill(HFmeson.Pt(), nTracks);
            
            h2_denom_efficiency_ptz->Fill(truth_z, jet_pt);
            h2_denom_efficiency_ptjt->Fill(truth_jt, jet_pt);
            h2_denom_efficiency_ptr->Fill(truth_r, jet_pt);
    
            h3_denom_efficiency_ptzjt->Fill(truth_z, truth_jt, jet_pt);
            h3_denom_efficiency_ptzr->Fill(truth_z, truth_r, jet_pt);
            h3_denom_efficiency_ptjtr->Fill(truth_jt, truth_r, jet_pt);
                      
            h1_denom_efficiency_z->Fill(truth_z);
            h1_denom_efficiency_jt->Fill(truth_jt);
            h1_denom_efficiency_r->Fill(truth_r);
  
            h2_denom_efficiency_zjt->Fill(truth_z, truth_jt);
            h2_denom_efficiency_zr->Fill(truth_z, truth_r);
            h2_denom_efficiency_jtr->Fill(truth_jt, truth_r);            
       
            
        if (Hasbbbar)
        {
//          h2_SVTag_eff_denom_z->Fill(HFmeson.Pt() / jet_pt, jet_pt);
        
//          h2_SVTag_eff_denom->Fill(HFmeson.Pt(), tr_jet_pt);
          h2_SVTag_eff_denom_z->Fill(truth_z, jet_pt);
        }

      }
        
      ///////////////////////////////////////////////////
      // Fill the numerator of the efficiency & response matrix
      // Numerator = all reco B jets in fiducial region
      ////////////////////////////////////////////////////  
      if (jet_pt > ptMin && jet_pt < ptMax && rap_cond && meas_pt_cond && hasRecoHF)
      {
          
          h1_num_efficiency_HFpt->Fill(HF_pt);
          h1_num_efficiency_HFeta->Fill(HFmeson.Eta());
          h1_num_efficiency_jetrap->Fill(jet_rap);
          
          h2_num_efficiency_HFpteta->Fill(HFmeson.Pt(), HFmeson.Rapidity());
          h1_num_efficiency_jetpt->Fill(jet_pt);
          h2_num_efficiency_jetpteta->Fill(jet_pt, jet_rap);
          
          h1_num_efficiency_z->Fill(truth_z);
          h1_num_efficiency_jt->Fill(truth_jt);
          h1_num_efficiency_r->Fill(truth_r);

          h2_num_efficiency_zjt->Fill(truth_z, truth_jt);
          h2_num_efficiency_zr->Fill(truth_z, truth_r);
          h2_num_efficiency_jtr->Fill(truth_jt, truth_r);
          
          h2_num_efficiency_ptz->Fill(truth_z, jet_pt);
          h2_num_efficiency_ptjt->Fill(truth_jt, jet_pt);
          h2_num_efficiency_ptr->Fill(truth_r, jet_pt);
  
          h3_num_efficiency_ptzjt->Fill(truth_z, truth_jt, jet_pt);
          h3_num_efficiency_ptzr->Fill(truth_z, truth_r, jet_pt);
          h3_num_efficiency_ptjtr->Fill( truth_jt, truth_r, jet_pt);
          
      }


        if (meas_pt_cond)
        {
            h1_num_purity_jetpt->Fill(jet_pt);
            h1_num_purity_rap->Fill(jet_rap);
            
//            h1_num_purity_z->Fill(meas_z);
//            h1_num_purity_jt->Fill(meas_jt);
//            h1_num_purity_r->Fill(meas_r);
// 
//            h2_num_purity_zjt->Fill(meas_z, meas_jt);
//            h2_num_purity_zr->Fill(meas_z, meas_r);
//            h2_num_purity_jtr->Fill(meas_jt, meas_r);
//
//            h2_num_purity_ptz->Fill(meas_z, meas_jet_pt);
//            h2_num_purity_ptjt->Fill(meas_jt, meas_jet_pt);
//            h2_num_purity_ptr->Fill(meas_r, meas_jet_pt);
//            
//            h3_num_purity_ptzjt->Fill(meas_z, meas_jt, meas_jet_pt );
//            h3_num_purity_ptzr->Fill(meas_z, meas_r, meas_jet_pt);
//            h3_num_purity_ptjtr->Fill(meas_jt, meas_r, meas_jet_pt);
//  
            
        }

        if (!pt_cond || !rap_cond)
          continue;
  
        h1_jetpt_truth->Fill(jet_pt);
        h1_z_truth->Fill(truth_z);
        h1_jt_truth->Fill(truth_jt);
        h1_r_truth->Fill(truth_r);
        
        h2_zjt_truth->Fill(truth_z, truth_jt);
        h2_zr_truth->Fill(truth_z, truth_r);
        h2_jtr_truth->Fill( truth_jt, truth_r);

        h3_truthptzjt->Fill(truth_z, truth_jt, jet_pt);
        h3_truthptzr->Fill(truth_z, truth_r, jet_pt);
        h3_truthptjtr->Fill( truth_jt, truth_r, jet_pt );
        
    
        if (jet_pt > 12.5 && jet_pt < 15.0 ) {
            
            h1_z_truth_12515->Fill(truth_z);
            h1_jt_truth_12515->Fill(truth_jt);
            h1_r_truth_12515->Fill(truth_r);
          
        }
        if (jet_pt > 15.0 && jet_pt < 20.0 ) {
          
            h1_z_truth_1520->Fill(truth_z);
            h1_jt_truth_1520->Fill(truth_jt);
            h1_r_truth_1520->Fill(truth_r);
            
        }
        if (jet_pt > 20.0 && jet_pt < 30.0 ) {
          
            h1_z_truth_2030->Fill(truth_z);
            h1_jt_truth_2030->Fill(truth_jt);
            h1_r_truth_2030->Fill(truth_r);
        }
        if (jet_pt > 30.0 && jet_pt < 50.0 ) {
            h1_z_truth_3050->Fill(truth_z);
            h1_jt_truth_3050->Fill(truth_jt);
            h1_r_truth_3050->Fill(truth_r);
        }
        if (jet_pt > 50.0 && jet_pt < 100.0 ) {
            h1_z_truth_50100->Fill(truth_z);
            h1_jt_truth_50100->Fill(truth_jt);
            h1_r_truth_50100->Fill(truth_r);
        }
        if (jet_pt > 100.0 && jet_pt < 150.0 ) {
            h1_z_truth_100150->Fill(truth_z);
            h1_jt_truth_100150->Fill(truth_jt);
            h1_r_truth_100150->Fill(truth_r);
        }
        
        if (!meas_pt_cond || !rap_cond || !hasRecoHF)
          continue;
        
        
        // if(meas_kts->size() < 1) continue;
        
        h1_z->Fill(meas_z);
        h1_jt->Fill(meas_jt);
        h1_r->Fill(meas_r);

        h2_truthreco_jetpt->Fill(meas_jet_pt, jet_pt);
        h2_truthreco_z->Fill(meas_z, truth_z);
        h2_truthreco_jt->Fill(meas_jt, truth_jt);
        h2_truthreco_r->Fill(meas_r, truth_r);
        
        h2_zjt->Fill(meas_z, meas_jt);
        h2_zr->Fill(meas_z, meas_r);
        h2_jtr->Fill( meas_jt, meas_r);
        
        h3_ptzjt->Fill(meas_z, meas_jt, meas_jet_pt);
        h3_ptzr->Fill(meas_z, meas_r, meas_jet_pt);
        h3_ptjtr->Fill( meas_jt,meas_r, meas_jet_pt);
        
//        if (jet_pt > ptMin && jet_pt < ptMax)
//        {
//            response_zjt->Fill(meas_z, meas_jt, truth_z, truth_jt);
//            response_zr->Fill(meas_z, meas_r, truth_z, truth_r);
//        }

        event_counter++;
    }

//  RooUnfoldResponse * response_jetpt = new RooUnfoldResponse(h1_jetpt, h1_jetpt_truth, h2_truthreco_jetpt, "response_jetpt");
    TH2 *h2_response_jetpt = response_jetpt->Hresponse();
    TH2 *h2_response_z = response_z->Hresponse();
    TH2 *h2_response_jt = response_jt->Hresponse();
    TH2 *h2_response_r = response_r->Hresponse();

    TH2 *h2_response_zjt = response_zjt->Hresponse() ;
    TH2 *h2_response_zr = response_zr->Hresponse();
    TH2 *h2_response_jtr = response_jtr->Hresponse();
    
    TH2 *h2_response_ptjt = response_ptjt->Hresponse();
    TH2 *h2_response_ptr = response_ptr->Hresponse();
    TH2 *h2_response_ptz =  response_ptz->Hresponse();
    
    TH2 *h2_response_ptzjt = response_ptzjt->Hresponse();
    TH2 *h2_response_ptzr = response_ptzr->Hresponse();
    TH2 *h2_response_ptjtr = response_ptjtr->Hresponse();
    

    response_jetpt->Write("response_jetpt");
    response_z->Write("response_z");
    response_jt->Write("response_jt");
    response_r->Write("response_r");
    response_zjt->Write("response_zjt") ;
    response_zr->Write("response_zr" );
    response_jtr->Write( "response_jtr");
    
    response_ptjt->Write("response_ptjt" );
    response_ptr->Write( "response_ptr");
    response_ptz->Write("response_ptz" );
    
    response_ptzjt->Write("response_ptzjt" );
    response_ptzr->Write( "response_ptzr");
    response_ptjtr->Write( "response_ptjtr") ;
    
    

    TH1D *h1_efficiency_jetpt = (TH1D *)h1_num_efficiency_jetpt->Clone("h1_efficiency_jetpt");
    h1_efficiency_jetpt->Divide(h1_num_efficiency_jetpt, h1_denom_efficiency_jetpt, 1, 1, "B");

    TH1D *h1_efficiency_jetrap = (TH1D *)h1_num_efficiency_jetrap->Clone("h1_efficiency_jetrap");
    h1_efficiency_jetrap->Divide(h1_num_efficiency_jetrap, h1_denom_efficiency_jetrap, 1, 1, "B");

    TH1D *h1_efficiency_HFpt = (TH1D *)h1_num_efficiency_HFpt->Clone("h1_efficiency_HFpt");
    h1_efficiency_HFpt->Divide(h1_num_efficiency_HFpt, h1_denom_efficiency_HFpt, 1, 1, "B");

    TH1D *h1_efficiency_HFeta = (TH1D *)h1_num_efficiency_HFeta->Clone("h1_efficiency_HFeta");
    h1_efficiency_HFeta->Divide(h1_num_efficiency_HFeta, h1_denom_efficiency_HFeta, 1, 1, "B");
    
    TH1D *h1_efficiency_z = (TH1D *)h1_num_efficiency_z->Clone("h1_efficiency_z");
    h1_efficiency_z->Divide(h1_num_efficiency_z, h1_denom_efficiency_z, 1, 1, "B");
    
    TH1D *h1_efficiency_jt = (TH1D *)h1_num_efficiency_jt->Clone("h1_efficiency_jt");
    h1_efficiency_jt->Divide(h1_num_efficiency_jt, h1_denom_efficiency_jt, 1, 1, "B");
    
    TH1D *h1_efficiency_r = (TH1D *)h1_num_efficiency_r->Clone("h1_efficiency_r");
    h1_efficiency_r->Divide(h1_num_efficiency_r, h1_denom_efficiency_r, 1, 1, "B");

    TH2D *h2_efficiency_jetpteta = (TH2D *)h2_num_efficiency_jetpteta->Clone("h2_efficiency_jetpteta");
    h2_efficiency_jetpteta->Divide(h2_num_efficiency_jetpteta, h2_denom_efficiency_jetpteta, 1, 1, "B");

    TH2D *h2_efficiency_HFpteta = (TH2D *)h2_num_efficiency_HFpteta->Clone("h2_efficiency_HFpteta");
    h2_efficiency_HFpteta->Divide(h2_num_efficiency_HFpteta, h2_denom_efficiency_HFpteta, 1, 1, "B");


    // TEfficiency * h1_efficiency_jetpt = new TEfficiency(*h1_num_efficiency_jetpt, *h1_denom_efficiency_jetpt);
    // TEfficiency * h1_efficiency_jeteta = new TEfficiency(*h1_num_efficiency_jeteta, *h1_denom_efficiency_jeteta);

    /////////////////////////////////

    TSpline3 sp(h1_efficiency_HFpt);
    TF1 *f_eff1 = new TF1("f_eff1", "[0]*tanh([1]*(x-[3]))+[2]", 2, 25);
    TF1 *f_eff2 = new TF1("f_eff2", "[0]*tanh([1]*([3] - x))+[2]", 25, 100);
    f_eff1->SetParameters(0.2, 0.1, 0.2, 6.4);
    f_eff2->SetParameters(0.2, 0.1, 0.2, 100);
    f_eff2->SetParLimits(3, 80, 200);
//  h1_efficiency_HFpt->Fit("f_eff1", "E", "SAME", 2, 25);
//  h1_efficiency_HFpt->Fit("f_eff2", "E", "SAME", 25, 100);

//  f_eff1->Write("f_eff1");
//  f_eff2->Write("f_eff2");
// TF1 *f_eff = new TF1(
//     "f_eff", [=](double *x, double *p)
//     { return (sp.Eval(x[0])); },
//     sp.GetXmin(), sp.GetXmax(), 2);

    /////////////////////////////////
    
    h1_efficiency_jetpt->Write("efficiency_jetpt");
    h1_denom_efficiency_HFpt->Write("h1_denom_efficiency_HFpt");
    h1_denom_efficiency_jetpt->Write("h1_denom_efficiency_jetpt");
    h1_efficiency_HFpt->Write("efficiency_HFpt");
    h1_efficiency_HFeta->Write("efficiency_HFeta");
    h1_efficiency_jetrap->Write("efficiency_jetrap");
    h1_efficiency_z->Write("efficiency_z");
    h1_efficiency_jt->Write("efficiency_jt");
    h1_efficiency_r->Write("efficiency_r");
    
    h2_denom_efficiency_HFptjetpt->Write("h2_denom_efficiency_HFptjetpt");
    h2_efficiency_jetpteta->Write("efficiency_jetpteta");
    h2_efficiency_HFpteta->Write("efficiency_HFpteta");

    cout << event_counter << " events processed" << endl;

    //
    int ican = -1, iframe = -1, itext = -1;
    TCanvas *ccan[1000];
    TH1F *frame[1000];
    TLatex *text[1000];

    for (int i = 0; i < 1000; i++) {
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
    Double_t lhcbWidth = 2.00; // Old LHCb style: 3.00;
    // Text size
    Double_t lhcbTSize = 0.06;

    // gROOT->SetStyle("Plain");
    gStyle->SetCanvasDefH(1414);
    gStyle->SetCanvasDefW(1000);
    gStyle->SetPadBottomMargin(0.08);
    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetTitleStyle(0);
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

    gPad->SetLogz();
    // Tl.drawLatex(0.7, 0.8, "No cuts");
    ccan[ican]->cd(2);
    // Tl.drawLatex(0.7, 0.8, "No cuts");
//  ccan[ican]->cd(3);
    h2_response_jetpt->SetStats(0);
    gPad->SetLogz();
    h2_response_jetpt->Draw("COLZ");
    h2_response_jetpt->SetTitle("; Reco Jet pT; Truth Jet pT");

    ccan[ican]->cd(3);
    NormalizeHist(h1_jetpt);
    h1_jetpt->Draw("same");

//    NormalizeHist(h1_jetpt_truth);
    h1_jetpt_truth->SetLineColor(kGreen);
    h1_jetpt_truth->Draw("same");

    ccan[ican]->cd();
    ccan[ican]->Update();

    if (ican == 0) {
        ccan[ican]->Print(plotfileO.Data());
    } else {
        ccan[ican]->Print(plotfilePDF.Data());
    }

    //
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

    h1_efficiency_jetpt->SetMarkerColor(kBlue);
    h1_efficiency_jetpt->SetMarkerStyle(21);
    h1_efficiency_jetpt->Draw("P E SAME");
    h1_efficiency_jetpt->SetTitle("; Jet pT; HF Jet Reco. Efficiency");
    h1_efficiency_jetpt->SetMinimum(0.);
    h1_efficiency_jetpt->SetMaximum(1.);
    // h1_efficiency_jetpt->SetXTitle("E_{radiator}");
    // h1_efficiency_jetpt->SetYTitle("Counts");
    // gPad->Update();
    // auto graphpt3 = h1_efficiency_jetpt->GetPaintedGraph();
    // graphpt3->SetMinimum(0);
    // graphpt3->SetMaximum(1);
    // gPad->Update();
    // h1_efficiency_jetpt->SetMinimum(0.0);
    ccan[ican]->cd(2);
    h1_efficiency_jetrap->SetMarkerColor(kBlue);
    h1_efficiency_jetrap->SetMarkerStyle(21);
    h1_efficiency_jetrap->Draw("P E same");
    h1_efficiency_jetrap->SetTitle("; Jet #eta; HF Jet Reco. Efficiency");
    h1_efficiency_jetrap->SetMinimum(0.);
    h1_efficiency_jetrap->SetMaximum(1.);
    // h1_efficiency_jeteta->SetXTitle("E_{radiator}");
    // h1_efficiency_jeteta->SetYTitle("Counts");

    ccan[ican]->cd(3);

    h1_efficiency_HFpt->SetMarkerColor(kBlue);
    h1_efficiency_HFpt->SetMarkerStyle(21);
    h1_efficiency_HFpt->Draw("P E SAME");
    h1_efficiency_HFpt->SetTitle("; HF pT; HF Jet Reco. Efficiency");
    h1_efficiency_HFpt->SetMinimum(0.);
    h1_efficiency_HFpt->SetMaximum(1.05);
//
//  f_eff1->Draw("SAME");
//  f_eff2->Draw("SAME");

    ccan[ican]->cd(4);

    h1_efficiency_HFeta->SetMarkerColor(kBlue);
    h1_efficiency_HFeta->SetMarkerStyle(21);
    h1_efficiency_HFeta->Draw("P E SAME");
    h1_efficiency_HFeta->SetTitle("; HF #eta; HF Jet Reco. Efficiency");
    h1_efficiency_HFeta->SetMinimum(0.);
    h1_efficiency_HFeta->SetMaximum(1.);

    ccan[ican]->cd(5);

    h2_efficiency_HFpteta->Draw("COLZ TEXT");
    h2_efficiency_HFpteta->SetTitle("; HF p_{T}; HF #eta; Jet Reco Efficiency");
    h2_efficiency_HFpteta->SetMinimum(0.);
    // h2_efficiency_HFpteta->SetMaximum(1.);

    ccan[ican]->cd(6);

    h2_efficiency_jetpteta->Draw("COLZ TEXT");
    h2_efficiency_jetpteta->SetTitle("; jet p_{T}; jet #eta; Jet Reco Efficiency");
    h2_efficiency_jetpteta->SetMinimum(0.);
    // h2_efficiency_jetpteta->SetMaximum(1.);
    
    ccan[ican]->cd(7);

    h1_efficiency_z->SetMarkerColor(kBlue);
    h1_efficiency_z->SetMarkerStyle(21);
    h1_efficiency_z->Draw("P E SAME");
    h1_efficiency_z->SetTitle("; z ; HF Jet Reco. Efficiency");
    h1_efficiency_z->SetMinimum(0.);
    h1_efficiency_z->SetMaximum(1.);
    
    ccan[ican]->cd(8);

    h1_efficiency_jt->SetMarkerColor(kBlue);
    h1_efficiency_jt->SetMarkerStyle(21);
    h1_efficiency_jt->Draw("P E SAME");
    h1_efficiency_jt->SetTitle("; j_{T} [GeV] ; HF Jet Reco. Efficiency");
    h1_efficiency_jt->SetMinimum(0.);
    h1_efficiency_jt->SetMaximum(1.);
    
    ccan[ican]->cd(9);

    h1_efficiency_r->SetMarkerColor(kBlue);
    h1_efficiency_r->SetMarkerStyle(21);
    h1_efficiency_r->Draw("P E SAME");
    h1_efficiency_r->SetTitle("; r ; HF Jet Reco. Efficiency");
    h1_efficiency_r->SetMinimum(0.);
    h1_efficiency_r->SetMaximum(1.);

//    ccan[ican]->cd(7);
    ccan[ican]->cd();
    ccan[ican]->Update();
    


    if (ican == 0) {
        ccan[ican]->Print(plotfileO.Data());
    } else {
        ccan[ican]->Print(plotfilePDF.Data());
    }

    f.Write();
    f.Close();

    if (ican > -1) {
        cout << " You plotted " << ican + 1 << " canvasses......." << endl;
        ccan[ican]->Print(plotfileC.Data());

        cout << "Num of True B jets = " << NumBjets << endl;
        cout << "Num of Reco B jets = " << NumRecoBjets << endl;
    }
}
