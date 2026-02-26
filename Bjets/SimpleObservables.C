#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <TF1.h>
#include <TLatex.h>
#include <THStack.h>
#include <TChain.h>
#include <TStyle.h>
#include "Settings.h"
#include "../Helpers.h"


using namespace std;
void SimpleObservables(int NumEvts = -1, int dataset = 91599,
                bool isData = true,
                int DoTrackEff = 2,
                int DoTrigEff = 0,
                int DoPIDEff = 0,
                bool DoRecSelEff = false,
                int DoMassFit = 0,
                bool DoSignalSys = false,
                bool DoJetID = false,                                
                bool SubtractGS = false,
                bool sPlotFit = false,
                bool L0MuonDiMuon = false,
                bool DoTrackEff_SysCrossCheck = false)
{

  bool MCflag = !isData;
  followHardest = false;


  int year = (dataset / 10000) % 10;
  int JetMeth = (dataset / 1000) % 10;
  int flavor = (dataset / 100) % 10;
  int ptRange = (dataset / 10) % 10;
  int Mag = (dataset / 1) % 10;
  int HF_pdgcode = -99;

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
  TString str_year = "";
  TString str_PID = "";
  TString str_GS = "";

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

  if (flavor == 1)
    str_flavor = "_udsg";
  else if (flavor == 4)
    str_flavor = "_c";
  else if (flavor == 5)
    str_flavor = "_b";


  if (PID_cut)
    str_PID = "_PID";

  if (isData)
    str_level = "data";
  else
  {
    if (truthLevel)
      str_level = "truth";
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
  

  if (ghostCut)
    str_ghost = Form("_ghost_%.1f", ghostProb);

  TString str_L0 = "";
  if (L0MuonDiMuon)
    str_L0 = "_L0MuonDiMuon";

  TString str_sPlot = "";
  if (sPlotFit)
    str_sPlot = "_splotfit";


  TString str_tree;
  TString extension_RootFilesMC, extension_RootFiles;
  TString extension_read, extension_eff;
  TString eff_path;

  TString extension_prefix, extension_trackeff;
  TString extension;

  extension = str_level + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_PID + str_GS + str_sPlot + str_L0 + Form("_%d", dataset);
    
  extension_RootFilesMC = TString("../../root_files/BjetsMC/");
  extension_RootFiles = isData ?  TString("../../root_files/Bjets/") : extension_RootFilesMC;
    
  eff_path = TString( "./efficiencies/");
    
  float minimum_dR = 0.1;

  extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_L0 + Form("_%d", dataset);

  if (DoJetID)
    extension_prefix = TString("jetid_");    
  if (DoRecSelEff)
    extension_prefix = TString("recselsys_");
  if (DoSignalSys)
    extension_prefix = TString("signalsys_");
    
  if (DoTrackEff == 1)
    extension = TString("trackingsysup_") + extension;
  if (DoTrackEff == 2)
    extension = TString("trackingsysdown_") + extension;
  if (DoTrackEff_SysCrossCheck && DoTrackEff != 0) 
    extension = TString("crosscheck_") + extension;
  if (DoPIDEff == 1)
    extension = TString("pidsysup_") + extension;
  if (DoPIDEff == 2)
    extension = TString("pidsysdown_") + extension;
  if (DoTrigEff == 1)
    extension = TString("trigsysup_") + extension;
  if (DoTrigEff == 2)
    extension = TString("trigsysdown_") + extension;  
  if (DoMassFit == 1)
    extension_prefix = TString("massfitsysnear_");
  if (DoMassFit == 2)
    extension_prefix = TString("massfitsysfar_");
  if (DoMassFit == 3)
    extension_prefix = TString("massfitsyslower_");
  if (DoMassFit == 4)
    extension_prefix = TString("massfitsysupper_");

  

  extension = extension_prefix + extension;
  extension_trackeff = "trackEff_" + str_year + "_Ratio_Full_Long_method";

  /////////////////// Mass Fit Parameters /////////////////////////////////
  TString massfit("");
  if (sPlotFit)
  {
    massfit = TString("splotfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.4_b" + str_PID + str_L0 + "_91599.root"; 
  }
  else
  {
    massfit = TString("massfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.4_b" + str_PID + str_L0 + "_91599.root"; 
  }

    
  if (DoRecSelEff)
    massfit = "recselsys_" + massfit;
  if (DoSignalSys)
  {
    massfit = "sys_" + massfit;
  }
  
  TString extension_mass = "../../root_files/Bjets/" + massfit;  
  std::cout << "extension_mass : " << extension_mass << std::endl;   
  TFile f_massfit( extension_mass, "READ");
  TH1D *h1_MassMin = (TH1D *)f_massfit.Get("h1_MassMin");
  TH1D *h1_MassMax = (TH1D *)f_massfit.Get("h1_MassMax");
  TH1D *h1_Sideband1Min = (TH1D *)f_massfit.Get("h1_Sideband1Min");
  TH1D *h1_Sideband1Max = (TH1D *)f_massfit.Get("h1_Sideband1Max");
  TH1D *h1_Sideband2Min = (TH1D *)f_massfit.Get("h1_Sideband2Min");
  TH1D *h1_Sideband2Max = (TH1D *)f_massfit.Get("h1_Sideband2Max");
  TH1D *h1_BkgScale = (TH1D *)f_massfit.Get("h1_BkgScale");
  TH1D *h1_BkgScale_forSysNear = (TH1D *)f_massfit.Get("h1_BkgScale_forSysNear");
  TH1D *h1_BkgScale_forSysFar = (TH1D *)f_massfit.Get("h1_BkgScale_forSysFar");
  TH1D *h1_BkgScale_forSysLower = (TH1D *)f_massfit.Get("h1_BkgScale_forSysLower");
  TH1D *h1_BkgScale_forSysUpper = (TH1D *)f_massfit.Get("h1_BkgScale_forSysUpper");

  TTree *sWeightTree;
  if (sPlotFit)
  {
    sWeightTree = (TTree *)f_massfit.Get("sWeightTree");
  }
    
  /////////////////// Trigger Ratio (TISTOS/MC True) /////////////////////////////////
  ////////////////////////////////////////////////////
  ///              Trigger
  //////////////////////////////////////////////////
  TString extension_trig_MC, extension_trig_Data;
  TH2D *h2_trigeff_Data;
  TH2D *h2_trigeff_MC;

  //extension_trig_MC = "PhotonHadronElectronTIS_jpsieff_reco_ev_-1_b_PID_91599";
  //extension_trig_Data = "PhotonHadronElectronTIS_jpsieff_data_ev_-1_b_PID_91599";
  extension_trig_MC = "jpsieff_reco_ev_-1_b_PID" + str_L0 + "_91599.root";
  extension_trig_Data = "jpsieff_data_ev_-1_b_PID" + str_L0 + "_91599.root";

  TFile file_trigeffMC("./efficiencies/TrigEff/MC/" + extension_trig_MC, "READ");
  TFile file_trigeffData("./efficiencies/TrigEff/data/" + extension_trig_Data, "READ");

  h2_trigeff_Data = (TH2D *)file_trigeffData.Get("efficiency_Jpsiptrap");
  h2_trigeff_MC = (TH2D *)file_trigeffMC.Get("efficiency_Jpsiptrap");
  TH2D *h2_trigeff_ratio = (TH2D *)h2_trigeff_Data->Clone("h2_trigeff_ratio");
  h2_trigeff_ratio->Divide(h2_trigeff_MC);

  TH2D *h2_trig_ratio = (TH2D *)file_trigeffMC.Get("h2_eff_ratio");

  TChain *BTree = new TChain("BTree", "B-jets Tree Variables");
  vector<int> vec_datasets;
  if (Mag == 0)
    vec_datasets = {61590, 71590, 81590};
  else if (Mag == 1)
    vec_datasets = {61591, 71591, 81591};
  else
    vec_datasets = {61590, 61591, 71590, 71591, 81590, 81591};

  if (year == 9 && JetMeth != 9)
  {
    for (int i = 0; i < vec_datasets.size(); i++)
    {
      Mag = (vec_datasets[i] / 1) % 10;
      if (Mag == 0)
        str_Mag = "_MD";
      else if (Mag == 1)
        str_Mag = "_MU";
      extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost  + str_Mag + str_flavor + str_L0 + Form("_%d", vec_datasets[i]);
      if (!DoRecSelEff && DoMassFit == 0 && DoSignalSys == 0)
      {
        extension_read = extension_prefix + extension_read;
        cout << "Loading Dataset(s) " << extension_read << endl;        
      }

     BTree->Add(extension_RootFiles + extension_read + ".root/BTree");
    }
  }
  else
  {
    extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_L0 + Form("_%d", dataset);
    if (!DoRecSelEff && DoMassFit == 0 && DoSignalSys == 0)
    {   
      extension_read = extension_prefix + extension_read;
      cout << "Loading Dataset(s) " << extension_read << endl;       
    }
      BTree->Add(extension_RootFiles + extension_read + ".root/BTree");
  }

  cout << BTree->GetTreeNumber() << endl;
  if (NumEvts > BTree->GetEntries())
    NumEvts = BTree->GetEntries();
  if (NumEvts == -1)
    NumEvts = BTree->GetEntries();
  cout << BTree->GetEntries() << endl;


  TFile f(extension_RootFiles  + extension + ".root", "RECREATE");

  TH2D *h2_zjt = new TH2D("zjt", ";z; j_{T} [GeV/c]", zbinsize, z_binedges, jtbinsize, jt_binedges);
  TH2D *h2_zr = new TH2D("zr", ";z; r", zbinsize, z_binedges, rbinsize, r_binedges);
  TH2D *h2_jtr = new TH2D("jtr", "; j_{T} [GeV/c]; r", jtbinsize, jt_binedges, rbinsize, r_binedges);
    
  TH1D *h1_z = new TH1D("z", ";z", zbinsize, z_binedges);
  TH1D *h1_jt = new TH1D("jt", ";j_{T} [GeV/c]", jtbinsize, jt_binedges);
  TH1D *h1_r = new TH1D("r", ";r", rbinsize, r_binedges);

  TH2D *h2_ptz = new TH2D("ptz", ";z;p_{T,jet} [GeV/c]", zbinsize, z_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptjt = new TH2D("ptjt", ";j_{T} [GeV/c];p_{T,jet} [GeV/c]", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptr = new TH2D("ptr", ";r;p_{T,jet} [GeV/c]", rbinsize, r_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptz_sweight = new TH2D("ptz_sweight", ";z;p_{T,jet} [GeV/c]", zbinsize, z_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptjt_sweight = new TH2D("ptjt_sweight", ";j_{T} [GeV/c];p_{T,jet} [GeV/c]", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptr_sweight = new TH2D("ptr_sweight", ";r;p_{T,jet} [GeV/c]", rbinsize, r_binedges, ptbinsize, pt_binedges);
    
  TH3D *h3_ptzjt = new TH3D("ptzjt", ";z;j_{T} [GeV/c];p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptzr = new TH3D("ptzr", ";z;r;p_{T} [GeV/c]",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptjtr = new TH3D("ptjtr", ";j_{T} [GeV/c];r;p_{T} [GeV/c]",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptzjt_sweight = new TH3D("ptzjt_sweight", ";z;j_{T} [GeV/c];p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptzr_sweight = new TH3D("ptzr_sweight", ";z;r;p_{T} [GeV/c]",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptjtr_sweight = new TH3D("ptjtr_sweight", ";j_{T} [GeV/c];r;p_{T} [GeV/c]",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
    
  TH2D *h2_ptz_uncorrected = new TH2D("ptz_uncorrected", ";z;p_{T,jet} [GeV/c]", zbinsize, z_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptjt_uncorrected = new TH2D("ptjt_uncorrected", ";j_{T} [GeV/c];p_{T,jet} [GeV/c]", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptr_uncorrected = new TH2D("ptr_uncorrected", ";r;p_{T,jet} [GeV/c]", rbinsize, r_binedges, ptbinsize, pt_binedges);

  TH2D *h2_ptz_uncorrected_nomasscond = new TH2D("ptz_uncorrected_nomasscond", ";z;p_{T,jet} [GeV/c]", zbinsize, z_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptjt_uncorrected_nomasscond = new TH2D("ptjt_uncorrected_nomasscond", ";j_{T} [GeV/c];p_{T,jet} [GeV/c]", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptr_uncorrected_nomasscond = new TH2D("ptr_uncorrected_nomasscond", ";r;p_{T,jet} [GeV/c]", rbinsize, r_binedges, ptbinsize, pt_binedges);
    
  TH3D *h3_ptzjt_uncorrected = new TH3D("ptzjt_uncorrected", ";z;j_{T} [GeV/c];p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptzr_uncorrected = new TH3D("ptzr_uncorrected", ";z;r;p_{T} [GeV/c]",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptjtr_uncorrected = new TH3D("ptjtr_uncorrected", ";j_{T} [GeV/c];r;p_{T} [GeV/c]",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );    
     
  TH1D *h1_z_ptHFcut_g5 = new TH1D("z_ptHFcut_g5", ";z;", zbinsize, z_binedges);
  TH1D *h1_z_ptHFcut_l5 = new TH1D("z_ptHFcut_l5", ";z;", zbinsize, z_binedges);    
    
  /// ------------------------------------------ COMBINATORIAL --------------------------------------------- ///
    
  TH1D *h1_jet_pt_comb = new TH1D("Jet_pT_comb", "", ptbinsize, pt_binedges); 
  TH1D *h1_z_comb = new TH1D("z_comb", "Combinatoric; z", zbinsize, z_binedges);
  TH1D *h1_jt_comb = new TH1D("jt_comb", "Combinatoric; j_{T}", jtbinsize, jt_binedges);
  TH1D *h1_r_comb = new TH1D("r_comb", "Combinatoric; r", rbinsize, r_binedges);
    
  TH2D *h2_zjt_comb = new TH2D("zjt_comb", "", zbinsize, z_binedges,  jtbinsize, jt_binedges);
  TH2D *h2_zr_comb = new TH2D("zr_comb", "", zbinsize, z_binedges,  rbinsize, r_binedges);
  TH2D *h2_jtr_comb = new TH2D("jtr_comb", "", jtbinsize, jt_binedges,  rbinsize, r_binedges);

  TH2D *h2_ptz_comb = new TH2D("ptz_comb", "", zbinsize, z_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptjt_comb = new TH2D("ptjt_comb", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptr_comb = new TH2D("ptr_comb", "", rbinsize, r_binedges, ptbinsize, pt_binedges);
    
  TH3D *h3_ptzjt_comb = new TH3D("ptzjt_comb", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptzr_comb = new TH3D("ptzr_comb ", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
  TH3D *h3_ptjtr_comb = new TH3D("ptjtr_comb", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );

  TH2D *h2_ptz_comb_nobkgweight = new TH2D("ptz_comb_nobkgweight", "", zbinsize, z_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptjt_comb_nobkgweight = new TH2D("ptjt_comb_nobkgweight", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
  TH2D *h2_ptr_comb_nobkgweight = new TH2D("ptr_comb_nobkgweight", "", rbinsize, r_binedges, ptbinsize, pt_binedges);

  /// ------------------------------------------ END COMBINATORIAL --------------------------------------------- ///

  /// ------------------------------------------ QA --------------------------------------------- ///


  TH2D *h2_jetpteta = new TH2D("h2_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
  TH2D *h2_jetpt_HFpt = new TH2D("jetpt_HFpt", "", 50, 0, 100, 50, 0, 100);

  TH1D *h1_jet_pt = new TH1D("Jet_pT", "", ptbinsize, pt_binedges);
  TH1D *h1_jet_eta = new TH1D("Jet_eta", "", 12, etaMin, etaMax);
  TH1D *h1_jet_rap = new TH1D("Jet_rap", "", 12, etaMin, etaMax);
  TH1D *h1_jet_phi = new TH1D("Jet_phi", "", 20, -3.14, 3.14);

  TH1D *h1_Jpsi_rap = new TH1D("Jpsi_rap", "", 12, etaMin - jetradius, etaMax + jetradius);
  TH1D *h1_Jpsi_pt = new TH1D("Jpsi_pT", "", 50, 0, 100);
  TH1D *h1_Jpsi_phi = new TH1D("Jpsi_phi", "", 20, -3.14, 3.14);
  TH1D *h1_Jpsi_mass = new TH1D("Jpsi_mass", "", 30, 3.1 - 0.2, 3.1 + 0.2);
  TH1D *h1_Jpsi_ipchi2 = new TH1D("Jpsi_ipchi2", "", 80, -3, 5);

  TH1D *h1_HF_rap = new TH1D("HF_rap", "", 12, etaMin - jetradius, etaMax + jetradius);
  TH1D *h1_HF_pt = new TH1D("HF_pT", "", 50, 0, 100);
  TH1D *h1_HF_phi = new TH1D("HF_phi", "", 20, -3.14, 3.14);
  TH1D *h1_HF_mass = new TH1D("HF_mass", "", 80, 5.15, 5.55);
  TH1D *h1_HF_mass_dtfcut = new TH1D("HF_mass_dtfcut", "", 80, 5.15, 5.55);
  TH1D *h1_HF_mass_sweight = new TH1D("HF_mass_sweight", "", 80, 5.15, 5.55);
  vector<TH1D *> h1_mass_HFpt, h1_mass_HFpt_sweight;

  for (int i = 0; i < ptHFbinsize; i++)
  {
      TH1D *h1_temp = new TH1D(Form("HF_mass%d", i), "", 80, 5.15, 5.55);
      TH1D *h1_temp_sweight = new TH1D(Form("HF_sweight_mass%d", i), "", 80, 5.15, 5.55);

      h1_mass_HFpt.push_back(h1_temp);
      h1_mass_HFpt_sweight.push_back(h1_temp_sweight);
  }

  TH1D *h1_HFpt = new TH1D("h1_HFpt", "", ptHFbinsize, ptHF_binedges);
  TH2D *h2_HFptjetpt = new TH2D("h2_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
  TH2D *h2_HFptrap = new TH2D("h2_HFptrap", "", ptHFbinsize, ptHF_binedges, etabinsize, eta_binedges);
  TH3D *h3_HFptjetptrap = new TH3D("h3_HFptjetptrap", "", ptHFbinsize, ptHF_binedges, fineptbinsize, finept_binedges, etabinsize, eta_binedges);

  TH2D *h2_mB_mJpsi = new TH2D("mB_mJpsi", ";m_{#mu #mu K}; m_{#mu #mu}", 40, 5.15, 5.55, 40, 2.997, 3.197);

  TH1D *h1_nSPDHits = new TH1D("nSPDHits", ";N_{SPD};", 100, 0., 1000.);
  
  TH1D *h1_z_nobgsub, *h1_jt_nobgsub, *h1_r_nobgsub;
  TH2D *h2_ptz_nobgsub, *h2_ptjt_nobgsub, *h2_ptr_nobgsub;
  h1_z_nobgsub = (TH1D*)h1_z->Clone("z_nobgsub");
  h1_jt_nobgsub = (TH1D*)h1_jt->Clone("jt_nobgsub");
  h1_r_nobgsub = (TH1D*)h1_r->Clone("r_nobgsub");
  h2_ptz_nobgsub = (TH2D*)h2_ptz->Clone("ptz_nobgsub");    
  h2_ptjt_nobgsub = (TH2D*)h2_ptjt->Clone("ptjt_nobgsub");            
  h2_ptr_nobgsub = (TH2D*)h2_ptr->Clone("ptr_nobgsub");     
   
        
  vector<float> *dtr_pt(0), *dtr_rap(0), *dtr_id(0), *dtr_3charge(0);
  vector<float> *sv_ipchi2(0), *sv_keys(0), *sv_pids(0);

  float jet_pt, jet_eta, tr_jet_pt, tr_jet_eta;
  float jet_rap, tr_jet_rap;
  float jet_px, jet_py, jet_pz, jet_e, jet_charge;
  float HF_px, HF_py, HF_pz, HF_e, HF_pt;
  float Z_px, Z_py, Z_pz, Z_e;
  float mup_px, mup_py, mup_pz, mup_e;
  float mum_px, mum_py, mum_pz, mum_e;
  float K_px, K_py, K_pz, K_e;
  float mup_CHI2NDOF, mup_GHOSTPROB, mup_IPCHI2;
  float mum_CHI2NDOF, mum_GHOSTPROB, mum_IPCHI2;
  float K_CHI2NDOF, K_GHOSTPROB, K_IPCHI2;
  float Jpsi_CHI2NDOF, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_IPCHI2, Jpsi_BPVDLS;
  float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
  float K_PIDK;
  int nTracks;

  float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
  float tr_mup_px, tr_mup_py, tr_mup_pz, tr_mup_e;
  float tr_mum_px, tr_mum_py, tr_mum_pz, tr_mum_e;
  float tr_K_px, tr_K_py, tr_K_pz, tr_K_e;
    
  double z, jt, r;
  double tr_z, tr_jt, tr_r;
  double zg, jtg, rg;

  bool isTrueBjet, Hasbbbar, TOS;
  float bmass_dtf;
  int nSV, NumBHads_tr;
  float jpsi_ipchi2;
  int eventNumber;
  float sv_mass, sv_chi2, sv_ntrks, sv_cosine;
  int SVTag;
  float chi2ndf_dtf, tau_dtf;
  int nSPDHits;
    
  float pideff_K(1.0), pideff_mup(1.0), pideff_mum(1.0);
  float pideff_K_err(1.0), pideff_mup_err(1.0), pideff_mum_err(1.0);
  float trkeff_ratio_K(1.0), trkeff_ratio_mup(1.0), trkeff_ratio_mum(1.0);
  float trkeff_ratio_K_errhi(1.0), trkeff_ratio_mup_errhi(1.0), trkeff_ratio_mum_errhi(1.0);
  float trkeff_ratio_K_errlo(1.0), trkeff_ratio_mup_errlo(1.0), trkeff_ratio_mum_errlo(1.0);
  float trigeff_Data(1.0), trigeff_MC(1.0), trigeff_ratio(1.0);
    
  BTree->SetBranchAddress("eventNumber", &eventNumber);

  BTree->SetBranchAddress("dtr_pt", &dtr_pt);
  BTree->SetBranchAddress("dtr_rap", &dtr_rap);
  BTree->SetBranchAddress("dtr_id", &dtr_id);
  BTree->SetBranchAddress("dtr_3charge", &dtr_3charge);

  BTree->SetBranchAddress("jet_pt", &jet_pt);
  BTree->SetBranchAddress("jet_eta", &jet_eta);
  BTree->SetBranchAddress("jet_rap", &jet_rap);

  BTree->SetBranchAddress("jet_px", &jet_px);
  BTree->SetBranchAddress("jet_py", &jet_py);
  BTree->SetBranchAddress("jet_pz", &jet_pz);
  BTree->SetBranchAddress("jet_e", &jet_e);

  BTree->SetBranchAddress("HF_px", &HF_px);
  BTree->SetBranchAddress("HF_py", &HF_py);
  BTree->SetBranchAddress("HF_pz", &HF_pz);
  BTree->SetBranchAddress("HF_e", &HF_e);
  BTree->SetBranchAddress("HF_pt", &HF_pt);

  BTree->SetBranchAddress("Bu_IPCHI2", &Bu_IPCHI2);
  BTree->SetBranchAddress("Bu_CHI2", &Bu_CHI2);
  BTree->SetBranchAddress("Bu_CHI2NDOF", &Bu_CHI2NDOF);
    
  BTree->SetBranchAddress("mum_px", &mum_px);
  BTree->SetBranchAddress("mum_py", &mum_py);
  BTree->SetBranchAddress("mum_pz", &mum_pz);
  BTree->SetBranchAddress("mum_e", &mum_e);
  BTree->SetBranchAddress("mum_IPCHI2", &mum_IPCHI2);
  BTree->SetBranchAddress("mum_CHI2NDOF", &mum_CHI2NDOF);
  BTree->SetBranchAddress("mum_GHOSTPROB", &mum_GHOSTPROB);
    
  BTree->SetBranchAddress("mup_px", &mup_px);
  BTree->SetBranchAddress("mup_py", &mup_py);
  BTree->SetBranchAddress("mup_pz", &mup_pz);
  BTree->SetBranchAddress("mup_e", &mup_e);
  BTree->SetBranchAddress("mup_IPCHI2", &mup_IPCHI2);
  BTree->SetBranchAddress("mup_CHI2NDOF", &mup_CHI2NDOF);
  BTree->SetBranchAddress("mup_GHOSTPROB", &mup_GHOSTPROB);

  BTree->SetBranchAddress("K_px", &K_px);
  BTree->SetBranchAddress("K_py", &K_py);
  BTree->SetBranchAddress("K_pz", &K_pz);
  BTree->SetBranchAddress("K_e", &K_e);
  BTree->SetBranchAddress("K_IPCHI2", &K_IPCHI2);
  BTree->SetBranchAddress("K_CHI2NDOF", &K_CHI2NDOF);
  BTree->SetBranchAddress("K_GHOSTPROB", &K_GHOSTPROB);
    
  BTree->SetBranchAddress("nTracks", &nTracks);
    
  BTree->SetBranchAddress("Jpsi_FDCHI2", &Jpsi_FDCHI2);
  BTree->SetBranchAddress("Jpsi_CHI2", &Jpsi_CHI2);
  BTree->SetBranchAddress("Jpsi_CHI2NDOF", &Jpsi_CHI2NDOF);
  BTree->SetBranchAddress("Jpsi_BPVDLS", &Jpsi_BPVDLS);  
    
  BTree->SetBranchAddress("jpsi_ipchi2", &jpsi_ipchi2);


  BTree->SetBranchAddress("tr_jet_pt", &tr_jet_pt);
  BTree->SetBranchAddress("tr_jet_eta", &tr_jet_eta);
  BTree->SetBranchAddress("tr_jet_rap", &tr_jet_rap);

  BTree->SetBranchAddress("tr_jet_px", &tr_jet_px);
  BTree->SetBranchAddress("tr_jet_py", &tr_jet_py);
  BTree->SetBranchAddress("tr_jet_pz", &tr_jet_pz);
  BTree->SetBranchAddress("tr_jet_e", &tr_jet_e);

  BTree->SetBranchAddress("tr_mum_px", &tr_mum_px);
  BTree->SetBranchAddress("tr_mum_py", &tr_mum_py);
  BTree->SetBranchAddress("tr_mum_pz", &tr_mum_pz);
  BTree->SetBranchAddress("tr_mum_e", &tr_mum_e);

  BTree->SetBranchAddress("tr_mup_px", &tr_mup_px);
  BTree->SetBranchAddress("tr_mup_py", &tr_mup_py);
  BTree->SetBranchAddress("tr_mup_pz", &tr_mup_pz);
  BTree->SetBranchAddress("tr_mup_e", &tr_mup_e);

  BTree->SetBranchAddress("tr_K_px", &tr_K_px);
  BTree->SetBranchAddress("tr_K_py", &tr_K_py);
  BTree->SetBranchAddress("tr_K_pz", &tr_K_pz);
  BTree->SetBranchAddress("tr_K_e", &tr_K_e);

  BTree->SetBranchAddress("isTrueBjet", &isTrueBjet);

  BTree->SetBranchAddress("nSV", &nSV);
  BTree->SetBranchAddress("sv_mass", &sv_mass);
  BTree->SetBranchAddress("sv_chi2", &sv_chi2);
  BTree->SetBranchAddress("sv_cosine", &sv_cosine);
  BTree->SetBranchAddress("sv_ntrks", &sv_ntrks);
  BTree->SetBranchAddress("sv_ipchi2", &sv_ipchi2);
  BTree->SetBranchAddress("sv_keys", &sv_keys);
  BTree->SetBranchAddress("sv_pids", &sv_pids);    
  BTree->SetBranchAddress("SVTag", &SVTag);
    
  BTree->SetBranchAddress("bmass_dtf", &bmass_dtf);
  BTree->SetBranchAddress("chi2ndf_dtf", &chi2ndf_dtf);
  BTree->SetBranchAddress("tau_dtf", &tau_dtf);
  BTree->SetBranchAddress("NumBHads_tr", &NumBHads_tr);
  BTree->SetBranchAddress("Hasbbbar", &Hasbbbar);
  BTree->SetBranchAddress("K_PIDK", &K_PIDK);
  BTree->SetBranchAddress("pideff_K", &pideff_K);
  BTree->SetBranchAddress("pideff_mum", &pideff_mum);
  BTree->SetBranchAddress("pideff_mup", &pideff_mup);
  BTree->SetBranchAddress("pideff_K_err", &pideff_K_err);
  BTree->SetBranchAddress("pideff_mum_err", &pideff_mum_err);
  BTree->SetBranchAddress("pideff_mup_err", &pideff_mup_err);

  BTree->SetBranchAddress("trkeff_ratio_K", &trkeff_ratio_K);
  BTree->SetBranchAddress("trkeff_ratio_mup", &trkeff_ratio_mup);
  BTree->SetBranchAddress("trkeff_ratio_mum", &trkeff_ratio_mum);

  BTree->SetBranchAddress("trkeff_ratio_K_errhi", &trkeff_ratio_K_errhi);
  BTree->SetBranchAddress("trkeff_ratio_mup_errhi", &trkeff_ratio_mup_errhi);
  BTree->SetBranchAddress("trkeff_ratio_mum_errhi", &trkeff_ratio_mum_errhi);
  BTree->SetBranchAddress("trkeff_ratio_K_errlo", &trkeff_ratio_K_errlo);
  BTree->SetBranchAddress("trkeff_ratio_mup_errlo", &trkeff_ratio_mup_errlo);
  BTree->SetBranchAddress("trkeff_ratio_mum_errlo", &trkeff_ratio_mum_errlo);

  BTree->SetBranchAddress("trigeff_Data", &trigeff_Data);
  BTree->SetBranchAddress("trigeff_MC", &trigeff_MC);
  //BTree->SetBranchAddress("trigeff_ratio", &trigeff_ratio);
  BTree->SetBranchAddress("TOS", &TOS);
  BTree->SetBranchAddress("nSPDHits", &nSPDHits);

    
  BTree->SetBranchAddress("z", &z);
  BTree->SetBranchAddress("jt", &jt);
  BTree->SetBranchAddress("r", &r);
  BTree->SetBranchAddress("zg", &zg);
  BTree->SetBranchAddress("jtg", &jtg);
  BTree->SetBranchAddress("rg", &rg);
  BTree->SetBranchAddress("tr_z", &tr_z);
  BTree->SetBranchAddress("tr_jt", &tr_jt);
  BTree->SetBranchAddress("tr_r", &tr_r);

  float sweight;
  if (sPlotFit)
  {
    sWeightTree->SetBranchAddress("sweight", &sweight);
  }

  int eventNum;

  cout << "Requested # of events : " << NumEvts << endl;
  if (sPlotFit)
  {
    cout << "NumEvents in the sweight tree : " << sWeightTree->GetEntries() << std::endl; 
    if (isData && !DoJetID && NumEvts != sWeightTree->GetEntries())
    {
      std::cout << "Different number of entries in sWeightTree vs. Data tree is different!" << std::endl;
      return;
    }
  }

  int nan_counter = 0;
  //
  // Event loop
  //
  for (int ev = 0; ev < NumEvts; ev++)
  {
    BTree->GetEntry(ev);
    if (isData && sPlotFit)
    {
      sWeightTree->GetEntry(ev);
    }
    else
    {
      sweight = 1.0;
    }
    if (ev % 100000 == 0)
      cout << "Executing event " << ev << endl;

    TLorentzVector HFmeson, HFjet, mum, mup, Kmeson, Jpsi;
    TLorentzVector tr_HFmeson, tr_HFjet, tr_mum, tr_mup, tr_Kmeson, tr_Jpsi;
    HFjet.SetPxPyPzE(jet_px, jet_py, jet_pz, jet_e);
    mup.SetPxPyPzE(mup_px, mup_py, mup_pz, mup_e);
    mum.SetPxPyPzE(mum_px, mum_py, mum_pz, mum_e);
    Kmeson.SetPxPyPzE(K_px, K_py, K_pz, K_e);

    bool pt_cond = (jet_pt > pTLow);
    bool rap_cond = (jet_rap > etaMin && jet_rap < etaMax);

    if (!TOS)
      continue;
    if (!pt_cond || !rap_cond)
      continue;

    bool StrippingCuts = false;
    if (StrippingCuts)
    {
      if (mup_CHI2NDOF > 4 || mum_CHI2NDOF > 4 || K_CHI2NDOF > 4)
        continue;
      if (mup_GHOSTPROB > 0.4 || mum_GHOSTPROB > 0.4 || K_GHOSTPROB > 0.4)
        continue;
      if (mup_IPCHI2 < 25 || mum_IPCHI2 < 25 || K_IPCHI2 < 25)
        continue;
    }


    if (sPlotFit)
    {
      float mass_low = 5.15;
      float mass_high = 5.55;
      if (bmass_dtf < mass_low || bmass_dtf > mass_high)
        continue;
    }

    //std:cout << "sweight : " << sweight << std::endl;


    tr_HFjet.SetPxPyPzE(tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e);
    tr_mup.SetPxPyPzE(tr_mup_px, tr_mup_py, tr_mup_pz, tr_mup_e);
    tr_mum.SetPxPyPzE(tr_mum_px, tr_mum_py, tr_mum_pz, tr_mum_e);
    tr_Kmeson.SetPxPyPzE(tr_K_px, tr_K_py, tr_K_pz, tr_K_e);
    tr_HFmeson = tr_mup + tr_mum + tr_Kmeson;
    tr_Jpsi = tr_mup + tr_mum;
    // HFmeson = mup + mum + Kmeson;
    HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
    Jpsi = mup + mum;
    float HF_eta = HFmeson.Eta();

    float event_weight = 1.0;

    float reweight = 1.0;
    float bkg_weight = h1_BkgScale != NULL ? h1_BkgScale->GetBinContent(h1_BkgScale->FindBin(HFmeson.Pt())) : 1.0;
    float MassHigh = h1_MassMax != NULL ? h1_MassMax->GetBinContent(h1_MassMax->FindBin(HFmeson.Pt())) : 5.31;
    float MassLow = h1_MassMin != NULL ? h1_MassMin->GetBinContent(h1_MassMin->FindBin(HFmeson.Pt())) : 5.24;

    float Sideband1Min = h1_Sideband1Min != NULL ? h1_Sideband1Min->GetBinContent(h1_Sideband1Min->FindBin(HFmeson.Pt())) : 5.24;
    float Sideband1Max = h1_Sideband1Max != NULL ? h1_Sideband1Max->GetBinContent(h1_Sideband1Max->FindBin(HFmeson.Pt())) : 5.31;
    float Sideband2Min = h1_Sideband2Min != NULL ? h1_Sideband2Min->GetBinContent(h1_Sideband2Min->FindBin(HFmeson.Pt())) : 5.24;
    float Sideband2Max = h1_Sideband2Max != NULL ? h1_Sideband2Max->GetBinContent(h1_Sideband2Max->FindBin(HFmeson.Pt())) : 5.31;
    
    //  Set Trigger Correction Weight
    if (isData && h2_trigeff_ratio != NULL)
    {
      double rap = (Jpsi.Rapidity() > 2.00) ? Jpsi.Rapidity() : 2.1;
      double pt = (Jpsi.Pt() > 2.00) ? Jpsi.Pt() : 2.1;
      trigeff_ratio = h2_trigeff_ratio->GetBinContent(h2_trigeff_ratio->GetXaxis()->FindBin(pt), h2_trigeff_ratio->GetYaxis()->FindBin(rap));
      if (std::isnan(trigeff_ratio) || std::isinf(trigeff_ratio) || trigeff_ratio == 0)
        trigeff_ratio = trigeff_Data = trigeff_MC = 1.0;
    }
    else
    {
      if (isData)
      {
        std::cout << "h2_trigeff_ratio is NULL! Setting ratio to 1.0, need to investigate" << std::endl;
      }
      trigeff_ratio = 1.0;
    }
    

    // cout << HFmeson.Pt() << ", " << MassLow << ", " << MassHigh << ", " << Sideband1Min << ", " << Sideband1Max << ", " << Sideband2Min << ", " << Sideband2Max << endl;
    if (DoMassFit == 1)
    {
      TH1D *h1_Sideband1Min_forSysNear = (TH1D *)f_massfit.Get("h1_Sideband1Min_forSysNear");
      TH1D *h1_Sideband1Max_forSysNear = (TH1D *)f_massfit.Get("h1_Sideband1Max_forSysNear");
      TH1D *h1_Sideband2Min_forSysNear = (TH1D *)f_massfit.Get("h1_Sideband2Min_forSysNear");
      TH1D *h1_Sideband2Max_forSysNear = (TH1D *)f_massfit.Get("h1_Sideband2Max_forSysNear");
      Sideband1Min = h1_Sideband1Min_forSysNear != NULL ? h1_Sideband1Min_forSysNear->GetBinContent(h1_Sideband1Min_forSysNear->FindBin(HFmeson.Pt())) : 5.24;
      Sideband1Max = h1_Sideband1Max_forSysNear != NULL ? h1_Sideband1Max_forSysNear->GetBinContent(h1_Sideband1Max_forSysNear->FindBin(HFmeson.Pt())) : 5.31;
      Sideband2Min = h1_Sideband2Min_forSysNear != NULL ? h1_Sideband2Min_forSysNear->GetBinContent(h1_Sideband2Min_forSysNear->FindBin(HFmeson.Pt())) : 5.24;
      Sideband2Max = h1_Sideband2Max_forSysNear != NULL ? h1_Sideband2Max_forSysNear->GetBinContent(h1_Sideband2Max_forSysNear->FindBin(HFmeson.Pt())) : 5.31;
      bkg_weight = h1_BkgScale_forSysNear != NULL ? h1_BkgScale_forSysNear->GetBinContent(h1_BkgScale_forSysNear->FindBin(HFmeson.Pt())) : 1.0;
    }
    if (DoMassFit == 2)
    {
      TH1D *h1_Sideband1Min_forSysFar = (TH1D *)f_massfit.Get("h1_Sideband1Min_forSysFar");
      TH1D *h1_Sideband1Max_forSysFar = (TH1D *)f_massfit.Get("h1_Sideband1Max_forSysFar");
      TH1D *h1_Sideband2Min_forSysFar = (TH1D *)f_massfit.Get("h1_Sideband2Min_forSysFar");
      TH1D *h1_Sideband2Max_forSysFar = (TH1D *)f_massfit.Get("h1_Sideband2Max_forSysFar");
      Sideband1Min = h1_Sideband1Min_forSysFar != NULL ? h1_Sideband1Min_forSysFar->GetBinContent(h1_Sideband1Min_forSysFar->FindBin(HFmeson.Pt())) : 5.24;
      Sideband1Max = h1_Sideband1Max_forSysFar != NULL ? h1_Sideband1Max_forSysFar->GetBinContent(h1_Sideband1Max_forSysFar->FindBin(HFmeson.Pt())) : 5.31;
      Sideband2Min = h1_Sideband2Min_forSysFar != NULL ? h1_Sideband2Min_forSysFar->GetBinContent(h1_Sideband2Min_forSysFar->FindBin(HFmeson.Pt())) : 5.24;
      Sideband2Max = h1_Sideband2Max_forSysFar != NULL ? h1_Sideband2Max_forSysFar->GetBinContent(h1_Sideband2Max_forSysFar->FindBin(HFmeson.Pt())) : 5.31;
      bkg_weight = h1_BkgScale_forSysFar != NULL ? h1_BkgScale_forSysFar->GetBinContent(h1_BkgScale_forSysFar->FindBin(HFmeson.Pt())) : 1.0;
    }
    if (DoMassFit == 3)
    {
      bkg_weight = h1_BkgScale_forSysLower != NULL ? h1_BkgScale_forSysLower->GetBinContent(h1_BkgScale_forSysLower->FindBin(HFmeson.Pt())) : 1.0;
    }    
    if (DoMassFit == 4)
    {
      bkg_weight = h1_BkgScale_forSysUpper != NULL ? h1_BkgScale_forSysUpper->GetBinContent(h1_BkgScale_forSysUpper->FindBin(HFmeson.Pt())) : 1.0;
    }    

    bool PID_cond = (K_PIDK > 0);
    bool mass_cond = (bmass_dtf > MassLow && bmass_dtf < MassHigh);
    bool bkg_cond = (bmass_dtf > Sideband1Min && bmass_dtf < Sideband1Max) || (bmass_dtf > Sideband2Min && bmass_dtf < Sideband2Max);
    if (DoMassFit == 3)
    {
      bkg_cond = (bmass_dtf > Sideband1Min && bmass_dtf < Sideband1Max);
    }
    else if (DoMassFit == 4)
    {
      bkg_cond = (bmass_dtf > Sideband2Min && bmass_dtf < Sideband2Max);
    }
    bool gluon_cond = mass_cond && Hasbbbar;
    bool SV_cond = (nSV > 0) && mass_cond && sv_mass > 0.4;
    bool signal_cond = mass_cond;
    if (PID_cut && isData)
    {
      signal_cond = signal_cond && PID_cond;
      bkg_cond = bkg_cond && PID_cond;
      SV_cond = SV_cond && PID_cond;
    }

    if (DoTrackEff == 1)
    {
      trkeff_ratio_K = trkeff_ratio_K + trkeff_ratio_K_errhi;
      trkeff_ratio_mum = trkeff_ratio_mum + trkeff_ratio_mum_errhi;
      trkeff_ratio_mup = trkeff_ratio_mup + trkeff_ratio_mup_errhi;
      if (DoTrackEff_SysCrossCheck) // Standard Run 2 track eff systematic uncertainty is 0.8%
      {
        trkeff_ratio_K = (1 + sqrt(0.0127*0.0127 + 0.008*0.008)) * trkeff_ratio_K; // Material Budget of 1.27% on Kaons
        trkeff_ratio_mum = (1 + 0.008) * trkeff_ratio_mum; 
        trkeff_ratio_mup = (1 + 0.008) * trkeff_ratio_mup; 
      }
      else
      {
        trkeff_ratio_K = (1 + 0.0127) * trkeff_ratio_K; // Material Budget of 1.27% on Kaons
      }
    }
    else if (DoTrackEff == 2)
    {
      trkeff_ratio_K = trkeff_ratio_K - trkeff_ratio_K_errlo;
      trkeff_ratio_K = (1 - 0.0127) * trkeff_ratio_K; // Material Budget of 1.27% on Kaons
      trkeff_ratio_mum = trkeff_ratio_mum - trkeff_ratio_mum_errlo;
      trkeff_ratio_mup = trkeff_ratio_mup - trkeff_ratio_mup_errlo;

      if (DoTrackEff_SysCrossCheck) // Standard Run 2 track eff systematic uncertainty is 0.8%
      {
        trkeff_ratio_K = (1 - sqrt(0.0127*0.0127 + 0.008*0.008)) * trkeff_ratio_K; // Material Budget of 1.27% on Kaons
        trkeff_ratio_mum = (1 - 0.008) * trkeff_ratio_mum; 
        trkeff_ratio_mup = (1 - 0.008) * trkeff_ratio_mup; 
      }
      else
      {
        trkeff_ratio_K = (1 - 0.0127) * trkeff_ratio_K; // Material Budget of 1.27% on Kaons
      }

      if (trkeff_ratio_K_errlo > 0.1 || trkeff_ratio_mup_errlo > 0.1 || trkeff_ratio_mum_errlo > 0.1)
      {
        cout << trkeff_ratio_K << endl;
        cout << trkeff_ratio_mum_errlo << ", " << trkeff_ratio_mup_errlo << ", " << trkeff_ratio_K_errlo << endl;
        cout << mup.Pt() << ", " << mup.Eta() << endl;
        cout << mum.Pt() << ", " << mum.Eta() << endl;
        cout << Kmeson.P() << ", " << Kmeson.Eta() << endl;
      }
    }
    if (DoPIDEff == 1)
    {
      pideff_K = pideff_K + pideff_K_err;
      pideff_mum = pideff_mum + pideff_mum_err;
      pideff_mup = pideff_mup + pideff_mup_err;
    }
    else if (DoPIDEff == 2)
    {
      pideff_K = pideff_K - pideff_K_err;
      pideff_mum = pideff_mum - pideff_mum_err;
      pideff_mup = pideff_mup - pideff_mup_err;
    }
    if (DoRecSelEff)
    {
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

    if (DoTrigEff == 1)
    {
      double trig_var = h2_trig_ratio->GetBinContent(h2_trig_ratio->GetXaxis()->FindBin(HFmeson.Pt()), h2_trig_ratio->GetYaxis()->FindBin(HFmeson.Rapidity()));
      double trig_err = h2_trig_ratio->GetBinError(h2_trig_ratio->GetXaxis()->FindBin(HFmeson.Pt()), h2_trig_ratio->GetYaxis()->FindBin(HFmeson.Rapidity()));
      trig_var = fabs(1.0 - trig_var);
      trig_var = (trig_err > trig_var) ? trig_err : trig_var;
      trigeff_ratio = trigeff_ratio * (1 + trig_var);
    }
    if (DoTrigEff == 2)
    {
      double trig_var = h2_trig_ratio->GetBinContent(h2_trig_ratio->GetXaxis()->FindBin(HFmeson.Pt()), h2_trig_ratio->GetYaxis()->FindBin(HFmeson.Rapidity()));
      double trig_err = h2_trig_ratio->GetBinError(h2_trig_ratio->GetXaxis()->FindBin(HFmeson.Pt()), h2_trig_ratio->GetYaxis()->FindBin(HFmeson.Rapidity()));
      trig_var = fabs(1.0 - trig_var);
      trig_var = (trig_err > trig_var) ? trig_err : trig_var;
      trigeff_ratio = trigeff_ratio * (1 - trig_var);
    }

      
         
    event_weight = (1./ (trkeff_ratio_K * trkeff_ratio_mum * trkeff_ratio_mup * pideff_mum * pideff_mup * trigeff_ratio));       
    
    if (PID_cut)
      event_weight *= (1. / (pideff_K));
    
    if (std::isinf(event_weight) || std::isnan(event_weight))
    {
      // This seems to happen less than .01% of the time, so it is negligible
      //std::cout << "trkeff_ratio_K : " << trkeff_ratio_K << " trkeff_ratio_mum : " << trkeff_ratio_mum << " trkeff_ratio_mup : " << trkeff_ratio_mup << " pideff_mum : " << pideff_mum << " pideff_mup : " << pideff_mup << " pideff_K : " << pideff_K << " trigeff_ratio : " << trigeff_ratio << std::endl; 
      //std::cout << "entry number : " << ev << " event number : " << eventNumber <<  std::endl;
      //std::cout << "event_weight : " << event_weight << " is_inf or is_nan -- setting to 1.0, need to investigate" << std::endl;
      //nan_counter++;
      event_weight = 1.0;
    }

    float dphi_HF_jet = checkphi(checkphi(HFmeson.Phi()) - checkphi(HFjet.Phi()));
    float dy_HF_jet = HFjet.Eta() - HFmeson.Rapidity();

    if (bkg_cond)
    {

      h1_z_comb->Fill(z, event_weight * bkg_weight);
      h1_jt_comb->Fill(jt, event_weight * bkg_weight);
      h1_r_comb->Fill(r, event_weight * bkg_weight);
        
      h2_zjt_comb->Fill(z, jt, event_weight * bkg_weight);
      h2_zr_comb->Fill(z, r, event_weight * bkg_weight);
      h2_jtr_comb->Fill(jt, r, event_weight * bkg_weight);

      h2_ptz_comb->Fill(z, jet_pt, event_weight * bkg_weight);
      h2_ptjt_comb->Fill(jt, jet_pt, event_weight * bkg_weight);
      h2_ptr_comb->Fill(r, jet_pt, event_weight * bkg_weight);
        
      h3_ptzjt_comb->Fill(z, jt, jet_pt, event_weight * bkg_weight);
      h3_ptzr_comb->Fill(z, jt, jet_pt, event_weight * bkg_weight);
      h3_ptjtr_comb->Fill(jt, r,  jet_pt, event_weight * bkg_weight);

      h2_ptz_comb_nobkgweight->Fill(z, jet_pt, event_weight);
      h2_ptjt_comb_nobkgweight->Fill(jt, jet_pt, event_weight);
      h2_ptr_comb_nobkgweight->Fill(r, jet_pt, event_weight);

      h1_jet_pt_comb->Fill(jet_pt, event_weight * bkg_weight);
        
    }

    if (signal_cond)
    {

      // QA HISTS -- note there is no event weight applied here... //
      h2_jetpt_HFpt->Fill(HFmeson.Pt(), HFjet.Pt());
      h2_jetpteta->Fill(jet_pt, jet_eta);

      h1_Jpsi_mass->Fill(Jpsi.M());
      h1_Jpsi_pt->Fill(Jpsi.Pt());
      h1_Jpsi_rap->Fill(Jpsi.Rapidity());
      h1_Jpsi_phi->Fill(Jpsi.Phi());
      h1_Jpsi_ipchi2->Fill(jpsi_ipchi2);

      h1_HF_rap->Fill(HFmeson.Rapidity());
      h1_HF_pt->Fill(HFmeson.Pt());
      h1_HF_phi->Fill(HFmeson.Phi());
      h1_jet_eta->Fill(jet_eta);
      h1_jet_rap->Fill(jet_rap);
      h1_jet_phi->Fill(HFjet.Phi());


      // HISTS FOR ANALYSIS //
      h1_jet_pt->Fill(jet_pt, event_weight);

      h3_ptzjt->Fill(z, jt, jet_pt, event_weight);
      h3_ptzr->Fill(z, r, jet_pt, event_weight);
      h3_ptjtr->Fill(jt, r, jet_pt, event_weight);
      
      h2_ptz->Fill(z, jet_pt, event_weight);
      h2_ptjt->Fill(jt, jet_pt, event_weight);
      h2_ptr->Fill(r, jet_pt, event_weight); 
      
      h1_z->Fill(z, event_weight);
      h1_jt->Fill(jt, event_weight);
      h1_r->Fill(r, event_weight);
      
      h2_zjt->Fill(z, jt, event_weight);
      h2_zr->Fill(z, r, event_weight);
      h2_jtr->Fill(jt, r, event_weight);
      
      
      h1_z_nobgsub->Fill(z, event_weight);
      h1_jt_nobgsub->Fill(jt, event_weight);
      h1_r_nobgsub->Fill(r, event_weight);
      
      h2_ptz_nobgsub->Fill(z, jet_pt, event_weight);
      h2_ptjt_nobgsub->Fill(jt, jet_pt, event_weight);
      h2_ptr_nobgsub->Fill(r, jet_pt, event_weight);
      
      h3_ptzjt_uncorrected->Fill(z, jt, jet_pt);
      h3_ptzr_uncorrected->Fill(z, r, jet_pt);
      h3_ptjtr_uncorrected->Fill(jt, r, jet_pt);
      
      h2_ptz_uncorrected->Fill(z, jet_pt);
      h2_ptjt_uncorrected->Fill(jt, jet_pt);
      h2_ptr_uncorrected->Fill(r, jet_pt);     
              
      if (HF_pt < 5.) { h1_z_ptHFcut_l5->Fill(z); }        
      else { h1_z_ptHFcut_g5->Fill(z); }

      h1_nSPDHits->Fill(nSPDHits);

    } // end signal cond 
    
    if (PID_cut && PID_cond)
    {
      h3_ptzjt_sweight->Fill(z, jt, jet_pt, event_weight*sweight);
      h3_ptzr_sweight->Fill(z, r, jet_pt, event_weight*sweight);
      h3_ptjtr_sweight->Fill(jt, r, jet_pt, event_weight*sweight);

      h2_ptz_sweight->Fill(z, jet_pt, event_weight*sweight);
      h2_ptjt_sweight->Fill(jt, jet_pt, event_weight*sweight);
      h2_ptr_sweight->Fill(r, jet_pt, event_weight*sweight);   

      h2_ptz_uncorrected_nomasscond->Fill(z, jet_pt);
      h2_ptjt_uncorrected_nomasscond->Fill(jt, jet_pt);
      h2_ptr_uncorrected_nomasscond->Fill(r, jet_pt);

      h1_HF_mass->Fill(bmass_dtf);
      h1_HF_mass_sweight->Fill(bmass_dtf, sweight);
      h2_mB_mJpsi->Fill(bmass_dtf, Jpsi.M()); 

      for (int i = 0; i < ptHFbinsize; i++)
      {
          if (HFmeson.Pt() > ptHF_binedges[i] && HFmeson.Pt() < ptHF_binedges[i + 1])
          {
              h1_mass_HFpt[i]->Fill(bmass_dtf);
              h1_mass_HFpt_sweight[i]->Fill(bmass_dtf, sweight);
              break;
          }
      }

    }
    
  } // end event loop

  //std::cout << " number of events (entries) : " << NumEvts << std::endl;
  //std::cout << " number of nan events (entries) : " << nan_counter << std::endl;
  //std::cout << " fraction of nan events (entries) : " << nan_counter / NumEvts << std::endl;

  ///////////////// ////////////////
  // Background Subtraction
  ////////////// ////////////////

  // NOTE!! Integrals do not yield the exact same numbers in all cases before sideband subtraction unless underflow and overflow bins are supplied as bounds (e.g. ptz->Integral(0, zbinsize+1, 0, ptbinsize+1))
  cout << "Before sub: " << endl;
  cout << "Num Jet pt = " << h1_jet_pt->Integral() << endl;
  cout << "In ptz = " << h2_ptz->Integral() << endl;
  cout << "In ptjt = " << h2_ptjt->Integral() << endl;
  cout << "In ptr = " << h2_ptr->Integral() << endl;
  cout << "In ptzjt = " << h3_ptzjt->Integral() << endl;
  cout << "In ptzr = " << h3_ptzr->Integral() << endl;
  cout << "In ptjtr = " << h3_ptjtr->Integral() << endl;
  if (isData)
  {
  
    // Why don't we call MakeHistPositive on the 1D distributions?
    h1_jet_pt->Add(h1_jet_pt_comb, -1);
                        
    h1_z->Add(h1_z_comb, -1);
    h1_jt->Add(h1_jt_comb, -1);
    h1_r->Add(h1_r_comb, -1);
    
    h2_zjt->Add(h2_zjt_comb, -1);
    MakeHistPositive(h2_zjt);
    h2_zr->Add(h2_zr_comb, -1);
    MakeHistPositive(h2_zr);
    h2_jtr->Add(h2_jtr_comb, -1);
    MakeHistPositive(h2_jtr);
 
    h2_ptz->Add(h2_ptz_comb, -1);
    MakeHistPositive(h2_ptz);    
    h2_ptjt->Add(h2_ptjt_comb, -1);
    MakeHistPositive(h2_ptjt);
    h2_ptr->Add(h2_ptr_comb, -1);
    MakeHistPositive(h2_ptr); 
      
    h3_ptzjt->Add(h3_ptzjt_comb, -1);
    MakeHistPositive(h3_ptzjt);    
    h3_ptzr->Add(h3_ptzr_comb, -1);
    MakeHistPositive(h3_ptzr);
    h3_ptjtr->Add(h3_ptjtr_comb, -1);
    MakeHistPositive(h3_ptjtr);
           
  }
  cout << "After sub: " << endl;
  cout << "Num Jet pt = " << h1_jet_pt->Integral() << endl;
  cout << "In ptz = " << h2_ptz->Integral() << endl;
  cout << "In ptjt = " << h2_ptjt->Integral() << endl;
  cout << "In ptr = " << h2_ptr->Integral() << endl;
  cout << "In ptzjt = " << h3_ptzjt->Integral() << endl;
  cout << "In ptzr = " << h3_ptzr->Integral() << endl;
  cout << "In ptjtr = " << h3_ptjtr->Integral() << endl;

  THStack *hs_ptz = new THStack("hs_ptz", ";z;#frac{1}{N_{jets}}#frac{dN}{dz}");
  THStack *hs_ptjt = new THStack("hs_ptjt", ";j_{T} [GeV/c];#frac{1}{N_{jets}}#frac{dN}{dj_{T}}");
  THStack *hs_ptr = new THStack("hs_ptr", ";r;#frac{1}{N_{jets}}#frac{dN}{dr}");    
  THStack *hs_ptz_sweight = new THStack("hs_ptz_sweight", ";z;#frac{1}{N_{jets}}#frac{dN}{dz}");
  THStack *hs_ptjt_sweight = new THStack("hs_ptjt_sweight", ";j_{T} [GeV/c];#frac{1}{N_{jets}}#frac{dN}{dj_{T}}");
  THStack *hs_ptr_sweight = new THStack("hs_ptr_sweight", ";r;#frac{1}{N_{jets}}#frac{dN}{dr}");        
  THStack *hs_ptz_uncorrected = new THStack("z_uncorrected_data_all", ";z;#frac{1}{N_{jets}}#frac{dN}{dz}");
  THStack *hs_ptjt_uncorrected = new THStack("jt_uncorrected_data_all", ";j_{T} [GeV/c];#frac{1}{N_{jets}}#frac{dN}{dj_{T}}");
  THStack *hs_ptr_uncorrected = new THStack("r_uncorrected_data_all", ";r;#frac{1}{N_{jets}}#frac{dN}{dr}");
  THStack *hs_ptz_nobgsub = new THStack("z_evtweights_nosbsub_data_all", ";z;#frac{1}{N_{jets}}#frac{dN}{dz}");
  THStack *hs_ptjt_nobgsub = new THStack("jt_evtweights_nosbsub_data_all", ";j_{T} [GeV/c];#frac{1}{N_{jets}}#frac{dN}{dj_{T}}");
  THStack *hs_ptr_nobgsub = new THStack("r_evtweights_nosbsub_data_all", ";r;#frac{1}{N_{jets}}#frac{dN}{dr}");  
  //vector<THStack> z_sbsub_inputs_pt;
  for (int j = 1; j < ptbinsize; j++)
  {
    THStack *z_sbsub_inputs_temp = new THStack(Form("z_sbsub_inputs_pt%d", j),"");
    TH1D *h1_temp_z_uncorr = (TH1D *)h2_ptz_uncorrected->ProjectionX(Form("z_uncorr_pt%d", j), j + 1, j + 1);
    TH1D *h1_temp_z_comb = (TH1D *)h2_ptz_comb->ProjectionX(Form("z_comb_pt%d", j), j + 1, j + 1);
    TH1D *h1_temp_z_comb_nobkgweight = (TH1D *)h2_ptz_comb_nobkgweight->ProjectionX(Form("z_comb_nobkgweight_pt%d", j), j + 1, j + 1);
    h1_temp_z_uncorr->SetLineColor(kBlack);
    h1_temp_z_uncorr->Sumw2();
    h1_temp_z_uncorr->SetTitle("S+B");
    h1_temp_z_comb->SetLineColor(kRed);
    h1_temp_z_comb->SetTitle("B");
    h1_temp_z_comb_nobkgweight->SetLineColor(kBlue);
    h1_temp_z_comb_nobkgweight->SetTitle("sidebands");
    z_sbsub_inputs_temp->Add(h1_temp_z_uncorr);
    z_sbsub_inputs_temp->Add(h1_temp_z_comb);
    z_sbsub_inputs_temp->Add(h1_temp_z_comb_nobkgweight);
    z_sbsub_inputs_temp->Write();
  }
  
  TH2D *h2_zjt_ptbinned_uncorrected[ptbinsize-1], *h2_zjt_ptbinned[ptbinsize-1];
  TH2D *h2_zr_ptbinned_uncorrected[ptbinsize-1], *h2_zr_ptbinned[ptbinsize-1];
  TH2D *h2_jtr_ptbinned_uncorrected[ptbinsize-1], *h2_jtr_ptbinned[ptbinsize-1];       
  for (int j = 1; j < ptbinsize; j++)
  {   
    
    TH1D *h1_ptz_temp = (TH1D *)h2_ptz->ProjectionX(Form("z_pt%d", j), j + 1, j + 1); 
    NormalizeHist(h1_ptz_temp);
    h1_ptz_temp->SetStats(0);
    h1_ptz_temp->SetMarkerStyle(j + 20);
    if (j!=5) 
    {
      h1_ptz_temp->SetMarkerColor(j);
      h1_ptz_temp->SetLineColor(j);
    }
    else
    {
      h1_ptz_temp->SetMarkerColor(j*j+3);
      h1_ptz_temp->SetLineColor(j*j+3);
    }
    h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
    //h1_ptz_temp->SetOption("PE HIST");
    hs_ptz->Add(h1_ptz_temp);

    TH1D *h1_ptjt_temp = (TH1D *)h2_ptjt->ProjectionX(Form("jt_pt%d", j), j + 1, j + 1); 
    NormalizeHist(h1_ptjt_temp);        
    h1_ptjt_temp->SetStats(0);
    h1_ptjt_temp->SetMarkerStyle(j + 20);
    if (j!=5) 
    {
      h1_ptjt_temp->SetMarkerColor(j);
      h1_ptjt_temp->SetLineColor(j);
    }
    else
    {
      h1_ptjt_temp->SetMarkerColor(j*j+3);
      h1_ptjt_temp->SetLineColor(j*j+3);
    }
    h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));        
    hs_ptjt->Add(h1_ptjt_temp);
    
    TH1D *h1_ptr_temp = (TH1D *)h2_ptr->ProjectionX(Form("r_pt%d", j), j + 1, j + 1); 
    NormalizeHist(h1_ptr_temp);        
    h1_ptr_temp->SetStats(0);
    h1_ptr_temp->SetMarkerStyle(j + 20);
    if (j!=5) 
    {
      h1_ptr_temp->SetMarkerColor(j);
      h1_ptr_temp->SetLineColor(j);
    }
    else
    {
      h1_ptr_temp->SetMarkerColor(j*j+3);
      h1_ptr_temp->SetLineColor(j*j+3);
    }
    h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));        
    hs_ptr->Add(h1_ptr_temp);   
    

    // Sweighted distributions
    h1_ptz_temp = (TH1D *)h2_ptz_sweight->ProjectionX(Form("z_sweight_pt%d", j), j + 1, j + 1); 
    NormalizeHist(h1_ptz_temp);
    h1_ptz_temp->SetStats(0);
    h1_ptz_temp->SetMarkerStyle(j + 20);
    if (j!=5) 
    {
      h1_ptz_temp->SetMarkerColor(j);
      h1_ptz_temp->SetLineColor(j);
    }
    else
    {
      h1_ptz_temp->SetMarkerColor(j*j+3);
      h1_ptz_temp->SetLineColor(j*j+3);
    }
    h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));
    //h1_ptz_temp->SetOption("PE HIST");
    hs_ptz_sweight->Add(h1_ptz_temp);

    h1_ptjt_temp = (TH1D *)h2_ptjt_sweight->ProjectionX(Form("jt_sweight_pt%d", j), j + 1, j + 1); 
    NormalizeHist(h1_ptjt_temp);        
    h1_ptjt_temp->SetStats(0);
    h1_ptjt_temp->SetMarkerStyle(j + 20);
    if (j!=5) 
    {
      h1_ptjt_temp->SetMarkerColor(j);
      h1_ptjt_temp->SetLineColor(j);
    }
    else
    {
      h1_ptjt_temp->SetMarkerColor(j*j+3);
      h1_ptjt_temp->SetLineColor(j*j+3);
    }
    h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));        
    hs_ptjt_sweight->Add(h1_ptjt_temp);
    
    h1_ptr_temp = (TH1D *)h2_ptr_sweight->ProjectionX(Form("r_sweight_pt%d", j), j + 1, j + 1); 
    NormalizeHist(h1_ptr_temp);        
    h1_ptr_temp->SetStats(0);
    h1_ptr_temp->SetMarkerStyle(j + 20);
    if (j!=5) 
    {
      h1_ptr_temp->SetMarkerColor(j);
      h1_ptr_temp->SetLineColor(j);
    }
    else
    {
      h1_ptr_temp->SetMarkerColor(j*j+3);
      h1_ptr_temp->SetLineColor(j*j+3);
    }
    h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));        
    hs_ptr_sweight->Add(h1_ptr_temp);  

    if (isData)
    {
      // Add uncorrected z histos (pt binned) to THStack
      h1_ptz_temp = (TH1D *)h2_ptz_uncorrected->ProjectionX(Form("z_pt%d_raw", j), j + 1, j + 1); 
      h1_ptz_temp->SetStats(0);
      NormalizeHist(h1_ptz_temp);
      h1_ptz_temp->SetMarkerStyle(j + 20);
      if (j!=5) 
      {
        h1_ptz_temp->SetMarkerColor(j);
        h1_ptz_temp->SetLineColor(j);
      }
      else
      {
        h1_ptz_temp->SetMarkerColor(j*j+3);
        h1_ptz_temp->SetLineColor(j*j+3);
      }
      h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));                
      hs_ptz_uncorrected->Add(h1_ptz_temp);
      
      // Add uncorrected jt histos (pt binned) to THStack        
      h1_ptjt_temp = (TH1D *)h2_ptjt_uncorrected->ProjectionX(Form("jt_pt%d_raw", j), j + 1, j + 1); 
      h1_ptjt_temp->SetStats(0);
      NormalizeHist(h1_ptjt_temp);
      h1_ptjt_temp->SetMarkerStyle(j + 20);
      if (j!=5) 
      {
        h1_ptjt_temp->SetMarkerColor(j);
        h1_ptjt_temp->SetLineColor(j);
      }
      else
      {
        h1_ptjt_temp->SetMarkerColor(j*j+3);
        h1_ptjt_temp->SetLineColor(j*j+3);
      }
      h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));                
      hs_ptjt_uncorrected->Add(h1_ptjt_temp);  
      
      // Add uncorrected r histos (pt binned) to THStack
      h1_ptr_temp = (TH1D *)h2_ptr_uncorrected->ProjectionX(Form("r_pt%d_raw", j), j + 1, j + 1); 
      h1_ptr_temp->SetStats(0);
      NormalizeHist(h1_ptr_temp);
      h1_ptr_temp->SetMarkerStyle(j + 20);
      if (j!=5) 
      {
        h1_ptr_temp->SetMarkerColor(j);
        h1_ptr_temp->SetLineColor(j);
      }
      else
      {
        h1_ptr_temp->SetMarkerColor(j*j+3);
        h1_ptr_temp->SetLineColor(j*j+3);
      }
      h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));                
      hs_ptr_uncorrected->Add(h1_ptr_temp);
      
      // Add no bg sub z histos (pt binned) to THStack
      h1_ptz_temp = (TH1D *)h2_ptz_nobgsub->ProjectionX(Form("z_pt%d_nobgsub", j), j + 1, j + 1); 
      h1_ptz_temp->SetStats(0);
      NormalizeHist(h1_ptz_temp);
      h1_ptz_temp->SetMarkerStyle(j + 20);
      if (j!=5) 
      {
        h1_ptz_temp->SetMarkerColor(j);
        h1_ptz_temp->SetLineColor(j);
      }
      else
      {
        h1_ptz_temp->SetMarkerColor(j*j+3);
        h1_ptz_temp->SetLineColor(j*j+3);
      }
      h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));                
      hs_ptz_nobgsub->Add(h1_ptz_temp);
      
      // Add no bg sub jt histos (pt binned) to THStack        
      h1_ptjt_temp = (TH1D *)h2_ptjt_nobgsub->ProjectionX(Form("jt_pt%d_nobgsub", j), j + 1, j + 1); 
      h1_ptjt_temp->SetStats(0);
      NormalizeHist(h1_ptjt_temp);
      h1_ptjt_temp->SetMarkerStyle(j + 20);
      if (j!=5) 
      {
        h1_ptjt_temp->SetMarkerColor(j);
        h1_ptjt_temp->SetLineColor(j);
      }
      else
      {
        h1_ptjt_temp->SetMarkerColor(j*j+3);
        h1_ptjt_temp->SetLineColor(j*j+3);
      }
      h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));                
      hs_ptjt_nobgsub->Add(h1_ptjt_temp);  
      
      // Add no bg sub r histos (pt binned) to THStack
      h1_ptr_temp = (TH1D *)h2_ptr_nobgsub->ProjectionX(Form("r_pt%d_nobgsub", j), j + 1, j + 1); 
      h1_ptr_temp->SetStats(0);
      NormalizeHist(h1_ptr_temp);
      h1_ptr_temp->SetMarkerStyle(j + 20);
      if (j!=5) 
      {
        h1_ptr_temp->SetMarkerColor(j);
        h1_ptr_temp->SetLineColor(j);
      }
      else
      {
        h1_ptr_temp->SetMarkerColor(j*j+3);
        h1_ptr_temp->SetLineColor(j*j+3);
      }
      h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[j], pt_binedges[j + 1]));                
      hs_ptr_nobgsub->Add(h1_ptr_temp);                                    
      
        
    
      h3_ptzjt_uncorrected->GetZaxis()->SetRange(j+1, j+1);      
      h2_zjt_ptbinned_uncorrected[j-1] = (TH2D *)h3_ptzjt_uncorrected->Project3D("yx");
      h2_zjt_ptbinned_uncorrected[j-1]->SetName(Form("zjt_uncorrected_pt%d",j)); 
      NormalizeHist(h2_zjt_ptbinned_uncorrected[j-1]);	

      h3_ptzr_uncorrected->GetZaxis()->SetRange(j+1, j+1);        
      h2_zr_ptbinned_uncorrected[j-1] = (TH2D *)h3_ptzr_uncorrected->Project3D("yx");
      h2_zr_ptbinned_uncorrected[j-1]->SetName(Form("zr_uncorrected_pt%d",j));
      NormalizeHist(h2_zr_ptbinned_uncorrected[j-1]);	
      
      h3_ptjtr_uncorrected->GetZaxis()->SetRange(j+1, j+1);        
      h2_jtr_ptbinned_uncorrected[j-1] = (TH2D *)h3_ptjtr_uncorrected->Project3D("yx");
      h2_jtr_ptbinned_uncorrected[j-1]->SetName(Form("jtr_uncorrected_pt%d",j));
      NormalizeHist(h2_jtr_ptbinned_uncorrected[j-1]);	
    }  
    
    h3_ptzjt->GetZaxis()->SetRange(j+1, j+1);             
    h2_zjt_ptbinned[j-1] = (TH2D *)h3_ptzjt->Project3D("yx");
    if (isData)
      h2_zjt_ptbinned[j-1]->SetName(Form("zjt_pt%d",j));
    else
      h2_zjt_ptbinned[j-1]->SetName(Form("zjt_reco_pt%d",j));                              
    NormalizeHist(h2_zjt_ptbinned[j-1]);

    h3_ptzr->GetZaxis()->SetRange(j+1, j+1);  
    h2_zr_ptbinned[j-1] = (TH2D *)h3_ptzr->Project3D("yx");     
    if (isData)
      h2_zr_ptbinned[j-1]->SetName(Form("zr_pt%d",j));
    else
      h2_zr_ptbinned[j-1]->SetName(Form("zr_reco_pt%d",j));                    
    NormalizeHist(h2_zr_ptbinned[j-1]);

    h3_ptjtr->GetZaxis()->SetRange(j+1, j+1);
    h2_jtr_ptbinned[j-1] = (TH2D *)h3_ptjtr->Project3D("yx");    
    if (isData)
      h2_jtr_ptbinned[j-1]->SetName(Form("jtr_pt%d",j));
    else
      h2_jtr_ptbinned[j-1]->SetName(Form("jtr_reco_pt%d",j));             
    NormalizeHist(h2_jtr_ptbinned[j-1]);        
  }  

  h1_nSPDHits->Write();
  if (isData)
  {
      hs_ptz->SetName("z_evtweights_sbsub_data_all");
      hs_ptjt->SetName("jt_evtweights_sbsub_data_all");
      hs_ptr->SetName("r_evtweights_sbsub_data_all");    
      //hs_ptz_sweight->SetName("z_evtweights_sbsub_data_sweight_all");
      //hs_ptjt_sweight->SetName("jt_evtweights_sbsub_data_sweight_all");
      //hs_ptr_sweight->SetName("r_evtweights_sbsub_data_sweight_all");   
      hs_ptz_uncorrected->Write();
      hs_ptjt_uncorrected->Write();   
      hs_ptr_uncorrected->Write();      
      hs_ptz_nobgsub->Write();
      hs_ptjt_nobgsub->Write();
      hs_ptr_nobgsub->Write();
  }
  else
  {
      hs_ptz->SetName("z_reco_all");
      hs_ptjt->SetName("jt_reco_all");
      hs_ptr->SetName("r_reco_all");
  }

  hs_ptz->Write();
  hs_ptjt->Write();
  hs_ptr->Write();
  hs_ptz_sweight->Write();
  hs_ptjt_sweight->Write();
  hs_ptr_sweight->Write();

  h3_ptzjt->GetZaxis()->SetRange(1, ptbinsize+1);
  h3_ptzr->GetZaxis()->SetRange(1, ptbinsize+1);
  h3_ptjtr->GetZaxis()->SetRange(1, ptbinsize+1);
  h3_ptzjt_uncorrected->GetZaxis()->SetRange(1, ptbinsize+1);
  h3_ptzr_uncorrected->GetZaxis()->SetRange(1, ptbinsize+1);
  h3_ptjtr_uncorrected->GetZaxis()->SetRange(1, ptbinsize+1);    
        
  SetRecoStyle(h1_jet_eta);
  SetDataStyle(h1_jet_rap);
  SetRecoStyle(h1_jet_pt);
  SetRecoStyle(h1_jet_phi);

  SetRecoStyle(h1_HF_rap);
  SetRecoStyle(h1_HF_pt);
  SetRecoStyle(h1_HF_mass);
  SetRecoStyle(h1_HF_phi);
    
  SetRecoStyle(h1_z);
  SetRecoStyle(h1_jt);
  SetRecoStyle(h1_r);
    
  SetRecoStyle(h1_Jpsi_rap);
  SetRecoStyle(h1_Jpsi_pt);
  SetRecoStyle(h1_Jpsi_mass);
  SetRecoStyle(h1_Jpsi_ipchi2);

  SetRecoStyle(h1_z_comb);
  SetRecoStyle(h1_jt_comb);
  SetRecoStyle(h1_r_comb);

  f.Write();
  f.Close();
}
