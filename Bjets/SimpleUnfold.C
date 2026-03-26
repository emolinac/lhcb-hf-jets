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
                  int dataset = 91599,
                  bool SubtractGS = false,                  
                  bool DoJESJER = false,
                  bool DoJetID = false,                  
                  bool DoRecSelEff = false,
                  bool DoSignalSys = false,
                  bool DoUnfoldPrior = false,
                  bool sPlotFit = false,
                  bool L0MuonDiMuon = false,
                  bool onlysim9 = false)
{
        followHardest = false;
        
        int ptRange = (dataset / 10) % 10;
        int Mag = (dataset / 1) % 10;
        
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
        TString str_year = "2016";
        TString str_GS = "";
        
        if (SubtractGS)
                str_GS = "_GSsub";
        
        if (Mag == 0)
                str_Mag = "_MD";
        else if (Mag == 1)
                str_Mag = "_MU";
        
        if (PID_cut)
                str_PID = "_PID";
        
        str_level = "reco";
        
        TString str_L0 = "";
        
        TString str_sPlot = "";
        
        TString str_simversion = "";     
        
        TString str_tree;   
        TString extension = TString("unfold_") + str_level + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_PID + str_GS + str_sPlot + str_L0 + str_simversion + Form("_%d", dataset);
        TString extension_prefix;

        if (DoJESJER)
                extension_prefix = TString("JESJER_");
        if (DoJetID)
                extension_prefix = TString("jetid_");
        if (DoRecSelEff)
                extension_prefix = TString("recselsys_");
        if (DoSignalSys)
                extension_prefix = TString("signalsys_");
        if (DoUnfoldPrior)
                extension_prefix = TString("priorsys_");
        
        extension = extension_prefix + extension;

        TString extension_RootFilesMC   = "../../root_files/BjetsMC/";
        TString extension_RootFilesData = "../../root_files/Bjets/";    
        
        TString extension_read = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_L0 + str_simversion + Form("_%d", dataset);
        TString extension_eff  = TString("efficiency_truth") + Form("_ev_%d", NumEvts) + Form("_ptj_%d%d", int(pTLow), int(ptMax)) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_Mag + str_flavor + str_GS + str_simversion + Form("_%d", dataset);
        
        TFile *file_eff = new TFile((output_folder + "bjets_efficiencies.root").c_str(), "READ");
        
        // TFile *file_decay = new TFile(extension_RootFilesMC + "HFeff_truth_ev_-1" + str_Mag + str_flavor + str_PID + Form("_%d.root", dataset), "READ");
        // TFile *file_decay = new TFile("HFeff_truth_ev_-1_b_PID_91599.root" , "READ");
        
        TChain *BTree = new TChain("BTree", "");
        BTree->Add((output_folder + "ntuple_bjets_mcreco.root/BTree").c_str());
        
        if (NumEvts > BTree->GetEntries())
                NumEvts = BTree->GetEntries();
        if (NumEvts == -1)
                NumEvts = BTree->GetEntries();
        
        TH1D *h1_denom_efficiency_HFpt = (TH1D*) file_eff->Get("denom_efficiency_HFpt");
        TH1D *h1_denom_efficiency_z    = (TH1D*) file_eff->Get("denom_efficiency_z");
        TH1D *h1_denom_efficiency_jt   = (TH1D*) file_eff->Get("denom_efficiency_jt");
        TH1D *h1_denom_efficiency_r    = (TH1D*) file_eff->Get("denom_efficiency_r");
        
        //    TH2D *h2_denom_efficiency_fragzjetpt = (TH2D *)file_eff->Get("h2_denom_efficiency_fragzjetpt");
        TH1D *h1_denom_efficiency_jetpt  = (TH1D*) file_eff->Get("denom_efficiency_jetpt");
        TH1D *h1_denom_efficiency_jetrap = (TH1D*) file_eff->Get("denom_efficiency_jetrap");
        
        TH2D *h2_denom_efficiency_HFptjetpt = (TH2D*) file_eff->Get("denom_efficiency_HFptjetpt");
        TH2D *h2_denom_efficiency_HFpteta   = (TH2D*) file_eff->Get("denom_efficiency_HFpteta");
        TH2D *h2_denom_efficiency_jetpteta  = (TH2D*) file_eff->Get("denom_efficiency_jetpteta");
        TH2D *h2_denom_efficiency_zjt       = (TH2D*) file_eff->Get("denom_efficiency_zjt");
        TH2D *h2_denom_efficiency_zr        = (TH2D*) file_eff->Get("denom_efficiency_zr");
        TH2D *h2_denom_efficiency_jtr       = (TH2D*) file_eff->Get("denom_efficiency_jtr");
        
        TH3D *h3_denom_efficiency_HFptetajetpt = (TH3D*) file_eff->Get("denom_efficiency_HFptetajetpt");
        
        TH2D *h2_denom_efficiency_ptz  = (TH2D*) file_eff->Get("denom_efficiency_ptz");
        TH2D *h2_denom_efficiency_ptjt = (TH2D*) file_eff->Get("denom_efficiency_ptjt");
        TH2D *h2_denom_efficiency_ptr  = (TH2D*) file_eff->Get("denom_efficiency_ptr");
        
        TH3D *h3_denom_efficiency_ptzjt = (TH3D*) file_eff->Get("denom_efficiency_ptzjt");
        TH3D *h3_denom_efficiency_ptzr  = (TH3D*) file_eff->Get("denom_efficiency_ptzr");
        TH3D *h3_denom_efficiency_ptjtr = (TH3D*) file_eff->Get("denom_efficiency_ptjtr");  
        
        TH2D *h2_SVTag_eff_denom_z = (TH2D*) file_eff->Get("h2_SVTag_eff_denom_z");
        TH2D *h2_SVTag_eff_denom   = (TH2D*) file_eff->Get("h2_SVTag_eff_denom");
        
        TH1D *h1_recoHF_HFpt         = (TH1D*) file_eff->Get("num_efficiency_HFpt");
        TH1D *h1_recoHF_jetpt        = (TH1D*) file_eff->Get("num_efficiency_jetpt");
        TH2D *h2_recoHF_HFpteta      = (TH2D*) file_eff->Get("num_efficiency_HFpteta");
        TH3D *h3_recoHF_HFptetajetpt = (TH3D*) file_eff->Get("num_efficiency_HFptetajetpt");
        
        TH1D *h1_HFeff_HFpt         = (TH1D*) file_eff->Get("efficiency_HFpt");
        TH1D *h1_HFeff_jetpt        = (TH1D*) file_eff->Get("efficiency_jetpt");
        TH2D *h2_HFeff_HFpteta      = (TH2D*) file_eff->Get("efficiency_HFpteta");
        TH3D *h3_HFeff_HFptetajetpt = (TH3D*) file_eff->Get("efficiency_HFptetajetpt");
        
        //    /////////////////// Mass Fit Parameters /////////////////////////////////
        TString extension_mass("");
        if (sPlotFit)
                extension_mass = TString("splotfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.4_b" + str_PID + str_L0 + "_91599.root";
        else
                extension_mass = TString("massfit_data_ev_-1") + Form("_ptj_%d%d", int(pTLow), int(250.)) + "_eta_2.54.0_ghost_0.4_b" + str_PID + str_L0 + "_91599.root";
        
        if (DoRecSelEff)
                extension_mass = "recselsys_" + extension_mass;

        if (DoSignalSys)
                extension_mass = "sys_" + extension_mass;

        extension_mass = extension_RootFilesData + extension_mass;
        
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

                h3_ptzjt_ratio = (TH3D *)file_reco_weights->Get("ptzjt_ratio");
        }
        
        TFile *f = TFile::Open((output_folder + "bjets_corrections.root").c_str(), "RECREATE");
        
        //Reco 1D Observables
        TH1D *h1_z  = new TH1D("z", "", zbinsize, z_binedges);
        TH1D *h1_jt = new TH1D("jt", "", jtbinsize, jt_binedges);
        TH1D *h1_r  = new TH1D("r", "", rbinsize, r_binedges);
        
        //Truth 1D Observables
        TH1D *h1_z_truth  = new TH1D("z_truth", "", zbinsize, z_binedges);
        TH1D *h1_jt_truth = new TH1D("jt_truth", "", jtbinsize, jt_binedges);
        TH1D *h1_r_truth  = new TH1D("r_truth", "", rbinsize, r_binedges);
        
        //Reco 2D Observables (208 - 217)
        TH2D *h2_zjt = new TH2D("zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        TH2D *h2_zr  = new TH2D("zr", "", zbinsize, z_binedges,  rbinsize, r_binedges);
        TH2D *h2_jtr = new TH2D("jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
        
        //Truth 2D Observables
        TH2D *h2_zjt_truth = new TH2D("zjt_truth", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        TH2D *h2_zr_truth  = new TH2D("zr_truth", "", zbinsize, z_binedges,  rbinsize, r_binedges);
        TH2D *h2_jtr_truth = new TH2D("jtr_truth", "", jtbinsize, jt_binedges,  rbinsize, r_binedges);
        
        // 2D Truth-Reco Correspondence (219 - 224)
        TH2D *h2_truthreco_z  = new TH2D("truthreco_z", ";Reco z; Truth z", zbinsize, z_binedges, zbinsize, z_binedges);
        TH2D *h2_truthreco_jt = new TH2D("truthreco_jt", ";Reco jT; Truth jT", jtbinsize, jt_binedges, jtbinsize, jt_binedges);
        TH2D *h2_truthreco_r  = new TH2D("truthreco_r", ";Reco r; Truth r", rbinsize, r_binedges, rbinsize, r_binedges);
        ////
        TH2D *h2_truthreco_jetpt = new TH2D("truthreco_jetpt", "", longptbinsize, longpt_binedges, longptbinsize, longpt_binedges);
        //////
        
        // 1D Num Efficiencies and Purities of Observables (226 - 235)
        TH1D *h1_num_purity_z  = new TH1D("num_purity_z", "", zbinsize, z_binedges);
        TH1D *h1_num_purity_jt = new TH1D("num_purity_jt", "", jtbinsize, jt_binedges );
        TH1D *h1_num_purity_r  = new TH1D("num_purity_r", "", rbinsize, r_binedges);
        
        TH1D *h1_num_efficiency_z  = new TH1D("num_efficiency_z", "", zbinsize, z_binedges);
        TH1D *h1_num_efficiency_jt = new TH1D("num_efficiency_jt", "", jtbinsize, jt_binedges );
        TH1D *h1_num_efficiency_r  = new TH1D("num_efficiency_r", "", rbinsize, r_binedges);
        
        // 1D Denom Efficiencies and Purities of Observables (237 - 246)
        TH1D *h1_denom_purity_z  = new TH1D("denom_purity_z", "", zbinsize, z_binedges);
        TH1D *h1_denom_purity_jt = new TH1D("denom_purity_jt", "", jtbinsize, jt_binedges );
        TH1D *h1_denom_purity_r  = new TH1D("denom_purity_r", "", rbinsize, r_binedges);
        //    TH1D *h1_denom_efficiency_z = new TH1D("denom_efficiency_z", "", zbinsize, z_binedges);
        //    TH1D *h1_denom_efficiency_jt = new TH1D("denom_efficiency_jt", "", jtbinsize, jt_binedges );
        //    TH1D *h1_denom_efficiency_r = new TH1D("denom_efficiency_r", "", rbinsize, r_binedges);
        
        TH1D *h1_num_purity_jetpt   = new TH1D("num_purity_jetpt", "", ptbinsize, pt_binedges);
        TH1D *h1_denom_purity_jetpt = new TH1D("denom_purity_jetpt", "", ptbinsize, pt_binedges);
        TH1D *h1_num_purity_HFpt    = new TH1D("num_purity_HFpt", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_denom_purity_HFpt  = new TH1D("denom_purity_HFpt", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_num_purity_rap     = new TH1D("num_purity_rap", "", etabinsize, eta_binedges);
        TH1D *h1_denom_purity_rap   = new TH1D("denom_purity_rap", "", etabinsize, eta_binedges);
        
        TH1D *h1_num_efficiency_HFpt      = new TH1D("num_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_num_efficiency_jetpt     = new TH1D("num_efficiency_jetpt", "", ptbinsize, pt_binedges);
        TH1D *h1_num_efficiency_jetrap    = new TH1D("num_efficiency_jetrap", "", ptbinsize, pt_binedges);
        TH2D *h2_num_efficiency_HFptjetpt = new TH2D("num_efficiency_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
        TH2D *h2_num_purity_HFptjetpt     = new TH2D("num_purity_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
        TH2D *h2_denom_purity_HFptjetpt   = new TH2D("denom_purity_HFptjetpt", "", ptHFbinsize, ptHF_binedges, customptbinsize, custompt_binedges);
        
        TH2D *h2_num_purity_HFpteta     = new TH2D("num_purity_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
        TH2D *h2_denom_purity_HFpteta   = new TH2D("denom_purity_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
        TH2D *h2_num_efficiency_HFpteta = new TH2D("num_efficiency_HFpteta", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges);
        
        TH3D *h3_num_purity_HFptetajetpt     = new TH3D("num_purity_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, customptbinsize, custompt_binedges);
        TH3D *h3_denom_purity_HFptetajetpt   = new TH3D("denom_purity_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, customptbinsize, custompt_binedges);
        TH3D *h3_num_efficiency_HFptetajetpt = new TH3D("num_efficiency_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, customptbinsize, custompt_binedges);
        
        TH2D *h2_denom_purity_jetpteta = new TH2D("denom_purity_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
        TH2D *h2_num_purity_jetpteta   = new TH2D("num_purity_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
        
        TH2D *h2_num_efficiency_jetpteta = new TH2D("num_efficiency_jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
        
        //  TEfficiency *t1_efficiency_dR = new TEfficiency("t1_efficiency_dR", "", dRbinsize, dR_binedges);
        
        // TH1D *h1_jetpt = new TH1D("jetpt", "", ptbinsize, pt_binedges);
        // TH1D *h1_jetpt_truth = new TH1D("jetpt_truth", "", ptbinsize_truth, pt_binedges_truth);
        // TH1D *h1_jetpt_test = new TH1D("jetpt_test", "", ptbinsize, pt_binedges);
        
        TH1D *h1_jetpt       = new TH1D("jetpt", "", ptbinsize, pt_binedges);
        TH1D *h1_jetpt_truth = new TH1D("jetpt_truth", "", ptbinsize, pt_binedges);
        TH1D *h1_jetpt_test  = new TH1D("jetpt_test", "", ptbinsize, pt_binedges);
        
        TH2D *h2_jetpteta = new TH2D("jetpteta", "", ptbinsize, pt_binedges, etabinsize, eta_binedges);
        
        RooUnfoldResponse *response_jetpt     = new RooUnfoldResponse(h1_jetpt, h1_jetpt_truth, "response_jetpt");
        RooUnfoldResponse *response_jeteta    = new RooUnfoldResponse(h1_num_purity_rap, h1_num_purity_rap, "response_jeteta");
        RooUnfoldResponse *response_HFptjetpt = new RooUnfoldResponse(h2_denom_efficiency_HFptjetpt, h2_denom_efficiency_HFptjetpt, "response_HFptjetpt");
        RooUnfoldResponse *response_jetpteta  = new RooUnfoldResponse(h2_jetpteta, h2_jetpteta, "response_jetpteta");
        
        // 2D Efficiencies and Purities Observables
        TH2D *h2_num_purity_zjt   = new TH2D("num_purity_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        TH2D *h2_denom_purity_zjt = new TH2D("denom_purity_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        
        TH2D *h2_num_purity_zr   = new TH2D("num_purity_zr", "", zbinsize, z_binedges, rbinsize, r_binedges);
        TH2D *h2_denom_purity_zr = new TH2D("denom_purity_zr", "", zbinsize, z_binedges, rbinsize, r_binedges);
        
        TH2D *h2_num_purity_jtr   = new TH2D("num_purity_jtr", "", jtbinsize, jt_binedges, rbinsize, r_binedges);
        TH2D *h2_denom_purity_jtr = new TH2D("denom_purity_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
        
        TH2D *h2_num_efficiency_zjt = new TH2D("num_efficiency_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        //    TH2D *h2_denom_efficiency_zjt = new TH2D("denom_efficiency_zjt", "", zbinsize, z_binedges, jtbinsize, jt_binedges);
        
        TH2D *h2_num_efficiency_zr = new TH2D("num_efficiency_zr", "",  zbinsize, z_binedges, rbinsize, r_binedges);
        //    TH2D *h2_denom_efficiency_zr = new TH2D("denom_efficiency_zr", "",  zbinsize, z_binedges, rbinsize, r_binedges);
        
        TH2D *h2_num_efficiency_jtr = new TH2D("num_efficiency_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
        //    TH2D *h2_denom_efficiency_jtr = new TH2D("denom_efficiency_jtr", "",  jtbinsize, jt_binedges, rbinsize, r_binedges);
        
        // 2D and 4D Observables response Matrices
        RooUnfoldResponse *response_zjt = new RooUnfoldResponse(h2_zjt, h2_zjt_truth, "response_zjt" );
        RooUnfoldResponse *response_zr  = new RooUnfoldResponse(h2_zr, h2_zr_truth, "respose_zr");
        RooUnfoldResponse *response_jtr = new RooUnfoldResponse(h2_jtr , h2_jtr_truth, "respose_jtr");
        RooUnfoldResponse *response_z   = new RooUnfoldResponse(h1_denom_purity_z, h1_denom_efficiency_z, "response_z");
        RooUnfoldResponse *response_jt  = new RooUnfoldResponse(h1_denom_purity_jt, h1_denom_efficiency_jt, "response_jt");
        RooUnfoldResponse *response_r   = new RooUnfoldResponse(h1_denom_purity_r, h1_denom_efficiency_r, "response_r");
        
        
        // 3D Efficiencies/Purities of 2D Observables vs jet_pt; 6D Response Matrices (324 - 357)
        /// Q: Is Ibrahim's "matched obsv." is my "truth obsv." ?
        
        TH3D *h3_ptzjt                = new TH3D("ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_truthptzjt           = new TH3D("truthptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_num_efficiency_ptzjt = new TH3D("num_efficiency_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_efficiency_ptzjt     = new TH3D("efficiency_ptzjt", "efficiency;z; j_{T} [GeV/c];p_{T,jet} [GeV/c];", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_num_purity_ptzjt     = new TH3D("num_purity_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_denom_purity_ptzjt   = new TH3D("denom_purity_ptzjt", "", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_purity_ptzjt         = new TH3D("purity_ptzjt", "purity;z; j_{T} [GeV/c];p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, jtbinsize2D, jt_binedges2D, ptbinsize, pt_binedges );
        
        RooUnfoldResponse *response_ptzjt = new RooUnfoldResponse(h3_ptzjt, h3_truthptzjt, "response_ptzjt");
        
        TH3D *h3_ptzr                = new TH3D("ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_truthptzr           = new TH3D("truthptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_num_efficiency_ptzr = new TH3D("num_efficiency_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_efficiency_ptzr     = new TH3D("efficiency_ptzr", "efficiency;z;r;p_{T,jet} [GeV/c]", zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_num_purity_ptzr     = new TH3D("num_purity_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_denom_purity_ptzr   = new TH3D("denom_purity_ptzr", "",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_purity_ptzr         = new TH3D("purity_ptzr", "purity;z;r;p_{T,jet} [GeV/c]",  zbinsize2D, z_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        
        RooUnfoldResponse *response_ptzr = new RooUnfoldResponse(h3_ptzr, h3_truthptzr, "response_ptzr");   
        
        TH3D *h3_ptjtr                = new TH3D("ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_truthptjtr           = new TH3D("truthptjtr", "", jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_num_efficiency_ptjtr = new TH3D("num_efficiency_ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_efficiency_ptjtr     = new TH3D("efficiency_ptjtr", "efficiency;j_{T} [GeV/c];r;p_{T,jet} [GeV/c];",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_num_purity_ptjtr     = new TH3D("num_purity_ptjtr", "",  jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_denom_purity_ptjtr   = new TH3D("denom_purity_ptjtr", "", jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D, ptbinsize, pt_binedges );
        TH3D *h3_purity_ptjtr         = new TH3D("purity_ptjtr", "purity;j_{T} [GeV/c];r;p_{T,jet} [GeV/c]", jtbinsize2D, jt_binedges2D, rbinsize2D, r_binedges2D,  ptbinsize, pt_binedges );
        
        RooUnfoldResponse *response_ptjtr = new RooUnfoldResponse(h3_ptjtr, h3_truthptjtr, "response_ptjtr");
        
        // 2D Jet_pt vs Observable histograms, Purities/Efficiencies, Responses 
        TH2D *h2_ptz                = new TH2D("ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_truthptz           = new TH2D("truthptz", "", zbinsize, z_binedges, ptbinsize, pt_binedges );
        TH2D *h2_num_efficiency_ptz = new TH2D("num_efficiency_ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_efficiency_ptz     = new TH2D("efficiency_ptz", "efficiency;z;p_{T,jet} [GeV/c]", zbinsize, z_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_num_purity_ptz     = new TH2D("num_purity_ptz", "", zbinsize, z_binedges, ptbinsize, pt_binedges );
        TH2D *h2_denom_purity_ptz   = new TH2D("denom_purity_ptz", "", zbinsize, z_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_purity_ptz         = new TH2D("purity_ptz", "purity;z;p_{T,jet} [GeV/c]", zbinsize, z_binedges,  ptbinsize, pt_binedges );
        
        RooUnfoldResponse *response_ptz = new RooUnfoldResponse(h2_ptz, h2_truthptz, "response_ptz");
        
        TH2D *h2_ptjt                = new TH2D("ptjt", "", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_truthptjt           = new TH2D("truthptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges );
        TH2D *h2_num_efficiency_ptjt = new TH2D("num_efficiency_ptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
        TH2D *h2_efficiency_ptjt     = new TH2D("efficiency_ptjt", "efficiency;j_{T} [GeV/c];p_{T,jet} [GeV/c]", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_num_purity_ptjt     = new TH2D("num_purity_ptjt", "", jtbinsize, jt_binedges, ptbinsize, pt_binedges );
        TH2D *h2_denom_purity_ptjt   = new TH2D("denom_purity_ptjt", "", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_purity_ptjt         = new TH2D("purity_ptjt", "purity;j_{T} [GeV/c];p_{T,jet} [GeV/c]", jtbinsize, jt_binedges,  ptbinsize, pt_binedges );
        
        RooUnfoldResponse *response_ptjt = new RooUnfoldResponse(h2_ptjt, h2_truthptjt, "response_ptjt");
        
        TH2D *h2_ptr                = new TH2D("ptr", "", rbinsize, r_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_truthptr           = new TH2D("truthptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges );
        TH2D *h2_num_efficiency_ptr = new TH2D("num_efficiency_ptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges);
        TH2D *h2_efficiency_ptr     = new TH2D("efficiency_ptr", "efficiency;r;p_{T,jet} [GeV/c]", rbinsize, r_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_num_purity_ptr     = new TH2D("num_purity_ptr", "", rbinsize, r_binedges, ptbinsize, pt_binedges );
        TH2D *h2_denom_purity_ptr   = new TH2D("denom_purity_ptr", "", rbinsize, r_binedges,  ptbinsize, pt_binedges );
        TH2D *h2_purity_ptr         = new TH2D("purity_ptr", "purity;r;p_{T,jet} [GeV/c]", rbinsize, r_binedges,  ptbinsize, pt_binedges );
        
        RooUnfoldResponse *response_ptr = new RooUnfoldResponse(h2_ptr, h2_truthptr, "response_ptr"); 
        
        TH1D *h1_z_true      = new TH1D("h1_z_true", "", 20, 0, 1.05);
        TH1D *h1_z_false     = new TH1D("h1_z_false", "", 20, 0, 1.05);
        TH1D *h1_z_truefalse = new TH1D("h1_z_truefalse", "", 20, 0, 1.05);
        
        TH1D *h1_r_true  = new TH1D("h1_r_true", "", 20, 0, 0.5);
        TH1D *h1_r_false = new TH1D("h1_r_false", "", 20, 0, 0.5);
        
        TH1D *h1_chi2_true  = new TH1D("h1_chi2_true", "", 20, -4., 4.);
        TH1D *h1_chi2_false = new TH1D("h1_chi2_false", "", 20, -4., 4.);
        
        TH2D *h2_sv_mass_z_truepos     = new TH2D("h2_sv_mass_z_truepos", "", 30, 0, 4., 20, 0, 1.05);
        TH2D *h2_sv_mass_ntrks_truepos = new TH2D("h2_sv_mass_ntrks_truepos", "", 30, 0, 4., 10, -0.5, 9.5);
        
        TH2D *h2_sv_mass_z_falsepos     = new TH2D("h2_sv_mass_z_falsepos", "", 30, 0, 4., 20, 0, 1.05);
        TH2D *h2_sv_mass_ntrks_falsepos = new TH2D("h2_sv_mass_ntrks_falsepos", "", 30, 0, 4., 10, -0.5, 9.5);
        
        TH2D *h2_ptz_weighted  = new TH2D("ptz_weighted", ";z;p_{T} (GeV/c)", zbinsize, z_binedges, ptbinsize, pt_binedges);
        TH2D *h2_ptjt_weighted = new TH2D("ptjt_weighted", ";j_{T} (GeV/c); p_{T} (GeV/c)", jtbinsize, jt_binedges, ptbinsize, pt_binedges);
        TH2D *h2_ptr_weighted  = new TH2D("ptr_weighted", ";r;p_{T} (GeV/c)", rbinsize, r_binedges, ptbinsize, pt_binedges);
        
        TH1D *h1_z_weighted  = new TH1D("z_weighted", ";z;", zbinsize, z_binedges);
        TH1D *h1_jt_weighted = new TH1D("jt_weighted", ";j_{T} (GeV/c);", jtbinsize, jt_binedges);
        TH1D *h1_r_weighted  = new TH1D("r_weighted", ";r;", rbinsize, r_binedges);
        
        //
        // Event loop
        //

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
        
        double z, jt, r;
        double tr_z, tr_jt, tr_r;
        double zg, jtg, rg;
        
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
        
        BTree->SetBranchAddress("jet_pt", &jet_pt);
        BTree->SetBranchAddress("jet_rap", &jet_rap);
        BTree->SetBranchAddress("jet_px", &jet_px);
        BTree->SetBranchAddress("jet_py", &jet_py);
        BTree->SetBranchAddress("jet_pz", &jet_pz);
        BTree->SetBranchAddress("jet_e", &jet_e);
        
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
        
        BTree->SetBranchAddress("z", &z);
        BTree->SetBranchAddress("jt", &jt);
        BTree->SetBranchAddress("r", &r);
        BTree->SetBranchAddress("zg", &zg);
        BTree->SetBranchAddress("jtg", &jtg);
        BTree->SetBranchAddress("rg", &rg);
        BTree->SetBranchAddress("tr_z", &tr_z);
        BTree->SetBranchAddress("tr_jt", &tr_jt);
        BTree->SetBranchAddress("tr_r", &tr_r);
        
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
                
                if (DoUnfoldPrior)
                        prior_weight = h3_ptzjt_ratio->GetBinContent(h3_ptzjt_ratio->GetXaxis()->FindBin(z),
                                                                     h3_ptzjt_ratio->GetYaxis()->FindBin(jt),
                                                                     h3_ptzjt_ratio->GetZaxis()->FindBin(jet_pt));
                
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
        
                ////////////////////////////////////////////////////////////////////
                // Two methods to compute efficiency:
                ///////////////////////////////////////////////////////////////////

                ///////////////////////////////////////////////////////////////////
                // Method 1:
                // Fill denominator with all truth B jets in fiducial region in MCSimpleEff
                // Fill numerator with all reconstructed B jets in fiducial region in MCSimpleEff
                // Take the ratio as the efficiency
                // Problem: This efficiency is limited, cannot apply trigger in the MCTruth Tree
                /////////////////////////////////////////////////////////////////////

                /////////////////////////////////////////////////////////////////////
                // Method 2:
                // Fill denominator with all truth B jets in fiducial region in MCSimpleEff
                // Fill numerator with all reconstructed B jets in fiducial region in SimpleUnfold
                // Advantage: more freedom with variables (like trigger) in MCReco Tree
                /////////////////////////////////////////////////////////////////////
                
                bool num_cond_eff = tr_pt_cond && tr_rap_cond && isTrueBjet && pt_cond && rap_cond;
                bool num_cond_pur = tr_pt_cond && tr_rap_cond && isTrueBjet && pt_cond && rap_cond;
                bool denom_cond   = pt_cond && rap_cond;
                
                
                if (num_cond_eff) {
                        if (eventNumber == 6176925)
                                cout << jet_pt << ", " << jet_rap << ", " << tr_jet_pt << ", " << tr_jet_rap << endl;
                        
                        h1_num_efficiency_jetpt->Fill(tr_jet_pt, eff_weight);
                                        
                        h2_num_efficiency_ptz->Fill(tr_z, tr_HFjet.Pt(), eff_weight);
                        h2_num_efficiency_ptjt->Fill(tr_jt, tr_HFjet.Pt(), eff_weight);
                        h2_num_efficiency_ptr->Fill(tr_r, tr_HFjet.Pt(), eff_weight);
                                
                        h3_num_efficiency_ptzjt->Fill(tr_z, tr_jt, tr_HFjet.Pt(), eff_weight);
                        h3_num_efficiency_ptzr->Fill(tr_z, tr_r, tr_HFjet.Pt(), eff_weight);
                        h3_num_efficiency_ptjtr->Fill( tr_jt, tr_r, tr_HFjet.Pt(), eff_weight);
                }
                
                if (denom_cond) {
                        h1_denom_purity_jetpt->Fill(jet_pt, eff_weight);
                        
                        h2_denom_purity_ptz->Fill(z, jet_pt, eff_weight);
                        h2_denom_purity_ptjt->Fill(jt, jet_pt, eff_weight);
                        h2_denom_purity_ptr->Fill(r, jet_pt, eff_weight);
                        
                        h3_denom_purity_ptzjt->Fill(z, jt, jet_pt, eff_weight);
                        h3_denom_purity_ptzr->Fill(z, r, jet_pt, eff_weight);
                        h3_denom_purity_ptjtr->Fill(jt, r, jet_pt, eff_weight);
                }    

                if (num_cond_pur) {
                        h1_num_purity_jetpt->Fill(jet_pt, eff_weight);
                        
                        h2_num_purity_ptz->Fill(z, jet_pt, eff_weight);
                        h2_num_purity_ptjt->Fill(jt, jet_pt, eff_weight);
                        h2_num_purity_ptr->Fill(r, jet_pt, eff_weight);
                        
                        h3_num_purity_ptzjt->Fill(z, jt, jet_pt, eff_weight);
                        h3_num_purity_ptzr->Fill(z, r, jet_pt, eff_weight);
                        h3_num_purity_ptjtr->Fill(jt, r, jet_pt, eff_weight);  
                        NumTrueBjets++;              
                }
                
                if (!pt_cond || !tr_pt_cond || !isTrueBjet || !rap_cond || !tr_rap_cond )
                        continue;
                
                response_jetpt->Fill(jet_pt, tr_jet_pt, prior_weight);

                response_ptz->Fill( z,  jet_pt, tr_z, tr_jet_pt, prior_weight);        
                response_ptjt->Fill( jt, jet_pt, tr_jt, tr_jet_pt, prior_weight);
                response_ptr->Fill( r,  jet_pt, tr_r, tr_jet_pt, prior_weight);
                
                response_ptzjt->Fill( z, jt, jet_pt, tr_z, tr_jt, tr_jet_pt, prior_weight);
                response_ptzr->Fill( z, r, jet_pt, tr_z, tr_r, tr_jet_pt, prior_weight);
                response_ptjtr->Fill( jt,  r, jet_pt, tr_jt, tr_r, tr_jet_pt, prior_weight);  
                
                // if(jet_pt > ptMax) continue;
                h2_truthreco_jetpt->Fill(jet_pt, tr_jet_pt, prior_weight);
                h2_truthreco_z->Fill(z, tr_z, prior_weight);
                h2_truthreco_jt->Fill(jt, tr_jt, prior_weight);
                h2_truthreco_r->Fill(r, tr_r, prior_weight);
                
                // For a cross check on the unfold prior systematic
                h2_ptz_weighted->Fill(z, jet_pt, prior_weight);
                h2_ptjt_weighted->Fill(jt, jet_pt, prior_weight);
                h2_ptr_weighted->Fill(r, jet_pt, prior_weight);  
                h1_z_weighted->Fill(z, prior_weight);
                h1_jt_weighted->Fill(jt, prior_weight);
                h1_r_weighted->Fill(r, prior_weight);       
                
                h2_ptz->Fill(z, jet_pt);
                h2_ptjt->Fill(jt, jet_pt);
                h2_ptr->Fill(r, jet_pt);  
                h1_z->Fill(z);
                h1_jt->Fill(jt);
                h1_r->Fill(r);   
                                                                
                event_counter++;
        }
                
        // RooUnfoldResponse * response_jetpt = new RooUnfoldResponse(h1_jetpt, h1_jetpt_truth, h2_truthreco_jetpt, "response_jetpt");
                
        TH2 *h2_response_jetpt = response_jetpt->Hresponse();
        
        TH2 *h2_response_ptz  =  response_ptz->Hresponse();        
        TH2 *h2_response_ptjt = response_ptjt->Hresponse();
        TH2 *h2_response_ptr  = response_ptr->Hresponse();
                
        TH2 *h2_response_ptzjt = response_ptzjt->Hresponse();
        TH2 *h2_response_ptzr  = response_ptzr->Hresponse();
        TH2 *h2_response_ptjtr = response_ptjtr->Hresponse();

        h2_response_jetpt->GetXaxis()->SetTitle("reco p_{T}^{j} [GeV/c]");
        h2_response_jetpt->GetYaxis()->SetTitle("truth p_{T}^{j} [GeV/c]");
        h2_response_ptz->GetXaxis()->SetTitle("reco bin # (flattened)");
        h2_response_ptz->GetYaxis()->SetTitle("truth bin # (flattened)");
        h2_response_ptjt->GetXaxis()->SetTitle("reco bin # (flattened)");
        h2_response_ptjt->GetYaxis()->SetTitle("truth bin # (flattened)");
        h2_response_ptr->GetXaxis()->SetTitle("reco bin # (flattened)");
        h2_response_ptr->GetYaxis()->SetTitle("truth bin # (flattened)");
        h2_response_ptzjt->GetXaxis()->SetTitle("reco bin # (flattened)");
        h2_response_ptzjt->GetYaxis()->SetTitle("truth bin # (flattened)");
        h2_response_ptzr->GetXaxis()->SetTitle("reco bin # (flattened)");
        h2_response_ptzr->GetYaxis()->SetTitle("truth bin # (flattened)");
        h2_response_ptjtr->GetXaxis()->SetTitle("reco bin # (flattened)");
        h2_response_ptjtr->GetYaxis()->SetTitle("truth bin # (flattened)");

        h2_response_jetpt->Write("response_jetpt");
        
        h2_response_ptz->Write("response_ptz");
        h2_response_ptjt->Write("response_ptjt");
        h2_response_ptr->Write("response_ptr");    

        h2_response_ptzjt->Write("response_ptzjt");
        h2_response_ptzr->Write("response_ptzr");
        h2_response_ptjtr->Write("response_ptjtr");     
                        
        response_jetpt->Write("Roo_response_jetpt");
        
        response_ptz->Write("Roo_response_ptz" );        
        response_ptjt->Write("Roo_response_ptjt" );
        response_ptr->Write( "Roo_response_ptr");
                
        response_ptzjt->Write("Roo_response_ptzjt");
        response_ptzr->Write( "Roo_response_ptzr");
        response_ptjtr->Write( "Roo_response_ptjtr");  
                
                
        TH1D *h1_purity_jetpt = (TH1D *)h1_num_purity_jetpt->Clone("h1_purity_jetpt");
        h1_purity_jetpt->Divide(h1_num_purity_jetpt, h1_denom_purity_jetpt, 1, 1, "B");
                
        TH1D *h1_efficiency_jetpt = (TH1D *)h1_num_efficiency_jetpt->Clone("h1_efficiency_jetpt");
        h1_efficiency_jetpt->Divide(h1_num_efficiency_jetpt, h1_denom_efficiency_jetpt, 1, 1, "B");
                                                
        h2_efficiency_ptz->Divide(h2_num_efficiency_ptz, h2_denom_efficiency_ptz, 1, 1, "B");    
        h2_efficiency_ptjt->Divide(h2_num_efficiency_ptjt, h2_denom_efficiency_ptjt, 1, 1, "B");    
        h2_efficiency_ptr->Divide(h2_num_efficiency_ptr, h2_denom_efficiency_ptr, 1, 1, "B");    
        
        h2_purity_ptz->Divide(h2_num_purity_ptz, h2_denom_purity_ptz, 1, 1, "B");     
        h2_purity_ptjt->Divide(h2_num_purity_ptjt, h2_denom_purity_ptjt, 1, 1, "B"); 
        h2_purity_ptr->Divide(h2_num_purity_ptr, h2_denom_purity_ptr, 1, 1, "B");
        
        h3_efficiency_ptzjt->Divide(h3_num_efficiency_ptzjt, h3_denom_efficiency_ptzjt, 1, 1, "B");
        h3_efficiency_ptzr->Divide(h3_num_efficiency_ptzr, h3_denom_efficiency_ptzr, 1, 1, "B");
        h3_efficiency_ptjtr->Divide(h3_num_efficiency_ptjtr, h3_denom_efficiency_ptjtr, 1, 1, "B");    
        
        h3_purity_ptzjt->Divide(h3_num_purity_ptzjt, h3_denom_purity_ptzjt, 1, 1, "B");    
        h3_purity_ptzr->Divide(h3_num_purity_ptzr, h3_denom_purity_ptzr, 1, 1, "B");
        h3_purity_ptjtr->Divide(h3_num_purity_ptjtr, h3_denom_purity_ptjtr, 1, 1, "B");    
                

        h1_efficiency_jetpt->Write("efficiency_jetpt");
        
        h1_purity_jetpt->Write("purity_jetpt");
        
        h2_efficiency_ptz->Write("efficiency_ptz");
        h2_efficiency_ptjt->Write("efficiency_ptjt");
        h2_efficiency_ptr->Write("efficiency_ptr");
        
        h2_purity_ptz->Write("purity_ptz");
        h2_purity_ptjt->Write("purity_ptjt");
        h2_num_purity_ptjt->Write("num_purity_ptjt");
        h2_denom_purity_ptjt->Write("denom_purity_ptjt");
        h2_purity_ptr->Write("purity_ptr");
        
        h3_efficiency_ptzjt->Write("efficiency_ptzjt");
        h3_efficiency_ptzr->Write("efficiency_ptzr");
        h3_efficiency_ptjtr->Write("efficiency_ptjtr");
        
        h3_purity_ptzjt->Write("purity_ptzjt");
        h3_purity_ptzr->Write("purity_ptzr");
        h3_purity_ptjtr->Write("purity_ptjtr"); 
        
        TH2D *h2_efficiency_zjt_ptbinned[ptbinsize-1];
        TH2D *h2_efficiency_zr_ptbinned[ptbinsize-1];
        TH2D *h2_efficiency_jtr_ptbinned[ptbinsize-1];  
        TH2D *h2_purity_zjt_ptbinned[ptbinsize-1];
        TH2D *h2_purity_zr_ptbinned[ptbinsize-1];
        TH2D *h2_purity_jtr_ptbinned[ptbinsize-1];      
        
        THStack *hs_efficiency_ptz  = new THStack("efficiency_z_all", ";z;efficiency");
        THStack *hs_efficiency_ptjt = new THStack("efficiency_jt_all", ";j_{T} [GeV/c];efficiency");
        THStack *hs_efficiency_ptr  = new THStack("efficiency_r_all", ";r;efficiency");
        THStack *hs_purity_ptz  = new THStack("purity_z_all", ";z;purity");
        THStack *hs_purity_ptjt = new THStack("purity_jt_all", ";j_{T} [GeV/c];purity");
        THStack *hs_purity_ptr  = new THStack("purity_r_all", ";r;purity");
                                
        for (int i = 1; i < ptbinsize; i++) {   
                TH1D *h1_ptz_temp = (TH1D *)h2_efficiency_ptz->ProjectionX(Form("efficiency_z_pt%d", i), i + 1, i + 1); 
                h1_ptz_temp->SetStats(0);
                h1_ptz_temp->SetMarkerStyle(i + 20);
                
                if (i!=5) {
                        h1_ptz_temp->SetMarkerColor(i);
                        h1_ptz_temp->SetLineColor(i);
                } else {
                        h1_ptz_temp->SetMarkerColor(i*i+3);
                        h1_ptz_temp->SetLineColor(i*i+3);
                }

                h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));
                hs_efficiency_ptz->Add(h1_ptz_temp);
                
                TH1D *h1_ptjt_temp = (TH1D *)h2_efficiency_ptjt->ProjectionX(Form("efficiency_jt_pt%d", i), i + 1, i + 1); 
                h1_ptjt_temp->SetStats(0);
                h1_ptjt_temp->SetMarkerStyle(i + 20);
                
                if (i!=5) {
                        h1_ptjt_temp->SetMarkerColor(i);
                        h1_ptjt_temp->SetLineColor(i);
                } else {
                        h1_ptjt_temp->SetMarkerColor(i*i+3);
                        h1_ptjt_temp->SetLineColor(i*i+3);
                }

                h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
                hs_efficiency_ptjt->Add(h1_ptjt_temp);
                
                TH1D *h1_ptr_temp = (TH1D *)h2_efficiency_ptr->ProjectionX(Form("efficiency_r_pt%d", i), i + 1, i + 1); 
                h1_ptr_temp->SetStats(0);
                h1_ptr_temp->SetMarkerStyle(i + 20);

                if (i!=5) {
                        h1_ptr_temp->SetMarkerColor(i);
                        h1_ptr_temp->SetLineColor(i);
                } else {
                        h1_ptr_temp->SetMarkerColor(i*i+3);
                        h1_ptr_temp->SetLineColor(i*i+3);
                }

                h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
                hs_efficiency_ptr->Add(h1_ptr_temp);
                
                
                h1_ptz_temp = (TH1D *)h2_purity_ptz->ProjectionX(Form("purity_z_pt%d", i), i + 1, i + 1); 
                h1_ptz_temp->SetStats(0);
                h1_ptz_temp->SetMarkerStyle(i + 20);
                
                if (i!=5) {
                        h1_ptz_temp->SetMarkerColor(i);
                        h1_ptz_temp->SetLineColor(i);
                } else {
                        h1_ptz_temp->SetMarkerColor(i*i+3);
                        h1_ptz_temp->SetLineColor(i*i+3);
                }

                h1_ptz_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
                hs_purity_ptz->Add(h1_ptz_temp);
                
                h1_ptjt_temp = (TH1D *)h2_purity_ptjt->ProjectionX(Form("purity_jt_pt%d", i), i + 1, i + 1); 
                h1_ptjt_temp->SetStats(0);
                h1_ptjt_temp->SetMarkerStyle(i + 20);

                if (i!=5) {
                        h1_ptjt_temp->SetMarkerColor(i);
                        h1_ptjt_temp->SetLineColor(i);
                } else {
                        h1_ptjt_temp->SetMarkerColor(i*i+3);
                        h1_ptjt_temp->SetLineColor(i*i+3);
                }

                h1_ptjt_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
                hs_purity_ptjt->Add(h1_ptjt_temp);
                
                h1_ptr_temp = (TH1D *)h2_purity_ptr->ProjectionX(Form("purity_r_pt%d", i), i + 1, i + 1); 
                h1_ptr_temp->SetStats(0);
                h1_ptr_temp->SetMarkerStyle(i + 20);

                if (i!=5) {
                        h1_ptr_temp->SetMarkerColor(i);
                        h1_ptr_temp->SetLineColor(i);
                } else {
                        h1_ptr_temp->SetMarkerColor(i*i+3);
                        h1_ptr_temp->SetLineColor(i*i+3);
                }

                h1_ptr_temp->SetTitle(Form("%.1f < p_{T, j} < %.1f GeV", pt_binedges[i], pt_binedges[i + 1]));        
                hs_purity_ptr->Add(h1_ptr_temp); 
                                        
                h3_efficiency_ptzjt->GetZaxis()->SetRange(i+1, i+1);      
                h2_efficiency_zjt_ptbinned[i-1] = (TH2D *)h3_efficiency_ptzjt->Project3D("yx");
                h2_efficiency_zjt_ptbinned[i-1]->SetStats(0);
                h2_efficiency_zjt_ptbinned[i-1]->SetName(Form("efficiency_zjt_pt%d", i));
                h2_efficiency_zjt_ptbinned[i-1]->Write(); 

                h3_efficiency_ptzr->GetZaxis()->SetRange(i+1, i+1); 
                h2_efficiency_zr_ptbinned[i-1] = (TH2D *)h3_efficiency_ptzr->Project3D("yx");
                h2_efficiency_zr_ptbinned[i-1]->SetStats(0);       
                h2_efficiency_zr_ptbinned[i-1]->SetName(Form("efficiency_zr_pt%d",i)); 
                h2_efficiency_zr_ptbinned[i-1]->Write();
                
                h3_efficiency_ptjtr->GetZaxis()->SetRange(i+1, i+1);        
                h2_efficiency_jtr_ptbinned[i-1] = (TH2D *)h3_efficiency_ptjtr->Project3D("yx");
                h2_efficiency_jtr_ptbinned[i-1]->SetStats(0);        
                h2_efficiency_jtr_ptbinned[i-1]->SetName(Form("efficiency_jtr_pt%d",i));        
                h2_efficiency_jtr_ptbinned[i-1]->Write();
                
                h3_purity_ptzjt->GetZaxis()->SetRange(i+1, i+1);      
                h2_purity_zjt_ptbinned[i-1] = (TH2D *)h3_purity_ptzjt->Project3D("yx");
                h2_purity_zjt_ptbinned[i-1]->SetStats(0);                       
                h2_purity_zjt_ptbinned[i-1]->SetName(Form("purity_zjt_pt%d",i)); 
                h2_purity_zjt_ptbinned[i-1]->Write();         

                h3_purity_ptzr->GetZaxis()->SetRange(i+1, i+1);        
                h2_purity_zr_ptbinned[i-1] = (TH2D *)h3_purity_ptzr->Project3D("yx");
                h2_purity_zr_ptbinned[i-1]->SetStats(0);
                h2_purity_zr_ptbinned[i-1]->SetName(Form("purity_zr_pt%d",i));
                h2_purity_zr_ptbinned[i-1]->Write();        
                
                h3_purity_ptjtr->GetZaxis()->SetRange(i+1, i+1);        
                h2_purity_jtr_ptbinned[i-1] = (TH2D *)h3_purity_ptjtr->Project3D("yx");
                h2_purity_jtr_ptbinned[i-1]->SetStats(0);
                h2_purity_jtr_ptbinned[i-1]->SetName(Form("purity_jtr_pt%d",i));
                h2_purity_jtr_ptbinned[i-1]->Write();                        
        }   
        
        //hs_efficiency_ptz->SetOption("nostack");
        //hs_efficiency_ptjt->SetOption("nostack");
        //hs_efficiency_ptr->SetOption("nostack");
        //hs_purity_ptz->SetOption("nostack");
        //hs_purity_ptjt->SetOption("nostack");
        //hs_purity_ptr->SetOption("nostack");
        hs_efficiency_ptz->Write();
        hs_efficiency_ptjt->Write();
        hs_efficiency_ptr->Write();
        hs_purity_ptz->Write();
        hs_purity_ptjt->Write();
        hs_purity_ptr->Write();       
                        
        NormalizeHist(h2_ptz_weighted);
        NormalizeHist(h2_ptjt_weighted);
        NormalizeHist(h2_ptr_weighted);   
        
        NormalizeHist(h1_z_weighted);
        NormalizeHist(h1_jt_weighted);
        NormalizeHist(h1_r_weighted);    
                
        NormalizeHist(h2_ptz);
        NormalizeHist(h2_ptjt);
        NormalizeHist(h2_ptr);   
        
        NormalizeHist(h1_z);
        NormalizeHist(h1_jt);
        NormalizeHist(h1_r); 
                        
        h2_ptz_weighted->Write();
        h2_ptjt_weighted->Write();
        h2_ptr_weighted->Write();
        h1_z_weighted->Write();
        h1_jt_weighted->Write();
        h1_r_weighted->Write();   
        
        h2_ptz->Write();
        h2_ptjt->Write();
        h2_ptr->Write();
        h1_z->Write();
        h1_jt->Write();
        h1_r->Write();       

        cout << event_counter << " events processed" << endl;
                
        /////////////////////////////////////////////////////
        TH1D *h1_jetpt_true = (TH1D *)h1_denom_efficiency_jetpt->Clone("h1_jetpt_true");
        TH1D *h1_jetpt_reco = (TH1D *)h1_denom_purity_jetpt->Clone("h1_jetpt_reco");
        
        h1_jetpt_reco->Multiply(h1_purity_jetpt);
        
        RooUnfoldBayes unfold_jetpt(response_jetpt, h1_jetpt_reco, 4);
        
        h1_jetpt_reco = (TH1D *)unfold_jetpt.Hreco();
        
        h1_jetpt_reco->Divide(h1_efficiency_jetpt);
        
        ////////////////////////////////////////////////////
                

//         int ican = -1, iframe = -1, itext = -1;
//         TCanvas *ccan[1000];
//         TH1F *frame[1000];
//         TLatex *text[1000];
//         for (int i = 0; i < 1000; i++)
//         {
//             text[i] = new TLatex();
//             text[i]->SetNDC(kTRUE);
//             text[i]->SetTextSize(0.06);
//         }
//         TLatex Tl;
//         Tl.SetNDC(kTRUE);
//         Tl.SetTextSize(0.04);
//         //
//         Int_t lhcbFont = 132; // Old LHCb style: 62;
//         // Line thickness
//         float_t lhcbWidth = 2.00; // Old LHCb style: 3.00;
//         // Text size
//         float_t lhcbTSize = 0.06;
        
//         // gROOT->SetStyle("Plain");
//         gStyle->SetCanvasDefH(1414);
//         gStyle->SetCanvasDefW(1000);
//         gStyle->SetPadBottomMargin(0.08);
//         gStyle->SetPadTopMargin(0.08);
//         gStyle->SetPadLeftMargin(0.12);
//         gStyle->SetPadRightMargin(0.12);
//         gStyle->SetTitleStyle(0);
//         gStyle->SetStatStyle(0);
//         // gStyle->SetLineWidth(3);
//         gStyle->SetOptStat(0);
//         //---- paint...
//         char buf[100];
//         char bufb[100];
        
//         TString rootfile;
//         TString plotfile;
//         TString plotfilePDF;
//         TString plotfileO;
//         TString plotfileC;
//         // TString OutputFileBase    = outbase+outinfo;
//         TString plotextension = TString("../../plots/BjetsMC/");
//         rootfile = extension_RootFilesMC + extension + TString(".root");
//         plotfile = plotextension + extension + TString(".ps");
        
//         plotfilePDF = "./plots/testing_simpleunfoldout.pdf";
//         plotfileO = plotfilePDF + TString("(");
//         plotfileC = plotfilePDF + TString("]");
        
//         // c->SaveAs("plots/"+extension+".pdf");
        
//         // c1->SaveAs(Form("plots_misc/Misc_%s"+extension+".pdf", file_name.c_str()));
//         // cout<<"...outbase   = "<<outbase.Data()<<endl;
//         // cout<<"...rootfile  = "<<rootfile.Data()<<endl;
//         // cout<<"...plotfile  = "<<plotfile.Data()<<endl;
        
//         // //
//         //
//           ++ican;
//           sprintf(buf, "ccan%d", ican);
//           ccan[ican] = new TCanvas(buf, buf, 1);
//           ccan[ican]->SetFillColor(10);
//           gPad->SetLeftMargin(0.16);
//           gPad->SetBottomMargin(0.06);
//           gPad->SetRightMargin(0.15);
//           ccan[ican]->cd();
//           ccan[ican]->Divide(2, 3);
//           ccan[ican]->cd(1);
    
//           h2_response_ptz->SetStats(0);
//           gPad->SetLogz();
//           h2_response_ptz->Draw("COLZ");
//           h2_response_ptz->SetYTitle("Truth z / p_{T}^{Jet}");
//           h2_response_ptz->SetXTitle("Reco z / p_{T}^{Jet}");
    
    
//           // Tl.drawLatex(0.7, 0.8, "No cuts");
//           ccan[ican]->cd(2);
//           h2_response_ptjt->SetStats(0);
//           gPad->SetLogz();
//           h2_response_ptjt->Draw("COLZ");
//           h2_response_ptjt->SetYTitle("Truth j_{T} / p_{T}^{Jet}");
//           h2_response_ptjt->SetXTitle("Reco j_{T} / p_{T}^{Jet}");
    
    
//           // Tl.drawLatex(0.7, 0.8, "No cuts");
//           ccan[ican]->cd(3);
//          h2_response_ptr->SetStats(0);
//          gPad->SetLogz();
//          h2_response_ptr->Draw("COLZ");
//          h2_response_ptr->SetYTitle("Truth r / p_{T}^{Jet}");
//          h2_response_ptr->SetXTitle("Reco r / p_{T}^{Jet}");
    
//         ccan[ican]->cd(4);
//         gPad->SetLogy();
//         TRandom3 *randomGenerator = new TRandom3();
//         for (int i = 0; i < ptbinsize; i++)
//         {
//             TH1D *h1_temp = (TH1D *)h2_response_jetpt->ProjectionX(Form("htemp_%djetpt", i), i + 1, i + 1);
//             h1_temp->Scale(1. / h1_temp->Integral());
//             h1_temp->SetMarkerStyle(i + 23);
//             h1_temp->SetMarkerColor(i + 4);
//             h1_temp->SetLineColor(i + 4);
//             h1_temp->Draw("P E SAME");
//             // int selectbin = SampleTH1(h1_temp, randomGenerator);
//             // h1_temp->SetMinimum(0.);
//             // h1_temp->SetMaximum(1.05);
//         }
        
//         ccan[ican]->cd(5);
//         gPad->SetLogz();
        
//         h2_response_jetpt->SetStats(0);
//         h2_response_jetpt->Draw("COLZ");
        
// //         ccan[ican]->cd(4);
// //         NormalizeHist(h1_jetpt);
// //         h1_jetpt->Draw("same");
// //        
// //         NormalizeHist(h1_jetpt_truth);
// //         h1_jetpt_truth->SetLineColor(kGreen);
// //         h1_jetpt_truth->Draw("same");
        
//         ccan[ican]->cd(6);
//         // h2_response_ptktdR->GetZaxis()->SetRange(2,5);
//         gPad->SetLogy();
//         // h2_response_ptktdR->Draw("COLZ");
//         h1_denom_efficiency_jetpt->SetXTitle("p_{T}^{jet}");
//         h1_denom_efficiency_jetpt->Draw("SAME");
//         h1_denom_purity_jetpt->Draw("SAME");
        
// //        ccan[ican]->cd(6);
// //        gPad->SetLogy();
// //        //  stack_theta_resolution->Draw("nostack");
// //        //  for (int i = 0; i < MaxEmissions; i++)
// //        //  {
// //        //    h1_theta_resolution_vec[i]->SetLineColor(1 + i);
// //        //    // h1_theta_resolution_vec[i]->Fit()
// //        //  }
// //        
//         ccan[ican]->cd();
//         ccan[ican]->Update();
//         if (ican == 0)
//         {
//             ccan[ican]->Print(plotfileO.Data());
//         }
//         else
//         {
//             ccan[ican]->Print(plotfilePDF.Data());
//         }
// //        //
// //        //
//         ++ican;
//         sprintf(buf, "ccan%d", ican);
//         ccan[ican] = new TCanvas(buf, buf, 1);
//         ccan[ican]->SetFillColor(10);
//         gPad->SetLeftMargin(0.16);
//         gPad->SetBottomMargin(0.06);
//         gPad->SetRightMargin(0.15);
        
//           ccan[ican]->cd();
//           ccan[ican]->Divide(2, 3);
//           ccan[ican]->cd(1);
//           gPad->SetLogz();
//           h2_truthreco_z->SetStats(0);
//           h2_truthreco_z->SetXTitle("MCTrue z");
//           h2_truthreco_z->SetYTitle("MCReco z");
//           // h2_lundplane->Sethflocimum(.01);
//           // h2_lundplane->SetMinimum(0.00001);
//           h2_truthreco_z->Draw("COLZ");
//           Tl.DrawLatex(.5, 0.75, Form("%.1f<p_{T,jet}<%.1f GeV", ptMin, ptMax));
//           // Tl.drawLatex(0.7, 0.8, "No cuts");
// //        //
//           ccan[ican]->cd(2);
//           gPad->SetLogz();
//           h2_truthreco_jt->SetStats(0);
//           h2_truthreco_jt->SetXTitle("MCTrue jt");
//           h2_truthreco_jt->SetYTitle("MCReco jt");
// //        //  // h2_lundplane->Sethflocimum(.01);
// //        //  // h2_lundplane->SetMinimum(0.00001);
//           h2_truthreco_jt->Draw("COLZ");
//           ccan[ican]->cd(3);
//           gPad->SetLogz();
//           h2_truthreco_r->SetStats(0);
//           h2_truthreco_r->SetXTitle("MCTrue r");
//           h2_truthreco_r->SetYTitle("MCReco r");
//           h2_truthreco_r->Draw("COLZ");


        
//           ccan[ican]->cd();
//           ccan[ican]->Update();
//           if (ican == 0)
//           {
//             ccan[ican]->Print(plotfileO.Data());
//           }
//           else
//           {
//             ccan[ican]->Print(plotfilePDF.Data());
//           }

//           ++ican;
 
    

//         ccan[ican]->cd();
//         ccan[ican]->Update();
//         if (ican == 0)
//         {
//             ccan[ican]->Print(plotfileO.Data());
//         }
//         else
//         {
//             ccan[ican]->Print(plotfilePDF.Data());
//         }
        
//         //
//         //
//         ++ican;
//         sprintf(buf, "ccan%d", ican);
//         ccan[ican] = new TCanvas(buf, buf, 1);
//         ccan[ican]->SetFillColor(10);
//         gPad->SetLeftMargin(0.16);
//         gPad->SetBottomMargin(0.06);
//         gPad->SetRightMargin(0.15);
//         ccan[ican]->cd();
//         ccan[ican]->Divide(2, 3);
//         ccan[ican]->cd(1);
        
//         h1_purity_jetpt->SetMarkerColor(kRed);
//         h1_purity_jetpt->SetMarkerStyle(21);
//         h1_purity_jetpt->SetMinimum(0.);
//         h1_purity_jetpt->SetMaximum(1.05);
//         h1_purity_jetpt->Draw("P E same");
//         h1_purity_jetpt->SetTitle("; p_{T}^{Jet};  Purity or Efficiency");
 
//         h1_efficiency_jetpt->SetMarkerColor(kBlue);
//         h1_efficiency_jetpt->SetMarkerStyle(22);
//         h1_efficiency_jetpt->SetMinimum(0.);
//         h1_efficiency_jetpt->SetMaximum(1.05);
//         h1_efficiency_jetpt->Draw("P E same");
    
//     auto legendeff_jetpt = new TLegend(0.5, 0.4, 0.7, 0.7);
//     legendeff_jetpt->SetTextSize(0.03);
//     legendeff_jetpt->SetBorderSize(0);
//     legendeff_jetpt->SetFillStyle(0);
//     legendeff_jetpt->SetFillColor(3);
//     legendeff_jetpt->AddEntry(h1_efficiency_jetpt, "HF Reco. Efficiency");
//     legendeff_jetpt->AddEntry(h1_purity_jetpt, "HF Reco. Purity");
//     legendeff_jetpt->Draw("SAME");


//         if (ican == 0)
//         {
//             ccan[ican]->Print(plotfileO.Data());
//         }
//         else
//         {
//             ccan[ican]->Print(plotfilePDF.Data());
//         }
        
//           ++ican;
//           sprintf(buf, "ccan%d", ican);
//           ccan[ican] = new TCanvas(buf, buf, 1);
//           ccan[ican]->SetFillColor(10);
//           gPad->SetLeftMargin(0.16);
//           gPad->SetBottomMargin(0.06);
//           gPad->SetRightMargin(0.15);
//           ccan[ican]->cd();
//           ccan[ican]->Divide(2, 3);
//           ccan[ican]->cd(1);
           
//         if (ican > -1)
//         {
//             cout << " You plotted " << ican + 1 << " canvasses......." << endl;
//             ccan[ican]->Print(plotfileC.Data());
//         }

        //f->Write();
        f->Close();     
        
}
