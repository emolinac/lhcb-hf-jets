

// #include "AnalyzeDijetData_withsubsets.h"
// #include "LundPlaneData.h"

#include "./BjetTree.C" 
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include "Settings.h"
#include "../Helpers.h"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/contrib/SoftDrop.hh"
//// #include "../RooHelpers.h"

using namespace std;

void MakeVarTree(int NumEvts_user = -1,
                 int dataset = 91599,
                 bool isData = true,
                 bool DoJESJER = false,
                 bool DoJetID = false,
                 bool L0MuonDiMuon = false)
{

    int NumEvts = NumEvts_user;
    int NumEvtsTruth = NumEvts_user;
    bool MCflag = !isData;
    
    int year = (dataset / 10000) % 10;
    int JetMeth = (dataset / 1000) % 10;
    int flavor = (dataset / 100) % 10;
    int ptRange = (dataset / 10) % 10;
    int Mag = (dataset / 1) % 10;
    int HF_pdgcode = -99;

    TString str_year = "2016";


    if (year == 6)
        str_year = "2016";
    else if (year == 7)
        str_year = "2017";
    else if (year == 8)
        str_year = "2018";

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


    if (isData)
        str_level = "data";
    else
    {
        str_level = "reco";
    }

    if (flavor == 1)
        str_followHard = "_hard";
    else
    {
        str_followHard = "_HF";
    }

    TString str_L0 = "";
    if (L0MuonDiMuon)
        str_L0 = "_L0MuonDiMuon";

    TString str_simversion = "";     
    if (onlysim9)
    {
        str_simversion = "_sim9";
    }    

    TString str_tree;
    TString extension_read, extension_RootFilesMC, extension_RootFilesData, extension_RootFiles;
    TString extension, extension_eff;
    
    extension = TString("tree_") + str_level + Form("_ev_%d", NumEvts) + Form("_eta_%.1f%.1f", etaMin, etaMax) + str_followHard + str_ghost + str_Mag + str_flavor + str_L0 + str_simversion + Form("_%d", dataset);
    
    //if (DoJES)
    //    extension = TString("JES_") + extension;
    //if (DoJER)
    //    extension = TString("JER_") + extension;
    if (DoJESJER)
        extension = TString("JESJER_") + extension;
    if (DoJetID)
        extension = TString("jetid_") + extension;
        
    ////////////////////////////////////////////////////
    ///              Tracking and PID
    //////////////////////////////////////////////////
    TString extension_pideffMu, extension_pideffK, extension_trackeff, extension_trackeff_Data, extension_trackeff_MC;
    TString extension_trackeff_Muon;
    TString eff_path;
    eff_path = TString( "../../Effs/");
   

    TH3D *h3_pideff_K_P_ETA_nTracks;
    TH3D *h3_pideff_Mu_P_ETA_nTracks;
    TH2D *h2_ratio_trkeff_P_ETA, *h2_ratio_trkeff_P_ETA_ERRHI, *h2_ratio_trkeff_P_ETA_ERRLO;
    TH2D *h2_ratio_trkeff_P_ETA_Muon, *h2_ratio_trkeff_P_ETA_Muon_ERRHI, *h2_ratio_trkeff_P_ETA_Muon_ERRLO;
    TH2D *h2_ratio_trkeff_P_ETA_DATA, *h2_ratio_trkeff_P_ETA_MC;

    extension_pideffK = "effhists-Turbo" + str_year + str_Mag + "-K-Brunel_DLLK>0-P.ETA.nTracks_Brunel";
    extension_pideffMu = "effhists-Turbo" + str_year + str_Mag + "-Mu-IsMuon&Brunel_DLLmu>0-P.ETA.nTracks_Brunel";
    extension_trackeff = "trackEff_" + str_year + "_Ratio_Full_Long_method";
    extension_trackeff_Muon = "trackEff_Muon_" + str_year + "_Ratio_Full_Long_method";
    extension_trackeff_Data = "trackEff_" + str_year + "_Data_Full_Long_method";
    extension_trackeff_MC = "trackEff_" + str_year + "_MC_Full_Long_method";
    
    TFile file_trackeff(eff_path+ "TrackEff/" + extension_trackeff + ".root", "READ");
    TFile file_trackeff_Muon(eff_path+ "TrackEff/" + extension_trackeff_Muon + ".root", "READ");
    TFile file_trackeff_Data(eff_path+ "TrackEff/" + extension_trackeff_Data + ".root", "READ");
    TFile file_trackeff_MC(eff_path+ "TrackEff/" + extension_trackeff_MC + ".root", "READ");
    TFile file_pideffK(eff_path + "PIDEff/" + extension_pideffK + ".root", "READ");
    TFile file_pideffMu(eff_path + "PIDEff/" + extension_pideffMu + ".root", "READ");

    cout << extension_pideffK << endl;
    cout << extension_pideffMu << endl;
    cout << extension_trackeff << endl;

    h2_ratio_trkeff_P_ETA_Muon = (TH2D *)file_trackeff_Muon.Get("hP_ETA");
    h2_ratio_trkeff_P_ETA_Muon_ERRHI = (TH2D *)file_trackeff_Muon.Get("hP_ETA_errhi");
    h2_ratio_trkeff_P_ETA_Muon_ERRLO = (TH2D *)file_trackeff_Muon.Get("hP_ETA_errlo");

    h2_ratio_trkeff_P_ETA = (TH2D *)file_trackeff.Get("hP_ETA");
    h2_ratio_trkeff_P_ETA_DATA = (TH2D *)file_trackeff_Data.Get("hP_ETA");
    h2_ratio_trkeff_P_ETA_MC = (TH2D *)file_trackeff_MC.Get("hP_ETA");
    h2_ratio_trkeff_P_ETA_ERRHI = (TH2D *)file_trackeff.Get("hP_ETA_errhi");
    h2_ratio_trkeff_P_ETA_ERRLO = (TH2D *)file_trackeff.Get("hP_ETA_errlo");
    h3_pideff_K_P_ETA_nTracks = (TH3D *)file_pideffK.Get("eff_Brunel_DLLK>0");
    h3_pideff_Mu_P_ETA_nTracks = (TH3D *)file_pideffMu.Get("eff_IsMuon&Brunel_DLLmu>0");

    //////////////////////////////////////////////////
    ///              Trigger
    ////////////////////////////////////////////////
    TString extension_trig_MC, extension_trig_Data, extension_RootFilesTrig;
    extension_RootFilesTrig = TString("../../Effs/TrigEff/");
    
    TH2D *h2_trigeff_Data;
    TH2D *h2_trigeff_MC;

    extension_trig_MC = "MC/jpsieff_reco_ev_-1_b_PID" + str_L0 +"_91599.root";
    extension_trig_Data = "data/jpsieff_data_ev_-1_b_PID" + str_L0 + "_91599.root";

    TFile file_trigeffMC(extension_RootFilesTrig + extension_trig_MC, "READ");
    TFile file_trigeffData(extension_RootFilesTrig + extension_trig_Data, "READ");

    h2_trigeff_Data = (TH2D *)file_trigeffData.Get("efficiency_Jpsiptrap");
    h2_trigeff_MC = (TH2D *)file_trigeffMC.Get("efficiency_Jpsiptrap");
    TH2D *h2_trigeff_ratio = (TH2D *)h2_trigeff_Data->Clone("h2_trigeff_ratio");
    h2_trigeff_ratio->Divide(h2_trigeff_MC);

    BjetTree Tree(0, dataset, isData, onlysim9);

    cout << "Total number of events = " << Tree.fChain->GetEntries() << endl;
    if (NumEvts == -1)
        NumEvts = Tree.fChain->GetEntries();
    // int NumEvts = 0;
    // cout <<"Choose number of events (-1: All Events, or enter integer): ";
    // cin>> NumEvts;

    cout << "Executing CAJetAlgo" << endl;

    extension_RootFilesMC = TString("../../root_files/BjetsMC/");
    extension_RootFilesData = TString("../../root_files/Bjets/");
    extension_RootFiles = isData ? "../../root_files/Bjets/" : "../../root_files/BjetsMC/"; 

    TFile *file_eff;
    TH2D *h2_trk_eff;
    TH2D *h2_neut_eff; 

    TFile f(extension_RootFiles + extension + ".root", "RECREATE");
    
    TH1F *h1_TIS = new TH1F("h1_TIS", "", ptJpsibinsize, ptJpsi_binedges);
    TH1F *h1_TISTOS = new TH1F("h1_TISTOS", "", ptJpsibinsize, ptJpsi_binedges);

    TH2F *h2_TIS_ptrap = new TH2F("h2_TIS_ptrap", "", ptJpsibinsize, ptJpsi_binedges, HFetabinsize, HFeta_binedges);
    TH2F *h2_TISTOS_ptrap = new TH2F("h2_TISTOS_ptrap", "", ptJpsibinsize, ptJpsi_binedges, HFetabinsize, HFeta_binedges);

    float jet_pt, jet_eta, tr_jet_pt, tr_jet_eta;
    float jet_rap, tr_jet_rap;
    float jet_px, jet_py, jet_pz, jet_e, jet_charge;
    float HF_px, HF_py, HF_pz, HF_e, HF_pt;
    float Jpsi_px, Jpsi_py, Jpsi_pz, Jpsi_e;
    float mup_px, mup_py, mup_pz, mup_e;
    float mum_px, mum_py, mum_pz, mum_e;
    float K_px, K_py, K_pz, K_e, K_p, K_eta;
    float mup_CHI2NDOF, mup_GHOSTPROB, mup_IPCHI2;
    float mum_CHI2NDOF, mum_GHOSTPROB, mum_IPCHI2;
    float K_CHI2NDOF, K_GHOSTPROB, K_IPCHI2;
    float Jpsi_CHI2NDOF, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_BPVDLS, Jpsi_IPCHI2;
    float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
    int nTracks, nSPDHits;
    float jpsi_ipchi2, k_ipchi2;
    float K_PIDK;

    float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
    float tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e, tr_HF_pt;
    float tr_mup_px, tr_mup_py, tr_mup_pz, tr_mup_e;
    float tr_mum_px, tr_mum_py, tr_mum_pz, tr_mum_e;
    float tr_K_px, tr_K_py, tr_K_pz, tr_K_e;
    
    double z, jt, r;
    double tr_z, tr_jt, tr_r;
    double zg, jtg, rg;

    int nSV;
    bool isTrueBjet, Hasbbbar;
    float bmass_dtf, chi2ndf_dtf, tau_dtf;
    int NumBHads_tr, eventNumber;
    
    float pideff_K(1.0), pideff_mup(1.0), pideff_mum(1.0);
    float pideff_K_err(1.0), pideff_mup_err(1.0), pideff_mum_err(1.0);
    float trkeff_ratio_K(1.0), trkeff_ratio_mup(1.0), trkeff_ratio_mum(1.0);
    float trkeff_ratio_K_errhi(1.0), trkeff_ratio_mup_errhi(1.0), trkeff_ratio_mum_errhi(1.0);
    float trkeff_ratio_K_errlo(1.0), trkeff_ratio_mup_errlo(1.0), trkeff_ratio_mum_errlo(1.0);
    float trigeff_Data(1.0), trigeff_MC(1.0), trigeff_ratio(1.0);
    
    vector<float> dtr_pt, dtr_rap, dtr_id, dtr_3charge;
    vector<float> sv_ipchi2, sv_keys, sv_pids;

    float sv_mass, sv_chi2, sv_cosine, sv_ntrks;
    int SVTag;

    bool mup_L0, mum_L0;
    bool jpsi_L0, jpsi_Hlt1, jpsi_Hlt2;
    bool Trig, TIS, TOS;
    // TLorentzVector
    TTree *BTree = new TTree("BTree", "B-jets Tree Variables");

    BTree->Branch("eventNumber", &eventNumber);

    BTree->Branch("dtr_pt", &dtr_pt);
    BTree->Branch("dtr_rap", &dtr_rap);
    BTree->Branch("dtr_id", &dtr_id);
    BTree->Branch("dtr_3charge", &dtr_3charge);
    
    BTree->Branch("jet_pt", &jet_pt);
    BTree->Branch("jet_eta", &jet_eta);
    BTree->Branch("jet_rap", &jet_rap);

    BTree->Branch("jet_px", &jet_px);
    BTree->Branch("jet_py", &jet_py);
    BTree->Branch("jet_pz", &jet_pz);
    BTree->Branch("jet_e", &jet_e);

    BTree->Branch("HF_px", &HF_px);
    BTree->Branch("HF_py", &HF_py);
    BTree->Branch("HF_pz", &HF_pz);
    BTree->Branch("HF_e", &HF_e);
    BTree->Branch("HF_pt", &HF_pt);

    BTree->Branch("Bu_IPCHI2", &Bu_IPCHI2);
    BTree->Branch("Bu_CHI2", &Bu_CHI2);
    BTree->Branch("Bu_CHI2NDOF", &Bu_CHI2NDOF);

    BTree->Branch("mum_px", &mum_px);
    BTree->Branch("mum_py", &mum_py);
    BTree->Branch("mum_pz", &mum_pz);
    BTree->Branch("mum_e", &mum_e);
    BTree->Branch("mum_IPCHI2", &mum_IPCHI2);
    BTree->Branch("mum_CHI2NDOF", &mum_CHI2NDOF);
    BTree->Branch("mum_GHOSTPROB", &mum_GHOSTPROB);

    BTree->Branch("mup_px", &mup_px);
    BTree->Branch("mup_py", &mup_py);
    BTree->Branch("mup_pz", &mup_pz);
    BTree->Branch("mup_e", &mup_e);
    BTree->Branch("mup_IPCHI2", &mup_IPCHI2);
    BTree->Branch("mup_CHI2NDOF", &mup_CHI2NDOF);
    BTree->Branch("mup_GHOSTPROB", &mup_GHOSTPROB);

    BTree->Branch("K_px", &K_px);
    BTree->Branch("K_py", &K_py);
    BTree->Branch("K_pz", &K_pz);
    BTree->Branch("K_e", &K_e);
    BTree->Branch("K_p", &K_p);
    BTree->Branch("K_eta", &K_eta);
    BTree->Branch("K_IPCHI2", &K_IPCHI2);
    BTree->Branch("K_CHI2NDOF", &K_CHI2NDOF);
    BTree->Branch("K_GHOSTPROB", &K_GHOSTPROB);
    BTree->Branch("k_ipchi2", &k_ipchi2);

    BTree->Branch("nTracks", &nTracks);
    BTree->Branch("nSPDHits", &nSPDHits);

    BTree->Branch("Jpsi_px", &Jpsi_px);
    BTree->Branch("Jpsi_py", &Jpsi_py);
    BTree->Branch("Jpsi_pz", &Jpsi_pz);
    BTree->Branch("Jpsi_e", &Jpsi_e);
    BTree->Branch("Jpsi_FDCHI2", &Jpsi_FDCHI2);
    BTree->Branch("Jpsi_BPVDLS", &Jpsi_BPVDLS);    
    BTree->Branch("Jpsi_CHI2", &Jpsi_CHI2);
    BTree->Branch("Jpsi_CHI2NDOF", &Jpsi_CHI2NDOF);
    BTree->Branch("jpsi_ipchi2", &jpsi_ipchi2);

    BTree->Branch("tr_jet_pt", &tr_jet_pt);
    BTree->Branch("tr_jet_eta", &tr_jet_eta);
    BTree->Branch("tr_jet_rap", &tr_jet_rap);

    BTree->Branch("tr_jet_px", &tr_jet_px);
    BTree->Branch("tr_jet_py", &tr_jet_py);
    BTree->Branch("tr_jet_pz", &tr_jet_pz);
    BTree->Branch("tr_jet_e", &tr_jet_e);

    BTree->Branch("tr_HF_px", &tr_HF_px);
    BTree->Branch("tr_HF_py", &tr_HF_py);
    BTree->Branch("tr_HF_pz", &tr_HF_pz);
    BTree->Branch("tr_HF_e", &tr_HF_e);
    BTree->Branch("tr_HF_pt", &tr_HF_pt);

    BTree->Branch("tr_mum_px", &tr_mum_px);
    BTree->Branch("tr_mum_py", &tr_mum_py);
    BTree->Branch("tr_mum_pz", &tr_mum_pz);
    BTree->Branch("tr_mum_e", &tr_mum_e);

    BTree->Branch("tr_mup_px", &tr_mup_px);
    BTree->Branch("tr_mup_py", &tr_mup_py);
    BTree->Branch("tr_mup_pz", &tr_mup_pz);
    BTree->Branch("tr_mup_e", &tr_mup_e);

    BTree->Branch("tr_K_px", &tr_K_px);
    BTree->Branch("tr_K_py", &tr_K_py);
    BTree->Branch("tr_K_pz", &tr_K_pz);
    BTree->Branch("tr_K_e", &tr_K_e);

    BTree->Branch("isTrueBjet", &isTrueBjet);

    BTree->Branch("nSV", &nSV);
    BTree->Branch("sv_mass", &sv_mass);
    BTree->Branch("sv_chi2", &sv_chi2);
    BTree->Branch("sv_ntrks", &sv_ntrks);
    BTree->Branch("sv_cosine", &sv_cosine);
    BTree->Branch("sv_ipchi2", &sv_ipchi2);
    BTree->Branch("sv_keys", &sv_keys);
    BTree->Branch("sv_pids", &sv_pids);

    BTree->Branch("SVTag", &SVTag);

    BTree->Branch("bmass_dtf", &bmass_dtf);
    BTree->Branch("chi2ndf_dtf", &chi2ndf_dtf);
    BTree->Branch("tau_dtf", &tau_dtf);
    BTree->Branch("NumBHads_tr", &NumBHads_tr);
    BTree->Branch("Hasbbbar", &Hasbbbar);
    BTree->Branch("K_PIDK", &K_PIDK);
    BTree->Branch("pideff_K", &pideff_K);
    BTree->Branch("pideff_mum", &pideff_mum);
    BTree->Branch("pideff_mup", &pideff_mup);
    BTree->Branch("pideff_K_err", &pideff_K_err);
    BTree->Branch("pideff_mum_err", &pideff_mum_err);
    BTree->Branch("pideff_mup_err", &pideff_mup_err);

    BTree->Branch("trkeff_ratio_K", &trkeff_ratio_K);
    BTree->Branch("trkeff_ratio_mup", &trkeff_ratio_mup);
    BTree->Branch("trkeff_ratio_mum", &trkeff_ratio_mum);

    BTree->Branch("trkeff_ratio_K_errhi", &trkeff_ratio_K_errhi);
    BTree->Branch("trkeff_ratio_mup_errhi", &trkeff_ratio_mup_errhi);
    BTree->Branch("trkeff_ratio_mum_errhi", &trkeff_ratio_mum_errhi);
    BTree->Branch("trkeff_ratio_K_errlo", &trkeff_ratio_K_errlo);
    BTree->Branch("trkeff_ratio_mup_errlo", &trkeff_ratio_mup_errlo);
    BTree->Branch("trkeff_ratio_mum_errlo", &trkeff_ratio_mum_errlo);

    BTree->Branch("trigeff_Data", &trigeff_Data);
    BTree->Branch("trigeff_MC", &trigeff_MC);
    BTree->Branch("trigeff_ratio", &trigeff_ratio);

    BTree->Branch("mup_L0", &mup_L0);
    BTree->Branch("mum_L0", &mum_L0);
    BTree->Branch("jpsi_L0", &jpsi_L0);
    BTree->Branch("jpsi_Hlt1", &jpsi_Hlt1);
    BTree->Branch("jpsi_Hlt2", &jpsi_Hlt2);
    BTree->Branch("Trig", &Trig);
    BTree->Branch("TIS", &TIS);
    BTree->Branch("TOS", &TOS);

    BTree->Branch("z", &z);
    BTree->Branch("jt", &jt);
    BTree->Branch("r", &r);
    BTree->Branch("zg", &zg);
    BTree->Branch("jtg", &jtg);
    BTree->Branch("rg", &rg);
    BTree->Branch("tr_z", &tr_z);
    BTree->Branch("tr_jt", &tr_jt);
    BTree->Branch("tr_r", &tr_r);
    
    //
    // Event loop
    //

    int eventNum;
    unsigned long long last_eventNum = 0;
    int events = 0;

    int ev_min = 0;
    TRandom3 *myRNG = new TRandom3();

    for (int ev = ev_min; ev < NumEvts + ev_min; ev++)
    {
        Tree.GetEntry(ev);
        if (ev % 10000 == 0)
            cout << "Executing event " << ev << endl;

        if (ev != 0)
        {
            if (Tree.eventNumber == last_eventNum)
            {
                continue;
            }
        }
        
        if (Tree.nPVs > 1)
            continue;
        
        
        bool WrongB = false;
        //std::cout << " Tree.totCandidates : " << Tree.totCandidates << std::endl;
        // if (Tree.totCandidates > 1)
        // {
        for (int dtrs0 = 0; dtrs0 < Tree.Jet_Dtr_nrecodtr; dtrs0++)
        {
            if (abs(Tree.Jet_Dtr_ID[dtrs0]) == HF_pdgcode)
            {
                if (fabs(Tree.Jet_Dtr_E[dtrs0] - Tree.Bu_PE) / (Tree.Bu_PE) > 0.001)
                {
                    WrongB = true;
                }    
                break;
            }
        }
        if (WrongB)
        {
            continue;
        }
        // }

        jpsi_L0 = L0MuonDiMuon ? Tree.Jpsi_L0MuonDecision_TOS || Tree.Jpsi_L0DiMuonDecision_TOS : Tree.Jpsi_L0DiMuonDecision_TOS;
        //mup_L0 = Tree.mup_L0MuonDecision_TOS || Tree.mup_L0DiMuonDecision_TOS;
        //mum_L0 = Tree.mum_L0MuonDecision_TOS || Tree.mum_L0DiMuonDecision_TOS;
        mup_L0 = Tree.mup_L0DiMuonDecision_TOS;
        mum_L0 = Tree.mum_L0DiMuonDecision_TOS;
        jpsi_Hlt1 = Tree.Jpsi_Hlt1DiMuonHighMassDecision_TOS;
        jpsi_Hlt2 = Tree.Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS || Tree.Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS || Tree.Jpsi_Hlt2DiMuonJPsiDecision_TOS;

        TIS = (Tree.Jpsi_L0Global_TIS && Tree.Jpsi_Hlt1Global_TIS && Tree.Jpsi_Hlt2Global_TIS);
        TOS = jpsi_L0 && jpsi_Hlt1 && jpsi_Hlt2;

        dtr_pt.clear();
        dtr_id.clear();
        dtr_rap.clear();
        dtr_3charge.clear();

        TLorentzVector HFjet, recojet, tr_truthjet, HFmeson, mup, mum, Jpsi, Kmeson;
        TLorentzVector tr_HFjet, tr_Kmeson, tr_mum, tr_mup, tr_HFmeson;

        mup.SetPxPyPzE(Tree.mup_PX / 1000., Tree.mup_PY / 1000., Tree.mup_PZ / 1000., Tree.mup_PE / 1000.);
        mum.SetPxPyPzE(Tree.mum_PX / 1000., Tree.mum_PY / 1000., Tree.mum_PZ / 1000., Tree.mum_PE / 1000.);
        Kmeson.SetPxPyPzE(Tree.K_PX / 1000., Tree.K_PY / 1000., Tree.K_PZ / 1000., Tree.K_PE / 1000.);
        HFmeson = mup + mum + Kmeson;

        HFjet.SetPxPyPzE(Tree.Jet_PX / 1000.,
                         Tree.Jet_PY / 1000.,
                         Tree.Jet_PZ / 1000.,
                         Tree.Jet_PE / 1000.);

        tr_HFjet.SetPxPyPzE(Tree.Jet_mcjet_PX / 1000.,
                            Tree.Jet_mcjet_PY / 1000.,
                            Tree.Jet_mcjet_PZ / 1000.,
                            Tree.Jet_mcjet_PE / 1000.);

        tr_mup.SetPxPyPzE(Tree.mup_TRUEP_X / 1000., Tree.mup_TRUEP_Y / 1000.,
                          Tree.mup_TRUEP_Z / 1000., Tree.mup_TRUEP_E / 1000.);
        tr_mum.SetPxPyPzE(Tree.mum_TRUEP_X / 1000., Tree.mum_TRUEP_Y / 1000.,
                          Tree.mum_TRUEP_Z / 1000., Tree.mum_TRUEP_E / 1000.);
        tr_Kmeson.SetPxPyPzE(Tree.K_TRUEP_X / 1000., Tree.K_TRUEP_Y / 1000.,
                             Tree.K_TRUEP_Z / 1000., Tree.K_TRUEP_E / 1000.);
        tr_HFmeson = tr_mup + tr_mum + tr_Kmeson;

        Jpsi.SetPxPyPzE(Tree.Jpsi_PX / 1000., Tree.Jpsi_PY / 1000., Tree.Jpsi_PZ / 1000., Tree.Jpsi_PE / 1000.);

        /*
        if (DoJES)
        {
            float rand = 0.976;
            HFjet -= HFmeson;
            HFjet.SetPx(HFjet.Px() * rand);
            HFjet.SetPy(HFjet.Py() * rand);
            HFjet.SetPz(HFjet.Pz() * rand);
            HFjet.SetE(HFjet.E() * rand);
            HFjet += HFmeson;
        }

        if (DoJER)
        {
            float rand = myRNG->Gaus(1, 0.114);
            HFjet -= HFmeson;
            HFjet.SetPx(HFjet.Px() * rand);
            HFjet.SetPy(HFjet.Py() * rand);
            HFjet.SetPz(HFjet.Pz() * rand);
            HFjet.SetE(HFjet.E() * rand);
            HFjet += HFmeson;
        }
        */

        if (DoJESJER)
        {
            const int n_iters = n_smearing_iter;
            for (int i_iter = 0; i_iter < n_iters; i_iter++)
            {

                // Reset HFjet using improved particle energy calculation
                //HFjet.SetPxPyPzE(mom.px(), mom.py(), mom.pz(), mom.E());

                double rand = 1;
                //std::cout << "Jet mass before: " << HFjet.M() << '\n';

                rand = get_JES_JER(HFjet.Pt(), myRNG);  // Standard inclusive Z+jet values
                //rand = get_JES_JER(HFjet.Pt(), myRNG, DoJESJER);  // Low-multiplicity Z+jet values

                // Testing for comparison to previous implementation
                //rand = 0.976 * myRNG->Gaus(1, 0.114);

                HFjet -= HFmeson;
                //double newE2 = HFjet.E()*HFjet.E() + (rand*rand - 1) * HFjet.Pt()*HFjet.Pt();
                double newE2 = HFjet.E()*HFjet.E() + (rand*rand - 1) * HFjet.P()*HFjet.P();
                double newE = (newE2 < 0) ? 0 : std::sqrt(newE2);

                //HFjet.SetPxPyPzE(HFjet.Px()*rand, HFjet.Py()*rand, HFjet.Pz(), newE);

                HFjet.SetPxPyPzE(HFjet.Px()*rand, HFjet.Py()*rand, HFjet.Pz()*rand, newE);
                // Testing for comparison to previous implementation
                //HFjet.SetPxPyPzE(HFjet.Px()*rand, HFjet.Py()*rand, HFjet.Pz()*rand, HFjet.E()*rand);
                HFjet += HFmeson;
                //std::cout << "Jet mass after: " << HFjet.M() << '\n';
            }
        }

        if (DoJetID)
        {
            double mpt = 0;
            int num_trk = 0;
            int num_neut = 0;
            int num_part = 0;
            double mtf = 0;
            for (int dtrs0 = 0; dtrs0 < Tree.Jet_Dtr_nrecodtr; dtrs0++)
            {

                if (fabs(Tree.Jet_Dtr_ThreeCharge[dtrs0]) == 0)
                {
                    num_neut++;
                }
                else
                {
                    num_trk++;
                }
            }

            if (num_trk < 2)
                continue;
        }


        bmass_dtf = Tree.Bu_ConsBu_M[0] / 1000.;
        chi2ndf_dtf = Tree.Bu_ConsBu_chi2[0] / Tree.Bu_ConsBu_nDOF[0];
        tau_dtf = Tree.Bu_ConsBu_ctau[0];


        nSV = Tree.Jet_SVTag_Nvertices;
        jpsi_ipchi2 = log10(Tree.Jpsi_IPCHI2_OWNPV);
        k_ipchi2 = (Tree.K_IPCHI2_OWNPV);

        float leading_pT = 0;
        int n_maxpT_cand = -999;
        int n_maxpT_entry = -999;
        int n_HFpt_entry = -999;
        int HF_counter = 0;

        float testjet_phi, testjet_px, testjet_py, testjet_pz, testjet_e;
        int ncand = Tree.nCandidate;


        bool hasHFhadron = false;

        int NumBHads = 0;

        TLorentzVector HFmatch;
        for (int dtrs0 = 0; dtrs0 < Tree.Jet_Dtr_nrecodtr; dtrs0++)
        {

            float trchi2ndf = Tree.Jet_Dtr_TrackChi2[dtrs0] / Tree.Jet_Dtr_TrackNDF[dtrs0];
            float dtr_charge = Tree.Jet_Dtr_ThreeCharge[dtrs0] / 3.;

            TLorentzVector dtr(Tree.Jet_Dtr_PX[dtrs0] / 1000.,
                               Tree.Jet_Dtr_PY[dtrs0] / 1000.,
                               Tree.Jet_Dtr_PZ[dtrs0] / 1000.,
                               Tree.Jet_Dtr_E[dtrs0] / 1000.);
            TLorentzVector tr_dtr(Tree.Jet_Dtr_TRUE_PX[dtrs0] / 1000.,
                                  Tree.Jet_Dtr_TRUE_PY[dtrs0] / 1000.,
                                  Tree.Jet_Dtr_TRUE_PZ[dtrs0] / 1000.,
                                  Tree.Jet_Dtr_TRUE_E[dtrs0] / 1000.);
            
            if (abs(Tree.Jet_Dtr_ID[dtrs0]) == HF_pdgcode)
            {
                HFmeson.SetPxPyPzE(dtr.Px(), dtr.Py(), dtr.Pz(), dtr.E());
                HFmatch.SetPxPyPzE(Tree.Jet_Dtr_TRUE_PX[dtrs0] / 1000.,
                                   Tree.Jet_Dtr_TRUE_PY[dtrs0] / 1000.,
                                   Tree.Jet_Dtr_TRUE_PZ[dtrs0] / 1000.,
                                   Tree.Jet_Dtr_TRUE_E[dtrs0] / 1000.);

                hasHFhadron = true;
                NumBHads++;
            }
            dtr_id.push_back(Tree.Jet_Dtr_ID[dtrs0]);
            dtr_pt.push_back(dtr.Pt());
            dtr_rap.push_back(dtr.Rapidity());
            dtr_3charge.push_back(Tree.Jet_Dtr_ThreeCharge[dtrs0]);
        }

        if (!hasHFhadron)
        {
            continue;
        }

        bool hasHFhadron_matched = false;
        NumBHads_tr = 0;
        bool hasb = false;
        bool hasbbar = false;
        if (!isData)
        {
            for (int dtrs0 = 0; dtrs0 < Tree.Jet_mcjet_nmcdtrs; dtrs0++)
            {
                float trchi2ndf = 0;
                TLorentzVector dtr(Tree.Jet_mcjet_dtrPX[dtrs0] / 1000.,
                                   Tree.Jet_mcjet_dtrPY[dtrs0] / 1000.,
                                   Tree.Jet_mcjet_dtrPZ[dtrs0] / 1000.,
                                   Tree.Jet_mcjet_dtrE[dtrs0] / 1000.);

                if (abs(Tree.Jet_mcjet_dtrID[dtrs0]) == HF_pdgcode)
                {
                    NumBHads_tr++;
                    if (fabs(dtr.Pt() - Tree.Bu_TRUEPT / 1000.) < 0.01)
                    {
                        tr_HFmeson.SetPxPyPzE(dtr.Px(), dtr.Py(), dtr.Pz(), dtr.E());
                        hasHFhadron_matched = true;
                    }
                }
            }
        }

        SVTag = 0;
        Hasbbbar = false;

        if (Tree.hasb && Tree.hasbbar)
        {
            Hasbbbar = true;
            if (Tree.Jet_SVTag_Tag)
            {
                SVTag = 1;
            }
        }
        if (Tree.Jet_SVTag_Tag && !(Tree.hasb && Tree.hasbbar)) // && !(Tree.tr_hasb && Tree.tr_hasbbar))
        {
            SVTag = 2;
        }

        TVector3 HF_meson = HFmeson.Vect();
        TVector3 HF_jet = HFjet.Vect();
        TVector3 tr_HF_meson = tr_HFmeson.Vect();
        TVector3 tr_HF_jet = tr_HFjet.Vect();

        z = (HF_meson.Dot(HF_jet) ) / (HF_jet.Mag2() );
        jt = (HF_jet.Cross(HF_meson).Mag()) / (HF_jet.Mag());
        r = static_cast < TLorentzVector > (HFmeson).DeltaR(HFjet, kTRUE);
        
        if(hasHFhadron_matched ){
            tr_z = (tr_HF_meson.Dot(tr_HF_jet) ) / (tr_HF_jet.Mag2() );
            tr_jt = (tr_HF_jet.Cross(tr_HF_meson).Mag()) / (tr_HF_jet.Mag());
            tr_r =  static_cast < TLorentzVector > (tr_HFmeson).DeltaR(tr_HFjet, kTRUE);

            if(Hasbbbar){
                zg = (tr_HF_meson.Dot(tr_HF_jet) ) / (tr_HF_jet.Mag2() );
                jtg = (tr_HF_jet.Cross(tr_HF_meson).Mag()) / (tr_HF_jet.Mag());
                rg = static_cast < TLorentzVector > (tr_HFmeson).DeltaR(tr_HFjet, kTRUE);
            }
            else {
                zg = -999.;
                jtg = -999.;
                rg = -999.;
            }
        }
        else {
            tr_z = -999.;
            tr_jt = -999.;
            tr_r = -999.;
        }

        
        jet_pt = HFjet.Pt();
        jet_eta = HFjet.Eta();
        jet_rap = HFjet.Rapidity();
        jet_px = HFjet.Px();
        jet_py = HFjet.Py();
        jet_pz = HFjet.Pz();
        jet_e = HFjet.E();
        HF_px = HFmeson.Px();
        HF_py = HFmeson.Py();
        HF_pz = HFmeson.Pz();
        HF_e = HFmeson.E();
        HF_pt = HFmeson.Pt();

        Bu_CHI2NDOF = Tree.Bu_ENDVERTEX_CHI2 / Tree.Bu_ENDVERTEX_NDOF;
        Bu_IPCHI2 = Tree.Bu_IPCHI2_OWNPV;
        // Bu_CHI2 = Tree.Bu_OWNPV_CHI2;
        Bu_CHI2 = Tree.Bu_ENDVERTEX_CHI2;

        mum_px = mum.Px();
        mum_py = mum.Py();
        mum_pz = mum.Pz();
        mum_e = mum.E();
        mum_CHI2NDOF = Tree.mum_TRACK_CHI2NDOF;
        mum_GHOSTPROB = Tree.mum_TRACK_GhostProb;
        mum_IPCHI2 = Tree.mum_IPCHI2_OWNPV;

        mup_px = mup.Px();
        mup_py = mup.Py();
        mup_pz = mup.Pz();
        mup_e = mup.E();
        mup_CHI2NDOF = Tree.mup_TRACK_CHI2NDOF;
        mup_GHOSTPROB = Tree.mup_TRACK_GhostProb;
        mup_IPCHI2 = Tree.mup_IPCHI2_OWNPV;

        K_px = Kmeson.Px();
        K_py = Kmeson.Py();
        K_pz = Kmeson.Pz();
        K_e = Kmeson.E();
        K_p = Kmeson.P() * 1000;
        K_eta = Kmeson.Eta();
        K_CHI2NDOF = Tree.K_TRACK_CHI2NDOF;
        K_GHOSTPROB = Tree.K_TRACK_GhostProb;
        K_IPCHI2 = Tree.K_IPCHI2_OWNPV;
        K_PIDK = Tree.K_PIDK;

        Jpsi_px = Jpsi.Px();
        Jpsi_py = Jpsi.Py();
        Jpsi_pz = Jpsi.Pz();
        Jpsi_e = Jpsi.E();
        Jpsi_CHI2NDOF = Tree.Jpsi_ENDVERTEX_CHI2 / Tree.Jpsi_ENDVERTEX_NDOF;
        Jpsi_IPCHI2 = Tree.Jpsi_IPCHI2_OWNPV;
        Jpsi_CHI2 = Tree.Jpsi_ENDVERTEX_CHI2;
        Jpsi_FDCHI2 = Tree.Jpsi_FDCHI2_OWNPV;
        Jpsi_BPVDLS = Tree.Jpsi_BPVDLS;
        
        nTracks = Tree.nTracks;
        nSPDHits = Tree.nSPDHits;

        tr_jet_eta = tr_HFjet.Eta();
        tr_jet_rap = tr_HFjet.Rapidity();
        tr_jet_pt = tr_HFjet.Pt();
        tr_jet_px = tr_HFjet.Px();
        tr_jet_py = tr_HFjet.Py();
        tr_jet_pz = tr_HFjet.Pz();
        tr_jet_e = tr_HFjet.E();

        tr_HF_pt = tr_HFmeson.Pt();
        tr_HF_px = tr_HFmeson.Px();
        tr_HF_py = tr_HFmeson.Py();
        tr_HF_pz = tr_HFmeson.Pz();
        tr_HF_e = tr_HFmeson.E();

        tr_mum_px = tr_mum.Px();
        tr_mum_py = tr_mum.Py();
        tr_mum_pz = tr_mum.Pz();
        tr_mum_e = tr_mum.E();
        tr_mup_px = tr_mup.Px();
        tr_mup_py = tr_mup.Py();
        tr_mup_pz = tr_mup.Pz();
        tr_mup_e = tr_mup.E();
        tr_K_px = tr_Kmeson.Px();
        tr_K_py = tr_Kmeson.Py();
        tr_K_pz = tr_Kmeson.Pz();
        tr_K_e = tr_Kmeson.E();

        sv_mass = Tree.Jet_SVTag_SigMaxMass / 1000.;
        sv_chi2 = log10(Tree.Jet_SVTag_SigMaxVtxChi2NDF);
        sv_ntrks = Tree.Jet_SVTag_SigMaxNtracks;
        sv_cosine = Tree.Jet_SVTag_SigMaxDirAngleS2S;

        if (isData && h3_pideff_K_P_ETA_nTracks != NULL && h2_ratio_trkeff_P_ETA != NULL)
        {
            float p_K = (Kmeson.P() >= 400) ? 399 * 1e3 : Kmeson.P() * 1e3;
            float p_mum = (mum.P() >= 400) ? 399 * 1e3 : mum.P() * 1e3;
            float p_mup = (mup.P() >= 400) ? 399 * 1e3 : mup.P() * 1e3;
            p_K = (p_K <= 4000.) ? 4000.1 : p_K;
            p_mum = (p_mum <= 4000.) ? 4000.1 : p_mum;
            p_mup = (p_mup <= 4000.) ? 4000.1 : p_mup;

            float eta_K = (Kmeson.Eta() >= 4.5) ? 4.4 : Kmeson.Eta();
            float eta_mum = (mum.Eta() >= 4.5) ? 4.4 : mum.Eta();
            float eta_mup = (mup.Eta() >= 4.5) ? 4.4 : mup.Eta();
            eta_K = (eta_K <= 2.) ? 2.1 : eta_K;
            eta_mum = (eta_mum <= 2.) ? 2.1 : eta_mum;
            eta_mup = (eta_mup <= 2.) ? 2.1 : eta_mup;

            int nTracks_new = (nTracks >= 500) ? 400 : nTracks;

            pideff_K = h3_pideff_K_P_ETA_nTracks->GetBinContent(h3_pideff_K_P_ETA_nTracks->GetXaxis()->FindBin(p_K), h3_pideff_K_P_ETA_nTracks->GetYaxis()->FindBin(eta_K), h3_pideff_K_P_ETA_nTracks->GetZaxis()->FindBin(nTracks_new));
            pideff_mup = h3_pideff_Mu_P_ETA_nTracks->GetBinContent(h3_pideff_Mu_P_ETA_nTracks->GetXaxis()->FindBin(p_mup), h3_pideff_Mu_P_ETA_nTracks->GetYaxis()->FindBin(eta_mup), h3_pideff_Mu_P_ETA_nTracks->GetZaxis()->FindBin(nTracks_new));
            pideff_mum = h3_pideff_Mu_P_ETA_nTracks->GetBinContent(h3_pideff_Mu_P_ETA_nTracks->GetXaxis()->FindBin(p_mum), h3_pideff_Mu_P_ETA_nTracks->GetYaxis()->FindBin(eta_mum), h3_pideff_Mu_P_ETA_nTracks->GetZaxis()->FindBin(nTracks_new));

            pideff_K_err = h3_pideff_K_P_ETA_nTracks->GetBinError(h3_pideff_K_P_ETA_nTracks->GetXaxis()->FindBin(p_K), h3_pideff_K_P_ETA_nTracks->GetYaxis()->FindBin(eta_K), h3_pideff_K_P_ETA_nTracks->GetZaxis()->FindBin(nTracks_new));
            pideff_mup_err = h3_pideff_Mu_P_ETA_nTracks->GetBinError(h3_pideff_Mu_P_ETA_nTracks->GetXaxis()->FindBin(p_mup), h3_pideff_Mu_P_ETA_nTracks->GetYaxis()->FindBin(eta_mup), h3_pideff_Mu_P_ETA_nTracks->GetZaxis()->FindBin(nTracks_new));
            pideff_mum_err = h3_pideff_Mu_P_ETA_nTracks->GetBinError(h3_pideff_Mu_P_ETA_nTracks->GetXaxis()->FindBin(p_mum), h3_pideff_Mu_P_ETA_nTracks->GetYaxis()->FindBin(eta_mum), h3_pideff_Mu_P_ETA_nTracks->GetZaxis()->FindBin(nTracks_new));

            // cout << pideff_mum << ", " << pideff_mup << ", " << pideff_mum_err << ", " << pideff_mup_err << endl;
            // if (pideff_mum == 0)
            //     cout << p_mum << ", " << eta_mum << ", " << nTracks_new << endl;
            trkeff_ratio_K = h2_ratio_trkeff_P_ETA->GetBinContent(h2_ratio_trkeff_P_ETA->GetXaxis()->FindBin(p_K), h2_ratio_trkeff_P_ETA->GetYaxis()->FindBin(eta_K));
            trkeff_ratio_mup = h2_ratio_trkeff_P_ETA_Muon->GetBinContent(h2_ratio_trkeff_P_ETA_Muon->GetXaxis()->FindBin(p_mup), h2_ratio_trkeff_P_ETA_Muon->GetYaxis()->FindBin(eta_mup));
            trkeff_ratio_mum = h2_ratio_trkeff_P_ETA_Muon->GetBinContent(h2_ratio_trkeff_P_ETA_Muon->GetXaxis()->FindBin(p_mum), h2_ratio_trkeff_P_ETA_Muon->GetYaxis()->FindBin(eta_mum));
            // double trkeff_K_DATA = h2_ratio_trkeff_P_ETA_DATA->GetBinContent(h2_ratio_trkeff_P_ETA_DATA->GetXaxis()->FindBin(p_K), h2_ratio_trkeff_P_ETA_DATA->GetYaxis()->FindBin(eta_K));
            // double trkeff_K_MC = h2_ratio_trkeff_P_ETA_MC->GetBinContent(h2_ratio_trkeff_P_ETA_MC->GetXaxis()->FindBin(p_K), h2_ratio_trkeff_P_ETA_MC->GetYaxis()->FindBin(eta_K));
            // cout << trkeff_ratio_K << ", " << trkeff_K_DATA / trkeff_K_MC << ", " << trkeff_K_MC / trkeff_K_DATA << endl;

            trkeff_ratio_K_errhi = h2_ratio_trkeff_P_ETA_ERRHI->GetBinContent(h2_ratio_trkeff_P_ETA_ERRHI->GetXaxis()->FindBin(p_K), h2_ratio_trkeff_P_ETA_ERRHI->GetYaxis()->FindBin(eta_K));
            trkeff_ratio_mup_errhi = h2_ratio_trkeff_P_ETA_Muon_ERRHI->GetBinContent(h2_ratio_trkeff_P_ETA_Muon_ERRHI->GetXaxis()->FindBin(p_mup), h2_ratio_trkeff_P_ETA_Muon_ERRHI->GetYaxis()->FindBin(eta_mup));
            trkeff_ratio_mum_errhi = h2_ratio_trkeff_P_ETA_Muon_ERRHI->GetBinContent(h2_ratio_trkeff_P_ETA_Muon_ERRHI->GetXaxis()->FindBin(p_mum), h2_ratio_trkeff_P_ETA_Muon_ERRHI->GetYaxis()->FindBin(eta_mum));

            trkeff_ratio_K_errlo = h2_ratio_trkeff_P_ETA_ERRLO->GetBinContent(h2_ratio_trkeff_P_ETA_ERRLO->GetXaxis()->FindBin(p_K), h2_ratio_trkeff_P_ETA_ERRLO->GetYaxis()->FindBin(eta_K));
            trkeff_ratio_mup_errlo = h2_ratio_trkeff_P_ETA_Muon_ERRLO->GetBinContent(h2_ratio_trkeff_P_ETA_Muon_ERRLO->GetXaxis()->FindBin(p_mup), h2_ratio_trkeff_P_ETA_Muon_ERRLO->GetYaxis()->FindBin(eta_mup));
            trkeff_ratio_mum_errlo = h2_ratio_trkeff_P_ETA_Muon_ERRLO->GetBinContent(h2_ratio_trkeff_P_ETA_Muon_ERRLO->GetXaxis()->FindBin(p_mum), h2_ratio_trkeff_P_ETA_Muon_ERRLO->GetYaxis()->FindBin(eta_mum));

            if (Kmeson.P() > 100 && Kmeson.P() < 200 && Kmeson.Eta() < 2.5)
            {
                trkeff_ratio_K_errhi = trkeff_ratio_K_errlo = 0.04;
            }
            if (Kmeson.P() < 10 && Kmeson.Eta() > 3.6)
            {
                trkeff_ratio_K_errhi = trkeff_ratio_K_errlo = 0.06;
            }

            if (pideff_K == 0 || trkeff_ratio_K == 0 || trkeff_ratio_mup == 0 || trkeff_ratio_mum == 0)
            {
                cout << pideff_K << ", " << trkeff_ratio_K << ", " << trkeff_ratio_mup << ", " << trkeff_ratio_mum << endl;
                cout << p_K << ", " << eta_K << endl;
                cout << p_mup << ", " << eta_mup << endl;
                cout << p_mum << ", " << eta_mum << endl;
                cout << nTracks << endl;
            }
        }
        else
        {
            if (isData)
            {
                std::cout << "input histograms not found! Tracking and pid efficiencies being set to 1.0." << std::endl;
            }
            pideff_K = pideff_mum = pideff_mup = trkeff_ratio_K = trkeff_ratio_mup = trkeff_ratio_mum = 1.0;
        }
        
        if (hasHFhadron_matched)
            isTrueBjet = true;
        else
            isTrueBjet = false;

        last_eventNum = Tree.eventNumber;
        eventNumber = Tree.eventNumber;

        // cout << endl;

        events++;
        BTree->Fill();
        // }
    }

    cout << "Number of saved events = " << events << endl;
    f.Write();
    f.Close();
}
//
