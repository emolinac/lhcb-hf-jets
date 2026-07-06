#include <iostream>
#include "./include/directories.h"

void macro_makeclasses_bmesonjets()
{
        TChain* sim_mc     = new TChain("MCJets/MCJetTree");  
        TChain* sim_mcreco = new TChain("Jets/DecayTree");
        TChain* data_2016_mu  = new TChain("Jets/DecayTree");
        TChain* data_2016_md  = new TChain("Jets/DecayTree");
        TChain* data_2017_mu  = new TChain("Jets/DecayTree");
        TChain* data_2017_md  = new TChain("Jets/DecayTree");
        TChain* data_2018_mu  = new TChain("Jets/DecayTree");
        TChain* data_2018_md  = new TChain("Jets/DecayTree");

        TChain* sim_mcreco_misid = new TChain("Jets/DecayTree");

        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root").c_str());
        
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root").c_str());

        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root").c_str());

        sim_mc->MakeClass("TBJetsMC");
        
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root").c_str());

        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root").c_str());
        
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root").c_str());

        sim_mcreco->MakeClass("TBJetsMCReco");

        data_2016_mu->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MU_02092025.root").c_str());
        data_2016_md->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MD_02092025.root").c_str());
        data_2017_mu->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MU_02092025.root").c_str());
        data_2017_md->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MD_02092025.root").c_str());
        data_2018_mu->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MU_02092025.root").c_str());
        data_2018_md->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MD_02092025.root").c_str());
        
        data_2016_mu->MakeClass("TBJetsData2016MU");
        data_2016_md->MakeClass("TBJetsData2016MD");
        data_2017_mu->MakeClass("TBJetsData2017MU");
        data_2017_md->MakeClass("TBJetsData2017MD");
        data_2018_mu->MakeClass("TBJetsData2018MU");
        data_2018_md->MakeClass("TBJetsData2018MD");

        sim_mcreco_misid->Add((input_folder + "BjetMisID_MC_nojetid_DTF_Sim09l_MD_05202024_full.root").c_str());
        sim_mcreco_misid->Add((input_folder + "BjetMisID_MC_nojetid_DTF_Sim09l_MU_05202024_full.root").c_str());
        sim_mcreco_misid->MakeClass("TBJetsMisID");
}