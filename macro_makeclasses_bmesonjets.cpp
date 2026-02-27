#include <iostream>
#include "./include/directories.h"

void macro_makeclasses_bmesonjets()
{
        TChain* sim_mc     = new TChain("MCJets/MCJetTree");  
        TChain* sim_mcreco = new TChain("Jets/DecayTree");
        TChain* data_2016  = new TChain("Jets/DecayTree");
        // TChain* data_2017 = new TChain("StdHltZJets/DecayTree");
        // TChain* data_2018 = new TChain("StdHltZJets/DecayTree");

        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root").c_str());
        sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root").c_str());
        
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root").c_str());

        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root").c_str());
        // sim_mc->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root").c_str());

        sim_mc->MakeClass("TBJetsMC");
        
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root").c_str());
        sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root").c_str());
        // sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root").c_str());
        // sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root").c_str());
        // sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root").c_str());
        // sim_mcreco->Add((input_folder + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root").c_str());
        sim_mcreco->MakeClass("TBJetsMCReco");

        data_2016->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MD_02092025.root").c_str());
        data_2016->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MU_02092025.root").c_str());
        // data_2017->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MD_02092025.root").c_str());
        // data_2017->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MU_02092025.root").c_str());
        // data_2018->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MD_02092025.root").c_str());
        // data_2018->Add((input_folder + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MU_02092025.root").c_str());
        data_2016->MakeClass("TBJetsData");
}