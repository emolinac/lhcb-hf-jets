#include "Settings.h"

#include "../Helpers_IC.h"
#include "../include/analysis-constants.h"
#include "../include/analysis-binning.h"
#include "../include/analysis-cuts.cpp"
#include "../include/analysis-cuts.h"
#include "../include/directories.h"
#include "../include/names.h"

#include "../include/utils.cpp"
#include "../include/utils.h"

using namespace std;

void MakeMCWeights(std::string variation = "nominal")
{
        if(gSystem->AccessPathName((output_folder + namef_simpleobservable_data[variation]).c_str())) {
                std::cout<<"Data file not found!"<<std::endl;
                std::cout<<"Check existence of input file or locations."<<std::endl;

                return;
        }

        if(gSystem->AccessPathName((output_folder + namef_simpleobservable_mcreco[variation]).c_str())) {
                std::cout<<"MCreco file not found!"<<std::endl;
                std::cout<<"Check existence of input file or locations."<<std::endl;

                return;
        }

        TFile f_data((output_folder + namef_simpleobservable_data[variation]).c_str());
        TFile f_mcreco((output_folder + namef_simpleobservable_mcreco[variation]).c_str());

        TFile f_output((output_folder + namef_data2mcreco_ratio).c_str(), "RECREATE");

        TH3D* h3_rl_jetpt_weight_data = (TH3D*) f_data.Get("h3_rl_jetpt_weight_uncorrected");
        TH3D* h3_rl_jetpt_weight_reco = (TH3D*) f_mcreco.Get("h3_rl_jetpt_weight_uncorrected");

        TH3D* h3_HFpt_eta_jetpt_data = (TH3D*) f_data.Get("h3_HFptetajetpt_uncorrected");
        TH3D* h3_HFpt_eta_jetpt_reco = (TH3D*) f_mcreco.Get("h3_HFptetajetpt_uncorrected");

        TH1D* h1_jetpt_data = (TH1D*) f_data.Get("Jet_pT_uncorrected");
        TH1D* h1_jetpt_reco = (TH1D*) f_mcreco.Get("Jet_pT_uncorrected");

        h3_rl_jetpt_weight_data->Scale(1./h3_rl_jetpt_weight_data->Integral());
        h3_rl_jetpt_weight_reco->Scale(1./h3_rl_jetpt_weight_reco->Integral());
        
        h3_HFpt_eta_jetpt_data->Scale(1./h3_HFpt_eta_jetpt_data->Integral());
        h3_HFpt_eta_jetpt_reco->Scale(1./h3_HFpt_eta_jetpt_reco->Integral());

        h1_jetpt_data->Scale(1./h1_jetpt_data->Integral());
        h1_jetpt_reco->Scale(1./h1_jetpt_reco->Integral());

        TH3D *h3_rl_jetpt_weight_ratio = (TH3D*) h3_rl_jetpt_weight_data->Clone(name_histo_data2mcreco_rl_jetpt_weight.c_str());
        TH3D *h3_HFpt_eta_jetpt_ratio  = (TH3D*) h3_HFpt_eta_jetpt_data->Clone(name_histo_data2mcreco_HFpt_eta_jetpt.c_str());
        TH3D *h1_jetpt_ratio           = (TH3D*) h1_jetpt_data->Clone(name_histo_data2mcreco_jetpt.c_str());
        
        h3_rl_jetpt_weight_ratio->Divide(h3_rl_jetpt_weight_reco);
        h3_rl_jetpt_weight_ratio->Write();

        h3_HFpt_eta_jetpt_ratio->Divide(h3_HFpt_eta_jetpt_reco);
        h3_HFpt_eta_jetpt_ratio->Write();

        h1_jetpt_ratio->Divide(h1_jetpt_reco);
        h1_jetpt_ratio->Write();

        f_data.Close();
        f_mcreco.Close();
}

