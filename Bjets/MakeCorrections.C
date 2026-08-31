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
#include "../include/names.h"

#include "../include/utils.cpp"
#include "../include/utils.h"

using namespace std;

void MakeCorrections(std::string variation = "nominal")
{
        if(gSystem->AccessPathName((output_folder + namef_mcreco_variations[variation]).c_str())) {
                std::cout<<"MCReco file not found. Check file or variation given as input."<<std::endl;

                return;
        }

        if(gSystem->AccessPathName((output_folder + "bjets_efficiencies.root").c_str())) {
                std::cout<<"Eff denominator file not found. Check file or variation given as input."<<std::endl;

                return;
        }

        if(gSystem->AccessPathName((output_folder_massfits + namef_massfits_results_mcreco[variation]).c_str())) {
                std::cout<<"Massfit file not found. Check file or variation given as input."<<std::endl;

                return;
        }

        bool SubtractGS    = false;
        bool DoJESJER      = false;
        bool DoJetID       = false;
        bool DoRecSelEff   = false;
        bool DoSignalSys   = false;
        bool DoUnfoldPrior = false;
        
        if (variation == "jesjer")
                DoJESJER = true;
        if (variation == "jetid")
                DoJetID = true;
        if (variation == "prior")
                DoUnfoldPrior = true;
        if (variation == "recseleff")
                DoRecSelEff = true;
        if (variation == "fitsignal")
                DoSignalSys = true;

        TFile* f_mcreco = new TFile((output_folder + namef_mcreco_variations[variation]).c_str());

        TTree* BTree = (TTree*) f_mcreco->Get("BTree");
        
        // Get denominators of the efficiencies
        TFile *file_eff = new TFile((output_folder + "bjets_efficiencies.root").c_str());
        
        TH1D *h1_denom_efficiency_HFpt   = (TH1D*) file_eff->Get("denom_efficiency_HFpt");
        TH1D *h1_denom_efficiency_jetpt  = (TH1D*) file_eff->Get("denom_efficiency_jetpt");
        
        TH2D *h2_denom_efficiency_HFptjetpt = (TH2D*) file_eff->Get("denom_efficiency_HFptjetpt");
        TH2D *h2_denom_efficiency_HFpteta   = (TH2D*) file_eff->Get("denom_efficiency_HFpteta");
        TH2D *h2_denom_efficiency_jetpteta  = (TH2D*) file_eff->Get("denom_efficiency_jetpteta");
        
        TH3D *h3_denom_efficiency_HFptetajetpt = (TH3D*) file_eff->Get("denom_efficiency_HFptetajetpt");
        
        //    /////////////////// Mass Fit Parameters /////////////////////////////////
        TFile f_massfit((output_folder_massfits + namef_massfits_results_mcreco[variation]).c_str(), "READ");
        
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
        TH3D*  h3_rl_jetpt_weight_data2mc_ratio;
        TH3D*  h3_HFpt_eta_jetpt_data2mc_ratio;
        TH1D*  h1_jetpt_data2mc_ratio;

        if (DoUnfoldPrior) {
                file_reco_weights = new TFile((output_folder + namef_data2mcreco_ratio).c_str(), "READ"); 

                h3_rl_jetpt_weight_data2mc_ratio = (TH3D*) file_reco_weights->Get(name_histo_data2mcreco_rl_jetpt_weight.c_str());
                h3_HFpt_eta_jetpt_data2mc_ratio  = (TH3D*) file_reco_weights->Get(name_histo_data2mcreco_HFpt_eta_jetpt.c_str());
                h1_jetpt_data2mc_ratio           = (TH1D*) file_reco_weights->Get(name_histo_data2mcreco_jetpt.c_str());
        }
        
        TFile *f = TFile::Open((output_folder + namef_corrections[variation]).c_str(), "RECREATE");
        gDirectory->cd();
        
        // 1D Denom Efficiencies and Purities of Observables (237 - 246)
        TH3D *h3_num_purity_HFptetajetpt     = new TH3D("num_purity_HFptetajetpt"    , "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, ptbinsize, pt_binedges);
        TH3D *h3_denom_purity_HFptetajetpt   = new TH3D("denom_purity_HFptetajetpt"  , "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, ptbinsize, pt_binedges);
        TH3D *h3_num_efficiency_HFptetajetpt = new TH3D("num_efficiency_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, ptbinsize, pt_binedges);
        
        TH2D *h2_num_purity_HFptjetpt     = new TH2D("num_purity_HFptjetpt"    , "", ptHFbinsize, ptHF_binedges, ptbinsize, pt_binedges);
        TH2D *h2_denom_purity_HFptjetpt   = new TH2D("denom_purity_HFptjetpt"  , "", ptHFbinsize, ptHF_binedges, ptbinsize, pt_binedges);
        TH2D *h2_num_efficiency_HFptjetpt = new TH2D("num_efficiency_HFptjetpt", "", ptHFbinsize, ptHF_binedges, ptbinsize, pt_binedges);
        
        TH1D *h1_num_purity_HFpt   = new TH1D("num_purity_HFpt", "", ptHFbinsize, ptHF_binedges);
        TH1D *h1_denom_purity_HFpt = new TH1D("denom_purity_HFpt", "", ptHFbinsize, ptHF_binedges);
        
        TH1D *h1_num_efficiency_HFpt = new TH1D("num_efficiency_HFpt", "", ptHFbinsize, ptHF_binedges);
        
        TH1D *h1_num_purity_jetpt   = new TH1D("num_purity_jetpt", "", ptbinsize, pt_binedges);
        TH1D *h1_denom_purity_jetpt = new TH1D("denom_purity_jetpt", "", ptbinsize, pt_binedges);
        
        TH1D *h1_num_efficiency_jetpt = new TH1D("num_efficiency_jetpt", "", ptbinsize, pt_binedges);
        
        TH1D *h1_jetpt          = new TH1D("jetpt"      , "", ptbinsize, pt_binedges);
        TH1D *h1_jetpt_truth    = new TH1D("jetpt_truth", "", ptbinsize, pt_binedges);
        TH2D *h2_response_jetpt = new TH2D("h2_response_jetpt", "", ptbinsize, pt_binedges, ptbinsize, pt_binedges);
        
        TH3D *h3_meas_HFptetajetpt = new TH3D("h3_meas_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, ptbinsize, pt_binedges);
        TH3D *h3_true_HFptetajetpt = new TH3D("h3_true_HFptetajetpt", "", ptHFbinsize, ptHF_binedges, HFetabinsize, HFeta_binedges, ptbinsize, pt_binedges);
        
        RooUnfoldResponse *response_HFptetajetpt = new RooUnfoldResponse(h3_meas_HFptetajetpt, h3_true_HFptetajetpt, "response_HFptetajetpt");
        
        // 1D RMs (for visualization purposes)
        TH1D *h1_form_rl     = new TH1D("h1_form_rl"    , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning);
        TH1D *h1_form_weight = new TH1D("h1_form_weight", "", nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_rl     = new RooUnfoldResponse(h1_form_rl    , h1_form_rl    , "response_rl");
        RooUnfoldResponse *response_weight = new RooUnfoldResponse(h1_form_weight, h1_form_weight, "response_weight");

        RooUnfoldResponse *response_rl_whf     = new RooUnfoldResponse(h1_form_rl    , h1_form_rl    , "response_rl_whf");
        RooUnfoldResponse *response_weight_whf = new RooUnfoldResponse(h1_form_weight, h1_form_weight, "response_weight_whf");

        RooUnfoldResponse *response_rl_wohf     = new RooUnfoldResponse(h1_form_rl    , h1_form_rl    , "response_rl_wohf");
        RooUnfoldResponse *response_weight_wohf = new RooUnfoldResponse(h1_form_weight, h1_form_weight, "response_weight_wohf");
        
        // 3D RM
        TH3D *h3_meas_rl_jetpt_weight = new TH3D("h3_meas_rl_jetpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight = new TH3D("h3_true_rl_jetpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair = new RooUnfoldResponse(h3_meas_rl_jetpt_weight, h3_true_rl_jetpt_weight, "response_npair");
        
        TH3D *h3_meas_rl_jetpt_weight_eqch = new TH3D("h3_meas_rl_jetpt_weight_eqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight_eqch = new TH3D("h3_true_rl_jetpt_weight_eqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair_eqch = new RooUnfoldResponse(h3_meas_rl_jetpt_weight_eqch, h3_true_rl_jetpt_weight_eqch, "response_npair_eqch");
        
        TH3D *h3_meas_rl_jetpt_weight_neqch = new TH3D("h3_meas_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight_neqch = new TH3D("h3_true_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair_neqch = new RooUnfoldResponse(h3_meas_rl_jetpt_weight_neqch, h3_true_rl_jetpt_weight_neqch, "response_npair_neqch");

        TH3D *h3_meas_rl_jetpt_weight_whf = new TH3D("h3_meas_rl_jetpt_weight_whf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight_whf = new TH3D("h3_true_rl_jetpt_weight_whf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair_whf = new RooUnfoldResponse(h3_meas_rl_jetpt_weight_whf, h3_true_rl_jetpt_weight_whf, "response_npair_whf");

        TH3D *h3_meas_rl_jetpt_weight_wohf = new TH3D("h3_meas_rl_jetpt_weight_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_true_rl_jetpt_weight_wohf = new TH3D("h3_true_rl_jetpt_weight_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        RooUnfoldResponse *response_npair_wohf = new RooUnfoldResponse(h3_meas_rl_jetpt_weight_wohf, h3_true_rl_jetpt_weight_wohf, "response_npair_wohf");
        
        // 3D RM alternative
        TH3D *h3_meas_rl_jetpt_HFpt = new TH3D("h3_meas_rl_jetpt_HFpt", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        TH3D *h3_true_rl_jetpt_HFpt = new TH3D("h3_true_rl_jetpt_HFpt", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        
        RooUnfoldResponse *response_npair_HFpt = new RooUnfoldResponse(h3_meas_rl_jetpt_HFpt, h3_true_rl_jetpt_HFpt, "response_npair_HFpt");
        
        // 3D Efficiencies
        TH3D *h3_denom_efficiency_rl_jetpt_weight       = new TH3D("h3_denom_efficiency_rl_jetpt_weight"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_denom_efficiency_rl_jetpt_weight_eqch  = new TH3D("h3_denom_efficiency_rl_jetpt_weight_eqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_denom_efficiency_rl_jetpt_weight_neqch = new TH3D("h3_denom_efficiency_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_denom_efficiency_rl_jetpt_weight_whf   = new TH3D("h3_denom_efficiency_rl_jetpt_weight_whf"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D *h3_denom_efficiency_rl_jetpt_weight_wohf  = new TH3D("h3_denom_efficiency_rl_jetpt_weight_wohf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);  

        TH3D* h3_num_efficiency_rl_jetpt_weight       = new TH3D("h3_num_efficiency_rl_jetpt_weight"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_num_efficiency_rl_jetpt_weight_eqch  = new TH3D("h3_num_efficiency_rl_jetpt_weight_eqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_num_efficiency_rl_jetpt_weight_neqch = new TH3D("h3_num_efficiency_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_num_efficiency_rl_jetpt_weight_whf   = new TH3D("h3_num_efficiency_rl_jetpt_weight_whf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_num_efficiency_rl_jetpt_weight_wohf  = new TH3D("h3_num_efficiency_rl_jetpt_weight_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_efficiency_rl_jetpt_weight       = new TH3D("h3_efficiency_rl_jetpt_weight"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight_eqch  = new TH3D("h3_efficiency_rl_jetpt_weight_eqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight_neqch = new TH3D("h3_efficiency_rl_jetpt_weight_neqch", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight_whf   = new TH3D("h3_efficiency_rl_jetpt_weight_whf"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_efficiency_rl_jetpt_weight_wohf  = new TH3D("h3_efficiency_rl_jetpt_weight_wohf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH2D* h2_num_efficiency_rl_jetpt   = new TH2D("h2_num_efficiency_rl_jetpt"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_denom_efficiency_rl_jetpt = new TH2D("h2_denom_efficiency_rl_jetpt", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_efficiency_rl_jetpt       = new TH2D("h2_efficiency_rl_jetpt"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        TH2D* h2_num_efficiency_rl_jetpt_whf   = new TH2D("h2_num_efficiency_rl_jetpt_whf"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_denom_efficiency_rl_jetpt_whf = new TH2D("h2_denom_efficiency_rl_jetpt_whf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_efficiency_rl_jetpt_whf       = new TH2D("h2_efficiency_rl_jetpt_whf"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);

        TH2D* h2_num_efficiency_rl_jetpt_wohf   = new TH2D("h2_num_efficiency_rl_jetpt_wohf"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_denom_efficiency_rl_jetpt_wohf = new TH2D("h2_denom_efficiency_rl_jetpt_wohf", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_efficiency_rl_jetpt_wohf       = new TH2D("h2_efficiency_rl_jetpt_wohf"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        
        TH3D *h3_denom_efficiency_rl_jetpt_HFpt = new TH3D("h3_denom_efficiency_rl_jetpt_HFpt", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        TH3D* h3_num_efficiency_rl_jetpt_HFpt   = new TH3D("h3_num_efficiency_rl_jetpt_HFpt"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        TH3D* h3_efficiency_rl_jetpt_HFpt       = new TH3D("h3_efficiency_rl_jetpt_HFpt"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);

        TH3F *h_denom_efficiency_rl_jetptHFpt_weight = new TH3F("h_denom_efficiency_rl_jetptHFpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jetpt_HFpt_1dim, jetpt_HFpt_1dim, nbin_weight, weight_binning);
        TH3F* h_num_efficiency_rl_jetptHFpt_weight   = new TH3F("h_num_efficiency_rl_jetptHFpt_weight"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jetpt_HFpt_1dim, jetpt_HFpt_1dim, nbin_weight, weight_binning);
        TH3F* h_efficiency_rl_jetptHFpt_weight       = new TH3F("h_efficiency_rl_jetptHFpt_weight"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jetpt_HFpt_1dim, jetpt_HFpt_1dim, nbin_weight, weight_binning);

        // 3D Purities
        TH3D* h3_num_purity_rl_jetpt_weight     = new TH3D("h3_num_purity_rl_jetpt_weight"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight   = new TH3D("h3_denom_purity_rl_jetpt_weight" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight         = new TH3D("h3_purity_rl_jetpt_weight"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_purity_rl_jetpt_weight_eqch   = new TH3D("h3_num_purity_rl_jetpt_weight_eqch"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight_eqch = new TH3D("h3_denom_purity_rl_jetpt_weight_eqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight_eqch       = new TH3D("h3_purity_rl_jetpt_weight_eqch"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_purity_rl_jetpt_weight_neqch   = new TH3D("h3_num_purity_rl_jetpt_weight_neqch"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight_neqch = new TH3D("h3_denom_purity_rl_jetpt_weight_neqch" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight_neqch       = new TH3D("h3_purity_rl_jetpt_weight_neqch"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_purity_rl_jetpt_weight_whf     = new TH3D("h3_num_purity_rl_jetpt_weight_whf"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight_whf   = new TH3D("h3_denom_purity_rl_jetpt_weight_whf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight_whf         = new TH3D("h3_purity_rl_jetpt_weight_whf"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH3D* h3_num_purity_rl_jetpt_weight_wohf     = new TH3D("h3_num_purity_rl_jetpt_weight_wohf"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_denom_purity_rl_jetpt_weight_wohf   = new TH3D("h3_denom_purity_rl_jetpt_weight_wohf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        TH3D* h3_purity_rl_jetpt_weight_wohf         = new TH3D("h3_purity_rl_jetpt_weight_wohf"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        TH2D* h2_num_purity_rl_jetpt   = new TH2D("h2_num_purity_rl_jetpt"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_denom_purity_rl_jetpt = new TH2D("h2_denom_purity_rl_jetpt" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_purity_rl_jetpt       = new TH2D("h2_purity_rl_jetpt"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);

        TH2D* h2_num_purity_rl_jetpt_whf   = new TH2D("h2_num_purity_rl_jetpt_whf"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_denom_purity_rl_jetpt_whf = new TH2D("h2_denom_purity_rl_jetpt_whf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_purity_rl_jetpt_whf       = new TH2D("h2_purity_rl_jetpt_whf"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);

        TH2D* h2_num_purity_rl_jetpt_wohf   = new TH2D("h2_num_purity_rl_jetpt_wohf"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_denom_purity_rl_jetpt_wohf = new TH2D("h2_denom_purity_rl_jetpt_wohf" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);
        TH2D* h2_purity_rl_jetpt_wohf       = new TH2D("h2_purity_rl_jetpt_wohf"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges);

        TH3D* h3_num_purity_rl_jetpt_HFpt     = new TH3D("h3_num_purity_rl_jetpt_HFpt"   , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        TH3D* h3_denom_purity_rl_jetpt_HFpt   = new TH3D("h3_denom_purity_rl_jetpt_HFpt" , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        TH3D* h3_purity_rl_jetpt_HFpt         = new TH3D("h3_purity_rl_jetpt_HFpt"       , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, ptHFbinsize, ptHF_binedges);
        
        TH3F *h_denom_purity_rl_jetptHFpt_weight = new TH3F("h_denom_purity_rl_jetptHFpt_weight", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jetpt_HFpt_1dim, jetpt_HFpt_1dim, nbin_weight, weight_binning);
        TH3F* h_num_purity_rl_jetptHFpt_weight   = new TH3F("h_num_purity_rl_jetptHFpt_weight"  , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jetpt_HFpt_1dim, jetpt_HFpt_1dim, nbin_weight, weight_binning);
        TH3F* h_purity_rl_jetptHFpt_weight       = new TH3F("h_purity_rl_jetptHFpt_weight"      , "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, nbin_jetpt_HFpt_1dim, jetpt_HFpt_1dim, nbin_weight, weight_binning);

        // Prior-weighted observables
        TH3D *h3_rl_jetpt_weight_weighted = new TH3D("h3_rl_jetpt_weight_weighted", "", nbin_rl_nominal_unfolding, unfolding_rl_nominal_binning, ptbinsize, pt_binedges, nbin_weight, weight_binning);
        
        // Event loop
        unsigned long long last_eventNum = 0;
        float last_jetpT = 0.;
        int event_counter = 0;
        
        cout << BTree->GetEntries() << endl;
        
        vector<float> *zs(0), *tr_zs(0);
        vector<float> *sv_ipchi2(0);
        
        float jet_pt, HF_pt, jet_rap, tr_jet_pt, tr_jet_rap, bmass_dtf;
        float tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e;
        float jet_px, jet_py, jet_pz, jet_e;
        float HF_px, HF_py, HF_pz, HF_e;
        float tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e, tr_HF_pt;
        bool isTrueBjet;
        bool truthmatched_jet_passed;
        int NumBHads_tr;
        int eventNumber, nSV;
        bool Hasbbbar, TOS;
        float Jpsi_CHI2NDOF, Jpsi_CHI2, Jpsi_FDCHI2, Jpsi_IPCHI2, Jpsi_BPVDLS;;
        float Bu_CHI2NDOF, Bu_CHI2, Bu_IPCHI2;
        float sv_mass, sv_chi2, sv_ntrks, sv_cosine;
        int SVTag;
        float K_PIDK, chi2ndf_dtf, tau_dtf;

        double WTA_reco_dist, WTA_true_dist;
        
        vector<float> *pair_rl = 0, *pair_weight = 0, *pair_chargeprod = 0, *pair_has_hf = 0;
        vector<float> *truthmatched_pair_rl = 0, *truthmatched_pair_weight = 0, *truthmatched_pair_has_hf = 0;
        vector<float> *true_pair_rl = 0, *true_pair_weight = 0, *true_pair_chargeprod = 0, *true_pair_has_hf = 0;
        
        BTree->SetBranchAddress("pair_rl"        , &pair_rl);
        BTree->SetBranchAddress("pair_weight"    , &pair_weight);
        BTree->SetBranchAddress("pair_chargeprod", &pair_chargeprod);
        BTree->SetBranchAddress("pair_has_hf"    , &pair_has_hf);
        
        BTree->SetBranchAddress("truthmatched_pair_rl"    , &truthmatched_pair_rl);
        BTree->SetBranchAddress("truthmatched_pair_weight", &truthmatched_pair_weight);
        BTree->SetBranchAddress("truthmatched_pair_has_hf", &truthmatched_pair_has_hf);

        BTree->SetBranchAddress("true_pair_rl"        , &true_pair_rl);
        BTree->SetBranchAddress("true_pair_weight"    , &true_pair_weight);
        BTree->SetBranchAddress("true_pair_chargeprod", &true_pair_chargeprod);
        BTree->SetBranchAddress("true_pair_has_hf"    , &true_pair_has_hf);

        BTree->SetBranchAddress("WTA_reco_dist", &WTA_reco_dist);
        BTree->SetBranchAddress("WTA_true_dist", &WTA_true_dist);
        
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
        
        BTree->SetBranchAddress("jet_pt" , &jet_pt);
        BTree->SetBranchAddress("jet_rap", &jet_rap);
        BTree->SetBranchAddress("jet_px" , &jet_px);
        BTree->SetBranchAddress("jet_py" , &jet_py);
        BTree->SetBranchAddress("jet_pz" , &jet_pz);
        BTree->SetBranchAddress("jet_e"  , &jet_e);

        BTree->SetBranchAddress("tr_HF_px", &tr_HF_px);
        BTree->SetBranchAddress("tr_HF_py", &tr_HF_py);
        BTree->SetBranchAddress("tr_HF_pz", &tr_HF_pz);
        BTree->SetBranchAddress("tr_HF_e", &tr_HF_e);
        BTree->SetBranchAddress("tr_HF_pt", &tr_HF_pt);
        
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
        
        BTree->SetBranchAddress("TOS", &TOS);

        BTree->SetBranchAddress("truthmatched_jet_passed", &truthmatched_jet_passed);
        
        int eventNum;
        int NumBjets = 0;
        int NumTrueBjets = 0;

        TLorentzVector HFmeson, HFjet, tr_HFjet, tr_HFmeson;
        
        for (int ev = 0; ev < BTree->GetEntries(); ev++) {
                BTree->GetEntry(ev);
                
                if (ev%10000 == 0) {
                        double percentage = 100.*ev/BTree->GetEntries();
                        std::cout<<"\r"<<percentage<<"\% jets processed."<< std::flush;
                }
                
                HFmeson.SetPxPyPzE(HF_px, HF_py, HF_pz, HF_e);
                HFjet.SetPxPyPzE(jet_px, jet_py, jet_pz, jet_e);
                tr_HFjet.SetPxPyPzE(tr_jet_px, tr_jet_py, tr_jet_pz, tr_jet_e);
                tr_HFmeson.SetPxPyPzE(tr_HF_px, tr_HF_py, tr_HF_pz, tr_HF_e);
                
                float prior_rescale_rl_jetpt_weight = 1.0;
                float prior_rescale_HFpt_eta_jetpt = 1.0;
                float prior_rescale_jetpt = 1.0;
                
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

                float jetptHFpt_mapped    = map_jetpt_HFpt(jet_pt, HF_pt);
                float tr_jetptHFpt_mapped = map_jetpt_HFpt(tr_jet_pt, tr_HF_pt);
        
                bool is_recojet_truthmatched = tr_pt_cond && tr_rap_cond && isTrueBjet && pt_cond && rap_cond && (WTA_reco_dist < WTA_dist_max && WTA_true_dist < WTA_dist_max);
                bool is_recojet_nominal      = pt_cond && rap_cond && (WTA_reco_dist < WTA_dist_max);
                
                if (is_recojet_truthmatched) {
                        h1_num_efficiency_jetpt->Fill(tr_jet_pt);
                        h1_num_efficiency_HFpt->Fill(tr_HF_pt);

                        h2_num_efficiency_HFptjetpt->Fill(tr_HF_pt, tr_jet_pt);

                        h3_num_efficiency_HFptetajetpt->Fill(tr_HF_pt, tr_HFmeson.Rapidity(), tr_jet_pt);

                        h1_num_purity_jetpt->Fill(jet_pt);
                        h1_num_purity_HFpt->Fill(HF_pt);

                        h2_num_purity_HFptjetpt->Fill(HF_pt, jet_pt);

                        h3_num_purity_HFptetajetpt->Fill(HF_pt, HFmeson.Rapidity(), jet_pt);

                        if (DoUnfoldPrior) {
                                prior_rescale_jetpt = h1_jetpt_data2mc_ratio->GetBinContent(h1_jetpt_data2mc_ratio->FindBin(jet_pt));
                                prior_rescale_HFpt_eta_jetpt = h3_HFpt_eta_jetpt_data2mc_ratio->GetBinContent(
                                        h3_HFpt_eta_jetpt_data2mc_ratio->FindBin(HF_pt, HFmeson.Rapidity(), jet_pt)
                                );
                        }

                        h2_response_jetpt->Fill(jet_pt, tr_jet_pt, prior_rescale_jetpt);
                        
                        response_HFptetajetpt->Fill(HF_pt, HFmeson.Rapidity(), jet_pt, tr_HF_pt, tr_HFmeson.Rapidity(), tr_jet_pt, prior_rescale_HFpt_eta_jetpt);

                        // Numerator of the in-jet pair efficiency
                        if (!truthmatched_pair_rl->empty()) {
                                ULong_t vector_size = truthmatched_pair_rl->size();

                                // Note: this vectors will have the same size always, dont worry
                                float *rl_info         = truthmatched_pair_rl->data();
                                float *weight_info     = truthmatched_pair_weight->data();
                                float *chargeprod_info = pair_chargeprod->data();
                                float *has_hf_info     = truthmatched_pair_has_hf->data();
                                
                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        if (rl_info[vector_index] == -999)
                                                continue;

                                        h3_num_efficiency_rl_jetpt_weight->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        h3_num_efficiency_rl_jetpt_HFpt->Fill(rl_info[vector_index], tr_jet_pt, tr_HF_pt);

                                        h_num_efficiency_rl_jetptHFpt_weight->Fill(rl_info[vector_index], tr_jetptHFpt_mapped, weight_info[vector_index]);

                                        h2_num_efficiency_rl_jetpt->Fill(rl_info[vector_index], tr_jet_pt);

                                        if (chargeprod_info[vector_index] > 0)
                                                h3_num_efficiency_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                        else if (chargeprod_info[vector_index] < 0)
                                                h3_num_efficiency_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        if (has_hf_info[vector_index] == 1) {
                                                h3_num_efficiency_rl_jetpt_weight_whf->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                                h2_num_efficiency_rl_jetpt_whf->Fill(rl_info[vector_index], tr_jet_pt);
                                        } else if (has_hf_info[vector_index] == 0 ) {
                                                h3_num_efficiency_rl_jetpt_weight_wohf->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                                h2_num_efficiency_rl_jetpt_wohf->Fill(rl_info[vector_index], tr_jet_pt);
                                        }
                                }
                        }

                        // Denominator of the in-jet pair effienciency
                        if (!true_pair_rl->empty()) {
                                ULong_t vector_size = true_pair_rl->size();

                                // Note: this vectors will have the same size always, dont worry
                                float *rl_info         = true_pair_rl->data();
                                float *weight_info     = true_pair_weight->data();
                                float *chargeprod_info = true_pair_chargeprod->data();
                                float *has_hf_info     = true_pair_has_hf->data();
                                
                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        h3_denom_efficiency_rl_jetpt_weight->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        h3_denom_efficiency_rl_jetpt_HFpt->Fill(rl_info[vector_index], tr_jet_pt, tr_HF_pt);

                                        h_denom_efficiency_rl_jetptHFpt_weight->Fill(rl_info[vector_index], tr_jetptHFpt_mapped, weight_info[vector_index]);

                                        h2_denom_efficiency_rl_jetpt->Fill(rl_info[vector_index], tr_jet_pt);

                                        if (chargeprod_info[vector_index] > 0)
                                                h3_denom_efficiency_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                        else if (chargeprod_info[vector_index] < 0)
                                                h3_denom_efficiency_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);

                                        if (has_hf_info[vector_index] == 1) {
                                                h3_denom_efficiency_rl_jetpt_weight_whf->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                                h2_denom_efficiency_rl_jetpt_whf->Fill(rl_info[vector_index], tr_jet_pt);
                                        } else if (has_hf_info[vector_index] == 0) {
                                                h3_denom_efficiency_rl_jetpt_weight_wohf->Fill(rl_info[vector_index], tr_jet_pt, weight_info[vector_index]);
                                                h2_denom_efficiency_rl_jetpt_wohf->Fill(rl_info[vector_index], tr_jet_pt);
                                        }
                                }
                        }
                
                        // In-jet pair purity
                        if (!pair_rl->empty()) {
                                ULong_t vector_size = pair_rl->size();

                                // Note: this vectors will have the same size always, dont worry
                                float *rl_info         = pair_rl->data();
                                float *weight_info     = pair_weight->data();
                                float *chargeprod_info = pair_chargeprod->data();
                                float *has_hf_info     = pair_has_hf->data();

                                float *truthmatched_rl_info     = truthmatched_pair_rl->data();
                                float *truthmatched_weight_info = truthmatched_pair_weight->data();
                                float *truthmatched_has_hf_info = truthmatched_pair_has_hf->data();

                                if (pair_rl->size() != truthmatched_pair_rl->size())
                                        std::cout<<"AAAAAAAAAAAAAAAAAAAAAA"<<std::endl;

                                for(int vector_index = 0 ; vector_index < vector_size ; vector_index++) {
                                        if (truthmatched_rl_info[vector_index] != -999){

                                                if (DoUnfoldPrior) {
                                                        prior_rescale_rl_jetpt_weight = h3_rl_jetpt_weight_data2mc_ratio->GetBinContent(
                                                                h3_rl_jetpt_weight_data2mc_ratio->FindBin(rl_info[vector_index], jet_pt, weight_info[vector_index])
                                                        );
                                                }
                                        
                                                h3_num_purity_rl_jetpt_weight->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);

                                                h3_num_purity_rl_jetpt_HFpt->Fill(rl_info[vector_index], jet_pt, HF_pt);

                                                h_num_purity_rl_jetptHFpt_weight->Fill(rl_info[vector_index], jetptHFpt_mapped, weight_info[vector_index]);

                                                h2_num_purity_rl_jetpt->Fill(rl_info[vector_index], jet_pt);

                                                response_npair->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], 
                                                        truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index], prior_rescale_rl_jetpt_weight);

                                                response_npair_HFpt->Fill(rl_info[vector_index], jet_pt, HF_pt, truthmatched_rl_info[vector_index], tr_jet_pt, tr_HF_pt);
                                      
                                                response_rl->Fill(rl_info[vector_index], truthmatched_rl_info[vector_index]);

                                                response_weight->Fill(weight_info[vector_index], truthmatched_weight_info[vector_index]);
                                                
                                                if (chargeprod_info[vector_index] > 0) {
                                                        h3_num_purity_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                                        
                                                        response_npair_eqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index]);

                                                } else if (chargeprod_info[vector_index] < 0) {
                                                        h3_num_purity_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);

                                                        response_npair_neqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index]);
                                                }

                                                if (has_hf_info[vector_index] == 1) {
                                                        h3_num_purity_rl_jetpt_weight_whf->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                                        h2_num_purity_rl_jetpt_whf->Fill(rl_info[vector_index], jet_pt);

                                                        response_npair_whf->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index],prior_rescale_rl_jetpt_weight);

                                                        response_rl_whf->Fill(rl_info[vector_index], truthmatched_rl_info[vector_index]);
                                                        response_weight_whf->Fill(weight_info[vector_index], truthmatched_weight_info[vector_index]);

                                                } else if (has_hf_info[vector_index] == 0) {
                                                        h3_num_purity_rl_jetpt_weight_wohf->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                                        h2_num_purity_rl_jetpt_wohf->Fill(rl_info[vector_index], jet_pt);

                                                        response_npair_wohf->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index], truthmatched_rl_info[vector_index], tr_jet_pt, truthmatched_weight_info[vector_index], prior_rescale_rl_jetpt_weight);
                                                        response_rl_wohf->Fill(rl_info[vector_index], truthmatched_rl_info[vector_index]);
                                                        response_weight_wohf->Fill(weight_info[vector_index], truthmatched_weight_info[vector_index]);
                                                }
                                        }

                                        h3_denom_purity_rl_jetpt_weight->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);

                                        h3_denom_purity_rl_jetpt_HFpt->Fill(rl_info[vector_index], jet_pt, HF_pt);

                                        h_denom_purity_rl_jetptHFpt_weight->Fill(rl_info[vector_index], jetptHFpt_mapped, weight_info[vector_index]);

                                        h2_denom_purity_rl_jetpt->Fill(rl_info[vector_index], jet_pt);

                                        if (chargeprod_info[vector_index] > 0)
                                                h3_denom_purity_rl_jetpt_weight_eqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                        else if (chargeprod_info[vector_index] < 0)
                                                h3_denom_purity_rl_jetpt_weight_neqch->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);

                                        if (has_hf_info[vector_index] == 1) {
                                                h3_denom_purity_rl_jetpt_weight_whf->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                                h2_denom_purity_rl_jetpt_whf->Fill(rl_info[vector_index], jet_pt);

                                        } else if (has_hf_info[vector_index] == 0) {
                                                h3_denom_purity_rl_jetpt_weight_wohf->Fill(rl_info[vector_index], jet_pt, weight_info[vector_index]);
                                                h2_denom_purity_rl_jetpt_wohf->Fill(rl_info[vector_index], jet_pt);
                                        }
                                }
                        }

                        NumTrueBjets++;              
                }
                
                if (is_recojet_nominal) {
                        h1_denom_purity_jetpt->Fill(jet_pt);
                        h1_denom_purity_HFpt->Fill(HF_pt);

                        h2_denom_purity_HFptjetpt->Fill(HF_pt, jet_pt);

                        h3_denom_purity_HFptetajetpt->Fill(HF_pt, HFmeson.Rapidity(), jet_pt);
                }    
                
                event_counter++;
        }

        RooUnfoldResponse *response_jetpt = new RooUnfoldResponse(h1_jetpt, h1_jetpt_truth, h2_response_jetpt, "response_jetpt");
        // TH2 *h2_response_jetpt  = response_jetpt->Hresponse();
        TH2 *h2_response_rl          = response_rl->Hresponse();
        TH2 *h2_response_weight      = response_weight->Hresponse();
        TH2 *h2_response_rl_whf      = response_rl_whf->Hresponse();
        TH2 *h2_response_weight_whf  = response_weight_whf->Hresponse();
        TH2 *h2_response_rl_wohf     = response_rl_wohf->Hresponse();
        TH2 *h2_response_weight_wohf = response_weight_wohf->Hresponse();
        
        TH2 *h3_response_npair       = response_npair->Hresponse();        
        TH2 *h3_response_npair_eqch  = response_npair_eqch->Hresponse();
        TH2 *h3_response_npair_neqch = response_npair_neqch->Hresponse();
        TH2 *h3_response_npair_whf   = response_npair_whf->Hresponse();
        TH2 *h3_response_npair_wohf  = response_npair_wohf->Hresponse();
        TH2 *h3_response_npair_HFpt  = response_npair_HFpt->Hresponse();
                
        h2_response_jetpt->GetXaxis()->SetTitle("reco p_{T, jet} [GeV/c]");
        h2_response_jetpt->GetYaxis()->SetTitle("truth p_{T, jet} [GeV/c]");
        h2_response_rl->GetXaxis()->SetTitle("reco R_{L}");
        h2_response_rl->GetYaxis()->SetTitle("truth R_{L}");
        h2_response_weight->GetXaxis()->SetTitle("reco w");
        h2_response_weight->GetYaxis()->SetTitle("truth w");
        h2_response_rl_whf->GetXaxis()->SetTitle("reco R_{L}");
        h2_response_rl_whf->GetYaxis()->SetTitle("truth R_{L}");
        h2_response_weight_whf->GetXaxis()->SetTitle("reco w");
        h2_response_weight_whf->GetYaxis()->SetTitle("truth w");
        h2_response_rl_wohf->GetXaxis()->SetTitle("reco R_{L}");
        h2_response_rl_wohf->GetYaxis()->SetTitle("truth R_{L}");
        h2_response_weight_wohf->GetXaxis()->SetTitle("reco w");
        h2_response_weight_wohf->GetYaxis()->SetTitle("truth w");
                
        f->cd();

        h2_response_jetpt->Write("h2_response_jetpt");
        h2_response_rl->Write("response_rl");
        h2_response_weight->Write("response_weight");
        h2_response_rl_whf->Write("response_rl_whf");
        h2_response_weight_whf->Write("response_weight_whf");
        h2_response_rl_wohf->Write("response_rl_wohf");
        h2_response_weight_wohf->Write("response_weight_wohf");
        
        h3_response_npair->Write("response_npair");
        h3_response_npair_eqch->Write("response_npair_eqch");
        h3_response_npair_neqch->Write("response_npair_neqch");    
        h3_response_npair_whf->Write("response_npair_whf");
        h3_response_npair_wohf->Write("response_npair_wohf");

        h3_response_npair_HFpt->Write("response_npair_HFpt");
        
        response_HFptetajetpt->Write("Roo_response_HFptetajetpt");
        response_jetpt->Write("Roo_response_jetpt");
        response_rl->Write("Roo_response_rl");
        response_weight->Write("Roo_response_weight");
        response_rl_whf->Write("Roo_response_rl_whf");
        response_weight_whf->Write("Roo_response_weight_whf");
        response_rl_wohf->Write("Roo_response_rl_wohf");
        response_weight_wohf->Write("Roo_response_weight_wohf");
        
        response_npair->Write("Roo_response_npair" );  
        response_npair_eqch->Write("Roo_response_npair_eqch" );
        response_npair_neqch->Write( "Roo_response_npair_neqch");
        response_npair_whf->Write("Roo_response_npair_whf" );
        response_npair_wohf->Write( "Roo_response_npair_wohf");
        
                
        response_npair_HFpt->Write("Roo_response_npair_HFpt");  
        
        TH3D *h3_purity_HFptetajetpt = (TH3D *)h3_num_purity_HFptetajetpt->Clone("h3_purity_HFptetajetpt");
        h3_purity_HFptetajetpt->Divide(h3_num_purity_HFptetajetpt, h3_denom_purity_HFptetajetpt, 1, 1, "B");
                
        TH3D *h3_efficiency_HFptetajetpt = (TH3D *)h3_num_efficiency_HFptetajetpt->Clone("h3_efficiency_HFptetajetpt");
        h3_efficiency_HFptetajetpt->Divide(h3_num_efficiency_HFptetajetpt, h3_denom_efficiency_HFptetajetpt, 1, 1, "B");
                                                
        TH2D *h2_purity_HFptjetpt = (TH2D *)h2_num_purity_HFptjetpt->Clone("h2_purity_HFptjetpt");
        h2_purity_HFptjetpt->Divide(h2_num_purity_HFptjetpt, h2_denom_purity_HFptjetpt, 1, 1, "B");
                
        TH2D *h2_efficiency_HFptjetpt = (TH2D *)h2_num_efficiency_HFptjetpt->Clone("h2_efficiency_HFptjetpt");
        h2_efficiency_HFptjetpt->Divide(h2_num_efficiency_HFptjetpt, h2_denom_efficiency_HFptjetpt, 1, 1, "B");
                                                
        TH1D *h1_purity_HFpt = (TH1D *)h1_num_purity_HFpt->Clone("h1_purity_HFpt");
        h1_purity_HFpt->Divide(h1_num_purity_HFpt, h1_denom_purity_HFpt, 1, 1, "B");
                
        TH1D *h1_efficiency_HFpt = (TH1D *)h1_num_efficiency_HFpt->Clone("h1_efficiency_HFpt");
        h1_efficiency_HFpt->Divide(h1_num_efficiency_HFpt, h1_denom_efficiency_HFpt, 1, 1, "B");
                                                
        TH1D *h1_purity_jetpt = (TH1D *)h1_num_purity_jetpt->Clone("h1_purity_jetpt");
        h1_purity_jetpt->Divide(h1_num_purity_jetpt, h1_denom_purity_jetpt, 1, 1, "B");
                
        TH1D *h1_efficiency_jetpt = (TH1D *)h1_num_efficiency_jetpt->Clone("h1_efficiency_jetpt");
        h1_efficiency_jetpt->Divide(h1_num_efficiency_jetpt, h1_denom_efficiency_jetpt, 1, 1, "B");
                                                
        h3_efficiency_rl_jetpt_weight->Divide(h3_num_efficiency_rl_jetpt_weight, h3_denom_efficiency_rl_jetpt_weight, 1, 1, "B");
        h3_efficiency_rl_jetpt_weight_eqch->Divide(h3_num_efficiency_rl_jetpt_weight_eqch, h3_denom_efficiency_rl_jetpt_weight_eqch, 1, 1, "B");
        h3_efficiency_rl_jetpt_weight_neqch->Divide(h3_num_efficiency_rl_jetpt_weight_neqch, h3_denom_efficiency_rl_jetpt_weight_neqch, 1, 1, "B");    
        h3_efficiency_rl_jetpt_weight_whf->Divide(h3_num_efficiency_rl_jetpt_weight_whf, h3_denom_efficiency_rl_jetpt_weight_whf, 1, 1, "B");
        h3_efficiency_rl_jetpt_weight_wohf->Divide(h3_num_efficiency_rl_jetpt_weight_wohf, h3_denom_efficiency_rl_jetpt_weight_wohf, 1, 1, "B");
        
        h3_efficiency_rl_jetpt_HFpt->Divide(h3_num_efficiency_rl_jetpt_HFpt, h3_denom_efficiency_rl_jetpt_HFpt, 1, 1, "B");
        
        h_efficiency_rl_jetptHFpt_weight->Divide(h_num_efficiency_rl_jetptHFpt_weight, h_denom_efficiency_rl_jetptHFpt_weight, 1, 1, "B");
        
        h2_efficiency_rl_jetpt->Divide(h2_num_efficiency_rl_jetpt, h2_denom_efficiency_rl_jetpt, 1, 1, "B");
        h2_efficiency_rl_jetpt_whf->Divide(h2_num_efficiency_rl_jetpt_whf, h2_denom_efficiency_rl_jetpt_whf, 1, 1, "B");
        h2_efficiency_rl_jetpt_wohf->Divide(h2_num_efficiency_rl_jetpt_wohf, h2_denom_efficiency_rl_jetpt_wohf, 1, 1, "B");

        h3_purity_rl_jetpt_weight->Divide(h3_num_purity_rl_jetpt_weight, h3_denom_purity_rl_jetpt_weight, 1, 1, "B");
        h3_purity_rl_jetpt_weight_eqch->Divide(h3_num_purity_rl_jetpt_weight_eqch, h3_denom_purity_rl_jetpt_weight_eqch, 1, 1, "B");
        h3_purity_rl_jetpt_weight_neqch->Divide(h3_num_purity_rl_jetpt_weight_neqch, h3_denom_purity_rl_jetpt_weight_neqch, 1, 1, "B");
        h3_purity_rl_jetpt_weight_whf->Divide(h3_num_purity_rl_jetpt_weight_whf, h3_denom_purity_rl_jetpt_weight_whf, 1, 1, "B");
        h3_purity_rl_jetpt_weight_wohf->Divide(h3_num_purity_rl_jetpt_weight_wohf, h3_denom_purity_rl_jetpt_weight_wohf, 1, 1, "B");
        
        h3_purity_rl_jetpt_HFpt->Divide(h3_num_purity_rl_jetpt_HFpt, h3_denom_purity_rl_jetpt_HFpt, 1, 1, "B");

        h_purity_rl_jetptHFpt_weight->Divide(h_num_purity_rl_jetptHFpt_weight, h_denom_purity_rl_jetptHFpt_weight, 1, 1, "B");
        
        h2_purity_rl_jetpt->Divide(h2_num_purity_rl_jetpt, h2_denom_purity_rl_jetpt, 1, 1, "B");
        h2_purity_rl_jetpt_whf->Divide(h2_num_purity_rl_jetpt_whf, h2_denom_purity_rl_jetpt_whf, 1, 1, "B");
        h2_purity_rl_jetpt_wohf->Divide(h2_num_purity_rl_jetpt_wohf, h2_denom_purity_rl_jetpt_wohf, 1, 1, "B");
               
        h3_efficiency_HFptetajetpt->Write("efficiency_HFptetajetpt");
        h3_purity_HFptetajetpt->Write("purity_HFptetajetpt");
        
        h2_efficiency_HFptjetpt->Write("efficiency_HFptjetpt");
        h2_purity_HFptjetpt->Write("purity_HFptjetpt");
        
        h1_efficiency_HFpt->Write("efficiency_HFpt");
        h1_purity_HFpt->Write("purity_HFpt");
        
        h1_efficiency_jetpt->Write("efficiency_jetpt");
        h1_purity_jetpt->Write("purity_jetpt");
        
        h2_efficiency_rl_jetpt->Write("efficiency_rl_jetpt");
        h2_efficiency_rl_jetpt_whf->Write("efficiency_rl_jetpt_whf");
        h2_efficiency_rl_jetpt_wohf->Write("efficiency_rl_jetpt_wohf");

        h3_efficiency_rl_jetpt_weight->Write("efficiency_rl_jetpt_weight");
        h3_efficiency_rl_jetpt_weight_whf->Write("efficiency_rl_jetpt_weight_whf");
        h3_efficiency_rl_jetpt_weight_wohf->Write("efficiency_rl_jetpt_weight_wohf");
        h3_efficiency_rl_jetpt_weight_eqch->Write("efficiency_rl_jetpt_weight_eqch");
        h3_efficiency_rl_jetpt_weight_neqch->Write("efficiency_rl_jetpt_weight_neqch");    
        
        h3_efficiency_rl_jetpt_HFpt->Write("efficiency_rl_jetpt_HFpt");

        h_efficiency_rl_jetptHFpt_weight->Write("efficiency_rl_jetptHFpt_weight");
        
        h2_purity_rl_jetpt->Write("purity_rl_jetpt");
        h2_purity_rl_jetpt_whf->Write("purity_rl_jetpt_whf");
        h2_purity_rl_jetpt_wohf->Write("purity_rl_jetpt_wohf");

        h3_purity_rl_jetpt_weight->Write("purity_rl_jetpt_weight");
        h3_purity_rl_jetpt_weight_whf->Write("purity_rl_jetpt_weight_whf");
        h3_purity_rl_jetpt_weight_wohf->Write("purity_rl_jetpt_weight_wohf");
        h3_purity_rl_jetpt_weight_eqch->Write("purity_rl_jetpt_weight_eqch");
        h3_purity_rl_jetpt_weight_neqch->Write("purity_rl_jetpt_weight_neqch");    
        
        h3_purity_rl_jetpt_HFpt->Write("purity_rl_jetpt_HFpt");

        h_purity_rl_jetptHFpt_weight->Write("purity_rl_jetptHFpt_weight");
        
        f->Close();
}
