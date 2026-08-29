#ifndef NAMES_H
#define NAMES_H

#include<map>

// Names of the files

// BJETS
std::string namef_ntuple_data          = "ntuple_bjets_data.root";
std::string namef_ntuple_data_jetid    = "ntuple_bjets_data_jetid.root";
std::string namef_ntuple_mcreco        = "ntuple_bjets_mcreco.root";
std::string namef_ntuple_mcreco_jetid  = "ntuple_bjets_mcreco_jetid.root";
std::string namef_ntuple_mcreco_jesjer = "ntuple_bjets_mcreco_jesjer.root";
std::string namef_ntuple_mc            = "ntuple_bjets_mc.root";
std::string namef_ntuple_bjets_misid   = "ntuple_bjets_misid.root";

std::string namef_data2mcreco_ratio = "data2mcreco_rl_jetpt_weight.root";
std::string name_histo_data2mcreco_rl_jetpt_weight = "data2mcreco_rl_jetpt_weight_ratio";
std::string name_histo_data2mcreco_HFpt_eta_jetpt  = "data2mcreco_HFpt_eta_jetpt_ratio";
std::string name_histo_data2mcreco_jetpt           = "data2mcreco_jetpt_ratio";


std::map<std::string, std::string> namef_data_variations = {
        {"nominal"      , namef_ntuple_data},
        {"jetid"        , namef_ntuple_data_jetid},
        {"jesjer"       , namef_ntuple_data},
        {"prior"        , namef_ntuple_data},
        {"recseleff"    , namef_ntuple_data},
        {"fitsignal"    , namef_ntuple_data},
        {"trackeff-hi"  , namef_ntuple_data},
        {"trackeff-lo"  , namef_ntuple_data},
        {"trigeff-hi"   , namef_ntuple_data},
        {"trigeff-lo"   , namef_ntuple_data},
        {"pideff-hi"    , namef_ntuple_data},
        {"pideff-lo"    , namef_ntuple_data},
        {"massfit-far"  , namef_ntuple_data},
        {"massfit-near" , namef_ntuple_data},
        {"massfit-upper", namef_ntuple_data},
        {"massfit-lower", namef_ntuple_data},
        {"trackeff_cc"  , namef_ntuple_data},
};

std::map<std::string, std::string> namef_mcreco_variations = {
        {"nominal"      , namef_ntuple_mcreco},
        {"jetid"        , namef_ntuple_mcreco_jetid},
        {"jesjer"       , namef_ntuple_mcreco_jesjer},
        {"prior"        , namef_ntuple_mcreco},
        {"recseleff"    , namef_ntuple_mcreco},
        {"fitsignal"    , namef_ntuple_mcreco},
        {"trackeff-hi"  , namef_ntuple_mcreco},
        {"trackeff-lo"  , namef_ntuple_mcreco},
        {"trigeff-hi"   , namef_ntuple_mcreco},
        {"trigeff-lo"   , namef_ntuple_mcreco},
        {"pideff-hi"    , namef_ntuple_mcreco},
        {"pideff-lo"    , namef_ntuple_mcreco},
        {"massfit-far"  , namef_ntuple_mcreco},
        {"massfit-near" , namef_ntuple_mcreco},
        {"massfit-upper", namef_ntuple_mcreco},
        {"massfit-lower", namef_ntuple_mcreco},
        {"trackeff_cc"  , namef_ntuple_mcreco},
};

std::map<std::string, std::string> namef_corrections = {
        {"nominal"      , "bjets_corrections.root"},
        {"jetid"        , "bjets_corrections_jetid.root"},
        {"jesjer"       , "bjets_corrections_jesjer.root"},
        {"prior"        , "bjets_corrections_prior.root"},
        {"recseleff"    , "bjets_corrections_recseleff.root"},
        {"fitsignal"    , "bjets_corrections_recseleff.root"},
        {"trackeff-hi"  , "bjets_corrections.root"},
        {"trackeff-lo"  , "bjets_corrections.root"},
        {"trigeff-hi"   , "bjets_corrections.root"},
        {"trigeff-lo"   , "bjets_corrections.root"},
        {"pideff-hi"    , "bjets_corrections.root"},
        {"pideff-lo"    , "bjets_corrections.root"},
        {"massfit-far"  , "bjets_corrections.root"},
        {"massfit-near" , "bjets_corrections.root"},
        {"massfit-upper", "bjets_corrections.root"},
        {"massfit-lower", "bjets_corrections.root"},
        {"trackeff_cc"  , "bjets_corrections.root"},
};

std::map<std::string, std::string> namef_simpleobservable_data = {
        {"nominal"      , "bjets_simpleobservable_data.root"},
        {"jetid"        , "bjets_simpleobservable_data_jetid.root"},
        {"jesjer"       , "bjets_simpleobservable_data_jesjer.root"},
        {"prior"        , "bjets_simpleobservable_data_prior.root"},
        {"recseleff"    , "bjets_simpleobservable_data_recseleff.root"},
        {"fitsignal"    , "bjets_simpleobservable_data_fitsignal.root"},
        {"trackeff-hi"  , "bjets_simpleobservable_data_trackeff_hi.root"},
        {"trackeff-lo"  , "bjets_simpleobservable_data_trackeff_lo.root"},
        {"trigeff-hi"   , "bjets_simpleobservable_data_trigeff_hi.root"},
        {"trigeff-lo"   , "bjets_simpleobservable_data_trigeff_lo.root"},
        {"pideff-hi"    , "bjets_simpleobservable_data_pideff_hi.root"},
        {"pideff-lo"    , "bjets_simpleobservable_data_pideff_lo.root"},
        {"massfit-far"  , "bjets_simpleobservable_data_massfit_far.root"},
        {"massfit-near" , "bjets_simpleobservable_data_massfit_near.root"},
        {"massfit-upper", "bjets_simpleobservable_data_massfit_upper.root"},
        {"massfit-lower", "bjets_simpleobservable_data_massfit_lower.root"},
        {"trackeff-cc"  , "bjets_simpleobservable_data_trackeff_cc.root"},
};

std::map<std::string, std::string> namef_simpleobservable_mcreco = {
        {"nominal"      , "bjets_simpleobservable_mcreco.root"},
        {"jetid"        , "bjets_simpleobservable_mcreco_jetid.root"},
        {"jesjer"       , "bjets_simpleobservable_mcreco_jesjer.root"},
        {"prior"        , "bjets_simpleobservable_mcreco_prior.root"},
        {"recseleff"    , "bjets_simpleobservable_mcreco_recseleff.root"},
        {"fitsignal"    , "bjets_simpleobservable_mcreco_fitsignal.root"},
        {"trackeff-hi"  , "bjets_simpleobservable_mcreco.root"},
        {"trackeff-lo"  , "bjets_simpleobservable_mcreco.root"},
        {"trigeff-hi"   , "bjets_simpleobservable_mcreco.root"},
        {"trigeff-lo"   , "bjets_simpleobservable_mcreco.root"},
        {"pideff-hi"    , "bjets_simpleobservable_mcreco.root"},
        {"pideff-lo"    , "bjets_simpleobservable_mcreco.root"},
        {"massfit-far"  , "bjets_simpleobservable_mcreco.root"},
        {"massfit-near" , "bjets_simpleobservable_mcreco.root"},
        {"massfit-upper", "bjets_simpleobservable_mcreco.root"},
        {"massfit-lower", "bjets_simpleobservable_mcreco.root"},
        {"trackeff-cc"  , "bjets_simpleobservable_mcreco.root"},
};

std::map<std::string, std::string> namef_massfits_results_data = {
        {"nominal"      , "results_mass_fit_data.root"},
        {"jetid"        , "results_mass_fit_data_jetid.root"},
        {"jesjer"       , "results_mass_fit_data_jesjer.root"},
        {"prior"        , "results_mass_fit_data.root"},
        {"recseleff"    , "results_mass_fit_data_recseleff.root"},
        {"fitsignal"    , "results_mass_fit_data_fitsignal.root"},
        {"trackeff-hi"  , "results_mass_fit_data.root"},
        {"trackeff-lo"  , "results_mass_fit_data.root"},
        {"trigeff-hi"   , "results_mass_fit_data.root"},
        {"trigeff-lo"   , "results_mass_fit_data.root"},
        {"pideff-hi"    , "results_mass_fit_data.root"},
        {"pideff-lo"    , "results_mass_fit_data.root"},
        {"massfit-far"  , "results_mass_fit_data.root"},
        {"massfit-near" , "results_mass_fit_data.root"},
        {"massfit-upper", "results_mass_fit_data.root"},
        {"massfit-lower", "results_mass_fit_data.root"},
        {"trackeff-cc"  , "results_mass_fit_data.root"},
};

std::map<std::string, std::string> namef_massfits_results_mcreco = {
        {"nominal"      , "results_mass_fit_mcreco.root"},
        {"jetid"        , "results_mass_fit_mcreco_jetid.root"},
        {"jesjer"       , "results_mass_fit_mcreco_jesjer.root"},
        {"prior"        , "results_mass_fit_mcreco.root"},
        {"recseleff"    , "results_mass_fit_mcreco_recseleff.root"},
        {"fitsignal"    , "results_mass_fit_mcreco_fitsignal.root"},
        {"trackeff-hi"  , "results_mass_fit_mcreco.root"},
        {"trackeff-lo"  , "results_mass_fit_mcreco.root"},
        {"trigeff-hi"   , "results_mass_fit_mcreco.root"},
        {"trigeff-lo"   , "results_mass_fit_mcreco.root"},
        {"pideff-hi"    , "results_mass_fit_mcreco.root"},
        {"pideff-lo"    , "results_mass_fit_mcreco.root"},
        {"massfit-far"  , "results_mass_fit_mcreco.root"},
        {"massfit-near" , "results_mass_fit_mcreco.root"},
        {"massfit-upper", "results_mass_fit_mcreco.root"},
        {"massfit-lower", "results_mass_fit_mcreco.root"},
        {"trackeff-cc"  , "results_mass_fit_mcreco.root"},
};

// About systematics
std::string available_systematics[] = {
        "nominal"      ,
        "jetid"        ,
        "jesjer"       ,
        "prior"        ,
        "recseleff"    ,
        "fitsignal"    ,
        "trackeff-hi"  ,
        "trackeff-lo"  ,
        "trigeff-hi"   ,
        "trigeff-lo"   ,
        "pideff-hi"    ,
        "pideff-lo"    ,
        "massfit-far"  ,
        "massfit-near" ,
        "massfit-upper",
        "massfit-lower",
        "trackeff-cc"  ,
};

// std::map<std::string, std::string> systematic_name  = {
//         {available_systematics[0],"Closure test n.c."},
//         {available_systematics[1],"Shape closure test n.c."},
//         {available_systematics[2],"JES-JER"},
//         {available_systematics[3],"Prior variation"},
//         {available_systematics[4],"Regularization parameter"},
//         {available_systematics[5],"Muon eff"},
//         {available_systematics[6],"Jet ID"},
//         {available_systematics[7],"ProbNNghost"}
// };

// std::map<std::string, std::string> systematic_namef = {
//         {available_systematics[0],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4_statct_niterct1.root"},
//         {available_systematics[1],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4_shapect.root"},
//         {available_systematics[2],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-jes-jer.root"},
//         {available_systematics[3],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-prior.root"},
//         {available_systematics[4],"histos_eec_3dcorr_rl_jetpt_weightpt_niterjet4--get-regpar.root"},
//         {available_systematics[5],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-muon.root"},
//         {available_systematics[6],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-jetid.root"},
//         {available_systematics[7],"histos_eec_3dcorr_rl_jetpt_weightpt_niter4_niterjet4--get-probnnghost.root"}
// };

// std::map<std::string, std::string> systematic_errtype = {
//         {available_systematics[0],"normal"},
//         {available_systematics[1],"normal"},
//         {available_systematics[2],"normal"},
//         {available_systematics[3],"normal"},
//         {available_systematics[4],"normal"},
//         {available_systematics[5],"normal"},
//         {available_systematics[6],"normal"},
//         {available_systematics[7],"normal"}
// };

#endif