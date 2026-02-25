
#include "analysis-binning.h"
#include "utils.h"
#include <iostream>

double get_hwhm(TH1F* h)
{
        // As written, this function works for single peak distributions
        double half_max = h->GetMaximum()/2.;
    
        int halfwidth_bin;
        for (int bin = 1 ; bin <= h->GetNbinsX() ; bin++)
                if (h->GetBinContent(bin) >= half_max) {
                        halfwidth_bin = bin; 
                        break;
                }
        
        
        // Return the Half width at half maximum value
        return std::abs(h->GetBinCenter(h->GetMaximumBin()) - h->GetBinCenter(halfwidth_bin));
}

double get_median_from_cumulative(TH1F* h)
{
        double total_entries = h->GetMaximum();
        double median = 0;
        std::cout<<"Maximum = "<<total_entries<<std::endl;

        for (int bin = 1 ; bin <= h->GetNbinsX() ; bin++) {
                if (h->GetBinContent(bin) > total_entries / 2.) {
                        median = h->GetBinCenter(bin);
                        break;
                }
        }
        
        return median;
}

void determine_log10binning(int Nbins, double x_i, double x_f, double* binning)
{
        double log_x_i = log10(x_i);
        double log_x_f = log10(x_f);
        double delta = (log_x_f-log_x_i)/Nbins;

        for (int i = 0 ; i <=Nbins ; i++) 
                binning[i] = pow(10,log_x_i+delta*i);
}

void determine_eqsizebinning(int Nbins, double x_i, double x_f, double* binning)
{
        double delta = (x_f - x_i)/Nbins;

        for (int i = 0;i <= Nbins;i++)
                binning[i] = x_i + delta*i;
}

void apply_unfolded_weights(TH3F* h_rl_jetpt_weight, TH2F* h_rl_jetpt)
{
        for (int i = 1 ; i <= h_rl_jetpt_weight->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h_rl_jetpt_weight->GetNbinsY() ; j++) {
                        double weighted_content = 0;
                        double weighted_content_err = 0;
                        
                        for (int k = 1 ; k <= h_rl_jetpt_weight->GetNbinsZ() ; k++) {
                                double weight = h_rl_jetpt_weight->GetZaxis()->GetBinCenter(k);
                                weighted_content += h_rl_jetpt_weight->GetBinContent(i,j,k)*weight;
                                weighted_content_err += pow(weight*h_rl_jetpt_weight->GetBinError(i,j,k),2);
                        }

                        h_rl_jetpt->SetBinContent(i, j, weighted_content);
                        h_rl_jetpt->SetBinError(i, j, sqrt(weighted_content_err));
                }
        }
}

void apply_unfolded_weights(TH3D* h_rl_jetpt_weight, TH2D* h_rl_jetpt)
{
        for (int i = 1 ; i <= h_rl_jetpt_weight->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h_rl_jetpt_weight->GetNbinsY() ; j++) {
                        double weighted_content = 0;
                        double weighted_content_err = 0;
                        
                        for (int k = 1 ; k <= h_rl_jetpt_weight->GetNbinsZ() ; k++) {
                                double weight = h_rl_jetpt_weight->GetZaxis()->GetBinCenter(k);
                                
                                weighted_content += h_rl_jetpt_weight->GetBinContent(i,j,k) * weight;
                                weighted_content_err += pow(weight*h_rl_jetpt_weight->GetBinError(i,j,k),2);
                        }

                        h_rl_jetpt->SetBinContent(i, j, weighted_content);
                        h_rl_jetpt->SetBinError(i, j, sqrt(weighted_content_err));
                }
        }
}

void apply_jet_weight_to_npairs(TH3D* h_npair, TH1F* h_purity_jet, TH1F* h_efficiency_jet)
{
        for (int i = 1 ; i <= h_npair->GetNbinsX(); i++) {
                for (int j = 1 ; j <= h_npair->GetNbinsY(); j++) {
                        for (int k = 1 ; k <= h_npair->GetNbinsZ() ; k++) {
                                double reweight = h_purity_jet->GetBinContent(j)/h_efficiency_jet->GetBinContent(j);
                                
                                h_npair->SetBinContent(i, j, k, h_npair->GetBinContent(i, j, k) * reweight);
                                h_npair->SetBinError(i, j, k, h_npair->GetBinError(i, j, k) * reweight);
                        }
                }
        }
}

void project_nominal_phase_space(TH2D* h_2d, TH1F* h_1d, int nominal_jet_pt_bin)
{
        // This function should be used only in the case of two jetpt underflow bins and one underflow/overflow for RL
        for (int i = 1 ; i <= h_1d->GetNbinsX(); i++) {
                if (i == 1 || i == h_1d->GetNbinsX()) {
                        h_1d->SetBinContent(i, 0);
                        h_1d->SetBinError(i, 0);
                        continue;
                }
                h_1d->SetBinContent(i, h_2d->GetBinContent(i, nominal_jet_pt_bin));
                h_1d->SetBinError(i, h_2d->GetBinError(i, nominal_jet_pt_bin));
        }
}

void substract_stat_error(TH1F* hnominal, TH1F* hsyst)
{
        for (int i = 1 ; i <= hnominal->GetNbinsX() ; i++) {
                double total_error = hsyst->GetBinError(i);
                double stat_error  = hnominal->GetBinError(i);
                double syst_error  = std::sqrt(total_error * total_error - stat_error * stat_error);

                hsyst->SetBinError(i, syst_error);
        }
}

void get_tau_from_uoflow_eec(TH1F* h_eec, TH1F* h_tau, double avge_pt2_jet)
{
        for (int i = 1 ; i <= h_eec->GetNbinsX() ; i++) {
                if (i == 1 || i == h_eec->GetNbinsX())
                        continue;
                
                h_tau->SetBinContent(i - 1, h_eec->GetBinContent(i));
                h_tau->SetBinError(i -1, h_eec->GetBinError(i));
        }

        h_tau->Scale(log(avge_pt2_jet) / avge_pt2_jet);
}

void get_tau_binning_from_eec_binning(double* tau_binning, const double* eec_binning, double average_pt2_jet)
{
        for(int i = 0 ; i <= nbin_rl_nominal ; i++)
                tau_binning[i] = eec_binning[i] * average_pt2_jet;
}

void normalize_by_njets(TH1F* h, double h_njet_content, double h_njet_error)
{
        if (h_njet_content < h_njet_error)
                std::cout<<"The jet error is bigge than the content. Something is wrong!"<<std::endl;
        
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                double new_content = h->GetBinContent(i) / h_njet_content;
                double new_error = sqrt(pow(h->GetBinError(i), 2) + pow(h_njet_error * h->GetBinContent(i) / h_njet_content, 2)) / h_njet_content;

                h->SetBinContent(i, new_content);
                h->SetBinError(i, new_error);
        }
}

void normalize_by_njets(TH1D* h, double h_njet_content, double h_njet_error)
{
        if (h_njet_content < h_njet_error)
                std::cout<<"The jet error is bigge than the content. Something is wrong!"<<std::endl;
        
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                double new_content = h->GetBinContent(i) / h_njet_content;
                double new_error = sqrt(pow(h->GetBinError(i), 2) + pow(h_njet_error * h->GetBinContent(i) / h_njet_content, 2)) / h_njet_content;

                h->SetBinContent(i, new_content);
                h->SetBinError(i, new_error);
        }
}

void square_root_bins(TH1F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                h->SetBinContent(i, std::sqrt(h->GetBinContent(i)));
                h->SetBinError(i, std::sqrt(h->GetBinError(i)));
        }
}

void regularize_correction_factors(TH2F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        double bin_content = h->GetBinContent(i, j);

                        if (bin_content > 1)
                                h->SetBinContent(i, j, 1.);
                }
        }
}

void regularize_correction_factors(TH2D* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        double bin_content = h->GetBinContent(i, j);

                        if (bin_content > 1)
                                h->SetBinContent(i, j, 1.);
                }
        }
}

void regularize_correction_factors(TH3F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        for (int k = 1 ; k <= h->GetNbinsZ() ; k++){
                                double bin_content = h->GetBinContent(i, j, k);

                                if (bin_content < 0)
                                        h->SetBinContent(i, j, k, 0.);
                                        
                                if (bin_content > 1)
                                        h->SetBinContent(i, j, k, 1.);
                        }
                }
        }
}

void regularize_correction_factors(TH3D* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                for (int j = 1 ; j <= h->GetNbinsY() ; j++) {
                        for (int k = 1 ; k <= h->GetNbinsZ() ; k++){
                                double bin_content = h->GetBinContent(i, j, k);

                                if (bin_content > 1)
                                        h->SetBinContent(i, j, k, 1.);
                        }
                }
        }
}

void set_histo_with_systematics(TH1F* hrelerror, TH1F* hnominal, TH1F* hsystematic, int syst_index, bool print_table)
{
        double total_err = 0, total = 0;
        double n_ct_sources = (syst_index != 2 && syst_index != 5 && syst_index != 6) ? 4. : 1.;

        for (int hbin = 1 ; hbin <= hrelerror->GetNbinsX() ; hbin++)
        {
                double dev     = hrelerror->GetBinContent(hbin);
                double dev_err = hrelerror->GetBinError(hbin);
                
                // Dev values are positive by construction
                total += hnominal->GetBinContent(hbin);

                // if (dev-dev_err < 0) 
                //         continue;

                double syst_error            = std::abs(dev)*hnominal->GetBinContent(hbin);
                double syst_error_percentage = std::abs(dev);
                
                hsystematic->SetBinError(hbin, std::sqrt(syst_error * syst_error / n_ct_sources + std::pow(hsystematic->GetBinError(hbin),2)));

                total_err += hnominal->GetBinContent(hbin)*syst_error_percentage*syst_error_percentage;
        }

        if (!print_table)
                return;

        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        std::cout.precision(2);

        if (total > 0)
                std::cout<<std::sqrt(total_err/total)*100.;
        else if (total == 0)
                std::cout<<0;
}

void set_histoa_errors_as_histob_content(TH1F* ha, TH1F* hb)
{
        for (int i = 1 ; i <= ha->GetNbinsX() ; i++) {
                double error = ha->GetBinError(i);
                double content = ha->GetBinContent(i);
                double nominal_relerror = error/content;

                if (std::isnan(nominal_relerror))
                        nominal_relerror = 0;

                hb->SetBinContent(i, nominal_relerror);
        }
}

void set_histo_sqrt_content(TH1F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++) {
                h->SetBinContent(i, sqrt(h->GetBinContent(i)));
                h->SetBinError(i, sqrt(h->GetBinError(i)));
        }
}

void set_histo_null_errors(TH1F* h)
{
        for (int i = 1 ; i <= h->GetNbinsX() ; i++)
                h->SetBinError(i, 0);
}

void set_shift_histo(TH2F* href, TH2F* hshift, TRandom3* rndm) // CHECK THIS FUNCTION!!
{
        for (int xbin = 1 ; xbin <= href->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= href->GetNbinsY() ; ybin++) {
                        double shift_window = href->GetBinError(href->GetBin(xbin, ybin));
                        double refdata      = href->GetBinContent(href->GetBin(xbin, ybin));
                        double shift        = rndm->Gaus(refdata, shift_window)/refdata;
                        hshift->SetBinContent(xbin, ybin, shift);
                }
        }
}

void set_shift_histo(TH2D* href, TH2D* hshift, TRandom3* rndm) // CHECK THIS FUNCTION!!
{
        for (int xbin = 1 ; xbin <= href->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= href->GetNbinsY() ; ybin++) {
                        double shift_window = href->GetBinError(href->GetBin(xbin, ybin));
                        double refdata      = href->GetBinContent(href->GetBin(xbin, ybin));
                        double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                        hshift->SetBinContent(xbin, ybin, shift);
                        hshift->SetBinError(xbin, ybin, shift_window/refdata);
                }
        }
}

void set_shift_histo(TH3D* href, TH3D* hshift, TRandom3* rndm) // CHECK THIS FUNCTION!!
{
        for (int xbin = 1 ; xbin <= href->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= href->GetNbinsY() ; ybin++) {
                        for (int zbin = 1 ; zbin <= href->GetNbinsZ() ; zbin++) {
                                double shift_window = href->GetBinError(href->GetBin(xbin, ybin, zbin));
                                double refdata      = href->GetBinContent(href->GetBin(xbin, ybin, zbin));
                                double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                                if (std::isnan(shift)){
                                        hshift->SetBinContent(xbin, ybin, zbin, 0);
                                        hshift->SetBinError(xbin, ybin, zbin, 0);

                                        continue;
                                }
                                
                                hshift->SetBinContent(xbin, ybin, zbin, shift);
                                hshift->SetBinError(xbin, ybin, zbin, shift_window/refdata);
                        }
                }
        }
}

void smear_pseudodata(TH1F* hpseudodata, TH1F* hrealdata, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= hpseudodata->GetNbinsX() ; xbin++) {
                double shift_window = hrealdata->GetBinError(hrealdata->GetBin(xbin));
                double refdata      = hrealdata->GetBinContent(hrealdata->GetBin(xbin));
                double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                if (std::isnan(shift)){
                        hpseudodata->SetBinContent(xbin, 0);
                        hpseudodata->SetBinError(xbin, 0);

                        continue;
                }

                double pseudodata_content = shift * hpseudodata->GetBinContent(hpseudodata->GetBin(xbin));
                
                hpseudodata->SetBinContent(xbin, pseudodata_content);
        }
}

void smear_pseudodata(TH3D* hpseudodata, TH3D* hrealdata, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= hpseudodata->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= hpseudodata->GetNbinsY() ; ybin++) {
                        for (int zbin = 1 ; zbin <= hpseudodata->GetNbinsZ() ; zbin++) {
                                double shift_window = hrealdata->GetBinError(hrealdata->GetBin(xbin, ybin, zbin));
                                double refdata      = hrealdata->GetBinContent(hrealdata->GetBin(xbin, ybin, zbin));
                                double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                                if (std::isnan(shift)){
                                        hpseudodata->SetBinContent(xbin, ybin, zbin, 0);
                                        hpseudodata->SetBinError(xbin, ybin, zbin, 0);

                                        continue;
                                }

                                double pseudodata_content = shift * hpseudodata->GetBinContent(hpseudodata->GetBin(xbin, ybin, zbin));
                                
                                hpseudodata->SetBinContent(xbin, ybin, zbin, pseudodata_content);
                        }
                }
        }
}

void smear_pseudodata(TH3F* hpseudodata, TH3F* hrealdata, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= hpseudodata->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= hpseudodata->GetNbinsY() ; ybin++) {
                        for (int zbin = 1 ; zbin <= hpseudodata->GetNbinsZ() ; zbin++) {
                                double shift_window = hrealdata->GetBinError(hrealdata->GetBin(xbin, ybin, zbin));
                                double refdata      = hrealdata->GetBinContent(hrealdata->GetBin(xbin, ybin, zbin));
                                double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                                if (std::isnan(shift)){
                                        hpseudodata->SetBinContent(xbin, ybin, zbin, 0);
                                        hpseudodata->SetBinError(xbin, ybin, zbin, 0);

                                        continue;
                                }

                                double pseudodata_content = shift * hpseudodata->GetBinContent(hpseudodata->GetBin(xbin, ybin, zbin));
                                
                                hpseudodata->SetBinContent(xbin, ybin, zbin, pseudodata_content);
                        }
                }
        }
}

void smear_pseudodata(TH2D* hpseudodata, TH2D* hrealdata, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= hpseudodata->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= hpseudodata->GetNbinsY() ; ybin++) {
                        double shift_window = hrealdata->GetBinError(hrealdata->GetBin(xbin, ybin));
                        double refdata      = hrealdata->GetBinContent(hrealdata->GetBin(xbin, ybin));
                        double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                        if (std::isnan(shift)){
                                hpseudodata->SetBinContent(xbin, ybin, 0);
                                hpseudodata->SetBinError(xbin, ybin, 0);

                                continue;
                        }

                        double pseudodata_content = shift * hpseudodata->GetBinContent(hpseudodata->GetBin(xbin, ybin));
                        
                        hpseudodata->SetBinContent(xbin, ybin, pseudodata_content);
                }
        }
}

void smear_pseudodata(TH2F* hpseudodata, TH2F* hrealdata, TRandom3* rndm)
{
        for (int xbin = 1 ; xbin <= hpseudodata->GetNbinsX() ; xbin++) {
                for (int ybin = 1 ; ybin <= hpseudodata->GetNbinsY() ; ybin++) {
                        double shift_window = hrealdata->GetBinError(hrealdata->GetBin(xbin, ybin));
                        double refdata      = hrealdata->GetBinContent(hrealdata->GetBin(xbin, ybin));
                        double shift        = rndm->Gaus(refdata, shift_window)/refdata;

                        if (std::isnan(shift)){
                                hpseudodata->SetBinContent(xbin, ybin, 0);
                                hpseudodata->SetBinError(xbin, ybin, 0);

                                continue;
                        }

                        double pseudodata_content = shift * hpseudodata->GetBinContent(hpseudodata->GetBin(xbin, ybin));
                        
                        hpseudodata->SetBinContent(xbin, ybin, pseudodata_content);
                }
        }
}

void set_data_ntuple_branches(TNtuple* ntuple, float* event_weight, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror)
{
        ntuple->SetBranchAddress("event_weight", event_weight);
        ntuple->SetBranchAddress("R_L", R_L);
        ntuple->SetBranchAddress("jet_pt", jet_pt);
        ntuple->SetBranchAddress("weight_pt", weight_pt);
        ntuple->SetBranchAddress("efficiency", efficiency);
        ntuple->SetBranchAddress("efficiency_relerror", efficiency_relerror);
        ntuple->SetBranchAddress("purity", purity);
        ntuple->SetBranchAddress("purity_relerror", purity_relerror);
}

void set_data_ntuple_branches(TNtuple* ntuple, float* event_weight, float* R_L, float* jet_pt, float* weight_pt, float* efficiency, float* purity, float* efficiency_relerror, float* purity_relerror, float* eq_charge)
{
        ntuple->SetBranchAddress("event_weight", event_weight);
        ntuple->SetBranchAddress("R_L", R_L);
        ntuple->SetBranchAddress("jet_pt", jet_pt);
        ntuple->SetBranchAddress("weight_pt", weight_pt);
        ntuple->SetBranchAddress("efficiency", efficiency);
        ntuple->SetBranchAddress("efficiency_relerror", efficiency_relerror);
        ntuple->SetBranchAddress("purity", purity);
        ntuple->SetBranchAddress("purity_relerror", purity_relerror);
        ntuple->SetBranchAddress("eq_charge", eq_charge);
}

void set_unfolding_ntuple_branches(TNtuple* ntuple, float* R_L_reco, float* R_L_truth, float* jet_pt_reco, float* jet_pt_truth, float* weight_pt_reco, float* weight_pt_truth)
{
        ntuple->SetBranchAddress("jet_pt", jet_pt_reco);
        ntuple->SetBranchAddress("jet_pt_truth", jet_pt_truth);
        ntuple->SetBranchAddress("R_L", R_L_reco);
        ntuple->SetBranchAddress("R_L_truth", R_L_truth);
        ntuple->SetBranchAddress("weight_pt", weight_pt_reco);
        ntuple->SetBranchAddress("weight_pt_truth", weight_pt_truth);
}

void set_unfolding_jet_ntuple_branches(TNtuple* ntuple, float* jet_pt_reco, float* jet_pt_truth)
{
        ntuple->SetBranchAddress("jet_pt", jet_pt_reco);
        ntuple->SetBranchAddress("jet_pt_truth", jet_pt_truth);
}

void set_unity_content(TH1F* h)
{
        for (int bin = 1 ; bin <= h->GetNbinsX() ; bin++) {
                h->SetBinContent(bin, 1.);
        }
}

double get_jes_jer_factor(const double jet_pt, TRandom3 *myRNG) 
{
        double jes_var = 1, jer_var = 1;

        if (jet_pt < 10) { 
                jes_var = 0.912; 
                jer_var = 0.006; 
        } else if (jet_pt < 12) { 
                jes_var = 0.948; 
                jer_var = 0.132; 
        } else if (jet_pt < 15) { 
                jes_var = 0.944; 
                jer_var = 0.192; 
        } else if (jet_pt < 17) {
                jes_var = 0.976;
                jer_var = 0.114; 
        } else if (jet_pt < 20) { 
                jes_var = 0.972; 
                jer_var = 0.096; 
        } else if (jet_pt < 25) { 
                jes_var = 0.976; 
                jer_var = 0.138; 
        } else if (jet_pt < 35) {
                jes_var = 0.984;
                jer_var = 0.054; 
        } else if (jet_pt < 50) { 
                jes_var = 1.000; 
                jer_var = 0.144; 
        } else { 
                jes_var = 1.008; 
                jer_var = 0.120; 
        }

        return jes_var * myRNG->Gaus(1, jer_var);
}

double weight(double h1_E, double h2_E, double jet_E)
{
        return h1_E*h2_E/jet_E/jet_E;
}

void smooth_nominal_phase_space(TH1F* h_to_smooth, TH1F* h_nominal_phase_space)
{
        for (int bin = 2 ; bin < h_to_smooth->GetNbinsX(); bin++) {
                h_nominal_phase_space->SetBinContent(bin - 1, h_to_smooth->GetBinContent(bin));
                h_nominal_phase_space->SetBinError(bin - 1, h_to_smooth->GetBinError(bin));
        }

        h_nominal_phase_space->Smooth();

        for (int bin = 2 ; bin < h_to_smooth->GetNbinsX(); bin++) {
                h_to_smooth->SetBinContent(bin, h_nominal_phase_space->GetBinContent(bin - 1));
                h_to_smooth->SetBinError(bin, h_nominal_phase_space->GetBinError(bin - 1));
        }
}
