#include "Settings.h"
#include "../Helpers.h"

void analyzeStatClosure()
{
    TH1D *z3_pt4 = new TH1D("z3_pt4", "0.2 < z < 0.3, 20 < p_{T} < 30;ratio;occurances", 20, 0.5, 1.5);
    TH1D *z10_pt4 = new TH1D("z10_pt4", "0.7 < z < 0.75, 20 < p_{T} < 30;ratio;occurances", 20, 0.9, 1.1);

    gStyle->SetImageScaling(1.5);
    int nRuns = 10;
    TFile *infile = new TFile("../../root_files/BjetsMC/closure_ev_-1_ptj_5100_eta_2.54.0_HF_ghost_0.4_b_PID_iters_4_91599_91599.root", "READ");
    TFile *outfile = new TFile("../../root_files/BjetsMC/stat_closure_analysis.root", "RECREATE");
    for (int i=0; i<nRuns; i++)
    {
        TH1D *htemp = (TH1D*)infile->Get(Form("z_ratio_smear%d_pt4", i));
        std::cout << "htemp Bin 3 : " << htemp->GetBinContent(3) << std::endl; 
        z3_pt4->Fill(htemp->GetBinContent(3));
        z10_pt4->Fill(htemp->GetBinContent(10));
        //z3_pt4->Draw();
    }
    TCanvas *c0 = new TCanvas("c0");
    z3_pt4->Draw();
    c0->SaveAs("../../plots/BjetsMC/stat_closure_analysis/z3_pt4.png");
    z3_pt4->Write();

    //TCanvas *c1 = new TCanvas("c1");
    z10_pt4->Draw();
    c0->SaveAs("../../plots/BjetsMC/stat_closure_analysis/z10_pt4.png");
    z10_pt4->Write();

    TH2D *h2_zjt_closure_error[ptbinsize-1], *h2_zr_closure_error[ptbinsize-1], *h2_jtr_closure_error[ptbinsize-1];
    TH2D *h2_zjt_closure_error_pulls[ptbinsize-1], *h2_zr_closure_error_pulls[ptbinsize-1], *h2_jtr_closure_error_pulls[ptbinsize-1];
    TH1D *h1_zjt_closure_error_pulldist[ptbinsize-1], *h1_zr_closure_error_pulldist[ptbinsize-1], *h1_jtr_closure_error_pulldist[ptbinsize-1];
    for (int i=1; i<ptbinsize; i++)
    {
        h2_zjt_closure_error[i-1] = (TH2D*)infile->Get(Form("zjt_pt%d_closure_error",i));
        h2_zr_closure_error[i-1] = (TH2D*)infile->Get(Form("zr_pt%d_closure_error",i));
        h2_jtr_closure_error[i-1] = (TH2D*)infile->Get(Form("jtr_pt%d_closure_error",i));

        h2_zjt_closure_error_pulls[i-1] = (TH2D*)h2_zjt_closure_error[i-1]->Clone(Form("zjt_pt%d_closure_error_pulls",i));
        h2_zr_closure_error_pulls[i-1] = (TH2D*)h2_zr_closure_error[i-1]->Clone(Form("zr_pt%d_closure_error_pulls",i));
        h2_jtr_closure_error_pulls[i-1] = (TH2D*)h2_jtr_closure_error[i-1]->Clone(Form("jtr_pt%d_closure_error_pulls",i));

        h1_zjt_closure_error_pulldist[i-1] = new TH1D(Form("zjt_pt%d_closure_error_pulldist", i), ";pulls", 20, -4.5, 4.5);
        h1_zr_closure_error_pulldist[i-1] = new TH1D(Form("zr_pt%d_closure_error_pulldist", i), ";pulls", 20, -4.5, 4.5);
        h1_jtr_closure_error_pulldist[i-1] = new TH1D(Form("jtr_pt%d_closure_error_pulldist", i), ";pulls", 20, -4.5, 4.5);

        GetClosurePulls(h2_zjt_closure_error[i-1], h2_zjt_closure_error_pulls[i-1], h1_zjt_closure_error_pulldist[i-1]);
        GetClosurePulls(h2_zr_closure_error[i-1], h2_zr_closure_error_pulls[i-1], h1_zr_closure_error_pulldist[i-1]);
        GetClosurePulls(h2_jtr_closure_error[i-1], h2_jtr_closure_error_pulls[i-1], h1_jtr_closure_error_pulldist[i-1]);

        h2_zjt_closure_error[i-1]->Write();
        h2_zr_closure_error[i-1]->Write();
        h2_jtr_closure_error[i-1]->Write();

        h2_zjt_closure_error_pulls[i-1]->Write();
        h2_zr_closure_error_pulls[i-1]->Write();
        h2_jtr_closure_error_pulls[i-1]->Write();

        h1_zjt_closure_error_pulldist[i-1]->Write();
        h1_zr_closure_error_pulldist[i-1]->Write();
        h1_jtr_closure_error_pulldist[i-1]->Write();
    }

    outfile->Close();

}
