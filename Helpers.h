#include <vector>
#include <string>
#include <iostream>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TTree.h>
#include <TH1D.h>
#include <TRandom3.h>

double checkphi(double phi);
double folddphi(double dphi);
double checkdR(TLorentzVector v1, TLorentzVector v2);

double checkphi(double phi)
{
  double returnphi = phi;

  if (phi > 1. * TMath::Pi())
    returnphi -= 2 * TMath::Pi();
  else if (phi < -1 * TMath::Pi())
    returnphi += 2 * TMath::Pi();

  return returnphi;
}

double folddphi(double dphi)
{
  double returndphi = dphi;
  if (dphi < 0)
    returndphi += 2. * TMath::Pi();
  if (dphi > 2. * TMath::Pi())
    returndphi -= 2. * TMath::Pi();
  if (dphi > TMath::Pi())
    returndphi = TMath::Pi() - (dphi - TMath::Pi());

  return returndphi;
}
double checkdR(TLorentzVector v1, TLorentzVector v2)
{

  double dr = -99;

  double dphi = checkphi(v1.Phi()) - checkphi(v2.Phi());
  double deta = v1.Rapidity() - v2.Rapidity();
  dphi = checkphi(dphi);

  dr = sqrt(dphi * dphi + deta * deta);

  return dr;
}
void BinLogX(TH1 *h)
{

  TAxis *axis = h->GetXaxis();
  int bins = axis->GetNbins();

  Axis_t from = TMath::Log10(axis->GetXmin());
  Axis_t to = TMath::Log10(axis->GetXmax());
  Axis_t width = (to - from) / bins;
  Axis_t *new_bins = new Axis_t[bins + 1];

  for (int i = 0; i <= bins; i++)
  {
    new_bins[i] = TMath::Power(10, from + i * width);
  }
  axis->Set(bins, new_bins);
  delete[] new_bins;
}

void NormalizeHist(TH1 *h)
{
  if (h == NULL)
    return;
  // h->Scale(1./h->Integral());
  h->Scale(1. / h->Integral(), "width");
}

void NormalizeHistScale(TH1 *h)
{
  if (h == NULL)
    return;
  // h->Scale(1./h->Integral());
  h->Scale(1.0 / h->Integral(), "width");
}

void NormalizeHist(TH2 *h)
{
  if (h == NULL)
    return;
  // h->Scale(1./h->Integral());
  h->Scale(1. / h->Integral(), "width");
}

void NormalizeHist(TH3 *h)
{
  if (h == NULL)
    return;
  // h->Scale(1./h->Integral());
  h->Scale(1. / h->Integral(), "width");
}

double eta_to_theta(double eta1, double eta2)
{
  return std::fabs(-2 * atan(exp(-eta1)) + 2 * atan(exp(-eta2)));
}

template <typename T>
double getAverage(std::vector<T> const &v)
{
  if (v.empty())
  {
    return 0;
  }
  return std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double getLargestAngle(double etaMin, double etaMax)
{
  double smallest_eta;
  if (etaMin < 0 && etaMax > 0)
    smallest_eta = 0;
  else if (std::fabs(etaMin) < std::fabs(etaMax))
    smallest_eta = etaMin;
  else
    smallest_eta = etaMax;
  double angle_high = eta_to_theta(jetradius, smallest_eta);
  return angle_high;
}

void GetPulls(TH2 *hvals, TH2 *hpull, TH1 *h1)
{
  for (int i = 1; i < hvals->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < hvals->GetNbinsY() + 1; j++)
    {
      double err = hvals->GetBinError(i, j);
      double content = hvals->GetBinContent(i, j);
      if (content == 0)
        continue;
      double pull = (content) / err;
      hpull->SetBinContent(i, j, pull);
      h1->Fill(pull);
    }
  }
}
void GetPullsRatio(TH2 *hvals, TH2 *hpull, TH1 *h1)
{
  for (int i = 1; i < hvals->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < hvals->GetNbinsY() + 1; j++)
    {
      double err = hvals->GetBinError(i, j);
      double content = hvals->GetBinContent(i, j);
      if (content == 0)
        continue;
      double pull = (content - 1) / err;
      hpull->SetBinContent(i, j, pull);
      h1->Fill(pull);
    }
  }
}
void GetPullsRatio(TH2 *hvals, TH2 *hpull, TH1 *h1, int xbinsize, int ybinsize)
{
  for (int i = 1; i < xbinsize + 1; i++)
  {
    for (int j = 1; j < ybinsize + 1; j++)
    {
      double err = hvals->GetBinError(i, j);
      double content = hvals->GetBinContent(i, j);
      if (content == 0)
        continue;
      double pull = (content - 1) / err;
      hpull->SetBinContent(i, j, pull);
      h1->Fill(pull);
    }
  }
}

void GetClosurePulls(TH2 *hvals, TH2 *hpull, TH1 *h1)
{
for (int i = 1; i < hvals->GetNbinsX() + 1; i++)
{
    for (int j = 1; j < hvals->GetNbinsY() + 1; j++)
    {
    double err = hvals->GetBinError(i, j);
    double content = hvals->GetBinContent(i, j);
    if (content == 0)
        continue;
    //double pull = (content - 1) / err;
    double pull = content / err;
    hpull->SetBinContent(i, j, pull);
    h1->Fill(pull);
    }
}
}

void GetMSE(TH2 *hunfold, TH2 *htrue, TH2 *hmse)
{
  for (int i = 1; i < hunfold->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < hunfold->GetNbinsY() + 1; j++)
    {
      // double err = hunfold->GetBinError(i, j);
      double content = hunfold->GetBinContent(i, j);
      if (content == 0)
        continue;
      double truecontent = htrue->GetBinContent(i, j);
      double mse = pow(content - truecontent, 2);
      hmse->SetBinContent(i, j, mse);
    }
  }
}
double GetPullPerformance(TH2 *hpull, double sigma)
{
  double percent;
  int num = 0;
  int denom = 0;
  for (int i = 1; i < hpull->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < hpull->GetNbinsY() + 1; j++)
    {
      // double err = hvals->GetBinError(i, j);
      double content = hpull->GetBinContent(i, j);
      if (content == 0)
        continue;
      if (std::fabs(content) < sigma)
        num++;
      denom++;
    }
  }
  percent = (double)num / (double)denom * 100;
  return percent;
}
double GetMSEPerformance(TH2 *hmse)
{
  double avg;
  double sum = 0;
  int denom = 0;
  for (int i = 1; i < hmse->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < hmse->GetNbinsY() + 1; j++)
    {
      // double err = hvals->GetBinError(i, j);
      double content = hmse->GetBinContent(i, j);
      if (content == 0)
        continue;
      sum += content;
      denom++;
    }
  }
  avg = sum / (double)denom;
  return sqrt(avg);
}

void SetDataStyle(TH1 *h, double markersize = 1.)
{
  if (h == NULL)
    return;
  h->SetMarkerStyle(20);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->SetMarkerSize(markersize);
}
void SetRecoStyle(TH1 *h, double markersize = 1.)
{
  if (h == NULL)
    return;
  h->SetMarkerStyle(kFullSquare);
  h->SetMarkerColor(kBlue);
  h->SetLineColor(kBlue);
  h->SetMarkerSize(markersize);
}
void SetTruthStyle(TH1 *h, double markersize = 1.)
{
  if (h == NULL)
    return;
  h->SetMarkerStyle(kFullDiamond);
  h->SetMarkerColor(kAzure-3);
  h->SetLineColor(kAzure-3);
  h->SetMarkerSize(markersize);
}
void SubtractUnity(TH2D *h2)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      // double err = hvals->GetBinError(i, j);
      if (h2->GetBinContent(i, j) != 0.)
      {
        h2->SetBinContent(i, j, h2->GetBinContent(i, j) - 1.);
      }
    }
  }
}

void SubtractUnity(TH1D *h1)
{
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    // double err = hvals->GetBinError(i, j);
    if (h1->GetBinContent(i) != 0.)
    {
      h1->SetBinContent(i, h1->GetBinContent(i) - 1.);
    }
  }
}

double GetWeightedAverage(TH2D *h2)
{
  double sumWeights = 0.;
  double num = 0.;
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      double weight = 1. / h2->GetBinError(i, j);
      double val = std::fabs(h2->GetBinContent(i, j) - 1.);
      if (val == 1.)
        continue;
      num += val * weight;
      sumWeights += weight;
    }
  }
  double result = num / sumWeights;
  return result;
}

double GetWeightedAverage(TH1D *h1)
{
  double sumWeights = 0.;
  double num = 0.;
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    double weight = 1. / h1->GetBinError(i);
    double val = std::fabs(h1->GetBinContent(i) - 1.);
    if (val == 1.)
      continue;
    num += val * weight;
    sumWeights += weight;
  }
  double result = num / sumWeights;
  return result;
}


double GetWeightedAverage(TH2D *h2, TH2D *h2_weight)
{
  double sumWeights = 0.;
  double num = 0.;
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      double val = std::fabs(h2->GetBinContent(i, j));
      double weight = h2_weight->GetBinContent(i, j) / h2_weight->GetBinError(i, j);
      if (val == 100.)
        continue;
      // cout<<"("<<val<<","<<weight<<")";
      num += val * weight;
      sumWeights += weight;
    }
  }
  double result = num / sumWeights;
  return result;
}

void GetErrorHist(TH2 *h2, TH2 *h2_err)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      double err = h2->GetBinError(i, j);
      h2_err->SetBinContent(i, j, err);
    }
  }
}

void SetHistConst(TH2 *h2, double val)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      h2->SetBinContent(i, j, val);
    }
  }
}
void SetHistConst(TH1 *h1, double val)
{
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    h1->SetBinContent(i, val);
  }
}
void SetHistError(TH2 *h2, TH2 *h2_err)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      h2->SetBinError(i, j, h2_err->GetBinContent(i, j) * h2->GetBinContent(i, j));
    }
  }
}
void SetHistError(TH1 *h1, TH1 *h1_err)
{
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    h1->SetBinError(i, h1_err->GetBinContent(i) * h1->GetBinContent(i));
  }
}
void SetHistErrorZero(TH2 *h2)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      h2->SetBinError(i, j, 0.);
    }
  }
}
void SetHistErrorZero(TH1 *h1)
{
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    h1->SetBinError(i, 0.);
  }
}

void SetHistErrCorr(TH2 *hist, TH2 *hist1, TH2 *hist2, double corr = 0.)
{
  for (int i = 1; i < hist->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < hist->GetNbinsY() + 1; j++)
    {
      double err1 = hist1->GetBinError(i, j);
      double err2 = hist2->GetBinError(i, j);
      hist->SetBinError(i, j, sqrt(err1 * err1 + err2 * err2 - 2 * corr * err1 * err2));
    }
  }
}

void SetHistErrCorr(TH1 *hist, TH1 *hist1, TH1 *hist2, double corr = 0.)
{
  for (int i = 1; i < hist->GetNbinsX() + 1; i++)
  {
    double err1 = hist1->GetBinError(i);
    double err2 = hist2->GetBinError(i);
    hist->SetBinError(i, sqrt(err1 * err1 + err2 * err2 - 2 * corr * err1 * err2));

  }
}

void GetSqrtHist(TH2 *h2)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      h2->SetBinContent(i, j, sqrt(h2->GetBinContent(i, j)));
    }
  }
}

void GetSqrtHist(TH1 *h1)
{
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    h1->SetBinContent(i, sqrt(h1->GetBinContent(i)));
  }
}

void FixEffHist(TH2 *h2eff, TH2 *h2purity)
{
  for (int i = 1; i < h2eff->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2eff->GetNbinsY() + 1; j++)
    {
      double pur = h2purity->GetBinContent(i, j);
      double eff = h2eff->GetBinContent(i, j);
      if (eff > 1.)
      {
        h2eff->SetBinContent(i, j, 1.);
        h2purity->SetBinContent(i, j, pur / eff);
      }
    }
  }
}

void SetFeldmanErr(TH2 *h2)
{
  TFeldmanCousins tfc(0.68);
  double low, high;
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      double err = h2->GetBinError(i, j);
      double val = h2->GetBinContent(i, j);
      if (val < 40)
      {
        low = tfc.CalculateLowerLimit(val, 0.);
        high = tfc.CalculateUpperLimit(val, 0.);
        h2->SetBinError(i, j, (low + high) / 2.);
      }
      else
      {
        h2->SetBinError(i, j, err);
      }
    }
  }
}

void PaintOverflow(TH1 *h)
{
  // This function paint the histogram h with an extra bin for overflows

  const char *name = h->GetName();
  const char *title = h->GetTitle();
  Int_t nx = h->GetNbinsX() + 1;
  Double_t x1 = h->GetBinLowEdge(1);
  Double_t bw = h->GetBinWidth(nx);
  Double_t x2 = h->GetBinLowEdge(nx) + bw;

  // Book a temporary histogram having ab extra bin for overflows
  TH1F *htmp = new TH1F(name, title, nx, x1, x2);

  // Fill the new hitogram including the extra bin for overflows
  for (Int_t i = 1; i <= nx; i++)
  {
    htmp->Fill(htmp->GetBinCenter(i), h->GetBinContent(i));
  }

  // Fill the underflows
  htmp->Fill(x1 - 1, h->GetBinContent(0));

  // Restore the number of entries
  htmp->SetEntries(h->GetEntries());

  // Draw the temporary histogram
  htmp->Draw("SAME");
  gPad->SetLogy();
  TText *t = new TText(x2 - bw / 2, h->GetBinContent(nx), "Overflow");
  t->SetTextAngle(90);
  t->SetTextAlign(12);
  t->SetTextSize(0.03);
  ;
  t->Draw("SAME");
}

void SmearObservables(TH3 *h3, TH3 *h3_data, TRandom3 *myRNG)
{
  cout << "Smearing Observables" << endl;
  for (int i = 1; i < h3->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h3->GetNbinsY() + 1; j++)
    {
      for (int k = 1; k < h3->GetNbinsZ() + 1; k++)
      {
        double sigma = h3_data->GetBinError(i, j, k);
        double data_val = h3_data->GetBinContent(i, j, k);
        double rand;
        if (sigma == data_val)
        {
          sigma = sqrt(data_val);
        }
        rand = myRNG->Gaus(1, sigma / data_val);
        rand = (rand < 0) ? 1. : rand;
        double val = h3->GetBinContent(i, j, k);
        if (data_val > 0)
        {
          // cout << sigma << ", " << data_val << ", " << rand << endl;
          h3->SetBinContent(i, j, k, val * rand);
        }
      }
    }
  }
}

void SmearObservables(TH2 *h2, TH2 *h2_data, TRandom3 *myRNG)
{
  cout << "Smearing Observables" << endl;
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {
      double sigma = h2_data->GetBinError(i, j);
      double data_val = h2_data->GetBinContent(i, j);
      double rand;
      if (sigma == data_val)
      {
        sigma = sqrt(data_val);
      }
      rand = myRNG->Gaus(1, sigma / data_val);
      rand = (rand < 0) ? 1. : rand;
      double val = h2->GetBinContent(i, j);
      if (data_val > 0)
      {
        // cout << sigma << ", " << data_val << ", " << rand << endl;
        h2->SetBinContent(i, j, val * rand);
      }
    }
  }
}

void SmearObservables(TH1 *h1, TH1 *h1_data, TRandom3 *myRNG)
{
  cout << "Smearing Observables" << endl;
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {

    double sigma = h1_data->GetBinError(i);
    double data_val = h1_data->GetBinContent(i);
    double rand;
    if (sigma == data_val)
    {
      sigma = sqrt(data_val);
    }
    rand = myRNG->Gaus(1, sigma / data_val);
    rand = (rand < 0) ? 1. : rand;
    double val = h1->GetBinContent(i);
    if (data_val > 0)
    {
      // cout << sigma << ", " << data_val << ", " << rand << endl;
      h1->SetBinContent(i, val * rand);
    }
  }
}

void SmearDtrPt(TLorentzVector &dtr, double charge, double id, TRandom3 *myRNG)
{
  if (std::fabs(charge) > 0)
  {
    double rand = myRNG->Gaus(1, 0.03);
    dtr.SetPxPyPzE(dtr.Px() * rand,
                   dtr.Py() * rand,
                   dtr.Pz() * rand,
                   dtr.E() * rand);
  }
  else
  {
    // dE/E = a/sqrt(E) + b + c/E
    // a = 10%
    // b = 1%
    // c is small
    double b, a, e;
    if (abs(id) == 22 || abs(id) == 111)
    {
      b = 0.01;
      a = 0.1;
    }
    else
    {
      b = 0.09;
      a = 0.67;
    }

    e = dtr.E();
    double sigma = sqrt(b * b + a / sqrt(e) * a / sqrt(e));
    double rand = myRNG->Gaus(1, sigma);
    dtr.SetPxPyPzE(dtr.Px() * rand,
                   dtr.Py() * rand,
                   dtr.Pz() * rand,
                   dtr.E() * rand);
  }
}
bool DropFSP(double prob, TRandom3 *myUniform)
{
  double rand = myUniform->Uniform(1.);
  if (rand < prob)
    return true;
  else
    return false;
}
bool DropFSPHist(TLorentzVector dtr, TH2 *h2, TRandom3 *myUniform)
{
  double rand = myUniform->Uniform(1.);
  int pbin = h2->GetXaxis()->FindBin(dtr.P());
  int etabin = h2->GetYaxis()->FindBin(dtr.Eta());
  double ratio = h2->GetBinContent(pbin, etabin);
  double prob_trk = std::fabs(ratio - 1.);
  if (rand < prob_trk)
    return true;
  else
    return false;
}

void MakeHistPositive(TH3 *h3)
{
  for (int i = 1; i < h3->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h3->GetNbinsY() + 1; j++)
    {
      for (int k = 1; k < h3->GetNbinsZ() + 1; k++)
      {
        if (h3->GetBinContent(i, j, k) < 0)
          h3->SetBinContent(i, j, k, 0.);
      }
    }
  }
}
void MakeHistPositive(TH2 *h2)
{
  for (int i = 1; i < h2->GetNbinsX() + 1; i++)
  {
    for (int j = 1; j < h2->GetNbinsY() + 1; j++)
    {

      if (h2->GetBinContent(i, j) < 0)
        h2->SetBinContent(i, j, 0.);
    }
  }
}
void MakeHistPositive(TH1 *h1)
{
  for (int i = 1; i < h1->GetNbinsX() + 1; i++)
  {
    if (h1->GetBinContent(i) < 0)
      h1->SetBinContent(i, 0.);
  }
}

int SampleTH1(TH1D *hist, TRandom3 *randomGenerator)
{
  const int numBins = hist->GetNbinsX();
  double cumulativeProbabilities[(const int) numBins + 1];

  // Fill the array with cumulative probabilities
  cumulativeProbabilities[0] = 0.0;
  for (int i = 1; i <= numBins; ++i)
  {
    cumulativeProbabilities[i] = cumulativeProbabilities[i - 1] + hist->GetBinContent(i);
  }

  // Generate a random number between 0 and 1
  // TRandom3 randomGenerator;

  double randomValue = randomGenerator->Uniform();

  // Find the bin corresponding to the generated random number
  int selectedBin = 0;
  for (int i = 1; i <= numBins; ++i)
  {
    if (randomValue <= cumulativeProbabilities[i])
    {
      selectedBin = i;
      break;
    }
  }

  // Get the interval corresponding to the selected bin
  double binLowEdge = hist->GetXaxis()->GetBinLowEdge(selectedBin);
  double binWidth = hist->GetXaxis()->GetBinWidth(selectedBin);

  double intervalStart = binLowEdge;
  double intervalEnd = binLowEdge + binWidth;

  // Output the result
  // std::cout << "Generated random number: " << randomValue << std::endl;
  // std::cout << "Selected bin: " << selectedBin << std::endl;
  // std::cout << "Corresponding interval: [" << intervalStart << ", " << intervalEnd << "]" << std::endl;

  return selectedBin;
}

// Returns JES + JER scale factor as a function of jet pT (in GeV)
// Values and function provided by Ezra Lesser (thank you!)
double get_JES_JER(const double jet_pt, TRandom3 *myRNG, bool do_sys = false) {
  double jes_var = 1, jer_var = 1;

  // Low-multiplicity (N_consit < 10) fits
  if (do_sys) {
      if      (jet_pt < 10) { jes_var = 0.904; jer_var = 0.006; }
      else if (jet_pt < 12) { jes_var = 0.936; jer_var = 0.180; }
      else if (jet_pt < 15) { jes_var = 0.924; jer_var = 0.264; }
      else if (jet_pt < 17) { jes_var = 0.972; jer_var = 0.192; }
      else if (jet_pt < 20) { jes_var = 0.970; jer_var = 0.100; }
      else if (jet_pt < 25) { jes_var = 0.976; jer_var = 0.114; }
      else if (jet_pt < 35) { jes_var = 0.984; jer_var = 0.036; }
      else if (jet_pt < 50) { jes_var = 1.000; jer_var = 0.144; }
      else                  { jes_var = 1.008; jer_var = 0.126; }
  }

  // Nominal (inclusive Z+jet) fits
  else {
      if      (jet_pt < 10) { jes_var = 0.912; jer_var = 0.006; }
      else if (jet_pt < 12) { jes_var = 0.948; jer_var = 0.132; }
      else if (jet_pt < 15) { jes_var = 0.944; jer_var = 0.192; }
      else if (jet_pt < 17) { jes_var = 0.976; jer_var = 0.114; }
      else if (jet_pt < 20) { jes_var = 0.972; jer_var = 0.096; }
      else if (jet_pt < 25) { jes_var = 0.976; jer_var = 0.138; }
      else if (jet_pt < 35) { jes_var = 0.984; jer_var = 0.054; }
      else if (jet_pt < 50) { jes_var = 1.000; jer_var = 0.144; }
      else                  { jes_var = 1.008; jer_var = 0.120; }
  }

  return jes_var * myRNG->Gaus(1, jer_var);
}



