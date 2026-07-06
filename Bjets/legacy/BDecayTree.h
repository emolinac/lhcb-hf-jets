//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 16 19:06:05 2024 by ROOT version 6.30/02
// from TChain Decay/DecayTree/
//////////////////////////////////////////////////////////

#ifndef BDecayTree_h
#define BDecayTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BDecayTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxBu_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxBu_OWNPV_COV = 1;
   static constexpr Int_t kMaxBplus0_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxBplus0_OWNPV_COV = 1;
   static constexpr Int_t kMaxJpsi_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxJpsi_OWNPV_COV = 1;
   static constexpr Int_t kMaxJpsi_ORIVX_COV = 1;
   static constexpr Int_t kMaxmup_OWNPV_COV = 1;
   static constexpr Int_t kMaxmup_ORIVX_COV = 1;
   static constexpr Int_t kMaxmum_OWNPV_COV = 1;
   static constexpr Int_t kMaxmum_ORIVX_COV = 1;
   static constexpr Int_t kMaxK_OWNPV_COV = 1;
   static constexpr Int_t kMaxK_ORIVX_COV = 1;

   // Declaration of leaf types
   Double_t        Bu_P;
   Double_t        Bu_PT;
   Double_t        Bu_PE;
   Double_t        Bu_PX;
   Double_t        Bu_PY;
   Double_t        Bu_PZ;
   Double_t        Bu_MM;
   Double_t        Bu_MMERR;
   Double_t        Bu_M;
   Int_t           Bu_ID;
   Bool_t          Bu_L0Global_Dec;
   Bool_t          Bu_L0Global_TIS;
   Bool_t          Bu_L0Global_TOS;
   Bool_t          Bu_Hlt1Global_Dec;
   Bool_t          Bu_Hlt1Global_TIS;
   Bool_t          Bu_Hlt1Global_TOS;
   Bool_t          Bu_Hlt1Phys_Dec;
   Bool_t          Bu_Hlt1Phys_TIS;
   Bool_t          Bu_Hlt1Phys_TOS;
   Bool_t          Bu_Hlt2Global_Dec;
   Bool_t          Bu_Hlt2Global_TIS;
   Bool_t          Bu_Hlt2Global_TOS;
   Bool_t          Bu_Hlt2Phys_Dec;
   Bool_t          Bu_Hlt2Phys_TIS;
   Bool_t          Bu_Hlt2Phys_TOS;
   Bool_t          Bu_L0MuonDecision_Dec;
   Bool_t          Bu_L0MuonDecision_TIS;
   Bool_t          Bu_L0MuonDecision_TOS;
   Bool_t          Bu_L0DiMuonDecision_Dec;
   Bool_t          Bu_L0DiMuonDecision_TIS;
   Bool_t          Bu_L0DiMuonDecision_TOS;
   Bool_t          Bu_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          Bu_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          Bu_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          Bu_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          Bu_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          Bu_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          Bu_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          Bu_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          Bu_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          Bu_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          Bu_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          Bu_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t          Bu_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t          Bu_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t          Bu_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t          Bu_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t          Bu_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t          Bu_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t          Bu_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t          Bu_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t          Bu_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Int_t           Bu_ConsBu_nPV;
   Float_t         Bu_ConsBu_J_psi_1S_M[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_MERR[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_P[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_PERR[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_ctau[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_ctauErr[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_decayLength[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_J_psi_1S_decayLengthErr[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_M[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_MERR[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_P[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_PERR[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_PV_X[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_PV_Y[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_PV_Z[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_PV_key[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_chi2[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_ctau[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_ctauErr[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_decayLength[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_decayLengthErr[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_nDOF[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_nIter[100];   //[Bu_ConsBu_nPV]
   Float_t         Bu_ConsBu_status[100];   //[Bu_ConsBu_nPV]
   Double_t        Bu_MCJet_Px;
   Double_t        Bu_MCJet_Py;
   Double_t        Bu_MCJet_Pz;
   Double_t        Bu_MCJet_E;
   Double_t        Bu_MCJet_PT;
   Double_t        Bu_MCJet_Rapidity;
   Int_t           Bu_MCJet_FLAV;
   Double_t        Bu_MCHF_Px;
   Double_t        Bu_MCHF_Py;
   Double_t        Bu_MCHF_Pz;
   Double_t        Bu_MCHF_E;
   Double_t        Bu_MCHF_PT;
   Double_t        Bu_MCHF_Rapidity;
   Int_t           Bu_MCJet_NumHFHads;
   Int_t           Bu_NumBHads;
   Int_t           Bu_NumDHads;
   Bool_t          Bu_hasc;
   Bool_t          Bu_hascbar;
   Bool_t          Bu_hasb;
   Bool_t          Bu_hasbbar;
   Bool_t          Bu_FromJet;
   Double_t        Bu_ENDVERTEX_X;
   Double_t        Bu_ENDVERTEX_Y;
   Double_t        Bu_ENDVERTEX_Z;
   Double_t        Bu_ENDVERTEX_XERR;
   Double_t        Bu_ENDVERTEX_YERR;
   Double_t        Bu_ENDVERTEX_ZERR;
   Double_t        Bu_ENDVERTEX_CHI2;
   Int_t           Bu_ENDVERTEX_NDOF;
   Float_t         Bu_ENDVERTEX_COV_[3][3];
   Double_t        Bu_OWNPV_X;
   Double_t        Bu_OWNPV_Y;
   Double_t        Bu_OWNPV_Z;
   Double_t        Bu_OWNPV_XERR;
   Double_t        Bu_OWNPV_YERR;
   Double_t        Bu_OWNPV_ZERR;
   Double_t        Bu_OWNPV_CHI2;
   Int_t           Bu_OWNPV_NDOF;
   Float_t         Bu_OWNPV_COV_[3][3];
   Double_t        Bu_IP_OWNPV;
   Double_t        Bu_IPCHI2_OWNPV;
   Double_t        Bu_FD_OWNPV;
   Double_t        Bu_FDCHI2_OWNPV;
   Double_t        Bu_DIRA_OWNPV;
   Int_t           Bu_TRUEID;
   Double_t        Bu_TRUEP_E;
   Double_t        Bu_TRUEP_X;
   Double_t        Bu_TRUEP_Y;
   Double_t        Bu_TRUEP_Z;
   Double_t        Bu_TRUEPT;
   Double_t        Bu_TRUEORIGINVERTEX_X;
   Double_t        Bu_TRUEORIGINVERTEX_Y;
   Double_t        Bu_TRUEORIGINVERTEX_Z;
   Double_t        Bu_TRUEENDVERTEX_X;
   Double_t        Bu_TRUEENDVERTEX_Y;
   Double_t        Bu_TRUEENDVERTEX_Z;
   Bool_t          Bu_TRUEISSTABLE;
   Double_t        Bu_TRUETAU;
   Double_t        Bplus0_ENDVERTEX_X;
   Double_t        Bplus0_ENDVERTEX_Y;
   Double_t        Bplus0_ENDVERTEX_Z;
   Double_t        Bplus0_ENDVERTEX_XERR;
   Double_t        Bplus0_ENDVERTEX_YERR;
   Double_t        Bplus0_ENDVERTEX_ZERR;
   Double_t        Bplus0_ENDVERTEX_CHI2;
   Int_t           Bplus0_ENDVERTEX_NDOF;
   Float_t         Bplus0_ENDVERTEX_COV_[3][3];
   Double_t        Bplus0_OWNPV_X;
   Double_t        Bplus0_OWNPV_Y;
   Double_t        Bplus0_OWNPV_Z;
   Double_t        Bplus0_OWNPV_XERR;
   Double_t        Bplus0_OWNPV_YERR;
   Double_t        Bplus0_OWNPV_ZERR;
   Double_t        Bplus0_OWNPV_CHI2;
   Int_t           Bplus0_OWNPV_NDOF;
   Float_t         Bplus0_OWNPV_COV_[3][3];
   Double_t        Bplus0_IP_OWNPV;
   Double_t        Bplus0_IPCHI2_OWNPV;
   Double_t        Bplus0_FD_OWNPV;
   Double_t        Bplus0_FDCHI2_OWNPV;
   Double_t        Bplus0_DIRA_OWNPV;
   Double_t        Bplus0_P;
   Double_t        Bplus0_PT;
   Double_t        Bplus0_PE;
   Double_t        Bplus0_PX;
   Double_t        Bplus0_PY;
   Double_t        Bplus0_PZ;
   Double_t        Bplus0_MM;
   Double_t        Bplus0_MMERR;
   Double_t        Bplus0_M;
   Int_t           Bplus0_ID;
   Bool_t          Bplus0_L0Global_Dec;
   Bool_t          Bplus0_L0Global_TIS;
   Bool_t          Bplus0_L0Global_TOS;
   Bool_t          Bplus0_Hlt1Global_Dec;
   Bool_t          Bplus0_Hlt1Global_TIS;
   Bool_t          Bplus0_Hlt1Global_TOS;
   Bool_t          Bplus0_Hlt1Phys_Dec;
   Bool_t          Bplus0_Hlt1Phys_TIS;
   Bool_t          Bplus0_Hlt1Phys_TOS;
   Bool_t          Bplus0_Hlt2Global_Dec;
   Bool_t          Bplus0_Hlt2Global_TIS;
   Bool_t          Bplus0_Hlt2Global_TOS;
   Bool_t          Bplus0_Hlt2Phys_Dec;
   Bool_t          Bplus0_Hlt2Phys_TIS;
   Bool_t          Bplus0_Hlt2Phys_TOS;
   Bool_t          Bplus0_L0MuonDecision_Dec;
   Bool_t          Bplus0_L0MuonDecision_TIS;
   Bool_t          Bplus0_L0MuonDecision_TOS;
   Bool_t          Bplus0_L0DiMuonDecision_Dec;
   Bool_t          Bplus0_L0DiMuonDecision_TIS;
   Bool_t          Bplus0_L0DiMuonDecision_TOS;
   Bool_t          Bplus0_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          Bplus0_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          Bplus0_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          Bplus0_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          Bplus0_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          Bplus0_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t          Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t          Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t          Bplus0_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t          Bplus0_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t          Bplus0_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t          Bplus0_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t          Bplus0_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t          Bplus0_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Double_t        Jpsi_P;
   Double_t        Jpsi_PT;
   Double_t        Jpsi_PE;
   Double_t        Jpsi_PX;
   Double_t        Jpsi_PY;
   Double_t        Jpsi_PZ;
   Double_t        Jpsi_MM;
   Double_t        Jpsi_MMERR;
   Double_t        Jpsi_M;
   Int_t           Jpsi_ID;
   Bool_t          Jpsi_L0Global_Dec;
   Bool_t          Jpsi_L0Global_TIS;
   Bool_t          Jpsi_L0Global_TOS;
   Bool_t          Jpsi_Hlt1Global_Dec;
   Bool_t          Jpsi_Hlt1Global_TIS;
   Bool_t          Jpsi_Hlt1Global_TOS;
   Bool_t          Jpsi_Hlt1Phys_Dec;
   Bool_t          Jpsi_Hlt1Phys_TIS;
   Bool_t          Jpsi_Hlt1Phys_TOS;
   Bool_t          Jpsi_Hlt2Global_Dec;
   Bool_t          Jpsi_Hlt2Global_TIS;
   Bool_t          Jpsi_Hlt2Global_TOS;
   Bool_t          Jpsi_Hlt2Phys_Dec;
   Bool_t          Jpsi_Hlt2Phys_TIS;
   Bool_t          Jpsi_Hlt2Phys_TOS;
   Bool_t          Jpsi_L0MuonDecision_Dec;
   Bool_t          Jpsi_L0MuonDecision_TIS;
   Bool_t          Jpsi_L0MuonDecision_TOS;
   Bool_t          Jpsi_L0DiMuonDecision_Dec;
   Bool_t          Jpsi_L0DiMuonDecision_TIS;
   Bool_t          Jpsi_L0DiMuonDecision_TOS;
   Bool_t          Jpsi_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          Jpsi_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          Jpsi_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Double_t        Jpsi_ENDVERTEX_X;
   Double_t        Jpsi_ENDVERTEX_Y;
   Double_t        Jpsi_ENDVERTEX_Z;
   Double_t        Jpsi_ENDVERTEX_XERR;
   Double_t        Jpsi_ENDVERTEX_YERR;
   Double_t        Jpsi_ENDVERTEX_ZERR;
   Double_t        Jpsi_ENDVERTEX_CHI2;
   Int_t           Jpsi_ENDVERTEX_NDOF;
   Float_t         Jpsi_ENDVERTEX_COV_[3][3];
   Double_t        Jpsi_OWNPV_X;
   Double_t        Jpsi_OWNPV_Y;
   Double_t        Jpsi_OWNPV_Z;
   Double_t        Jpsi_OWNPV_XERR;
   Double_t        Jpsi_OWNPV_YERR;
   Double_t        Jpsi_OWNPV_ZERR;
   Double_t        Jpsi_OWNPV_CHI2;
   Int_t           Jpsi_OWNPV_NDOF;
   Float_t         Jpsi_OWNPV_COV_[3][3];
   Double_t        Jpsi_IP_OWNPV;
   Double_t        Jpsi_IPCHI2_OWNPV;
   Double_t        Jpsi_FD_OWNPV;
   Double_t        Jpsi_FDCHI2_OWNPV;
   Double_t        Jpsi_DIRA_OWNPV;
   Double_t        Jpsi_ORIVX_X;
   Double_t        Jpsi_ORIVX_Y;
   Double_t        Jpsi_ORIVX_Z;
   Double_t        Jpsi_ORIVX_XERR;
   Double_t        Jpsi_ORIVX_YERR;
   Double_t        Jpsi_ORIVX_ZERR;
   Double_t        Jpsi_ORIVX_CHI2;
   Int_t           Jpsi_ORIVX_NDOF;
   Float_t         Jpsi_ORIVX_COV_[3][3];
   Double_t        Jpsi_FD_ORIVX;
   Double_t        Jpsi_FDCHI2_ORIVX;
   Double_t        Jpsi_DIRA_ORIVX;
   Int_t           Jpsi_TRUEID;
   Double_t        Jpsi_TRUEP_E;
   Double_t        Jpsi_TRUEP_X;
   Double_t        Jpsi_TRUEP_Y;
   Double_t        Jpsi_TRUEP_Z;
   Double_t        Jpsi_TRUEPT;
   Double_t        Jpsi_TRUEORIGINVERTEX_X;
   Double_t        Jpsi_TRUEORIGINVERTEX_Y;
   Double_t        Jpsi_TRUEORIGINVERTEX_Z;
   Double_t        Jpsi_TRUEENDVERTEX_X;
   Double_t        Jpsi_TRUEENDVERTEX_Y;
   Double_t        Jpsi_TRUEENDVERTEX_Z;
   Bool_t          Jpsi_TRUEISSTABLE;
   Double_t        Jpsi_TRUETAU;
    
   Double_t        mup_MC12TuneV2_ProbNNe;
   Double_t        mup_MC12TuneV2_ProbNNmu;
   Double_t        mup_MC12TuneV2_ProbNNpi;
   Double_t        mup_MC12TuneV2_ProbNNk;
   Double_t        mup_MC12TuneV2_ProbNNp;
   Double_t        mup_MC12TuneV2_ProbNNghost;
   Double_t        mup_MC12TuneV3_ProbNNe;
   Double_t        mup_MC12TuneV3_ProbNNmu;
   Double_t        mup_MC12TuneV3_ProbNNpi;
   Double_t        mup_MC12TuneV3_ProbNNk;
   Double_t        mup_MC12TuneV3_ProbNNp;
   Double_t        mup_MC12TuneV3_ProbNNghost;
   Double_t        mup_MC12TuneV4_ProbNNe;
   Double_t        mup_MC12TuneV4_ProbNNmu;
   Double_t        mup_MC12TuneV4_ProbNNpi;
   Double_t        mup_MC12TuneV4_ProbNNk;
   Double_t        mup_MC12TuneV4_ProbNNp;
   Double_t        mup_MC12TuneV4_ProbNNghost;
   Double_t        mup_MC15TuneV1_ProbNNe;
   Double_t        mup_MC15TuneV1_ProbNNmu;
   Double_t        mup_MC15TuneV1_ProbNNpi;
   Double_t        mup_MC15TuneV1_ProbNNk;
   Double_t        mup_MC15TuneV1_ProbNNp;
   Double_t        mup_MC15TuneV1_ProbNNghost;
   Double_t        mup_P;
   Double_t        mup_PT;
   Double_t        mup_PE;
   Double_t        mup_PX;
   Double_t        mup_PY;
   Double_t        mup_PZ;
   Double_t        mup_M;
   Int_t           mup_ID;
   Double_t        mup_PIDe;
   Double_t        mup_PIDmu;
   Double_t        mup_PIDK;
   Double_t        mup_PIDp;
   Double_t        mup_PIDd;
   Double_t        mup_ProbNNe;
   Double_t        mup_ProbNNk;
   Double_t        mup_ProbNNp;
   Double_t        mup_ProbNNpi;
   Double_t        mup_ProbNNmu;
   Double_t        mup_ProbNNd;
   Double_t        mup_ProbNNghost;
   Bool_t          mup_hasMuon;
   Bool_t          mup_isMuon;
   Bool_t          mup_hasRich;
   Bool_t          mup_UsedRichAerogel;
   Bool_t          mup_UsedRich1Gas;
   Bool_t          mup_UsedRich2Gas;
   Bool_t          mup_RichAboveElThres;
   Bool_t          mup_RichAboveMuThres;
   Bool_t          mup_RichAbovePiThres;
   Bool_t          mup_RichAboveKaThres;
   Bool_t          mup_RichAbovePrThres;
   Bool_t          mup_hasCalo;
   Bool_t          mup_L0Global_Dec;
   Bool_t          mup_L0Global_TIS;
   Bool_t          mup_L0Global_TOS;
   Bool_t          mup_Hlt1Global_Dec;
   Bool_t          mup_Hlt1Global_TIS;
   Bool_t          mup_Hlt1Global_TOS;
   Bool_t          mup_Hlt1Phys_Dec;
   Bool_t          mup_Hlt1Phys_TIS;
   Bool_t          mup_Hlt1Phys_TOS;
   Bool_t          mup_Hlt2Global_Dec;
   Bool_t          mup_Hlt2Global_TIS;
   Bool_t          mup_Hlt2Global_TOS;
   Bool_t          mup_Hlt2Phys_Dec;
   Bool_t          mup_Hlt2Phys_TIS;
   Bool_t          mup_Hlt2Phys_TOS;
   Bool_t          mup_L0MuonDecision_Dec;
   Bool_t          mup_L0MuonDecision_TIS;
   Bool_t          mup_L0MuonDecision_TOS;
   Bool_t          mup_L0DiMuonDecision_Dec;
   Bool_t          mup_L0DiMuonDecision_TIS;
   Bool_t          mup_L0DiMuonDecision_TOS;
   Bool_t          mup_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          mup_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          mup_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          mup_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          mup_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          mup_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          mup_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          mup_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          mup_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          mup_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          mup_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          mup_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t          mup_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t          mup_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t          mup_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t          mup_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t          mup_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t          mup_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t          mup_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t          mup_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t          mup_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Double_t        mup_OWNPV_X;
   Double_t        mup_OWNPV_Y;
   Double_t        mup_OWNPV_Z;
   Double_t        mup_OWNPV_XERR;
   Double_t        mup_OWNPV_YERR;
   Double_t        mup_OWNPV_ZERR;
   Double_t        mup_OWNPV_CHI2;
   Int_t           mup_OWNPV_NDOF;
   Float_t         mup_OWNPV_COV_[3][3];
   Double_t        mup_IP_OWNPV;
   Double_t        mup_IPCHI2_OWNPV;
   Double_t        mup_ORIVX_X;
   Double_t        mup_ORIVX_Y;
   Double_t        mup_ORIVX_Z;
   Double_t        mup_ORIVX_XERR;
   Double_t        mup_ORIVX_YERR;
   Double_t        mup_ORIVX_ZERR;
   Double_t        mup_ORIVX_CHI2;
   Int_t           mup_ORIVX_NDOF;
   Float_t         mup_ORIVX_COV_[3][3];
   Int_t           mup_TRUEID;
   Double_t        mup_TRUEP_E;
   Double_t        mup_TRUEP_X;
   Double_t        mup_TRUEP_Y;
   Double_t        mup_TRUEP_Z;
   Double_t        mup_TRUEPT;
   Double_t        mup_TRUEORIGINVERTEX_X;
   Double_t        mup_TRUEORIGINVERTEX_Y;
   Double_t        mup_TRUEORIGINVERTEX_Z;
   Double_t        mup_TRUEENDVERTEX_X;
   Double_t        mup_TRUEENDVERTEX_Y;
   Double_t        mup_TRUEENDVERTEX_Z;
   Bool_t          mup_TRUEISSTABLE;
   Double_t        mup_TRUETAU;
   Int_t           mup_TRACK_Type;
   Int_t           mup_TRACK_Key;
   Double_t        mup_TRACK_CHI2NDOF;
   Double_t        mup_TRACK_PCHI2;
   Double_t        mup_TRACK_MatchCHI2;
   Double_t        mup_TRACK_GhostProb;
   Double_t        mup_TRACK_CloneDist;
   Double_t        mup_TRACK_Likelihood;
    
   Double_t        mum_MC12TuneV2_ProbNNe;
   Double_t        mum_MC12TuneV2_ProbNNmu;
   Double_t        mum_MC12TuneV2_ProbNNpi;
   Double_t        mum_MC12TuneV2_ProbNNk;
   Double_t        mum_MC12TuneV2_ProbNNp;
   Double_t        mum_MC12TuneV2_ProbNNghost;
   Double_t        mum_MC12TuneV3_ProbNNe;
   Double_t        mum_MC12TuneV3_ProbNNmu;
   Double_t        mum_MC12TuneV3_ProbNNpi;
   Double_t        mum_MC12TuneV3_ProbNNk;
   Double_t        mum_MC12TuneV3_ProbNNp;
   Double_t        mum_MC12TuneV3_ProbNNghost;
   Double_t        mum_MC12TuneV4_ProbNNe;
   Double_t        mum_MC12TuneV4_ProbNNmu;
   Double_t        mum_MC12TuneV4_ProbNNpi;
   Double_t        mum_MC12TuneV4_ProbNNk;
   Double_t        mum_MC12TuneV4_ProbNNp;
   Double_t        mum_MC12TuneV4_ProbNNghost;
   Double_t        mum_MC15TuneV1_ProbNNe;
   Double_t        mum_MC15TuneV1_ProbNNmu;
   Double_t        mum_MC15TuneV1_ProbNNpi;
   Double_t        mum_MC15TuneV1_ProbNNk;
   Double_t        mum_MC15TuneV1_ProbNNp;
   Double_t        mum_MC15TuneV1_ProbNNghost;
   Double_t        mum_P;
   Double_t        mum_PT;
   Double_t        mum_PE;
   Double_t        mum_PX;
   Double_t        mum_PY;
   Double_t        mum_PZ;
   Double_t        mum_M;
   Int_t           mum_ID;
   Double_t        mum_PIDe;
   Double_t        mum_PIDmu;
   Double_t        mum_PIDK;
   Double_t        mum_PIDp;
   Double_t        mum_PIDd;
   Double_t        mum_ProbNNe;
   Double_t        mum_ProbNNk;
   Double_t        mum_ProbNNp;
   Double_t        mum_ProbNNpi;
   Double_t        mum_ProbNNmu;
   Double_t        mum_ProbNNd;
   Double_t        mum_ProbNNghost;
   Bool_t          mum_hasMuon;
   Bool_t          mum_isMuon;
   Bool_t          mum_hasRich;
   Bool_t          mum_UsedRichAerogel;
   Bool_t          mum_UsedRich1Gas;
   Bool_t          mum_UsedRich2Gas;
   Bool_t          mum_RichAboveElThres;
   Bool_t          mum_RichAboveMuThres;
   Bool_t          mum_RichAbovePiThres;
   Bool_t          mum_RichAboveKaThres;
   Bool_t          mum_RichAbovePrThres;
   Bool_t          mum_hasCalo;
   Bool_t          mum_L0Global_Dec;
   Bool_t          mum_L0Global_TIS;
   Bool_t          mum_L0Global_TOS;
   Bool_t          mum_Hlt1Global_Dec;
   Bool_t          mum_Hlt1Global_TIS;
   Bool_t          mum_Hlt1Global_TOS;
   Bool_t          mum_Hlt1Phys_Dec;
   Bool_t          mum_Hlt1Phys_TIS;
   Bool_t          mum_Hlt1Phys_TOS;
   Bool_t          mum_Hlt2Global_Dec;
   Bool_t          mum_Hlt2Global_TIS;
   Bool_t          mum_Hlt2Global_TOS;
   Bool_t          mum_Hlt2Phys_Dec;
   Bool_t          mum_Hlt2Phys_TIS;
   Bool_t          mum_Hlt2Phys_TOS;
   Bool_t          mum_L0MuonDecision_Dec;
   Bool_t          mum_L0MuonDecision_TIS;
   Bool_t          mum_L0MuonDecision_TOS;
   Bool_t          mum_L0DiMuonDecision_Dec;
   Bool_t          mum_L0DiMuonDecision_TIS;
   Bool_t          mum_L0DiMuonDecision_TOS;
   Bool_t          mum_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          mum_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          mum_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          mum_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          mum_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          mum_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          mum_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          mum_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          mum_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          mum_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          mum_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          mum_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t          mum_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t          mum_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t          mum_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t          mum_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t          mum_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t          mum_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t          mum_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t          mum_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t          mum_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Double_t        mum_OWNPV_X;
   Double_t        mum_OWNPV_Y;
   Double_t        mum_OWNPV_Z;
   Double_t        mum_OWNPV_XERR;
   Double_t        mum_OWNPV_YERR;
   Double_t        mum_OWNPV_ZERR;
   Double_t        mum_OWNPV_CHI2;
   Int_t           mum_OWNPV_NDOF;
   Float_t         mum_OWNPV_COV_[3][3];
   Double_t        mum_IP_OWNPV;
   Double_t        mum_IPCHI2_OWNPV;
   Double_t        mum_ORIVX_X;
   Double_t        mum_ORIVX_Y;
   Double_t        mum_ORIVX_Z;
   Double_t        mum_ORIVX_XERR;
   Double_t        mum_ORIVX_YERR;
   Double_t        mum_ORIVX_ZERR;
   Double_t        mum_ORIVX_CHI2;
   Int_t           mum_ORIVX_NDOF;
   Float_t         mum_ORIVX_COV_[3][3];
   Int_t           mum_TRUEID;
   Double_t        mum_TRUEP_E;
   Double_t        mum_TRUEP_X;
   Double_t        mum_TRUEP_Y;
   Double_t        mum_TRUEP_Z;
   Double_t        mum_TRUEPT;
   Double_t        mum_TRUEORIGINVERTEX_X;
   Double_t        mum_TRUEORIGINVERTEX_Y;
   Double_t        mum_TRUEORIGINVERTEX_Z;
   Double_t        mum_TRUEENDVERTEX_X;
   Double_t        mum_TRUEENDVERTEX_Y;
   Double_t        mum_TRUEENDVERTEX_Z;
   Bool_t          mum_TRUEISSTABLE;
   Double_t        mum_TRUETAU;
   Int_t           mum_TRACK_Type;
   Int_t           mum_TRACK_Key;
   Double_t        mum_TRACK_CHI2NDOF;
   Double_t        mum_TRACK_PCHI2;
   Double_t        mum_TRACK_MatchCHI2;
   Double_t        mum_TRACK_GhostProb;
   Double_t        mum_TRACK_CloneDist;
   Double_t        mum_TRACK_Likelihood;
   Double_t        K_MC12TuneV2_ProbNNe;
   Double_t        K_MC12TuneV2_ProbNNmu;
   Double_t        K_MC12TuneV2_ProbNNpi;
   Double_t        K_MC12TuneV2_ProbNNk;
   Double_t        K_MC12TuneV2_ProbNNp;
   Double_t        K_MC12TuneV2_ProbNNghost;
   Double_t        K_MC12TuneV3_ProbNNe;
   Double_t        K_MC12TuneV3_ProbNNmu;
   Double_t        K_MC12TuneV3_ProbNNpi;
   Double_t        K_MC12TuneV3_ProbNNk;
   Double_t        K_MC12TuneV3_ProbNNp;
   Double_t        K_MC12TuneV3_ProbNNghost;
   Double_t        K_MC12TuneV4_ProbNNe;
   Double_t        K_MC12TuneV4_ProbNNmu;
   Double_t        K_MC12TuneV4_ProbNNpi;
   Double_t        K_MC12TuneV4_ProbNNk;
   Double_t        K_MC12TuneV4_ProbNNp;
   Double_t        K_MC12TuneV4_ProbNNghost;
   Double_t        K_MC15TuneV1_ProbNNe;
   Double_t        K_MC15TuneV1_ProbNNmu;
   Double_t        K_MC15TuneV1_ProbNNpi;
   Double_t        K_MC15TuneV1_ProbNNk;
   Double_t        K_MC15TuneV1_ProbNNp;
   Double_t        K_MC15TuneV1_ProbNNghost;
   Double_t        K_P;
   Double_t        K_PT;
   Double_t        K_PE;
   Double_t        K_PX;
   Double_t        K_PY;
   Double_t        K_PZ;
   Double_t        K_M;
   Int_t           K_ID;
   Double_t        K_PIDe;
   Double_t        K_PIDmu;
   Double_t        K_PIDK;
   Double_t        K_PIDp;
   Double_t        K_PIDd;
   Double_t        K_ProbNNe;
   Double_t        K_ProbNNk;
   Double_t        K_ProbNNp;
   Double_t        K_ProbNNpi;
   Double_t        K_ProbNNmu;
   Double_t        K_ProbNNd;
   Double_t        K_ProbNNghost;
   Bool_t          K_hasMuon;
   Bool_t          K_isMuon;
   Bool_t          K_hasRich;
   Bool_t          K_UsedRichAerogel;
   Bool_t          K_UsedRich1Gas;
   Bool_t          K_UsedRich2Gas;
   Bool_t          K_RichAboveElThres;
   Bool_t          K_RichAboveMuThres;
   Bool_t          K_RichAbovePiThres;
   Bool_t          K_RichAboveKaThres;
   Bool_t          K_RichAbovePrThres;
   Bool_t          K_hasCalo;
   Bool_t          K_L0Global_Dec;
   Bool_t          K_L0Global_TIS;
   Bool_t          K_L0Global_TOS;
   Bool_t          K_Hlt1Global_Dec;
   Bool_t          K_Hlt1Global_TIS;
   Bool_t          K_Hlt1Global_TOS;
   Bool_t          K_Hlt1Phys_Dec;
   Bool_t          K_Hlt1Phys_TIS;
   Bool_t          K_Hlt1Phys_TOS;
   Bool_t          K_Hlt2Global_Dec;
   Bool_t          K_Hlt2Global_TIS;
   Bool_t          K_Hlt2Global_TOS;
   Bool_t          K_Hlt2Phys_Dec;
   Bool_t          K_Hlt2Phys_TIS;
   Bool_t          K_Hlt2Phys_TOS;
   Bool_t          K_L0MuonDecision_Dec;
   Bool_t          K_L0MuonDecision_TIS;
   Bool_t          K_L0MuonDecision_TOS;
   Bool_t          K_L0DiMuonDecision_Dec;
   Bool_t          K_L0DiMuonDecision_TIS;
   Bool_t          K_L0DiMuonDecision_TOS;
   Bool_t          K_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          K_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          K_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          K_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t          K_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t          K_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t          K_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          K_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          K_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          K_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          K_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          K_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t          K_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t          K_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t          K_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t          K_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t          K_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t          K_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t          K_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t          K_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t          K_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Double_t        K_OWNPV_X;
   Double_t        K_OWNPV_Y;
   Double_t        K_OWNPV_Z;
   Double_t        K_OWNPV_XERR;
   Double_t        K_OWNPV_YERR;
   Double_t        K_OWNPV_ZERR;
   Double_t        K_OWNPV_CHI2;
   Int_t           K_OWNPV_NDOF;
   Float_t         K_OWNPV_COV_[3][3];
   Double_t        K_IP_OWNPV;
   Double_t        K_IPCHI2_OWNPV;
   Double_t        K_ORIVX_X;
   Double_t        K_ORIVX_Y;
   Double_t        K_ORIVX_Z;
   Double_t        K_ORIVX_XERR;
   Double_t        K_ORIVX_YERR;
   Double_t        K_ORIVX_ZERR;
   Double_t        K_ORIVX_CHI2;
   Int_t           K_ORIVX_NDOF;
   Float_t         K_ORIVX_COV_[3][3];
   Int_t           K_TRUEID;
   Double_t        K_TRUEP_E;
   Double_t        K_TRUEP_X;
   Double_t        K_TRUEP_Y;
   Double_t        K_TRUEP_Z;
   Double_t        K_TRUEPT;
   Double_t        K_TRUEORIGINVERTEX_X;
   Double_t        K_TRUEORIGINVERTEX_Y;
   Double_t        K_TRUEORIGINVERTEX_Z;
   Double_t        K_TRUEENDVERTEX_X;
   Double_t        K_TRUEENDVERTEX_Y;
   Double_t        K_TRUEENDVERTEX_Z;
   Bool_t          K_TRUEISSTABLE;
   Double_t        K_TRUETAU;
   Int_t           K_TRACK_Type;
   Int_t           K_TRACK_Key;
   Double_t        K_TRACK_CHI2NDOF;
   Double_t        K_TRACK_PCHI2;
   Double_t        K_TRACK_MatchCHI2;
   Double_t        K_TRACK_GhostProb;
   Double_t        K_TRACK_CloneDist;
   Double_t        K_TRACK_Likelihood;
   UInt_t          nCandidate;
   ULong64_t       totCandidates;
   ULong64_t       EventInSequence;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          BCID;
   Int_t           BCType;
   UInt_t          OdinTCK;
   UInt_t          L0DUTCK;
   UInt_t          HLT1TCK;
   UInt_t          HLT2TCK;
   ULong64_t       GpsTime;
   Short_t         Polarity;
   Int_t           nPVs;
   Int_t           nTracks;
   Int_t           nLongTracks;
   Int_t           nDownstreamTracks;
   Int_t           nUpstreamTracks;
   Int_t           nVeloTracks;
   Int_t           nTTracks;
   Int_t           nBackTracks;
   Int_t           nRich1Hits;
   Int_t           nRich2Hits;
   Int_t           nVeloClusters;
   Int_t           nITClusters;
   Int_t           nTTClusters;
   Int_t           nOTClusters;
   Int_t           nSPDHits;
   Int_t           nMuonCoordsS0;
   Int_t           nMuonCoordsS1;
   Int_t           nMuonCoordsS2;
   Int_t           nMuonCoordsS3;
   Int_t           nMuonCoordsS4;
   Int_t           nMuonTracks;

   // List of branches
   TBranch        *b_Bu_P;   //!
   TBranch        *b_Bu_PT;   //!
   TBranch        *b_Bu_PE;   //!
   TBranch        *b_Bu_PX;   //!
   TBranch        *b_Bu_PY;   //!
   TBranch        *b_Bu_PZ;   //!
   TBranch        *b_Bu_MM;   //!
   TBranch        *b_Bu_MMERR;   //!
   TBranch        *b_Bu_M;   //!
   TBranch        *b_Bu_ID;   //!
   TBranch        *b_Bu_L0Global_Dec;   //!
   TBranch        *b_Bu_L0Global_TIS;   //!
   TBranch        *b_Bu_L0Global_TOS;   //!
   TBranch        *b_Bu_Hlt1Global_Dec;   //!
   TBranch        *b_Bu_Hlt1Global_TIS;   //!
   TBranch        *b_Bu_Hlt1Global_TOS;   //!
   TBranch        *b_Bu_Hlt1Phys_Dec;   //!
   TBranch        *b_Bu_Hlt1Phys_TIS;   //!
   TBranch        *b_Bu_Hlt1Phys_TOS;   //!
   TBranch        *b_Bu_Hlt2Global_Dec;   //!
   TBranch        *b_Bu_Hlt2Global_TIS;   //!
   TBranch        *b_Bu_Hlt2Global_TOS;   //!
   TBranch        *b_Bu_Hlt2Phys_Dec;   //!
   TBranch        *b_Bu_Hlt2Phys_TIS;   //!
   TBranch        *b_Bu_Hlt2Phys_TOS;   //!
   TBranch        *b_Bu_L0MuonDecision_Dec;   //!
   TBranch        *b_Bu_L0MuonDecision_TIS;   //!
   TBranch        *b_Bu_L0MuonDecision_TOS;   //!
   TBranch        *b_Bu_L0DiMuonDecision_Dec;   //!
   TBranch        *b_Bu_L0DiMuonDecision_TIS;   //!
   TBranch        *b_Bu_L0DiMuonDecision_TOS;   //!
   TBranch        *b_Bu_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_Bu_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_Bu_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedPsi2SDecision_Dec;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedPsi2SDecision_TIS;   //!
   TBranch        *b_Bu_Hlt2DiMuonDetachedPsi2SDecision_TOS;   //!
   TBranch        *b_Bu_Hlt2DiMuonJPsiDecision_Dec;   //!
   TBranch        *b_Bu_Hlt2DiMuonJPsiDecision_TIS;   //!
   TBranch        *b_Bu_Hlt2DiMuonJPsiDecision_TOS;   //!
   TBranch        *b_Bu_Hlt2DiMuonJPsiHighPTDecision_Dec;   //!
   TBranch        *b_Bu_Hlt2DiMuonJPsiHighPTDecision_TIS;   //!
   TBranch        *b_Bu_Hlt2DiMuonJPsiHighPTDecision_TOS;   //!
   TBranch        *b_Bu_ConsBu_nPV;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_M;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_MERR;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_P;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_PERR;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_ctau;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_ctauErr;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_decayLength;   //!
   TBranch        *b_Bu_ConsBu_J_psi_1S_decayLengthErr;   //!
   TBranch        *b_Bu_ConsBu_M;   //!
   TBranch        *b_Bu_ConsBu_MERR;   //!
   TBranch        *b_Bu_ConsBu_P;   //!
   TBranch        *b_Bu_ConsBu_PERR;   //!
   TBranch        *b_Bu_ConsBu_PV_X;   //!
   TBranch        *b_Bu_ConsBu_PV_Y;   //!
   TBranch        *b_Bu_ConsBu_PV_Z;   //!
   TBranch        *b_Bu_ConsBu_PV_key;   //!
   TBranch        *b_Bu_ConsBu_chi2;   //!
   TBranch        *b_Bu_ConsBu_ctau;   //!
   TBranch        *b_Bu_ConsBu_ctauErr;   //!
   TBranch        *b_Bu_ConsBu_decayLength;   //!
   TBranch        *b_Bu_ConsBu_decayLengthErr;   //!
   TBranch        *b_Bu_ConsBu_nDOF;   //!
   TBranch        *b_Bu_ConsBu_nIter;   //!
   TBranch        *b_Bu_ConsBu_status;   //!
   TBranch        *b_Bu_MCJet_Px;   //!
   TBranch        *b_Bu_MCJet_Py;   //!
   TBranch        *b_Bu_MCJet_Pz;   //!
   TBranch        *b_Bu_MCJet_E;   //!
   TBranch        *b_Bu_MCJet_PT;   //!
   TBranch        *b_Bu_MCJet_Rapidity;   //!
   TBranch        *b_Bu_MCJet_FLAV;   //!
   TBranch        *b_Bu_MCHF_Px;   //!
   TBranch        *b_Bu_MCHF_Py;   //!
   TBranch        *b_Bu_MCHF_Pz;   //!
   TBranch        *b_Bu_MCHF_E;   //!
   TBranch        *b_Bu_MCHF_PT;   //!
   TBranch        *b_Bu_MCHF_Rapidity;   //!
   TBranch        *b_Bu_MCJet_NumHFHads;   //!
   TBranch        *b_Bu_NumBHads;   //!
   TBranch        *b_Bu_NumDHads;   //!
   TBranch        *b_Bu_hasc;   //!
   TBranch        *b_Bu_hascbar;   //!
   TBranch        *b_Bu_hasb;   //!
   TBranch        *b_Bu_hasbbar;   //!
   TBranch        *b_Bu_FromJet;   //!
   TBranch        *b_Bu_ENDVERTEX_X;   //!
   TBranch        *b_Bu_ENDVERTEX_Y;   //!
   TBranch        *b_Bu_ENDVERTEX_Z;   //!
   TBranch        *b_Bu_ENDVERTEX_XERR;   //!
   TBranch        *b_Bu_ENDVERTEX_YERR;   //!
   TBranch        *b_Bu_ENDVERTEX_ZERR;   //!
   TBranch        *b_Bu_ENDVERTEX_CHI2;   //!
   TBranch        *b_Bu_ENDVERTEX_NDOF;   //!
   TBranch        *b_Bu_ENDVERTEX_COV_;   //!
   TBranch        *b_Bu_OWNPV_X;   //!
   TBranch        *b_Bu_OWNPV_Y;   //!
   TBranch        *b_Bu_OWNPV_Z;   //!
   TBranch        *b_Bu_OWNPV_XERR;   //!
   TBranch        *b_Bu_OWNPV_YERR;   //!
   TBranch        *b_Bu_OWNPV_ZERR;   //!
   TBranch        *b_Bu_OWNPV_CHI2;   //!
   TBranch        *b_Bu_OWNPV_NDOF;   //!
   TBranch        *b_Bu_OWNPV_COV_;   //!
   TBranch        *b_Bu_IP_OWNPV;   //!
   TBranch        *b_Bu_IPCHI2_OWNPV;   //!
   TBranch        *b_Bu_FD_OWNPV;   //!
   TBranch        *b_Bu_FDCHI2_OWNPV;   //!
   TBranch        *b_Bu_DIRA_OWNPV;   //!
   TBranch        *b_Bu_TRUEID;   //!
   TBranch        *b_Bu_TRUEP_E;   //!
   TBranch        *b_Bu_TRUEP_X;   //!
   TBranch        *b_Bu_TRUEP_Y;   //!
   TBranch        *b_Bu_TRUEP_Z;   //!
   TBranch        *b_Bu_TRUEPT;   //!
   TBranch        *b_Bu_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_Bu_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_Bu_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_Bu_TRUEENDVERTEX_X;   //!
   TBranch        *b_Bu_TRUEENDVERTEX_Y;   //!
   TBranch        *b_Bu_TRUEENDVERTEX_Z;   //!
   TBranch        *b_Bu_TRUEISSTABLE;   //!
   TBranch        *b_Bu_TRUETAU;   //!
   TBranch        *b_Bplus0_ENDVERTEX_X;   //!
   TBranch        *b_Bplus0_ENDVERTEX_Y;   //!
   TBranch        *b_Bplus0_ENDVERTEX_Z;   //!
   TBranch        *b_Bplus0_ENDVERTEX_XERR;   //!
   TBranch        *b_Bplus0_ENDVERTEX_YERR;   //!
   TBranch        *b_Bplus0_ENDVERTEX_ZERR;   //!
   TBranch        *b_Bplus0_ENDVERTEX_CHI2;   //!
   TBranch        *b_Bplus0_ENDVERTEX_NDOF;   //!
   TBranch        *b_Bplus0_ENDVERTEX_COV_;   //!
   TBranch        *b_Bplus0_OWNPV_X;   //!
   TBranch        *b_Bplus0_OWNPV_Y;   //!
   TBranch        *b_Bplus0_OWNPV_Z;   //!
   TBranch        *b_Bplus0_OWNPV_XERR;   //!
   TBranch        *b_Bplus0_OWNPV_YERR;   //!
   TBranch        *b_Bplus0_OWNPV_ZERR;   //!
   TBranch        *b_Bplus0_OWNPV_CHI2;   //!
   TBranch        *b_Bplus0_OWNPV_NDOF;   //!
   TBranch        *b_Bplus0_OWNPV_COV_;   //!
   TBranch        *b_Bplus0_IP_OWNPV;   //!
   TBranch        *b_Bplus0_IPCHI2_OWNPV;   //!
   TBranch        *b_Bplus0_FD_OWNPV;   //!
   TBranch        *b_Bplus0_FDCHI2_OWNPV;   //!
   TBranch        *b_Bplus0_DIRA_OWNPV;   //!
   TBranch        *b_Bplus0_P;   //!
   TBranch        *b_Bplus0_PT;   //!
   TBranch        *b_Bplus0_PE;   //!
   TBranch        *b_Bplus0_PX;   //!
   TBranch        *b_Bplus0_PY;   //!
   TBranch        *b_Bplus0_PZ;   //!
   TBranch        *b_Bplus0_MM;   //!
   TBranch        *b_Bplus0_MMERR;   //!
   TBranch        *b_Bplus0_M;   //!
   TBranch        *b_Bplus0_ID;   //!
   TBranch        *b_Bplus0_L0Global_Dec;   //!
   TBranch        *b_Bplus0_L0Global_TIS;   //!
   TBranch        *b_Bplus0_L0Global_TOS;   //!
   TBranch        *b_Bplus0_Hlt1Global_Dec;   //!
   TBranch        *b_Bplus0_Hlt1Global_TIS;   //!
   TBranch        *b_Bplus0_Hlt1Global_TOS;   //!
   TBranch        *b_Bplus0_Hlt1Phys_Dec;   //!
   TBranch        *b_Bplus0_Hlt1Phys_TIS;   //!
   TBranch        *b_Bplus0_Hlt1Phys_TOS;   //!
   TBranch        *b_Bplus0_Hlt2Global_Dec;   //!
   TBranch        *b_Bplus0_Hlt2Global_TIS;   //!
   TBranch        *b_Bplus0_Hlt2Global_TOS;   //!
   TBranch        *b_Bplus0_Hlt2Phys_Dec;   //!
   TBranch        *b_Bplus0_Hlt2Phys_TIS;   //!
   TBranch        *b_Bplus0_Hlt2Phys_TOS;   //!
   TBranch        *b_Bplus0_L0MuonDecision_Dec;   //!
   TBranch        *b_Bplus0_L0MuonDecision_TIS;   //!
   TBranch        *b_Bplus0_L0MuonDecision_TOS;   //!
   TBranch        *b_Bplus0_L0DiMuonDecision_Dec;   //!
   TBranch        *b_Bplus0_L0DiMuonDecision_TIS;   //!
   TBranch        *b_Bplus0_L0DiMuonDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedPsi2SDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonJPsiDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonJPsiDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonJPsiDecision_TOS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonJPsiHighPTDecision_Dec;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonJPsiHighPTDecision_TIS;   //!
   TBranch        *b_Bplus0_Hlt2DiMuonJPsiHighPTDecision_TOS;   //!
   TBranch        *b_Jpsi_P;   //!
   TBranch        *b_Jpsi_PT;   //!
   TBranch        *b_Jpsi_PE;   //!
   TBranch        *b_Jpsi_PX;   //!
   TBranch        *b_Jpsi_PY;   //!
   TBranch        *b_Jpsi_PZ;   //!
   TBranch        *b_Jpsi_MM;   //!
   TBranch        *b_Jpsi_MMERR;   //!
   TBranch        *b_Jpsi_M;   //!
   TBranch        *b_Jpsi_ID;   //!
   TBranch        *b_Jpsi_L0Global_Dec;   //!
   TBranch        *b_Jpsi_L0Global_TIS;   //!
   TBranch        *b_Jpsi_L0Global_TOS;   //!
   TBranch        *b_Jpsi_Hlt1Global_Dec;   //!
   TBranch        *b_Jpsi_Hlt1Global_TIS;   //!
   TBranch        *b_Jpsi_Hlt1Global_TOS;   //!
   TBranch        *b_Jpsi_Hlt1Phys_Dec;   //!
   TBranch        *b_Jpsi_Hlt1Phys_TIS;   //!
   TBranch        *b_Jpsi_Hlt1Phys_TOS;   //!
   TBranch        *b_Jpsi_Hlt2Global_Dec;   //!
   TBranch        *b_Jpsi_Hlt2Global_TIS;   //!
   TBranch        *b_Jpsi_Hlt2Global_TOS;   //!
   TBranch        *b_Jpsi_Hlt2Phys_Dec;   //!
   TBranch        *b_Jpsi_Hlt2Phys_TIS;   //!
   TBranch        *b_Jpsi_Hlt2Phys_TOS;   //!
   TBranch        *b_Jpsi_L0MuonDecision_Dec;   //!
   TBranch        *b_Jpsi_L0MuonDecision_TIS;   //!
   TBranch        *b_Jpsi_L0MuonDecision_TOS;   //!
   TBranch        *b_Jpsi_L0DiMuonDecision_Dec;   //!
   TBranch        *b_Jpsi_L0DiMuonDecision_TIS;   //!
   TBranch        *b_Jpsi_L0DiMuonDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedPsi2SDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonJPsiDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonJPsiDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonJPsiDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonJPsiHighPTDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonJPsiHighPTDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS;   //!
   TBranch        *b_Jpsi_ENDVERTEX_X;   //!
   TBranch        *b_Jpsi_ENDVERTEX_Y;   //!
   TBranch        *b_Jpsi_ENDVERTEX_Z;   //!
   TBranch        *b_Jpsi_ENDVERTEX_XERR;   //!
   TBranch        *b_Jpsi_ENDVERTEX_YERR;   //!
   TBranch        *b_Jpsi_ENDVERTEX_ZERR;   //!
   TBranch        *b_Jpsi_ENDVERTEX_CHI2;   //!
   TBranch        *b_Jpsi_ENDVERTEX_NDOF;   //!
   TBranch        *b_Jpsi_ENDVERTEX_COV_;   //!
   TBranch        *b_Jpsi_OWNPV_X;   //!
   TBranch        *b_Jpsi_OWNPV_Y;   //!
   TBranch        *b_Jpsi_OWNPV_Z;   //!
   TBranch        *b_Jpsi_OWNPV_XERR;   //!
   TBranch        *b_Jpsi_OWNPV_YERR;   //!
   TBranch        *b_Jpsi_OWNPV_ZERR;   //!
   TBranch        *b_Jpsi_OWNPV_CHI2;   //!
   TBranch        *b_Jpsi_OWNPV_NDOF;   //!
   TBranch        *b_Jpsi_OWNPV_COV_;   //!
   TBranch        *b_Jpsi_IP_OWNPV;   //!
   TBranch        *b_Jpsi_IPCHI2_OWNPV;   //!
   TBranch        *b_Jpsi_FD_OWNPV;   //!
   TBranch        *b_Jpsi_FDCHI2_OWNPV;   //!
   TBranch        *b_Jpsi_DIRA_OWNPV;   //!
   TBranch        *b_Jpsi_ORIVX_X;   //!
   TBranch        *b_Jpsi_ORIVX_Y;   //!
   TBranch        *b_Jpsi_ORIVX_Z;   //!
   TBranch        *b_Jpsi_ORIVX_XERR;   //!
   TBranch        *b_Jpsi_ORIVX_YERR;   //!
   TBranch        *b_Jpsi_ORIVX_ZERR;   //!
   TBranch        *b_Jpsi_ORIVX_CHI2;   //!
   TBranch        *b_Jpsi_ORIVX_NDOF;   //!
   TBranch        *b_Jpsi_ORIVX_COV_;   //!
   TBranch        *b_Jpsi_FD_ORIVX;   //!
   TBranch        *b_Jpsi_FDCHI2_ORIVX;   //!
   TBranch        *b_Jpsi_DIRA_ORIVX;   //!
   TBranch        *b_Jpsi_TRUEID;   //!
   TBranch        *b_Jpsi_TRUEP_E;   //!
   TBranch        *b_Jpsi_TRUEP_X;   //!
   TBranch        *b_Jpsi_TRUEP_Y;   //!
   TBranch        *b_Jpsi_TRUEP_Z;   //!
   TBranch        *b_Jpsi_TRUEPT;   //!
   TBranch        *b_Jpsi_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_Jpsi_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_Jpsi_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_Jpsi_TRUEENDVERTEX_X;   //!
   TBranch        *b_Jpsi_TRUEENDVERTEX_Y;   //!
   TBranch        *b_Jpsi_TRUEENDVERTEX_Z;   //!
   TBranch        *b_Jpsi_TRUEISSTABLE;   //!
   TBranch        *b_Jpsi_TRUETAU;   //!
   TBranch        *b_mup_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_mup_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_mup_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_mup_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_mup_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_mup_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_mup_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_mup_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_mup_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_mup_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_mup_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_mup_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_mup_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_mup_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_mup_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_mup_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_mup_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_mup_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_mup_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_mup_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_mup_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_mup_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_mup_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_mup_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_mup_P;   //!
   TBranch        *b_mup_PT;   //!
   TBranch        *b_mup_PE;   //!
   TBranch        *b_mup_PX;   //!
   TBranch        *b_mup_PY;   //!
   TBranch        *b_mup_PZ;   //!
   TBranch        *b_mup_M;   //!
   TBranch        *b_mup_ID;   //!
   TBranch        *b_mup_PIDe;   //!
   TBranch        *b_mup_PIDmu;   //!
   TBranch        *b_mup_PIDK;   //!
   TBranch        *b_mup_PIDp;   //!
   TBranch        *b_mup_PIDd;   //!
   TBranch        *b_mup_ProbNNe;   //!
   TBranch        *b_mup_ProbNNk;   //!
   TBranch        *b_mup_ProbNNp;   //!
   TBranch        *b_mup_ProbNNpi;   //!
   TBranch        *b_mup_ProbNNmu;   //!
   TBranch        *b_mup_ProbNNd;   //!
   TBranch        *b_mup_ProbNNghost;   //!
   TBranch        *b_mup_hasMuon;   //!
   TBranch        *b_mup_isMuon;   //!
   TBranch        *b_mup_hasRich;   //!
   TBranch        *b_mup_UsedRichAerogel;   //!
   TBranch        *b_mup_UsedRich1Gas;   //!
   TBranch        *b_mup_UsedRich2Gas;   //!
   TBranch        *b_mup_RichAboveElThres;   //!
   TBranch        *b_mup_RichAboveMuThres;   //!
   TBranch        *b_mup_RichAbovePiThres;   //!
   TBranch        *b_mup_RichAboveKaThres;   //!
   TBranch        *b_mup_RichAbovePrThres;   //!
   TBranch        *b_mup_hasCalo;   //!
   TBranch        *b_mup_L0Global_Dec;   //!
   TBranch        *b_mup_L0Global_TIS;   //!
   TBranch        *b_mup_L0Global_TOS;   //!
   TBranch        *b_mup_Hlt1Global_Dec;   //!
   TBranch        *b_mup_Hlt1Global_TIS;   //!
   TBranch        *b_mup_Hlt1Global_TOS;   //!
   TBranch        *b_mup_Hlt1Phys_Dec;   //!
   TBranch        *b_mup_Hlt1Phys_TIS;   //!
   TBranch        *b_mup_Hlt1Phys_TOS;   //!
   TBranch        *b_mup_Hlt2Global_Dec;   //!
   TBranch        *b_mup_Hlt2Global_TIS;   //!
   TBranch        *b_mup_Hlt2Global_TOS;   //!
   TBranch        *b_mup_Hlt2Phys_Dec;   //!
   TBranch        *b_mup_Hlt2Phys_TIS;   //!
   TBranch        *b_mup_Hlt2Phys_TOS;   //!
   TBranch        *b_mup_L0MuonDecision_Dec;   //!
   TBranch        *b_mup_L0MuonDecision_TIS;   //!
   TBranch        *b_mup_L0MuonDecision_TOS;   //!
   TBranch        *b_mup_L0DiMuonDecision_Dec;   //!
   TBranch        *b_mup_L0DiMuonDecision_TIS;   //!
   TBranch        *b_mup_L0DiMuonDecision_TOS;   //!
   TBranch        *b_mup_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_mup_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_mup_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedPsi2SDecision_Dec;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedPsi2SDecision_TIS;   //!
   TBranch        *b_mup_Hlt2DiMuonDetachedPsi2SDecision_TOS;   //!
   TBranch        *b_mup_Hlt2DiMuonJPsiDecision_Dec;   //!
   TBranch        *b_mup_Hlt2DiMuonJPsiDecision_TIS;   //!
   TBranch        *b_mup_Hlt2DiMuonJPsiDecision_TOS;   //!
   TBranch        *b_mup_Hlt2DiMuonJPsiHighPTDecision_Dec;   //!
   TBranch        *b_mup_Hlt2DiMuonJPsiHighPTDecision_TIS;   //!
   TBranch        *b_mup_Hlt2DiMuonJPsiHighPTDecision_TOS;   //!
   TBranch        *b_mup_OWNPV_X;   //!
   TBranch        *b_mup_OWNPV_Y;   //!
   TBranch        *b_mup_OWNPV_Z;   //!
   TBranch        *b_mup_OWNPV_XERR;   //!
   TBranch        *b_mup_OWNPV_YERR;   //!
   TBranch        *b_mup_OWNPV_ZERR;   //!
   TBranch        *b_mup_OWNPV_CHI2;   //!
   TBranch        *b_mup_OWNPV_NDOF;   //!
   TBranch        *b_mup_OWNPV_COV_;   //!
   TBranch        *b_mup_IP_OWNPV;   //!
   TBranch        *b_mup_IPCHI2_OWNPV;   //!
   TBranch        *b_mup_ORIVX_X;   //!
   TBranch        *b_mup_ORIVX_Y;   //!
   TBranch        *b_mup_ORIVX_Z;   //!
   TBranch        *b_mup_ORIVX_XERR;   //!
   TBranch        *b_mup_ORIVX_YERR;   //!
   TBranch        *b_mup_ORIVX_ZERR;   //!
   TBranch        *b_mup_ORIVX_CHI2;   //!
   TBranch        *b_mup_ORIVX_NDOF;   //!
   TBranch        *b_mup_ORIVX_COV_;   //!
   TBranch        *b_mup_TRUEID;   //!
   TBranch        *b_mup_TRUEP_E;   //!
   TBranch        *b_mup_TRUEP_X;   //!
   TBranch        *b_mup_TRUEP_Y;   //!
   TBranch        *b_mup_TRUEP_Z;   //!
   TBranch        *b_mup_TRUEPT;   //!
   TBranch        *b_mup_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_mup_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_mup_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_mup_TRUEENDVERTEX_X;   //!
   TBranch        *b_mup_TRUEENDVERTEX_Y;   //!
   TBranch        *b_mup_TRUEENDVERTEX_Z;   //!
   TBranch        *b_mup_TRUEISSTABLE;   //!
   TBranch        *b_mup_TRUETAU;   //!
   TBranch        *b_mup_TRACK_Type;   //!
   TBranch        *b_mup_TRACK_Key;   //!
   TBranch        *b_mup_TRACK_CHI2NDOF;   //!
   TBranch        *b_mup_TRACK_PCHI2;   //!
   TBranch        *b_mup_TRACK_MatchCHI2;   //!
   TBranch        *b_mup_TRACK_GhostProb;   //!
   TBranch        *b_mup_TRACK_CloneDist;   //!
   TBranch        *b_mup_TRACK_Likelihood;   //!
   TBranch        *b_mum_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_mum_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_mum_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_mum_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_mum_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_mum_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_mum_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_mum_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_mum_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_mum_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_mum_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_mum_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_mum_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_mum_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_mum_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_mum_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_mum_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_mum_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_mum_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_mum_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_mum_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_mum_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_mum_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_mum_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_mum_P;   //!
   TBranch        *b_mum_PT;   //!
   TBranch        *b_mum_PE;   //!
   TBranch        *b_mum_PX;   //!
   TBranch        *b_mum_PY;   //!
   TBranch        *b_mum_PZ;   //!
   TBranch        *b_mum_M;   //!
   TBranch        *b_mum_ID;   //!
   TBranch        *b_mum_PIDe;   //!
   TBranch        *b_mum_PIDmu;   //!
   TBranch        *b_mum_PIDK;   //!
   TBranch        *b_mum_PIDp;   //!
   TBranch        *b_mum_PIDd;   //!
   TBranch        *b_mum_ProbNNe;   //!
   TBranch        *b_mum_ProbNNk;   //!
   TBranch        *b_mum_ProbNNp;   //!
   TBranch        *b_mum_ProbNNpi;   //!
   TBranch        *b_mum_ProbNNmu;   //!
   TBranch        *b_mum_ProbNNd;   //!
   TBranch        *b_mum_ProbNNghost;   //!
   TBranch        *b_mum_hasMuon;   //!
   TBranch        *b_mum_isMuon;   //!
   TBranch        *b_mum_hasRich;   //!
   TBranch        *b_mum_UsedRichAerogel;   //!
   TBranch        *b_mum_UsedRich1Gas;   //!
   TBranch        *b_mum_UsedRich2Gas;   //!
   TBranch        *b_mum_RichAboveElThres;   //!
   TBranch        *b_mum_RichAboveMuThres;   //!
   TBranch        *b_mum_RichAbovePiThres;   //!
   TBranch        *b_mum_RichAboveKaThres;   //!
   TBranch        *b_mum_RichAbovePrThres;   //!
   TBranch        *b_mum_hasCalo;   //!
   TBranch        *b_mum_L0Global_Dec;   //!
   TBranch        *b_mum_L0Global_TIS;   //!
   TBranch        *b_mum_L0Global_TOS;   //!
   TBranch        *b_mum_Hlt1Global_Dec;   //!
   TBranch        *b_mum_Hlt1Global_TIS;   //!
   TBranch        *b_mum_Hlt1Global_TOS;   //!
   TBranch        *b_mum_Hlt1Phys_Dec;   //!
   TBranch        *b_mum_Hlt1Phys_TIS;   //!
   TBranch        *b_mum_Hlt1Phys_TOS;   //!
   TBranch        *b_mum_Hlt2Global_Dec;   //!
   TBranch        *b_mum_Hlt2Global_TIS;   //!
   TBranch        *b_mum_Hlt2Global_TOS;   //!
   TBranch        *b_mum_Hlt2Phys_Dec;   //!
   TBranch        *b_mum_Hlt2Phys_TIS;   //!
   TBranch        *b_mum_Hlt2Phys_TOS;   //!
   TBranch        *b_mum_L0MuonDecision_Dec;   //!
   TBranch        *b_mum_L0MuonDecision_TIS;   //!
   TBranch        *b_mum_L0MuonDecision_TOS;   //!
   TBranch        *b_mum_L0DiMuonDecision_Dec;   //!
   TBranch        *b_mum_L0DiMuonDecision_TIS;   //!
   TBranch        *b_mum_L0DiMuonDecision_TOS;   //!
   TBranch        *b_mum_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_mum_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_mum_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedPsi2SDecision_Dec;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedPsi2SDecision_TIS;   //!
   TBranch        *b_mum_Hlt2DiMuonDetachedPsi2SDecision_TOS;   //!
   TBranch        *b_mum_Hlt2DiMuonJPsiDecision_Dec;   //!
   TBranch        *b_mum_Hlt2DiMuonJPsiDecision_TIS;   //!
   TBranch        *b_mum_Hlt2DiMuonJPsiDecision_TOS;   //!
   TBranch        *b_mum_Hlt2DiMuonJPsiHighPTDecision_Dec;   //!
   TBranch        *b_mum_Hlt2DiMuonJPsiHighPTDecision_TIS;   //!
   TBranch        *b_mum_Hlt2DiMuonJPsiHighPTDecision_TOS;   //!
   TBranch        *b_mum_OWNPV_X;   //!
   TBranch        *b_mum_OWNPV_Y;   //!
   TBranch        *b_mum_OWNPV_Z;   //!
   TBranch        *b_mum_OWNPV_XERR;   //!
   TBranch        *b_mum_OWNPV_YERR;   //!
   TBranch        *b_mum_OWNPV_ZERR;   //!
   TBranch        *b_mum_OWNPV_CHI2;   //!
   TBranch        *b_mum_OWNPV_NDOF;   //!
   TBranch        *b_mum_OWNPV_COV_;   //!
   TBranch        *b_mum_IP_OWNPV;   //!
   TBranch        *b_mum_IPCHI2_OWNPV;   //!
   TBranch        *b_mum_ORIVX_X;   //!
   TBranch        *b_mum_ORIVX_Y;   //!
   TBranch        *b_mum_ORIVX_Z;   //!
   TBranch        *b_mum_ORIVX_XERR;   //!
   TBranch        *b_mum_ORIVX_YERR;   //!
   TBranch        *b_mum_ORIVX_ZERR;   //!
   TBranch        *b_mum_ORIVX_CHI2;   //!
   TBranch        *b_mum_ORIVX_NDOF;   //!
   TBranch        *b_mum_ORIVX_COV_;   //!
   TBranch        *b_mum_TRUEID;   //!
   TBranch        *b_mum_TRUEP_E;   //!
   TBranch        *b_mum_TRUEP_X;   //!
   TBranch        *b_mum_TRUEP_Y;   //!
   TBranch        *b_mum_TRUEP_Z;   //!
   TBranch        *b_mum_TRUEPT;   //!
   TBranch        *b_mum_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_mum_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_mum_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_mum_TRUEENDVERTEX_X;   //!
   TBranch        *b_mum_TRUEENDVERTEX_Y;   //!
   TBranch        *b_mum_TRUEENDVERTEX_Z;   //!
   TBranch        *b_mum_TRUEISSTABLE;   //!
   TBranch        *b_mum_TRUETAU;   //!
   TBranch        *b_mum_TRACK_Type;   //!
   TBranch        *b_mum_TRACK_Key;   //!
   TBranch        *b_mum_TRACK_CHI2NDOF;   //!
   TBranch        *b_mum_TRACK_PCHI2;   //!
   TBranch        *b_mum_TRACK_MatchCHI2;   //!
   TBranch        *b_mum_TRACK_GhostProb;   //!
   TBranch        *b_mum_TRACK_CloneDist;   //!
   TBranch        *b_mum_TRACK_Likelihood;   //!
   TBranch        *b_K_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_K_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_K_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_K_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_K_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_K_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_K_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_K_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_K_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_K_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_K_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_K_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_K_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_K_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_K_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_K_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_K_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_K_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_K_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_K_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_K_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_K_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_K_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_K_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_K_P;   //!
   TBranch        *b_K_PT;   //!
   TBranch        *b_K_PE;   //!
   TBranch        *b_K_PX;   //!
   TBranch        *b_K_PY;   //!
   TBranch        *b_K_PZ;   //!
   TBranch        *b_K_M;   //!
   TBranch        *b_K_ID;   //!
   TBranch        *b_K_PIDe;   //!
   TBranch        *b_K_PIDmu;   //!
   TBranch        *b_K_PIDK;   //!
   TBranch        *b_K_PIDp;   //!
   TBranch        *b_K_PIDd;   //!
   TBranch        *b_K_ProbNNe;   //!
   TBranch        *b_K_ProbNNk;   //!
   TBranch        *b_K_ProbNNp;   //!
   TBranch        *b_K_ProbNNpi;   //!
   TBranch        *b_K_ProbNNmu;   //!
   TBranch        *b_K_ProbNNd;   //!
   TBranch        *b_K_ProbNNghost;   //!
   TBranch        *b_K_hasMuon;   //!
   TBranch        *b_K_isMuon;   //!
   TBranch        *b_K_hasRich;   //!
   TBranch        *b_K_UsedRichAerogel;   //!
   TBranch        *b_K_UsedRich1Gas;   //!
   TBranch        *b_K_UsedRich2Gas;   //!
   TBranch        *b_K_RichAboveElThres;   //!
   TBranch        *b_K_RichAboveMuThres;   //!
   TBranch        *b_K_RichAbovePiThres;   //!
   TBranch        *b_K_RichAboveKaThres;   //!
   TBranch        *b_K_RichAbovePrThres;   //!
   TBranch        *b_K_hasCalo;   //!
   TBranch        *b_K_L0Global_Dec;   //!
   TBranch        *b_K_L0Global_TIS;   //!
   TBranch        *b_K_L0Global_TOS;   //!
   TBranch        *b_K_Hlt1Global_Dec;   //!
   TBranch        *b_K_Hlt1Global_TIS;   //!
   TBranch        *b_K_Hlt1Global_TOS;   //!
   TBranch        *b_K_Hlt1Phys_Dec;   //!
   TBranch        *b_K_Hlt1Phys_TIS;   //!
   TBranch        *b_K_Hlt1Phys_TOS;   //!
   TBranch        *b_K_Hlt2Global_Dec;   //!
   TBranch        *b_K_Hlt2Global_TIS;   //!
   TBranch        *b_K_Hlt2Global_TOS;   //!
   TBranch        *b_K_Hlt2Phys_Dec;   //!
   TBranch        *b_K_Hlt2Phys_TIS;   //!
   TBranch        *b_K_Hlt2Phys_TOS;   //!
   TBranch        *b_K_L0MuonDecision_Dec;   //!
   TBranch        *b_K_L0MuonDecision_TIS;   //!
   TBranch        *b_K_L0MuonDecision_TOS;   //!
   TBranch        *b_K_L0DiMuonDecision_Dec;   //!
   TBranch        *b_K_L0DiMuonDecision_TIS;   //!
   TBranch        *b_K_L0DiMuonDecision_TOS;   //!
   TBranch        *b_K_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_K_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_K_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedDecision_Dec;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedDecision_TIS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedDecision_TOS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedPsi2SDecision_Dec;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedPsi2SDecision_TIS;   //!
   TBranch        *b_K_Hlt2DiMuonDetachedPsi2SDecision_TOS;   //!
   TBranch        *b_K_Hlt2DiMuonJPsiDecision_Dec;   //!
   TBranch        *b_K_Hlt2DiMuonJPsiDecision_TIS;   //!
   TBranch        *b_K_Hlt2DiMuonJPsiDecision_TOS;   //!
   TBranch        *b_K_Hlt2DiMuonJPsiHighPTDecision_Dec;   //!
   TBranch        *b_K_Hlt2DiMuonJPsiHighPTDecision_TIS;   //!
   TBranch        *b_K_Hlt2DiMuonJPsiHighPTDecision_TOS;   //!
   TBranch        *b_K_OWNPV_X;   //!
   TBranch        *b_K_OWNPV_Y;   //!
   TBranch        *b_K_OWNPV_Z;   //!
   TBranch        *b_K_OWNPV_XERR;   //!
   TBranch        *b_K_OWNPV_YERR;   //!
   TBranch        *b_K_OWNPV_ZERR;   //!
   TBranch        *b_K_OWNPV_CHI2;   //!
   TBranch        *b_K_OWNPV_NDOF;   //!
   TBranch        *b_K_OWNPV_COV_;   //!
   TBranch        *b_K_IP_OWNPV;   //!
   TBranch        *b_K_IPCHI2_OWNPV;   //!
   TBranch        *b_K_ORIVX_X;   //!
   TBranch        *b_K_ORIVX_Y;   //!
   TBranch        *b_K_ORIVX_Z;   //!
   TBranch        *b_K_ORIVX_XERR;   //!
   TBranch        *b_K_ORIVX_YERR;   //!
   TBranch        *b_K_ORIVX_ZERR;   //!
   TBranch        *b_K_ORIVX_CHI2;   //!
   TBranch        *b_K_ORIVX_NDOF;   //!
   TBranch        *b_K_ORIVX_COV_;   //!
   TBranch        *b_K_TRUEID;   //!
   TBranch        *b_K_TRUEP_E;   //!
   TBranch        *b_K_TRUEP_X;   //!
   TBranch        *b_K_TRUEP_Y;   //!
   TBranch        *b_K_TRUEP_Z;   //!
   TBranch        *b_K_TRUEPT;   //!
   TBranch        *b_K_TRUEORIGINVERTEX_X;   //!
   TBranch        *b_K_TRUEORIGINVERTEX_Y;   //!
   TBranch        *b_K_TRUEORIGINVERTEX_Z;   //!
   TBranch        *b_K_TRUEENDVERTEX_X;   //!
   TBranch        *b_K_TRUEENDVERTEX_Y;   //!
   TBranch        *b_K_TRUEENDVERTEX_Z;   //!
   TBranch        *b_K_TRUEISSTABLE;   //!
   TBranch        *b_K_TRUETAU;   //!
   TBranch        *b_K_TRACK_Type;   //!
   TBranch        *b_K_TRACK_Key;   //!
   TBranch        *b_K_TRACK_CHI2NDOF;   //!
   TBranch        *b_K_TRACK_PCHI2;   //!
   TBranch        *b_K_TRACK_MatchCHI2;   //!
   TBranch        *b_K_TRACK_GhostProb;   //!
   TBranch        *b_K_TRACK_CloneDist;   //!
   TBranch        *b_K_TRACK_Likelihood;   //!
   TBranch        *b_nCandidate;   //!
   TBranch        *b_totCandidates;   //!
   TBranch        *b_EventInSequence;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_BCID;   //!
   TBranch        *b_BCType;   //!
   TBranch        *b_OdinTCK;   //!
   TBranch        *b_L0DUTCK;   //!
   TBranch        *b_HLT1TCK;   //!
   TBranch        *b_HLT2TCK;   //!
   TBranch        *b_GpsTime;   //!
   TBranch        *b_Polarity;   //!
   TBranch        *b_nPVs;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nLongTracks;   //!
   TBranch        *b_nDownstreamTracks;   //!
   TBranch        *b_nUpstreamTracks;   //!
   TBranch        *b_nVeloTracks;   //!
   TBranch        *b_nTTracks;   //!
   TBranch        *b_nBackTracks;   //!
   TBranch        *b_nRich1Hits;   //!
   TBranch        *b_nRich2Hits;   //!
   TBranch        *b_nVeloClusters;   //!
   TBranch        *b_nITClusters;   //!
   TBranch        *b_nTTClusters;   //!
   TBranch        *b_nOTClusters;   //!
   TBranch        *b_nSPDHits;   //!
   TBranch        *b_nMuonCoordsS0;   //!
   TBranch        *b_nMuonCoordsS1;   //!
   TBranch        *b_nMuonCoordsS2;   //!
   TBranch        *b_nMuonCoordsS3;   //!
   TBranch        *b_nMuonCoordsS4;   //!
   TBranch        *b_nMuonTracks;   //!

   BDecayTree(TTree *tree = 0, int dataset = 93539, bool isData = true);
   virtual ~BDecayTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BDecayTree_cxx
BDecayTree::BDecayTree(TTree *tree = 0, int dataset = 93539, bool isData = true) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) 
   {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("Decay/DecayTree",tree);

#else // SINGLE_TREE

        // The following code should be used if you want this class to access a chain
        // of trees.
        TChain * chain = new TChain("Decay/DecayTree","");
        TString data_dir = "../../data/";
       
        if (dataset == 91599)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MU_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MU_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MU_05182024.root/Decay/DecayTree");
            }
            else
            {

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/Decay/DecayTree");

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/Decay/DecayTree");  

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/Decay/DecayTree");   
            }
        }
        else if (dataset == 91590)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MD_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/Decay/DecayTree");

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/Decay/DecayTree");

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 91591)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MU_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MU_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MU_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/Decay/DecayTree");

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/Decay/DecayTree");

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/Decay/DecayTree");
            }
        }
        if (dataset == 61590)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MD_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 61591)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MU_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 61599)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2016_MU_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 71590)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MD_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 71591)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MU_05182024.root/Decay/DecayTree");
            }
            else
            {

                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 71599)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2017_MU_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/Decay/DecayTree");              
            }
        }
        else if (dataset == 81590)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MD_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 81591)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MU_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/Decay/DecayTree");
            }
        }
        else if (dataset == 81599)
        {
            if (isData)
            {
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MD_05182024.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_HighPT_2018_MU_05182024.root/Decay/DecayTree");
            }
            else
            {
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/Decay/DecayTree");
                chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/Decay/DecayTree");              
            }
        }
       

      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

BDecayTree::~BDecayTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BDecayTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BDecayTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BDecayTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Bu_P", &Bu_P, &b_Bu_P);
   fChain->SetBranchAddress("Bu_PT", &Bu_PT, &b_Bu_PT);
   fChain->SetBranchAddress("Bu_PE", &Bu_PE, &b_Bu_PE);
   fChain->SetBranchAddress("Bu_PX", &Bu_PX, &b_Bu_PX);
   fChain->SetBranchAddress("Bu_PY", &Bu_PY, &b_Bu_PY);
   fChain->SetBranchAddress("Bu_PZ", &Bu_PZ, &b_Bu_PZ);
   fChain->SetBranchAddress("Bu_MM", &Bu_MM, &b_Bu_MM);
   fChain->SetBranchAddress("Bu_MMERR", &Bu_MMERR, &b_Bu_MMERR);
   fChain->SetBranchAddress("Bu_M", &Bu_M, &b_Bu_M);
   fChain->SetBranchAddress("Bu_ID", &Bu_ID, &b_Bu_ID);
   fChain->SetBranchAddress("Bu_L0Global_Dec", &Bu_L0Global_Dec, &b_Bu_L0Global_Dec);
   fChain->SetBranchAddress("Bu_L0Global_TIS", &Bu_L0Global_TIS, &b_Bu_L0Global_TIS);
   fChain->SetBranchAddress("Bu_L0Global_TOS", &Bu_L0Global_TOS, &b_Bu_L0Global_TOS);
   fChain->SetBranchAddress("Bu_Hlt1Global_Dec", &Bu_Hlt1Global_Dec, &b_Bu_Hlt1Global_Dec);
   fChain->SetBranchAddress("Bu_Hlt1Global_TIS", &Bu_Hlt1Global_TIS, &b_Bu_Hlt1Global_TIS);
   fChain->SetBranchAddress("Bu_Hlt1Global_TOS", &Bu_Hlt1Global_TOS, &b_Bu_Hlt1Global_TOS);
   fChain->SetBranchAddress("Bu_Hlt1Phys_Dec", &Bu_Hlt1Phys_Dec, &b_Bu_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Bu_Hlt1Phys_TIS", &Bu_Hlt1Phys_TIS, &b_Bu_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Bu_Hlt1Phys_TOS", &Bu_Hlt1Phys_TOS, &b_Bu_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Bu_Hlt2Global_Dec", &Bu_Hlt2Global_Dec, &b_Bu_Hlt2Global_Dec);
   fChain->SetBranchAddress("Bu_Hlt2Global_TIS", &Bu_Hlt2Global_TIS, &b_Bu_Hlt2Global_TIS);
   fChain->SetBranchAddress("Bu_Hlt2Global_TOS", &Bu_Hlt2Global_TOS, &b_Bu_Hlt2Global_TOS);
   fChain->SetBranchAddress("Bu_Hlt2Phys_Dec", &Bu_Hlt2Phys_Dec, &b_Bu_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Bu_Hlt2Phys_TIS", &Bu_Hlt2Phys_TIS, &b_Bu_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Bu_Hlt2Phys_TOS", &Bu_Hlt2Phys_TOS, &b_Bu_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Bu_L0MuonDecision_Dec", &Bu_L0MuonDecision_Dec, &b_Bu_L0MuonDecision_Dec);
   fChain->SetBranchAddress("Bu_L0MuonDecision_TIS", &Bu_L0MuonDecision_TIS, &b_Bu_L0MuonDecision_TIS);
   fChain->SetBranchAddress("Bu_L0MuonDecision_TOS", &Bu_L0MuonDecision_TOS, &b_Bu_L0MuonDecision_TOS);
   fChain->SetBranchAddress("Bu_L0DiMuonDecision_Dec", &Bu_L0DiMuonDecision_Dec, &b_Bu_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("Bu_L0DiMuonDecision_TIS", &Bu_L0DiMuonDecision_TIS, &b_Bu_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("Bu_L0DiMuonDecision_TOS", &Bu_L0DiMuonDecision_TOS, &b_Bu_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt1DiMuonHighMassDecision_Dec", &Bu_Hlt1DiMuonHighMassDecision_Dec, &b_Bu_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt1DiMuonHighMassDecision_TIS", &Bu_Hlt1DiMuonHighMassDecision_TIS, &b_Bu_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt1DiMuonHighMassDecision_TOS", &Bu_Hlt1DiMuonHighMassDecision_TOS, &b_Bu_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedDecision_Dec", &Bu_Hlt2DiMuonDetachedDecision_Dec, &b_Bu_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedDecision_TIS", &Bu_Hlt2DiMuonDetachedDecision_TIS, &b_Bu_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedDecision_TOS", &Bu_Hlt2DiMuonDetachedDecision_TOS, &b_Bu_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedHeavyDecision_Dec", &Bu_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_Bu_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedHeavyDecision_TIS", &Bu_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_Bu_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedHeavyDecision_TOS", &Bu_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_Bu_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedJPsiDecision_Dec", &Bu_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_Bu_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedJPsiDecision_TIS", &Bu_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_Bu_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedJPsiDecision_TOS", &Bu_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_Bu_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedPsi2SDecision_Dec", &Bu_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_Bu_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedPsi2SDecision_TIS", &Bu_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_Bu_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonDetachedPsi2SDecision_TOS", &Bu_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_Bu_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonJPsiDecision_Dec", &Bu_Hlt2DiMuonJPsiDecision_Dec, &b_Bu_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonJPsiDecision_TIS", &Bu_Hlt2DiMuonJPsiDecision_TIS, &b_Bu_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonJPsiDecision_TOS", &Bu_Hlt2DiMuonJPsiDecision_TOS, &b_Bu_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonJPsiHighPTDecision_Dec", &Bu_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_Bu_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonJPsiHighPTDecision_TIS", &Bu_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_Bu_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("Bu_Hlt2DiMuonJPsiHighPTDecision_TOS", &Bu_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_Bu_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("Bu_ConsBu_nPV", &Bu_ConsBu_nPV, &b_Bu_ConsBu_nPV);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_M", Bu_ConsBu_J_psi_1S_M, &b_Bu_ConsBu_J_psi_1S_M);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_MERR", Bu_ConsBu_J_psi_1S_MERR, &b_Bu_ConsBu_J_psi_1S_MERR);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_P", Bu_ConsBu_J_psi_1S_P, &b_Bu_ConsBu_J_psi_1S_P);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_PERR", Bu_ConsBu_J_psi_1S_PERR, &b_Bu_ConsBu_J_psi_1S_PERR);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_ctau", Bu_ConsBu_J_psi_1S_ctau, &b_Bu_ConsBu_J_psi_1S_ctau);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_ctauErr", Bu_ConsBu_J_psi_1S_ctauErr, &b_Bu_ConsBu_J_psi_1S_ctauErr);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_decayLength", Bu_ConsBu_J_psi_1S_decayLength, &b_Bu_ConsBu_J_psi_1S_decayLength);
   fChain->SetBranchAddress("Bu_ConsBu_J_psi_1S_decayLengthErr", Bu_ConsBu_J_psi_1S_decayLengthErr, &b_Bu_ConsBu_J_psi_1S_decayLengthErr);
   fChain->SetBranchAddress("Bu_ConsBu_M", Bu_ConsBu_M, &b_Bu_ConsBu_M);
   fChain->SetBranchAddress("Bu_ConsBu_MERR", Bu_ConsBu_MERR, &b_Bu_ConsBu_MERR);
   fChain->SetBranchAddress("Bu_ConsBu_P", Bu_ConsBu_P, &b_Bu_ConsBu_P);
   fChain->SetBranchAddress("Bu_ConsBu_PERR", Bu_ConsBu_PERR, &b_Bu_ConsBu_PERR);
   fChain->SetBranchAddress("Bu_ConsBu_PV_X", Bu_ConsBu_PV_X, &b_Bu_ConsBu_PV_X);
   fChain->SetBranchAddress("Bu_ConsBu_PV_Y", Bu_ConsBu_PV_Y, &b_Bu_ConsBu_PV_Y);
   fChain->SetBranchAddress("Bu_ConsBu_PV_Z", Bu_ConsBu_PV_Z, &b_Bu_ConsBu_PV_Z);
   fChain->SetBranchAddress("Bu_ConsBu_PV_key", Bu_ConsBu_PV_key, &b_Bu_ConsBu_PV_key);
   fChain->SetBranchAddress("Bu_ConsBu_chi2", Bu_ConsBu_chi2, &b_Bu_ConsBu_chi2);
   fChain->SetBranchAddress("Bu_ConsBu_ctau", Bu_ConsBu_ctau, &b_Bu_ConsBu_ctau);
   fChain->SetBranchAddress("Bu_ConsBu_ctauErr", Bu_ConsBu_ctauErr, &b_Bu_ConsBu_ctauErr);
   fChain->SetBranchAddress("Bu_ConsBu_decayLength", Bu_ConsBu_decayLength, &b_Bu_ConsBu_decayLength);
   fChain->SetBranchAddress("Bu_ConsBu_decayLengthErr", Bu_ConsBu_decayLengthErr, &b_Bu_ConsBu_decayLengthErr);
   fChain->SetBranchAddress("Bu_ConsBu_nDOF", Bu_ConsBu_nDOF, &b_Bu_ConsBu_nDOF);
   fChain->SetBranchAddress("Bu_ConsBu_nIter", Bu_ConsBu_nIter, &b_Bu_ConsBu_nIter);
   fChain->SetBranchAddress("Bu_ConsBu_status", Bu_ConsBu_status, &b_Bu_ConsBu_status);
   fChain->SetBranchAddress("Bu_MCJet_Px", &Bu_MCJet_Px, &b_Bu_MCJet_Px);
   fChain->SetBranchAddress("Bu_MCJet_Py", &Bu_MCJet_Py, &b_Bu_MCJet_Py);
   fChain->SetBranchAddress("Bu_MCJet_Pz", &Bu_MCJet_Pz, &b_Bu_MCJet_Pz);
   fChain->SetBranchAddress("Bu_MCJet_E", &Bu_MCJet_E, &b_Bu_MCJet_E);
   fChain->SetBranchAddress("Bu_MCJet_PT", &Bu_MCJet_PT, &b_Bu_MCJet_PT);
   fChain->SetBranchAddress("Bu_MCJet_Rapidity", &Bu_MCJet_Rapidity, &b_Bu_MCJet_Rapidity);
   fChain->SetBranchAddress("Bu_MCJet_FLAV", &Bu_MCJet_FLAV, &b_Bu_MCJet_FLAV);
   fChain->SetBranchAddress("Bu_MCHF_Px", &Bu_MCHF_Px, &b_Bu_MCHF_Px);
   fChain->SetBranchAddress("Bu_MCHF_Py", &Bu_MCHF_Py, &b_Bu_MCHF_Py);
   fChain->SetBranchAddress("Bu_MCHF_Pz", &Bu_MCHF_Pz, &b_Bu_MCHF_Pz);
   fChain->SetBranchAddress("Bu_MCHF_E", &Bu_MCHF_E, &b_Bu_MCHF_E);
   fChain->SetBranchAddress("Bu_MCHF_PT", &Bu_MCHF_PT, &b_Bu_MCHF_PT);
   fChain->SetBranchAddress("Bu_MCHF_Rapidity", &Bu_MCHF_Rapidity, &b_Bu_MCHF_Rapidity);
   fChain->SetBranchAddress("Bu_MCJet_NumHFHads", &Bu_MCJet_NumHFHads, &b_Bu_MCJet_NumHFHads);
   fChain->SetBranchAddress("Bu_NumBHads", &Bu_NumBHads, &b_Bu_NumBHads);
   fChain->SetBranchAddress("Bu_NumDHads", &Bu_NumDHads, &b_Bu_NumDHads);
   fChain->SetBranchAddress("Bu_hasc", &Bu_hasc, &b_Bu_hasc);
   fChain->SetBranchAddress("Bu_hascbar", &Bu_hascbar, &b_Bu_hascbar);
   fChain->SetBranchAddress("Bu_hasb", &Bu_hasb, &b_Bu_hasb);
   fChain->SetBranchAddress("Bu_hasbbar", &Bu_hasbbar, &b_Bu_hasbbar);
   fChain->SetBranchAddress("Bu_FromJet", &Bu_FromJet, &b_Bu_FromJet);
   fChain->SetBranchAddress("Bu_ENDVERTEX_X", &Bu_ENDVERTEX_X, &b_Bu_ENDVERTEX_X);
   fChain->SetBranchAddress("Bu_ENDVERTEX_Y", &Bu_ENDVERTEX_Y, &b_Bu_ENDVERTEX_Y);
   fChain->SetBranchAddress("Bu_ENDVERTEX_Z", &Bu_ENDVERTEX_Z, &b_Bu_ENDVERTEX_Z);
   fChain->SetBranchAddress("Bu_ENDVERTEX_XERR", &Bu_ENDVERTEX_XERR, &b_Bu_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Bu_ENDVERTEX_YERR", &Bu_ENDVERTEX_YERR, &b_Bu_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Bu_ENDVERTEX_ZERR", &Bu_ENDVERTEX_ZERR, &b_Bu_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Bu_ENDVERTEX_CHI2", &Bu_ENDVERTEX_CHI2, &b_Bu_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Bu_ENDVERTEX_NDOF", &Bu_ENDVERTEX_NDOF, &b_Bu_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Bu_ENDVERTEX_COV_", Bu_ENDVERTEX_COV_, &b_Bu_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Bu_OWNPV_X", &Bu_OWNPV_X, &b_Bu_OWNPV_X);
   fChain->SetBranchAddress("Bu_OWNPV_Y", &Bu_OWNPV_Y, &b_Bu_OWNPV_Y);
   fChain->SetBranchAddress("Bu_OWNPV_Z", &Bu_OWNPV_Z, &b_Bu_OWNPV_Z);
   fChain->SetBranchAddress("Bu_OWNPV_XERR", &Bu_OWNPV_XERR, &b_Bu_OWNPV_XERR);
   fChain->SetBranchAddress("Bu_OWNPV_YERR", &Bu_OWNPV_YERR, &b_Bu_OWNPV_YERR);
   fChain->SetBranchAddress("Bu_OWNPV_ZERR", &Bu_OWNPV_ZERR, &b_Bu_OWNPV_ZERR);
   fChain->SetBranchAddress("Bu_OWNPV_CHI2", &Bu_OWNPV_CHI2, &b_Bu_OWNPV_CHI2);
   fChain->SetBranchAddress("Bu_OWNPV_NDOF", &Bu_OWNPV_NDOF, &b_Bu_OWNPV_NDOF);
   fChain->SetBranchAddress("Bu_OWNPV_COV_", Bu_OWNPV_COV_, &b_Bu_OWNPV_COV_);
   fChain->SetBranchAddress("Bu_IP_OWNPV", &Bu_IP_OWNPV, &b_Bu_IP_OWNPV);
   fChain->SetBranchAddress("Bu_IPCHI2_OWNPV", &Bu_IPCHI2_OWNPV, &b_Bu_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Bu_FD_OWNPV", &Bu_FD_OWNPV, &b_Bu_FD_OWNPV);
   fChain->SetBranchAddress("Bu_FDCHI2_OWNPV", &Bu_FDCHI2_OWNPV, &b_Bu_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Bu_DIRA_OWNPV", &Bu_DIRA_OWNPV, &b_Bu_DIRA_OWNPV);
   fChain->SetBranchAddress("Bu_TRUEID", &Bu_TRUEID, &b_Bu_TRUEID);
   fChain->SetBranchAddress("Bu_TRUEP_E", &Bu_TRUEP_E, &b_Bu_TRUEP_E);
   fChain->SetBranchAddress("Bu_TRUEP_X", &Bu_TRUEP_X, &b_Bu_TRUEP_X);
   fChain->SetBranchAddress("Bu_TRUEP_Y", &Bu_TRUEP_Y, &b_Bu_TRUEP_Y);
   fChain->SetBranchAddress("Bu_TRUEP_Z", &Bu_TRUEP_Z, &b_Bu_TRUEP_Z);
   fChain->SetBranchAddress("Bu_TRUEPT", &Bu_TRUEPT, &b_Bu_TRUEPT);
   fChain->SetBranchAddress("Bu_TRUEORIGINVERTEX_X", &Bu_TRUEORIGINVERTEX_X, &b_Bu_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Bu_TRUEORIGINVERTEX_Y", &Bu_TRUEORIGINVERTEX_Y, &b_Bu_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Bu_TRUEORIGINVERTEX_Z", &Bu_TRUEORIGINVERTEX_Z, &b_Bu_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Bu_TRUEENDVERTEX_X", &Bu_TRUEENDVERTEX_X, &b_Bu_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Bu_TRUEENDVERTEX_Y", &Bu_TRUEENDVERTEX_Y, &b_Bu_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Bu_TRUEENDVERTEX_Z", &Bu_TRUEENDVERTEX_Z, &b_Bu_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Bu_TRUEISSTABLE", &Bu_TRUEISSTABLE, &b_Bu_TRUEISSTABLE);
   fChain->SetBranchAddress("Bu_TRUETAU", &Bu_TRUETAU, &b_Bu_TRUETAU);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_X", &Bplus0_ENDVERTEX_X, &b_Bplus0_ENDVERTEX_X);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_Y", &Bplus0_ENDVERTEX_Y, &b_Bplus0_ENDVERTEX_Y);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_Z", &Bplus0_ENDVERTEX_Z, &b_Bplus0_ENDVERTEX_Z);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_XERR", &Bplus0_ENDVERTEX_XERR, &b_Bplus0_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_YERR", &Bplus0_ENDVERTEX_YERR, &b_Bplus0_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_ZERR", &Bplus0_ENDVERTEX_ZERR, &b_Bplus0_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_CHI2", &Bplus0_ENDVERTEX_CHI2, &b_Bplus0_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_NDOF", &Bplus0_ENDVERTEX_NDOF, &b_Bplus0_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Bplus0_ENDVERTEX_COV_", Bplus0_ENDVERTEX_COV_, &b_Bplus0_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Bplus0_OWNPV_X", &Bplus0_OWNPV_X, &b_Bplus0_OWNPV_X);
   fChain->SetBranchAddress("Bplus0_OWNPV_Y", &Bplus0_OWNPV_Y, &b_Bplus0_OWNPV_Y);
   fChain->SetBranchAddress("Bplus0_OWNPV_Z", &Bplus0_OWNPV_Z, &b_Bplus0_OWNPV_Z);
   fChain->SetBranchAddress("Bplus0_OWNPV_XERR", &Bplus0_OWNPV_XERR, &b_Bplus0_OWNPV_XERR);
   fChain->SetBranchAddress("Bplus0_OWNPV_YERR", &Bplus0_OWNPV_YERR, &b_Bplus0_OWNPV_YERR);
   fChain->SetBranchAddress("Bplus0_OWNPV_ZERR", &Bplus0_OWNPV_ZERR, &b_Bplus0_OWNPV_ZERR);
   fChain->SetBranchAddress("Bplus0_OWNPV_CHI2", &Bplus0_OWNPV_CHI2, &b_Bplus0_OWNPV_CHI2);
   fChain->SetBranchAddress("Bplus0_OWNPV_NDOF", &Bplus0_OWNPV_NDOF, &b_Bplus0_OWNPV_NDOF);
   fChain->SetBranchAddress("Bplus0_OWNPV_COV_", Bplus0_OWNPV_COV_, &b_Bplus0_OWNPV_COV_);
   fChain->SetBranchAddress("Bplus0_IP_OWNPV", &Bplus0_IP_OWNPV, &b_Bplus0_IP_OWNPV);
   fChain->SetBranchAddress("Bplus0_IPCHI2_OWNPV", &Bplus0_IPCHI2_OWNPV, &b_Bplus0_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Bplus0_FD_OWNPV", &Bplus0_FD_OWNPV, &b_Bplus0_FD_OWNPV);
   fChain->SetBranchAddress("Bplus0_FDCHI2_OWNPV", &Bplus0_FDCHI2_OWNPV, &b_Bplus0_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Bplus0_DIRA_OWNPV", &Bplus0_DIRA_OWNPV, &b_Bplus0_DIRA_OWNPV);
   fChain->SetBranchAddress("Bplus0_P", &Bplus0_P, &b_Bplus0_P);
   fChain->SetBranchAddress("Bplus0_PT", &Bplus0_PT, &b_Bplus0_PT);
   fChain->SetBranchAddress("Bplus0_PE", &Bplus0_PE, &b_Bplus0_PE);
   fChain->SetBranchAddress("Bplus0_PX", &Bplus0_PX, &b_Bplus0_PX);
   fChain->SetBranchAddress("Bplus0_PY", &Bplus0_PY, &b_Bplus0_PY);
   fChain->SetBranchAddress("Bplus0_PZ", &Bplus0_PZ, &b_Bplus0_PZ);
   fChain->SetBranchAddress("Bplus0_MM", &Bplus0_MM, &b_Bplus0_MM);
   fChain->SetBranchAddress("Bplus0_MMERR", &Bplus0_MMERR, &b_Bplus0_MMERR);
   fChain->SetBranchAddress("Bplus0_M", &Bplus0_M, &b_Bplus0_M);
   fChain->SetBranchAddress("Bplus0_ID", &Bplus0_ID, &b_Bplus0_ID);
   fChain->SetBranchAddress("Bplus0_L0Global_Dec", &Bplus0_L0Global_Dec, &b_Bplus0_L0Global_Dec);
   fChain->SetBranchAddress("Bplus0_L0Global_TIS", &Bplus0_L0Global_TIS, &b_Bplus0_L0Global_TIS);
   fChain->SetBranchAddress("Bplus0_L0Global_TOS", &Bplus0_L0Global_TOS, &b_Bplus0_L0Global_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt1Global_Dec", &Bplus0_Hlt1Global_Dec, &b_Bplus0_Hlt1Global_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt1Global_TIS", &Bplus0_Hlt1Global_TIS, &b_Bplus0_Hlt1Global_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt1Global_TOS", &Bplus0_Hlt1Global_TOS, &b_Bplus0_Hlt1Global_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt1Phys_Dec", &Bplus0_Hlt1Phys_Dec, &b_Bplus0_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt1Phys_TIS", &Bplus0_Hlt1Phys_TIS, &b_Bplus0_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt1Phys_TOS", &Bplus0_Hlt1Phys_TOS, &b_Bplus0_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2Global_Dec", &Bplus0_Hlt2Global_Dec, &b_Bplus0_Hlt2Global_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2Global_TIS", &Bplus0_Hlt2Global_TIS, &b_Bplus0_Hlt2Global_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2Global_TOS", &Bplus0_Hlt2Global_TOS, &b_Bplus0_Hlt2Global_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2Phys_Dec", &Bplus0_Hlt2Phys_Dec, &b_Bplus0_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2Phys_TIS", &Bplus0_Hlt2Phys_TIS, &b_Bplus0_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2Phys_TOS", &Bplus0_Hlt2Phys_TOS, &b_Bplus0_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Bplus0_L0MuonDecision_Dec", &Bplus0_L0MuonDecision_Dec, &b_Bplus0_L0MuonDecision_Dec);
   fChain->SetBranchAddress("Bplus0_L0MuonDecision_TIS", &Bplus0_L0MuonDecision_TIS, &b_Bplus0_L0MuonDecision_TIS);
   fChain->SetBranchAddress("Bplus0_L0MuonDecision_TOS", &Bplus0_L0MuonDecision_TOS, &b_Bplus0_L0MuonDecision_TOS);
   fChain->SetBranchAddress("Bplus0_L0DiMuonDecision_Dec", &Bplus0_L0DiMuonDecision_Dec, &b_Bplus0_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("Bplus0_L0DiMuonDecision_TIS", &Bplus0_L0DiMuonDecision_TIS, &b_Bplus0_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("Bplus0_L0DiMuonDecision_TOS", &Bplus0_L0DiMuonDecision_TOS, &b_Bplus0_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt1DiMuonHighMassDecision_Dec", &Bplus0_Hlt1DiMuonHighMassDecision_Dec, &b_Bplus0_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt1DiMuonHighMassDecision_TIS", &Bplus0_Hlt1DiMuonHighMassDecision_TIS, &b_Bplus0_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt1DiMuonHighMassDecision_TOS", &Bplus0_Hlt1DiMuonHighMassDecision_TOS, &b_Bplus0_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedDecision_Dec", &Bplus0_Hlt2DiMuonDetachedDecision_Dec, &b_Bplus0_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedDecision_TIS", &Bplus0_Hlt2DiMuonDetachedDecision_TIS, &b_Bplus0_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedDecision_TOS", &Bplus0_Hlt2DiMuonDetachedDecision_TOS, &b_Bplus0_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedHeavyDecision_Dec", &Bplus0_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_Bplus0_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedHeavyDecision_TIS", &Bplus0_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_Bplus0_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedHeavyDecision_TOS", &Bplus0_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_Bplus0_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedJPsiDecision_Dec", &Bplus0_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_Bplus0_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedJPsiDecision_TIS", &Bplus0_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_Bplus0_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedJPsiDecision_TOS", &Bplus0_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_Bplus0_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedPsi2SDecision_Dec", &Bplus0_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_Bplus0_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TIS", &Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TOS", &Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_Bplus0_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonJPsiDecision_Dec", &Bplus0_Hlt2DiMuonJPsiDecision_Dec, &b_Bplus0_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonJPsiDecision_TIS", &Bplus0_Hlt2DiMuonJPsiDecision_TIS, &b_Bplus0_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonJPsiDecision_TOS", &Bplus0_Hlt2DiMuonJPsiDecision_TOS, &b_Bplus0_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonJPsiHighPTDecision_Dec", &Bplus0_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_Bplus0_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonJPsiHighPTDecision_TIS", &Bplus0_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_Bplus0_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("Bplus0_Hlt2DiMuonJPsiHighPTDecision_TOS", &Bplus0_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_Bplus0_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("Jpsi_P", &Jpsi_P, &b_Jpsi_P);
   fChain->SetBranchAddress("Jpsi_PT", &Jpsi_PT, &b_Jpsi_PT);
   fChain->SetBranchAddress("Jpsi_PE", &Jpsi_PE, &b_Jpsi_PE);
   fChain->SetBranchAddress("Jpsi_PX", &Jpsi_PX, &b_Jpsi_PX);
   fChain->SetBranchAddress("Jpsi_PY", &Jpsi_PY, &b_Jpsi_PY);
   fChain->SetBranchAddress("Jpsi_PZ", &Jpsi_PZ, &b_Jpsi_PZ);
   fChain->SetBranchAddress("Jpsi_MM", &Jpsi_MM, &b_Jpsi_MM);
   fChain->SetBranchAddress("Jpsi_MMERR", &Jpsi_MMERR, &b_Jpsi_MMERR);
   fChain->SetBranchAddress("Jpsi_M", &Jpsi_M, &b_Jpsi_M);
   fChain->SetBranchAddress("Jpsi_ID", &Jpsi_ID, &b_Jpsi_ID);
   fChain->SetBranchAddress("Jpsi_L0Global_Dec", &Jpsi_L0Global_Dec, &b_Jpsi_L0Global_Dec);
   fChain->SetBranchAddress("Jpsi_L0Global_TIS", &Jpsi_L0Global_TIS, &b_Jpsi_L0Global_TIS);
   fChain->SetBranchAddress("Jpsi_L0Global_TOS", &Jpsi_L0Global_TOS, &b_Jpsi_L0Global_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1Global_Dec", &Jpsi_Hlt1Global_Dec, &b_Jpsi_Hlt1Global_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1Global_TIS", &Jpsi_Hlt1Global_TIS, &b_Jpsi_Hlt1Global_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1Global_TOS", &Jpsi_Hlt1Global_TOS, &b_Jpsi_Hlt1Global_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1Phys_Dec", &Jpsi_Hlt1Phys_Dec, &b_Jpsi_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1Phys_TIS", &Jpsi_Hlt1Phys_TIS, &b_Jpsi_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1Phys_TOS", &Jpsi_Hlt1Phys_TOS, &b_Jpsi_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2Global_Dec", &Jpsi_Hlt2Global_Dec, &b_Jpsi_Hlt2Global_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2Global_TIS", &Jpsi_Hlt2Global_TIS, &b_Jpsi_Hlt2Global_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2Global_TOS", &Jpsi_Hlt2Global_TOS, &b_Jpsi_Hlt2Global_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2Phys_Dec", &Jpsi_Hlt2Phys_Dec, &b_Jpsi_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2Phys_TIS", &Jpsi_Hlt2Phys_TIS, &b_Jpsi_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2Phys_TOS", &Jpsi_Hlt2Phys_TOS, &b_Jpsi_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Jpsi_L0MuonDecision_Dec", &Jpsi_L0MuonDecision_Dec, &b_Jpsi_L0MuonDecision_Dec);
   fChain->SetBranchAddress("Jpsi_L0MuonDecision_TIS", &Jpsi_L0MuonDecision_TIS, &b_Jpsi_L0MuonDecision_TIS);
   fChain->SetBranchAddress("Jpsi_L0MuonDecision_TOS", &Jpsi_L0MuonDecision_TOS, &b_Jpsi_L0MuonDecision_TOS);
   fChain->SetBranchAddress("Jpsi_L0DiMuonDecision_Dec", &Jpsi_L0DiMuonDecision_Dec, &b_Jpsi_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("Jpsi_L0DiMuonDecision_TIS", &Jpsi_L0DiMuonDecision_TIS, &b_Jpsi_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("Jpsi_L0DiMuonDecision_TOS", &Jpsi_L0DiMuonDecision_TOS, &b_Jpsi_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_Dec", &Jpsi_Hlt1DiMuonHighMassDecision_Dec, &b_Jpsi_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_TIS", &Jpsi_Hlt1DiMuonHighMassDecision_TIS, &b_Jpsi_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_TOS", &Jpsi_Hlt1DiMuonHighMassDecision_TOS, &b_Jpsi_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedDecision_Dec", &Jpsi_Hlt2DiMuonDetachedDecision_Dec, &b_Jpsi_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedDecision_TIS", &Jpsi_Hlt2DiMuonDetachedDecision_TIS, &b_Jpsi_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedDecision_TOS", &Jpsi_Hlt2DiMuonDetachedDecision_TOS, &b_Jpsi_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec", &Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS", &Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS", &Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec", &Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS", &Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS", &Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedPsi2SDecision_Dec", &Jpsi_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_Jpsi_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TIS", &Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TOS", &Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_Jpsi_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonJPsiDecision_Dec", &Jpsi_Hlt2DiMuonJPsiDecision_Dec, &b_Jpsi_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonJPsiDecision_TIS", &Jpsi_Hlt2DiMuonJPsiDecision_TIS, &b_Jpsi_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonJPsiDecision_TOS", &Jpsi_Hlt2DiMuonJPsiDecision_TOS, &b_Jpsi_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonJPsiHighPTDecision_Dec", &Jpsi_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_Jpsi_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonJPsiHighPTDecision_TIS", &Jpsi_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_Jpsi_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS", &Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_Jpsi_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_X", &Jpsi_ENDVERTEX_X, &b_Jpsi_ENDVERTEX_X);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_Y", &Jpsi_ENDVERTEX_Y, &b_Jpsi_ENDVERTEX_Y);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_Z", &Jpsi_ENDVERTEX_Z, &b_Jpsi_ENDVERTEX_Z);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_XERR", &Jpsi_ENDVERTEX_XERR, &b_Jpsi_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_YERR", &Jpsi_ENDVERTEX_YERR, &b_Jpsi_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_ZERR", &Jpsi_ENDVERTEX_ZERR, &b_Jpsi_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_CHI2", &Jpsi_ENDVERTEX_CHI2, &b_Jpsi_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_NDOF", &Jpsi_ENDVERTEX_NDOF, &b_Jpsi_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_COV_", Jpsi_ENDVERTEX_COV_, &b_Jpsi_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Jpsi_OWNPV_X", &Jpsi_OWNPV_X, &b_Jpsi_OWNPV_X);
   fChain->SetBranchAddress("Jpsi_OWNPV_Y", &Jpsi_OWNPV_Y, &b_Jpsi_OWNPV_Y);
   fChain->SetBranchAddress("Jpsi_OWNPV_Z", &Jpsi_OWNPV_Z, &b_Jpsi_OWNPV_Z);
   fChain->SetBranchAddress("Jpsi_OWNPV_XERR", &Jpsi_OWNPV_XERR, &b_Jpsi_OWNPV_XERR);
   fChain->SetBranchAddress("Jpsi_OWNPV_YERR", &Jpsi_OWNPV_YERR, &b_Jpsi_OWNPV_YERR);
   fChain->SetBranchAddress("Jpsi_OWNPV_ZERR", &Jpsi_OWNPV_ZERR, &b_Jpsi_OWNPV_ZERR);
   fChain->SetBranchAddress("Jpsi_OWNPV_CHI2", &Jpsi_OWNPV_CHI2, &b_Jpsi_OWNPV_CHI2);
   fChain->SetBranchAddress("Jpsi_OWNPV_NDOF", &Jpsi_OWNPV_NDOF, &b_Jpsi_OWNPV_NDOF);
   fChain->SetBranchAddress("Jpsi_OWNPV_COV_", Jpsi_OWNPV_COV_, &b_Jpsi_OWNPV_COV_);
   fChain->SetBranchAddress("Jpsi_IP_OWNPV", &Jpsi_IP_OWNPV, &b_Jpsi_IP_OWNPV);
   fChain->SetBranchAddress("Jpsi_IPCHI2_OWNPV", &Jpsi_IPCHI2_OWNPV, &b_Jpsi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Jpsi_FD_OWNPV", &Jpsi_FD_OWNPV, &b_Jpsi_FD_OWNPV);
   fChain->SetBranchAddress("Jpsi_FDCHI2_OWNPV", &Jpsi_FDCHI2_OWNPV, &b_Jpsi_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Jpsi_DIRA_OWNPV", &Jpsi_DIRA_OWNPV, &b_Jpsi_DIRA_OWNPV);
   fChain->SetBranchAddress("Jpsi_ORIVX_X", &Jpsi_ORIVX_X, &b_Jpsi_ORIVX_X);
   fChain->SetBranchAddress("Jpsi_ORIVX_Y", &Jpsi_ORIVX_Y, &b_Jpsi_ORIVX_Y);
   fChain->SetBranchAddress("Jpsi_ORIVX_Z", &Jpsi_ORIVX_Z, &b_Jpsi_ORIVX_Z);
   fChain->SetBranchAddress("Jpsi_ORIVX_XERR", &Jpsi_ORIVX_XERR, &b_Jpsi_ORIVX_XERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_YERR", &Jpsi_ORIVX_YERR, &b_Jpsi_ORIVX_YERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_ZERR", &Jpsi_ORIVX_ZERR, &b_Jpsi_ORIVX_ZERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_CHI2", &Jpsi_ORIVX_CHI2, &b_Jpsi_ORIVX_CHI2);
   fChain->SetBranchAddress("Jpsi_ORIVX_NDOF", &Jpsi_ORIVX_NDOF, &b_Jpsi_ORIVX_NDOF);
   fChain->SetBranchAddress("Jpsi_ORIVX_COV_", Jpsi_ORIVX_COV_, &b_Jpsi_ORIVX_COV_);
   fChain->SetBranchAddress("Jpsi_FD_ORIVX", &Jpsi_FD_ORIVX, &b_Jpsi_FD_ORIVX);
   fChain->SetBranchAddress("Jpsi_FDCHI2_ORIVX", &Jpsi_FDCHI2_ORIVX, &b_Jpsi_FDCHI2_ORIVX);
   fChain->SetBranchAddress("Jpsi_DIRA_ORIVX", &Jpsi_DIRA_ORIVX, &b_Jpsi_DIRA_ORIVX);
   fChain->SetBranchAddress("Jpsi_TRUEID", &Jpsi_TRUEID, &b_Jpsi_TRUEID);
   fChain->SetBranchAddress("Jpsi_TRUEP_E", &Jpsi_TRUEP_E, &b_Jpsi_TRUEP_E);
   fChain->SetBranchAddress("Jpsi_TRUEP_X", &Jpsi_TRUEP_X, &b_Jpsi_TRUEP_X);
   fChain->SetBranchAddress("Jpsi_TRUEP_Y", &Jpsi_TRUEP_Y, &b_Jpsi_TRUEP_Y);
   fChain->SetBranchAddress("Jpsi_TRUEP_Z", &Jpsi_TRUEP_Z, &b_Jpsi_TRUEP_Z);
   fChain->SetBranchAddress("Jpsi_TRUEPT", &Jpsi_TRUEPT, &b_Jpsi_TRUEPT);
   fChain->SetBranchAddress("Jpsi_TRUEORIGINVERTEX_X", &Jpsi_TRUEORIGINVERTEX_X, &b_Jpsi_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Jpsi_TRUEORIGINVERTEX_Y", &Jpsi_TRUEORIGINVERTEX_Y, &b_Jpsi_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Jpsi_TRUEORIGINVERTEX_Z", &Jpsi_TRUEORIGINVERTEX_Z, &b_Jpsi_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Jpsi_TRUEENDVERTEX_X", &Jpsi_TRUEENDVERTEX_X, &b_Jpsi_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Jpsi_TRUEENDVERTEX_Y", &Jpsi_TRUEENDVERTEX_Y, &b_Jpsi_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Jpsi_TRUEENDVERTEX_Z", &Jpsi_TRUEENDVERTEX_Z, &b_Jpsi_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Jpsi_TRUEISSTABLE", &Jpsi_TRUEISSTABLE, &b_Jpsi_TRUEISSTABLE);
   fChain->SetBranchAddress("Jpsi_TRUETAU", &Jpsi_TRUETAU, &b_Jpsi_TRUETAU);
   fChain->SetBranchAddress("mup_MC12TuneV2_ProbNNe", &mup_MC12TuneV2_ProbNNe, &b_mup_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("mup_MC12TuneV2_ProbNNmu", &mup_MC12TuneV2_ProbNNmu, &b_mup_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("mup_MC12TuneV2_ProbNNpi", &mup_MC12TuneV2_ProbNNpi, &b_mup_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("mup_MC12TuneV2_ProbNNk", &mup_MC12TuneV2_ProbNNk, &b_mup_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("mup_MC12TuneV2_ProbNNp", &mup_MC12TuneV2_ProbNNp, &b_mup_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("mup_MC12TuneV2_ProbNNghost", &mup_MC12TuneV2_ProbNNghost, &b_mup_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("mup_MC12TuneV3_ProbNNe", &mup_MC12TuneV3_ProbNNe, &b_mup_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("mup_MC12TuneV3_ProbNNmu", &mup_MC12TuneV3_ProbNNmu, &b_mup_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("mup_MC12TuneV3_ProbNNpi", &mup_MC12TuneV3_ProbNNpi, &b_mup_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("mup_MC12TuneV3_ProbNNk", &mup_MC12TuneV3_ProbNNk, &b_mup_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("mup_MC12TuneV3_ProbNNp", &mup_MC12TuneV3_ProbNNp, &b_mup_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("mup_MC12TuneV3_ProbNNghost", &mup_MC12TuneV3_ProbNNghost, &b_mup_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("mup_MC12TuneV4_ProbNNe", &mup_MC12TuneV4_ProbNNe, &b_mup_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("mup_MC12TuneV4_ProbNNmu", &mup_MC12TuneV4_ProbNNmu, &b_mup_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("mup_MC12TuneV4_ProbNNpi", &mup_MC12TuneV4_ProbNNpi, &b_mup_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("mup_MC12TuneV4_ProbNNk", &mup_MC12TuneV4_ProbNNk, &b_mup_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("mup_MC12TuneV4_ProbNNp", &mup_MC12TuneV4_ProbNNp, &b_mup_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("mup_MC12TuneV4_ProbNNghost", &mup_MC12TuneV4_ProbNNghost, &b_mup_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("mup_MC15TuneV1_ProbNNe", &mup_MC15TuneV1_ProbNNe, &b_mup_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("mup_MC15TuneV1_ProbNNmu", &mup_MC15TuneV1_ProbNNmu, &b_mup_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("mup_MC15TuneV1_ProbNNpi", &mup_MC15TuneV1_ProbNNpi, &b_mup_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("mup_MC15TuneV1_ProbNNk", &mup_MC15TuneV1_ProbNNk, &b_mup_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("mup_MC15TuneV1_ProbNNp", &mup_MC15TuneV1_ProbNNp, &b_mup_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("mup_MC15TuneV1_ProbNNghost", &mup_MC15TuneV1_ProbNNghost, &b_mup_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("mup_P", &mup_P, &b_mup_P);
   fChain->SetBranchAddress("mup_PT", &mup_PT, &b_mup_PT);
   fChain->SetBranchAddress("mup_PE", &mup_PE, &b_mup_PE);
   fChain->SetBranchAddress("mup_PX", &mup_PX, &b_mup_PX);
   fChain->SetBranchAddress("mup_PY", &mup_PY, &b_mup_PY);
   fChain->SetBranchAddress("mup_PZ", &mup_PZ, &b_mup_PZ);
   fChain->SetBranchAddress("mup_M", &mup_M, &b_mup_M);
   fChain->SetBranchAddress("mup_ID", &mup_ID, &b_mup_ID);
   fChain->SetBranchAddress("mup_PIDe", &mup_PIDe, &b_mup_PIDe);
   fChain->SetBranchAddress("mup_PIDmu", &mup_PIDmu, &b_mup_PIDmu);
   fChain->SetBranchAddress("mup_PIDK", &mup_PIDK, &b_mup_PIDK);
   fChain->SetBranchAddress("mup_PIDp", &mup_PIDp, &b_mup_PIDp);
   fChain->SetBranchAddress("mup_PIDd", &mup_PIDd, &b_mup_PIDd);
   fChain->SetBranchAddress("mup_ProbNNe", &mup_ProbNNe, &b_mup_ProbNNe);
   fChain->SetBranchAddress("mup_ProbNNk", &mup_ProbNNk, &b_mup_ProbNNk);
   fChain->SetBranchAddress("mup_ProbNNp", &mup_ProbNNp, &b_mup_ProbNNp);
   fChain->SetBranchAddress("mup_ProbNNpi", &mup_ProbNNpi, &b_mup_ProbNNpi);
   fChain->SetBranchAddress("mup_ProbNNmu", &mup_ProbNNmu, &b_mup_ProbNNmu);
   fChain->SetBranchAddress("mup_ProbNNd", &mup_ProbNNd, &b_mup_ProbNNd);
   fChain->SetBranchAddress("mup_ProbNNghost", &mup_ProbNNghost, &b_mup_ProbNNghost);
   fChain->SetBranchAddress("mup_hasMuon", &mup_hasMuon, &b_mup_hasMuon);
   fChain->SetBranchAddress("mup_isMuon", &mup_isMuon, &b_mup_isMuon);
   fChain->SetBranchAddress("mup_hasRich", &mup_hasRich, &b_mup_hasRich);
   fChain->SetBranchAddress("mup_UsedRichAerogel", &mup_UsedRichAerogel, &b_mup_UsedRichAerogel);
   fChain->SetBranchAddress("mup_UsedRich1Gas", &mup_UsedRich1Gas, &b_mup_UsedRich1Gas);
   fChain->SetBranchAddress("mup_UsedRich2Gas", &mup_UsedRich2Gas, &b_mup_UsedRich2Gas);
   fChain->SetBranchAddress("mup_RichAboveElThres", &mup_RichAboveElThres, &b_mup_RichAboveElThres);
   fChain->SetBranchAddress("mup_RichAboveMuThres", &mup_RichAboveMuThres, &b_mup_RichAboveMuThres);
   fChain->SetBranchAddress("mup_RichAbovePiThres", &mup_RichAbovePiThres, &b_mup_RichAbovePiThres);
   fChain->SetBranchAddress("mup_RichAboveKaThres", &mup_RichAboveKaThres, &b_mup_RichAboveKaThres);
   fChain->SetBranchAddress("mup_RichAbovePrThres", &mup_RichAbovePrThres, &b_mup_RichAbovePrThres);
   fChain->SetBranchAddress("mup_hasCalo", &mup_hasCalo, &b_mup_hasCalo);
   fChain->SetBranchAddress("mup_L0Global_Dec", &mup_L0Global_Dec, &b_mup_L0Global_Dec);
   fChain->SetBranchAddress("mup_L0Global_TIS", &mup_L0Global_TIS, &b_mup_L0Global_TIS);
   fChain->SetBranchAddress("mup_L0Global_TOS", &mup_L0Global_TOS, &b_mup_L0Global_TOS);
   fChain->SetBranchAddress("mup_Hlt1Global_Dec", &mup_Hlt1Global_Dec, &b_mup_Hlt1Global_Dec);
   fChain->SetBranchAddress("mup_Hlt1Global_TIS", &mup_Hlt1Global_TIS, &b_mup_Hlt1Global_TIS);
   fChain->SetBranchAddress("mup_Hlt1Global_TOS", &mup_Hlt1Global_TOS, &b_mup_Hlt1Global_TOS);
   fChain->SetBranchAddress("mup_Hlt1Phys_Dec", &mup_Hlt1Phys_Dec, &b_mup_Hlt1Phys_Dec);
   fChain->SetBranchAddress("mup_Hlt1Phys_TIS", &mup_Hlt1Phys_TIS, &b_mup_Hlt1Phys_TIS);
   fChain->SetBranchAddress("mup_Hlt1Phys_TOS", &mup_Hlt1Phys_TOS, &b_mup_Hlt1Phys_TOS);
   fChain->SetBranchAddress("mup_Hlt2Global_Dec", &mup_Hlt2Global_Dec, &b_mup_Hlt2Global_Dec);
   fChain->SetBranchAddress("mup_Hlt2Global_TIS", &mup_Hlt2Global_TIS, &b_mup_Hlt2Global_TIS);
   fChain->SetBranchAddress("mup_Hlt2Global_TOS", &mup_Hlt2Global_TOS, &b_mup_Hlt2Global_TOS);
   fChain->SetBranchAddress("mup_Hlt2Phys_Dec", &mup_Hlt2Phys_Dec, &b_mup_Hlt2Phys_Dec);
   fChain->SetBranchAddress("mup_Hlt2Phys_TIS", &mup_Hlt2Phys_TIS, &b_mup_Hlt2Phys_TIS);
   fChain->SetBranchAddress("mup_Hlt2Phys_TOS", &mup_Hlt2Phys_TOS, &b_mup_Hlt2Phys_TOS);
   fChain->SetBranchAddress("mup_L0MuonDecision_Dec", &mup_L0MuonDecision_Dec, &b_mup_L0MuonDecision_Dec);
   fChain->SetBranchAddress("mup_L0MuonDecision_TIS", &mup_L0MuonDecision_TIS, &b_mup_L0MuonDecision_TIS);
   fChain->SetBranchAddress("mup_L0MuonDecision_TOS", &mup_L0MuonDecision_TOS, &b_mup_L0MuonDecision_TOS);
   fChain->SetBranchAddress("mup_L0DiMuonDecision_Dec", &mup_L0DiMuonDecision_Dec, &b_mup_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("mup_L0DiMuonDecision_TIS", &mup_L0DiMuonDecision_TIS, &b_mup_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("mup_L0DiMuonDecision_TOS", &mup_L0DiMuonDecision_TOS, &b_mup_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt1DiMuonHighMassDecision_Dec", &mup_Hlt1DiMuonHighMassDecision_Dec, &b_mup_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt1DiMuonHighMassDecision_TIS", &mup_Hlt1DiMuonHighMassDecision_TIS, &b_mup_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt1DiMuonHighMassDecision_TOS", &mup_Hlt1DiMuonHighMassDecision_TOS, &b_mup_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedDecision_Dec", &mup_Hlt2DiMuonDetachedDecision_Dec, &b_mup_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedDecision_TIS", &mup_Hlt2DiMuonDetachedDecision_TIS, &b_mup_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedDecision_TOS", &mup_Hlt2DiMuonDetachedDecision_TOS, &b_mup_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedHeavyDecision_Dec", &mup_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_mup_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedHeavyDecision_TIS", &mup_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_mup_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedHeavyDecision_TOS", &mup_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_mup_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedJPsiDecision_Dec", &mup_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_mup_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedJPsiDecision_TIS", &mup_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_mup_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedJPsiDecision_TOS", &mup_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_mup_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedPsi2SDecision_Dec", &mup_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_mup_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedPsi2SDecision_TIS", &mup_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_mup_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonDetachedPsi2SDecision_TOS", &mup_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_mup_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonJPsiDecision_Dec", &mup_Hlt2DiMuonJPsiDecision_Dec, &b_mup_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt2DiMuonJPsiDecision_TIS", &mup_Hlt2DiMuonJPsiDecision_TIS, &b_mup_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonJPsiDecision_TOS", &mup_Hlt2DiMuonJPsiDecision_TOS, &b_mup_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonJPsiHighPTDecision_Dec", &mup_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_mup_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("mup_Hlt2DiMuonJPsiHighPTDecision_TIS", &mup_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_mup_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("mup_Hlt2DiMuonJPsiHighPTDecision_TOS", &mup_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_mup_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("mup_OWNPV_X", &mup_OWNPV_X, &b_mup_OWNPV_X);
   fChain->SetBranchAddress("mup_OWNPV_Y", &mup_OWNPV_Y, &b_mup_OWNPV_Y);
   fChain->SetBranchAddress("mup_OWNPV_Z", &mup_OWNPV_Z, &b_mup_OWNPV_Z);
   fChain->SetBranchAddress("mup_OWNPV_XERR", &mup_OWNPV_XERR, &b_mup_OWNPV_XERR);
   fChain->SetBranchAddress("mup_OWNPV_YERR", &mup_OWNPV_YERR, &b_mup_OWNPV_YERR);
   fChain->SetBranchAddress("mup_OWNPV_ZERR", &mup_OWNPV_ZERR, &b_mup_OWNPV_ZERR);
   fChain->SetBranchAddress("mup_OWNPV_CHI2", &mup_OWNPV_CHI2, &b_mup_OWNPV_CHI2);
   fChain->SetBranchAddress("mup_OWNPV_NDOF", &mup_OWNPV_NDOF, &b_mup_OWNPV_NDOF);
   fChain->SetBranchAddress("mup_OWNPV_COV_", mup_OWNPV_COV_, &b_mup_OWNPV_COV_);
   fChain->SetBranchAddress("mup_IP_OWNPV", &mup_IP_OWNPV, &b_mup_IP_OWNPV);
   fChain->SetBranchAddress("mup_IPCHI2_OWNPV", &mup_IPCHI2_OWNPV, &b_mup_IPCHI2_OWNPV);
   fChain->SetBranchAddress("mup_ORIVX_X", &mup_ORIVX_X, &b_mup_ORIVX_X);
   fChain->SetBranchAddress("mup_ORIVX_Y", &mup_ORIVX_Y, &b_mup_ORIVX_Y);
   fChain->SetBranchAddress("mup_ORIVX_Z", &mup_ORIVX_Z, &b_mup_ORIVX_Z);
   fChain->SetBranchAddress("mup_ORIVX_XERR", &mup_ORIVX_XERR, &b_mup_ORIVX_XERR);
   fChain->SetBranchAddress("mup_ORIVX_YERR", &mup_ORIVX_YERR, &b_mup_ORIVX_YERR);
   fChain->SetBranchAddress("mup_ORIVX_ZERR", &mup_ORIVX_ZERR, &b_mup_ORIVX_ZERR);
   fChain->SetBranchAddress("mup_ORIVX_CHI2", &mup_ORIVX_CHI2, &b_mup_ORIVX_CHI2);
   fChain->SetBranchAddress("mup_ORIVX_NDOF", &mup_ORIVX_NDOF, &b_mup_ORIVX_NDOF);
   fChain->SetBranchAddress("mup_ORIVX_COV_", mup_ORIVX_COV_, &b_mup_ORIVX_COV_);
   fChain->SetBranchAddress("mup_TRUEID", &mup_TRUEID, &b_mup_TRUEID);
   fChain->SetBranchAddress("mup_TRUEP_E", &mup_TRUEP_E, &b_mup_TRUEP_E);
   fChain->SetBranchAddress("mup_TRUEP_X", &mup_TRUEP_X, &b_mup_TRUEP_X);
   fChain->SetBranchAddress("mup_TRUEP_Y", &mup_TRUEP_Y, &b_mup_TRUEP_Y);
   fChain->SetBranchAddress("mup_TRUEP_Z", &mup_TRUEP_Z, &b_mup_TRUEP_Z);
   fChain->SetBranchAddress("mup_TRUEPT", &mup_TRUEPT, &b_mup_TRUEPT);
   fChain->SetBranchAddress("mup_TRUEORIGINVERTEX_X", &mup_TRUEORIGINVERTEX_X, &b_mup_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("mup_TRUEORIGINVERTEX_Y", &mup_TRUEORIGINVERTEX_Y, &b_mup_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("mup_TRUEORIGINVERTEX_Z", &mup_TRUEORIGINVERTEX_Z, &b_mup_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("mup_TRUEENDVERTEX_X", &mup_TRUEENDVERTEX_X, &b_mup_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("mup_TRUEENDVERTEX_Y", &mup_TRUEENDVERTEX_Y, &b_mup_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("mup_TRUEENDVERTEX_Z", &mup_TRUEENDVERTEX_Z, &b_mup_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("mup_TRUEISSTABLE", &mup_TRUEISSTABLE, &b_mup_TRUEISSTABLE);
   fChain->SetBranchAddress("mup_TRUETAU", &mup_TRUETAU, &b_mup_TRUETAU);
   fChain->SetBranchAddress("mup_TRACK_Type", &mup_TRACK_Type, &b_mup_TRACK_Type);
   fChain->SetBranchAddress("mup_TRACK_Key", &mup_TRACK_Key, &b_mup_TRACK_Key);
   fChain->SetBranchAddress("mup_TRACK_CHI2NDOF", &mup_TRACK_CHI2NDOF, &b_mup_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("mup_TRACK_PCHI2", &mup_TRACK_PCHI2, &b_mup_TRACK_PCHI2);
   fChain->SetBranchAddress("mup_TRACK_MatchCHI2", &mup_TRACK_MatchCHI2, &b_mup_TRACK_MatchCHI2);
   fChain->SetBranchAddress("mup_TRACK_GhostProb", &mup_TRACK_GhostProb, &b_mup_TRACK_GhostProb);
   fChain->SetBranchAddress("mup_TRACK_CloneDist", &mup_TRACK_CloneDist, &b_mup_TRACK_CloneDist);
   fChain->SetBranchAddress("mup_TRACK_Likelihood", &mup_TRACK_Likelihood, &b_mup_TRACK_Likelihood);
   fChain->SetBranchAddress("mum_MC12TuneV2_ProbNNe", &mum_MC12TuneV2_ProbNNe, &b_mum_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("mum_MC12TuneV2_ProbNNmu", &mum_MC12TuneV2_ProbNNmu, &b_mum_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("mum_MC12TuneV2_ProbNNpi", &mum_MC12TuneV2_ProbNNpi, &b_mum_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("mum_MC12TuneV2_ProbNNk", &mum_MC12TuneV2_ProbNNk, &b_mum_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("mum_MC12TuneV2_ProbNNp", &mum_MC12TuneV2_ProbNNp, &b_mum_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("mum_MC12TuneV2_ProbNNghost", &mum_MC12TuneV2_ProbNNghost, &b_mum_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("mum_MC12TuneV3_ProbNNe", &mum_MC12TuneV3_ProbNNe, &b_mum_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("mum_MC12TuneV3_ProbNNmu", &mum_MC12TuneV3_ProbNNmu, &b_mum_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("mum_MC12TuneV3_ProbNNpi", &mum_MC12TuneV3_ProbNNpi, &b_mum_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("mum_MC12TuneV3_ProbNNk", &mum_MC12TuneV3_ProbNNk, &b_mum_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("mum_MC12TuneV3_ProbNNp", &mum_MC12TuneV3_ProbNNp, &b_mum_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("mum_MC12TuneV3_ProbNNghost", &mum_MC12TuneV3_ProbNNghost, &b_mum_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("mum_MC12TuneV4_ProbNNe", &mum_MC12TuneV4_ProbNNe, &b_mum_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("mum_MC12TuneV4_ProbNNmu", &mum_MC12TuneV4_ProbNNmu, &b_mum_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("mum_MC12TuneV4_ProbNNpi", &mum_MC12TuneV4_ProbNNpi, &b_mum_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("mum_MC12TuneV4_ProbNNk", &mum_MC12TuneV4_ProbNNk, &b_mum_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("mum_MC12TuneV4_ProbNNp", &mum_MC12TuneV4_ProbNNp, &b_mum_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("mum_MC12TuneV4_ProbNNghost", &mum_MC12TuneV4_ProbNNghost, &b_mum_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("mum_MC15TuneV1_ProbNNe", &mum_MC15TuneV1_ProbNNe, &b_mum_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("mum_MC15TuneV1_ProbNNmu", &mum_MC15TuneV1_ProbNNmu, &b_mum_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("mum_MC15TuneV1_ProbNNpi", &mum_MC15TuneV1_ProbNNpi, &b_mum_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("mum_MC15TuneV1_ProbNNk", &mum_MC15TuneV1_ProbNNk, &b_mum_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("mum_MC15TuneV1_ProbNNp", &mum_MC15TuneV1_ProbNNp, &b_mum_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("mum_MC15TuneV1_ProbNNghost", &mum_MC15TuneV1_ProbNNghost, &b_mum_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("mum_P", &mum_P, &b_mum_P);
   fChain->SetBranchAddress("mum_PT", &mum_PT, &b_mum_PT);
   fChain->SetBranchAddress("mum_PE", &mum_PE, &b_mum_PE);
   fChain->SetBranchAddress("mum_PX", &mum_PX, &b_mum_PX);
   fChain->SetBranchAddress("mum_PY", &mum_PY, &b_mum_PY);
   fChain->SetBranchAddress("mum_PZ", &mum_PZ, &b_mum_PZ);
   fChain->SetBranchAddress("mum_M", &mum_M, &b_mum_M);
   fChain->SetBranchAddress("mum_ID", &mum_ID, &b_mum_ID);
   fChain->SetBranchAddress("mum_PIDe", &mum_PIDe, &b_mum_PIDe);
   fChain->SetBranchAddress("mum_PIDmu", &mum_PIDmu, &b_mum_PIDmu);
   fChain->SetBranchAddress("mum_PIDK", &mum_PIDK, &b_mum_PIDK);
   fChain->SetBranchAddress("mum_PIDp", &mum_PIDp, &b_mum_PIDp);
   fChain->SetBranchAddress("mum_PIDd", &mum_PIDd, &b_mum_PIDd);
   fChain->SetBranchAddress("mum_ProbNNe", &mum_ProbNNe, &b_mum_ProbNNe);
   fChain->SetBranchAddress("mum_ProbNNk", &mum_ProbNNk, &b_mum_ProbNNk);
   fChain->SetBranchAddress("mum_ProbNNp", &mum_ProbNNp, &b_mum_ProbNNp);
   fChain->SetBranchAddress("mum_ProbNNpi", &mum_ProbNNpi, &b_mum_ProbNNpi);
   fChain->SetBranchAddress("mum_ProbNNmu", &mum_ProbNNmu, &b_mum_ProbNNmu);
   fChain->SetBranchAddress("mum_ProbNNd", &mum_ProbNNd, &b_mum_ProbNNd);
   fChain->SetBranchAddress("mum_ProbNNghost", &mum_ProbNNghost, &b_mum_ProbNNghost);
   fChain->SetBranchAddress("mum_hasMuon", &mum_hasMuon, &b_mum_hasMuon);
   fChain->SetBranchAddress("mum_isMuon", &mum_isMuon, &b_mum_isMuon);
   fChain->SetBranchAddress("mum_hasRich", &mum_hasRich, &b_mum_hasRich);
   fChain->SetBranchAddress("mum_UsedRichAerogel", &mum_UsedRichAerogel, &b_mum_UsedRichAerogel);
   fChain->SetBranchAddress("mum_UsedRich1Gas", &mum_UsedRich1Gas, &b_mum_UsedRich1Gas);
   fChain->SetBranchAddress("mum_UsedRich2Gas", &mum_UsedRich2Gas, &b_mum_UsedRich2Gas);
   fChain->SetBranchAddress("mum_RichAboveElThres", &mum_RichAboveElThres, &b_mum_RichAboveElThres);
   fChain->SetBranchAddress("mum_RichAboveMuThres", &mum_RichAboveMuThres, &b_mum_RichAboveMuThres);
   fChain->SetBranchAddress("mum_RichAbovePiThres", &mum_RichAbovePiThres, &b_mum_RichAbovePiThres);
   fChain->SetBranchAddress("mum_RichAboveKaThres", &mum_RichAboveKaThres, &b_mum_RichAboveKaThres);
   fChain->SetBranchAddress("mum_RichAbovePrThres", &mum_RichAbovePrThres, &b_mum_RichAbovePrThres);
   fChain->SetBranchAddress("mum_hasCalo", &mum_hasCalo, &b_mum_hasCalo);
   fChain->SetBranchAddress("mum_L0Global_Dec", &mum_L0Global_Dec, &b_mum_L0Global_Dec);
   fChain->SetBranchAddress("mum_L0Global_TIS", &mum_L0Global_TIS, &b_mum_L0Global_TIS);
   fChain->SetBranchAddress("mum_L0Global_TOS", &mum_L0Global_TOS, &b_mum_L0Global_TOS);
   fChain->SetBranchAddress("mum_Hlt1Global_Dec", &mum_Hlt1Global_Dec, &b_mum_Hlt1Global_Dec);
   fChain->SetBranchAddress("mum_Hlt1Global_TIS", &mum_Hlt1Global_TIS, &b_mum_Hlt1Global_TIS);
   fChain->SetBranchAddress("mum_Hlt1Global_TOS", &mum_Hlt1Global_TOS, &b_mum_Hlt1Global_TOS);
   fChain->SetBranchAddress("mum_Hlt1Phys_Dec", &mum_Hlt1Phys_Dec, &b_mum_Hlt1Phys_Dec);
   fChain->SetBranchAddress("mum_Hlt1Phys_TIS", &mum_Hlt1Phys_TIS, &b_mum_Hlt1Phys_TIS);
   fChain->SetBranchAddress("mum_Hlt1Phys_TOS", &mum_Hlt1Phys_TOS, &b_mum_Hlt1Phys_TOS);
   fChain->SetBranchAddress("mum_Hlt2Global_Dec", &mum_Hlt2Global_Dec, &b_mum_Hlt2Global_Dec);
   fChain->SetBranchAddress("mum_Hlt2Global_TIS", &mum_Hlt2Global_TIS, &b_mum_Hlt2Global_TIS);
   fChain->SetBranchAddress("mum_Hlt2Global_TOS", &mum_Hlt2Global_TOS, &b_mum_Hlt2Global_TOS);
   fChain->SetBranchAddress("mum_Hlt2Phys_Dec", &mum_Hlt2Phys_Dec, &b_mum_Hlt2Phys_Dec);
   fChain->SetBranchAddress("mum_Hlt2Phys_TIS", &mum_Hlt2Phys_TIS, &b_mum_Hlt2Phys_TIS);
   fChain->SetBranchAddress("mum_Hlt2Phys_TOS", &mum_Hlt2Phys_TOS, &b_mum_Hlt2Phys_TOS);
   fChain->SetBranchAddress("mum_L0MuonDecision_Dec", &mum_L0MuonDecision_Dec, &b_mum_L0MuonDecision_Dec);
   fChain->SetBranchAddress("mum_L0MuonDecision_TIS", &mum_L0MuonDecision_TIS, &b_mum_L0MuonDecision_TIS);
   fChain->SetBranchAddress("mum_L0MuonDecision_TOS", &mum_L0MuonDecision_TOS, &b_mum_L0MuonDecision_TOS);
   fChain->SetBranchAddress("mum_L0DiMuonDecision_Dec", &mum_L0DiMuonDecision_Dec, &b_mum_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("mum_L0DiMuonDecision_TIS", &mum_L0DiMuonDecision_TIS, &b_mum_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("mum_L0DiMuonDecision_TOS", &mum_L0DiMuonDecision_TOS, &b_mum_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt1DiMuonHighMassDecision_Dec", &mum_Hlt1DiMuonHighMassDecision_Dec, &b_mum_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt1DiMuonHighMassDecision_TIS", &mum_Hlt1DiMuonHighMassDecision_TIS, &b_mum_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt1DiMuonHighMassDecision_TOS", &mum_Hlt1DiMuonHighMassDecision_TOS, &b_mum_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedDecision_Dec", &mum_Hlt2DiMuonDetachedDecision_Dec, &b_mum_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedDecision_TIS", &mum_Hlt2DiMuonDetachedDecision_TIS, &b_mum_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedDecision_TOS", &mum_Hlt2DiMuonDetachedDecision_TOS, &b_mum_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedHeavyDecision_Dec", &mum_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_mum_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedHeavyDecision_TIS", &mum_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_mum_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedHeavyDecision_TOS", &mum_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_mum_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedJPsiDecision_Dec", &mum_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_mum_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedJPsiDecision_TIS", &mum_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_mum_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedJPsiDecision_TOS", &mum_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_mum_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedPsi2SDecision_Dec", &mum_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_mum_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedPsi2SDecision_TIS", &mum_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_mum_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonDetachedPsi2SDecision_TOS", &mum_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_mum_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonJPsiDecision_Dec", &mum_Hlt2DiMuonJPsiDecision_Dec, &b_mum_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt2DiMuonJPsiDecision_TIS", &mum_Hlt2DiMuonJPsiDecision_TIS, &b_mum_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonJPsiDecision_TOS", &mum_Hlt2DiMuonJPsiDecision_TOS, &b_mum_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonJPsiHighPTDecision_Dec", &mum_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_mum_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("mum_Hlt2DiMuonJPsiHighPTDecision_TIS", &mum_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_mum_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("mum_Hlt2DiMuonJPsiHighPTDecision_TOS", &mum_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_mum_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("mum_OWNPV_X", &mum_OWNPV_X, &b_mum_OWNPV_X);
   fChain->SetBranchAddress("mum_OWNPV_Y", &mum_OWNPV_Y, &b_mum_OWNPV_Y);
   fChain->SetBranchAddress("mum_OWNPV_Z", &mum_OWNPV_Z, &b_mum_OWNPV_Z);
   fChain->SetBranchAddress("mum_OWNPV_XERR", &mum_OWNPV_XERR, &b_mum_OWNPV_XERR);
   fChain->SetBranchAddress("mum_OWNPV_YERR", &mum_OWNPV_YERR, &b_mum_OWNPV_YERR);
   fChain->SetBranchAddress("mum_OWNPV_ZERR", &mum_OWNPV_ZERR, &b_mum_OWNPV_ZERR);
   fChain->SetBranchAddress("mum_OWNPV_CHI2", &mum_OWNPV_CHI2, &b_mum_OWNPV_CHI2);
   fChain->SetBranchAddress("mum_OWNPV_NDOF", &mum_OWNPV_NDOF, &b_mum_OWNPV_NDOF);
   fChain->SetBranchAddress("mum_OWNPV_COV_", mum_OWNPV_COV_, &b_mum_OWNPV_COV_);
   fChain->SetBranchAddress("mum_IP_OWNPV", &mum_IP_OWNPV, &b_mum_IP_OWNPV);
   fChain->SetBranchAddress("mum_IPCHI2_OWNPV", &mum_IPCHI2_OWNPV, &b_mum_IPCHI2_OWNPV);
   fChain->SetBranchAddress("mum_ORIVX_X", &mum_ORIVX_X, &b_mum_ORIVX_X);
   fChain->SetBranchAddress("mum_ORIVX_Y", &mum_ORIVX_Y, &b_mum_ORIVX_Y);
   fChain->SetBranchAddress("mum_ORIVX_Z", &mum_ORIVX_Z, &b_mum_ORIVX_Z);
   fChain->SetBranchAddress("mum_ORIVX_XERR", &mum_ORIVX_XERR, &b_mum_ORIVX_XERR);
   fChain->SetBranchAddress("mum_ORIVX_YERR", &mum_ORIVX_YERR, &b_mum_ORIVX_YERR);
   fChain->SetBranchAddress("mum_ORIVX_ZERR", &mum_ORIVX_ZERR, &b_mum_ORIVX_ZERR);
   fChain->SetBranchAddress("mum_ORIVX_CHI2", &mum_ORIVX_CHI2, &b_mum_ORIVX_CHI2);
   fChain->SetBranchAddress("mum_ORIVX_NDOF", &mum_ORIVX_NDOF, &b_mum_ORIVX_NDOF);
   fChain->SetBranchAddress("mum_ORIVX_COV_", mum_ORIVX_COV_, &b_mum_ORIVX_COV_);
   fChain->SetBranchAddress("mum_TRUEID", &mum_TRUEID, &b_mum_TRUEID);
   fChain->SetBranchAddress("mum_TRUEP_E", &mum_TRUEP_E, &b_mum_TRUEP_E);
   fChain->SetBranchAddress("mum_TRUEP_X", &mum_TRUEP_X, &b_mum_TRUEP_X);
   fChain->SetBranchAddress("mum_TRUEP_Y", &mum_TRUEP_Y, &b_mum_TRUEP_Y);
   fChain->SetBranchAddress("mum_TRUEP_Z", &mum_TRUEP_Z, &b_mum_TRUEP_Z);
   fChain->SetBranchAddress("mum_TRUEPT", &mum_TRUEPT, &b_mum_TRUEPT);
   fChain->SetBranchAddress("mum_TRUEORIGINVERTEX_X", &mum_TRUEORIGINVERTEX_X, &b_mum_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("mum_TRUEORIGINVERTEX_Y", &mum_TRUEORIGINVERTEX_Y, &b_mum_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("mum_TRUEORIGINVERTEX_Z", &mum_TRUEORIGINVERTEX_Z, &b_mum_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("mum_TRUEENDVERTEX_X", &mum_TRUEENDVERTEX_X, &b_mum_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("mum_TRUEENDVERTEX_Y", &mum_TRUEENDVERTEX_Y, &b_mum_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("mum_TRUEENDVERTEX_Z", &mum_TRUEENDVERTEX_Z, &b_mum_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("mum_TRUEISSTABLE", &mum_TRUEISSTABLE, &b_mum_TRUEISSTABLE);
   fChain->SetBranchAddress("mum_TRUETAU", &mum_TRUETAU, &b_mum_TRUETAU);
   fChain->SetBranchAddress("mum_TRACK_Type", &mum_TRACK_Type, &b_mum_TRACK_Type);
   fChain->SetBranchAddress("mum_TRACK_Key", &mum_TRACK_Key, &b_mum_TRACK_Key);
   fChain->SetBranchAddress("mum_TRACK_CHI2NDOF", &mum_TRACK_CHI2NDOF, &b_mum_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("mum_TRACK_PCHI2", &mum_TRACK_PCHI2, &b_mum_TRACK_PCHI2);
   fChain->SetBranchAddress("mum_TRACK_MatchCHI2", &mum_TRACK_MatchCHI2, &b_mum_TRACK_MatchCHI2);
   fChain->SetBranchAddress("mum_TRACK_GhostProb", &mum_TRACK_GhostProb, &b_mum_TRACK_GhostProb);
   fChain->SetBranchAddress("mum_TRACK_CloneDist", &mum_TRACK_CloneDist, &b_mum_TRACK_CloneDist);
   fChain->SetBranchAddress("mum_TRACK_Likelihood", &mum_TRACK_Likelihood, &b_mum_TRACK_Likelihood);
   fChain->SetBranchAddress("K_MC12TuneV2_ProbNNe", &K_MC12TuneV2_ProbNNe, &b_K_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("K_MC12TuneV2_ProbNNmu", &K_MC12TuneV2_ProbNNmu, &b_K_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("K_MC12TuneV2_ProbNNpi", &K_MC12TuneV2_ProbNNpi, &b_K_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("K_MC12TuneV2_ProbNNk", &K_MC12TuneV2_ProbNNk, &b_K_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("K_MC12TuneV2_ProbNNp", &K_MC12TuneV2_ProbNNp, &b_K_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("K_MC12TuneV2_ProbNNghost", &K_MC12TuneV2_ProbNNghost, &b_K_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("K_MC12TuneV3_ProbNNe", &K_MC12TuneV3_ProbNNe, &b_K_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("K_MC12TuneV3_ProbNNmu", &K_MC12TuneV3_ProbNNmu, &b_K_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("K_MC12TuneV3_ProbNNpi", &K_MC12TuneV3_ProbNNpi, &b_K_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("K_MC12TuneV3_ProbNNk", &K_MC12TuneV3_ProbNNk, &b_K_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("K_MC12TuneV3_ProbNNp", &K_MC12TuneV3_ProbNNp, &b_K_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("K_MC12TuneV3_ProbNNghost", &K_MC12TuneV3_ProbNNghost, &b_K_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("K_MC12TuneV4_ProbNNe", &K_MC12TuneV4_ProbNNe, &b_K_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("K_MC12TuneV4_ProbNNmu", &K_MC12TuneV4_ProbNNmu, &b_K_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("K_MC12TuneV4_ProbNNpi", &K_MC12TuneV4_ProbNNpi, &b_K_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("K_MC12TuneV4_ProbNNk", &K_MC12TuneV4_ProbNNk, &b_K_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("K_MC12TuneV4_ProbNNp", &K_MC12TuneV4_ProbNNp, &b_K_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("K_MC12TuneV4_ProbNNghost", &K_MC12TuneV4_ProbNNghost, &b_K_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("K_MC15TuneV1_ProbNNe", &K_MC15TuneV1_ProbNNe, &b_K_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("K_MC15TuneV1_ProbNNmu", &K_MC15TuneV1_ProbNNmu, &b_K_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("K_MC15TuneV1_ProbNNpi", &K_MC15TuneV1_ProbNNpi, &b_K_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("K_MC15TuneV1_ProbNNk", &K_MC15TuneV1_ProbNNk, &b_K_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("K_MC15TuneV1_ProbNNp", &K_MC15TuneV1_ProbNNp, &b_K_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("K_MC15TuneV1_ProbNNghost", &K_MC15TuneV1_ProbNNghost, &b_K_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("K_P", &K_P, &b_K_P);
   fChain->SetBranchAddress("K_PT", &K_PT, &b_K_PT);
   fChain->SetBranchAddress("K_PE", &K_PE, &b_K_PE);
   fChain->SetBranchAddress("K_PX", &K_PX, &b_K_PX);
   fChain->SetBranchAddress("K_PY", &K_PY, &b_K_PY);
   fChain->SetBranchAddress("K_PZ", &K_PZ, &b_K_PZ);
   fChain->SetBranchAddress("K_M", &K_M, &b_K_M);
   fChain->SetBranchAddress("K_ID", &K_ID, &b_K_ID);
   fChain->SetBranchAddress("K_PIDe", &K_PIDe, &b_K_PIDe);
   fChain->SetBranchAddress("K_PIDmu", &K_PIDmu, &b_K_PIDmu);
   fChain->SetBranchAddress("K_PIDK", &K_PIDK, &b_K_PIDK);
   fChain->SetBranchAddress("K_PIDp", &K_PIDp, &b_K_PIDp);
   fChain->SetBranchAddress("K_PIDd", &K_PIDd, &b_K_PIDd);
   fChain->SetBranchAddress("K_ProbNNe", &K_ProbNNe, &b_K_ProbNNe);
   fChain->SetBranchAddress("K_ProbNNk", &K_ProbNNk, &b_K_ProbNNk);
   fChain->SetBranchAddress("K_ProbNNp", &K_ProbNNp, &b_K_ProbNNp);
   fChain->SetBranchAddress("K_ProbNNpi", &K_ProbNNpi, &b_K_ProbNNpi);
   fChain->SetBranchAddress("K_ProbNNmu", &K_ProbNNmu, &b_K_ProbNNmu);
   fChain->SetBranchAddress("K_ProbNNd", &K_ProbNNd, &b_K_ProbNNd);
   fChain->SetBranchAddress("K_ProbNNghost", &K_ProbNNghost, &b_K_ProbNNghost);
   fChain->SetBranchAddress("K_hasMuon", &K_hasMuon, &b_K_hasMuon);
   fChain->SetBranchAddress("K_isMuon", &K_isMuon, &b_K_isMuon);
   fChain->SetBranchAddress("K_hasRich", &K_hasRich, &b_K_hasRich);
   fChain->SetBranchAddress("K_UsedRichAerogel", &K_UsedRichAerogel, &b_K_UsedRichAerogel);
   fChain->SetBranchAddress("K_UsedRich1Gas", &K_UsedRich1Gas, &b_K_UsedRich1Gas);
   fChain->SetBranchAddress("K_UsedRich2Gas", &K_UsedRich2Gas, &b_K_UsedRich2Gas);
   fChain->SetBranchAddress("K_RichAboveElThres", &K_RichAboveElThres, &b_K_RichAboveElThres);
   fChain->SetBranchAddress("K_RichAboveMuThres", &K_RichAboveMuThres, &b_K_RichAboveMuThres);
   fChain->SetBranchAddress("K_RichAbovePiThres", &K_RichAbovePiThres, &b_K_RichAbovePiThres);
   fChain->SetBranchAddress("K_RichAboveKaThres", &K_RichAboveKaThres, &b_K_RichAboveKaThres);
   fChain->SetBranchAddress("K_RichAbovePrThres", &K_RichAbovePrThres, &b_K_RichAbovePrThres);
   fChain->SetBranchAddress("K_hasCalo", &K_hasCalo, &b_K_hasCalo);
   fChain->SetBranchAddress("K_L0Global_Dec", &K_L0Global_Dec, &b_K_L0Global_Dec);
   fChain->SetBranchAddress("K_L0Global_TIS", &K_L0Global_TIS, &b_K_L0Global_TIS);
   fChain->SetBranchAddress("K_L0Global_TOS", &K_L0Global_TOS, &b_K_L0Global_TOS);
   fChain->SetBranchAddress("K_Hlt1Global_Dec", &K_Hlt1Global_Dec, &b_K_Hlt1Global_Dec);
   fChain->SetBranchAddress("K_Hlt1Global_TIS", &K_Hlt1Global_TIS, &b_K_Hlt1Global_TIS);
   fChain->SetBranchAddress("K_Hlt1Global_TOS", &K_Hlt1Global_TOS, &b_K_Hlt1Global_TOS);
   fChain->SetBranchAddress("K_Hlt1Phys_Dec", &K_Hlt1Phys_Dec, &b_K_Hlt1Phys_Dec);
   fChain->SetBranchAddress("K_Hlt1Phys_TIS", &K_Hlt1Phys_TIS, &b_K_Hlt1Phys_TIS);
   fChain->SetBranchAddress("K_Hlt1Phys_TOS", &K_Hlt1Phys_TOS, &b_K_Hlt1Phys_TOS);
   fChain->SetBranchAddress("K_Hlt2Global_Dec", &K_Hlt2Global_Dec, &b_K_Hlt2Global_Dec);
   fChain->SetBranchAddress("K_Hlt2Global_TIS", &K_Hlt2Global_TIS, &b_K_Hlt2Global_TIS);
   fChain->SetBranchAddress("K_Hlt2Global_TOS", &K_Hlt2Global_TOS, &b_K_Hlt2Global_TOS);
   fChain->SetBranchAddress("K_Hlt2Phys_Dec", &K_Hlt2Phys_Dec, &b_K_Hlt2Phys_Dec);
   fChain->SetBranchAddress("K_Hlt2Phys_TIS", &K_Hlt2Phys_TIS, &b_K_Hlt2Phys_TIS);
   fChain->SetBranchAddress("K_Hlt2Phys_TOS", &K_Hlt2Phys_TOS, &b_K_Hlt2Phys_TOS);
   fChain->SetBranchAddress("K_L0MuonDecision_Dec", &K_L0MuonDecision_Dec, &b_K_L0MuonDecision_Dec);
   fChain->SetBranchAddress("K_L0MuonDecision_TIS", &K_L0MuonDecision_TIS, &b_K_L0MuonDecision_TIS);
   fChain->SetBranchAddress("K_L0MuonDecision_TOS", &K_L0MuonDecision_TOS, &b_K_L0MuonDecision_TOS);
   fChain->SetBranchAddress("K_L0DiMuonDecision_Dec", &K_L0DiMuonDecision_Dec, &b_K_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("K_L0DiMuonDecision_TIS", &K_L0DiMuonDecision_TIS, &b_K_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("K_L0DiMuonDecision_TOS", &K_L0DiMuonDecision_TOS, &b_K_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("K_Hlt1DiMuonHighMassDecision_Dec", &K_Hlt1DiMuonHighMassDecision_Dec, &b_K_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("K_Hlt1DiMuonHighMassDecision_TIS", &K_Hlt1DiMuonHighMassDecision_TIS, &b_K_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("K_Hlt1DiMuonHighMassDecision_TOS", &K_Hlt1DiMuonHighMassDecision_TOS, &b_K_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedDecision_Dec", &K_Hlt2DiMuonDetachedDecision_Dec, &b_K_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedDecision_TIS", &K_Hlt2DiMuonDetachedDecision_TIS, &b_K_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedDecision_TOS", &K_Hlt2DiMuonDetachedDecision_TOS, &b_K_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedHeavyDecision_Dec", &K_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_K_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedHeavyDecision_TIS", &K_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_K_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedHeavyDecision_TOS", &K_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_K_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedJPsiDecision_Dec", &K_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_K_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedJPsiDecision_TIS", &K_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_K_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedJPsiDecision_TOS", &K_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_K_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedPsi2SDecision_Dec", &K_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_K_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedPsi2SDecision_TIS", &K_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_K_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("K_Hlt2DiMuonDetachedPsi2SDecision_TOS", &K_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_K_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("K_Hlt2DiMuonJPsiDecision_Dec", &K_Hlt2DiMuonJPsiDecision_Dec, &b_K_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("K_Hlt2DiMuonJPsiDecision_TIS", &K_Hlt2DiMuonJPsiDecision_TIS, &b_K_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("K_Hlt2DiMuonJPsiDecision_TOS", &K_Hlt2DiMuonJPsiDecision_TOS, &b_K_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("K_Hlt2DiMuonJPsiHighPTDecision_Dec", &K_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_K_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("K_Hlt2DiMuonJPsiHighPTDecision_TIS", &K_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_K_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("K_Hlt2DiMuonJPsiHighPTDecision_TOS", &K_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_K_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("K_OWNPV_X", &K_OWNPV_X, &b_K_OWNPV_X);
   fChain->SetBranchAddress("K_OWNPV_Y", &K_OWNPV_Y, &b_K_OWNPV_Y);
   fChain->SetBranchAddress("K_OWNPV_Z", &K_OWNPV_Z, &b_K_OWNPV_Z);
   fChain->SetBranchAddress("K_OWNPV_XERR", &K_OWNPV_XERR, &b_K_OWNPV_XERR);
   fChain->SetBranchAddress("K_OWNPV_YERR", &K_OWNPV_YERR, &b_K_OWNPV_YERR);
   fChain->SetBranchAddress("K_OWNPV_ZERR", &K_OWNPV_ZERR, &b_K_OWNPV_ZERR);
   fChain->SetBranchAddress("K_OWNPV_CHI2", &K_OWNPV_CHI2, &b_K_OWNPV_CHI2);
   fChain->SetBranchAddress("K_OWNPV_NDOF", &K_OWNPV_NDOF, &b_K_OWNPV_NDOF);
   fChain->SetBranchAddress("K_OWNPV_COV_", K_OWNPV_COV_, &b_K_OWNPV_COV_);
   fChain->SetBranchAddress("K_IP_OWNPV", &K_IP_OWNPV, &b_K_IP_OWNPV);
   fChain->SetBranchAddress("K_IPCHI2_OWNPV", &K_IPCHI2_OWNPV, &b_K_IPCHI2_OWNPV);
   fChain->SetBranchAddress("K_ORIVX_X", &K_ORIVX_X, &b_K_ORIVX_X);
   fChain->SetBranchAddress("K_ORIVX_Y", &K_ORIVX_Y, &b_K_ORIVX_Y);
   fChain->SetBranchAddress("K_ORIVX_Z", &K_ORIVX_Z, &b_K_ORIVX_Z);
   fChain->SetBranchAddress("K_ORIVX_XERR", &K_ORIVX_XERR, &b_K_ORIVX_XERR);
   fChain->SetBranchAddress("K_ORIVX_YERR", &K_ORIVX_YERR, &b_K_ORIVX_YERR);
   fChain->SetBranchAddress("K_ORIVX_ZERR", &K_ORIVX_ZERR, &b_K_ORIVX_ZERR);
   fChain->SetBranchAddress("K_ORIVX_CHI2", &K_ORIVX_CHI2, &b_K_ORIVX_CHI2);
   fChain->SetBranchAddress("K_ORIVX_NDOF", &K_ORIVX_NDOF, &b_K_ORIVX_NDOF);
   fChain->SetBranchAddress("K_ORIVX_COV_", K_ORIVX_COV_, &b_K_ORIVX_COV_);
   fChain->SetBranchAddress("K_TRUEID", &K_TRUEID, &b_K_TRUEID);
   fChain->SetBranchAddress("K_TRUEP_E", &K_TRUEP_E, &b_K_TRUEP_E);
   fChain->SetBranchAddress("K_TRUEP_X", &K_TRUEP_X, &b_K_TRUEP_X);
   fChain->SetBranchAddress("K_TRUEP_Y", &K_TRUEP_Y, &b_K_TRUEP_Y);
   fChain->SetBranchAddress("K_TRUEP_Z", &K_TRUEP_Z, &b_K_TRUEP_Z);
   fChain->SetBranchAddress("K_TRUEPT", &K_TRUEPT, &b_K_TRUEPT);
   fChain->SetBranchAddress("K_TRUEORIGINVERTEX_X", &K_TRUEORIGINVERTEX_X, &b_K_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("K_TRUEORIGINVERTEX_Y", &K_TRUEORIGINVERTEX_Y, &b_K_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("K_TRUEORIGINVERTEX_Z", &K_TRUEORIGINVERTEX_Z, &b_K_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("K_TRUEENDVERTEX_X", &K_TRUEENDVERTEX_X, &b_K_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("K_TRUEENDVERTEX_Y", &K_TRUEENDVERTEX_Y, &b_K_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("K_TRUEENDVERTEX_Z", &K_TRUEENDVERTEX_Z, &b_K_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("K_TRUEISSTABLE", &K_TRUEISSTABLE, &b_K_TRUEISSTABLE);
   fChain->SetBranchAddress("K_TRUETAU", &K_TRUETAU, &b_K_TRUETAU);
   fChain->SetBranchAddress("K_TRACK_Type", &K_TRACK_Type, &b_K_TRACK_Type);
   fChain->SetBranchAddress("K_TRACK_Key", &K_TRACK_Key, &b_K_TRACK_Key);
   fChain->SetBranchAddress("K_TRACK_CHI2NDOF", &K_TRACK_CHI2NDOF, &b_K_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("K_TRACK_PCHI2", &K_TRACK_PCHI2, &b_K_TRACK_PCHI2);
   fChain->SetBranchAddress("K_TRACK_MatchCHI2", &K_TRACK_MatchCHI2, &b_K_TRACK_MatchCHI2);
   fChain->SetBranchAddress("K_TRACK_GhostProb", &K_TRACK_GhostProb, &b_K_TRACK_GhostProb);
   fChain->SetBranchAddress("K_TRACK_CloneDist", &K_TRACK_CloneDist, &b_K_TRACK_CloneDist);
   fChain->SetBranchAddress("K_TRACK_Likelihood", &K_TRACK_Likelihood, &b_K_TRACK_Likelihood);
   fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
   fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
   fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
   fChain->SetBranchAddress("BCType", &BCType, &b_BCType);
   fChain->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
   fChain->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
   fChain->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
   fChain->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
   fChain->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
   fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
   fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("nLongTracks", &nLongTracks, &b_nLongTracks);
   fChain->SetBranchAddress("nDownstreamTracks", &nDownstreamTracks, &b_nDownstreamTracks);
   fChain->SetBranchAddress("nUpstreamTracks", &nUpstreamTracks, &b_nUpstreamTracks);
   fChain->SetBranchAddress("nVeloTracks", &nVeloTracks, &b_nVeloTracks);
   fChain->SetBranchAddress("nTTracks", &nTTracks, &b_nTTracks);
   fChain->SetBranchAddress("nBackTracks", &nBackTracks, &b_nBackTracks);
   fChain->SetBranchAddress("nRich1Hits", &nRich1Hits, &b_nRich1Hits);
   fChain->SetBranchAddress("nRich2Hits", &nRich2Hits, &b_nRich2Hits);
   fChain->SetBranchAddress("nVeloClusters", &nVeloClusters, &b_nVeloClusters);
   fChain->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
   fChain->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
   fChain->SetBranchAddress("nOTClusters", &nOTClusters, &b_nOTClusters);
   fChain->SetBranchAddress("nSPDHits", &nSPDHits, &b_nSPDHits);
   fChain->SetBranchAddress("nMuonCoordsS0", &nMuonCoordsS0, &b_nMuonCoordsS0);
   fChain->SetBranchAddress("nMuonCoordsS1", &nMuonCoordsS1, &b_nMuonCoordsS1);
   fChain->SetBranchAddress("nMuonCoordsS2", &nMuonCoordsS2, &b_nMuonCoordsS2);
   fChain->SetBranchAddress("nMuonCoordsS3", &nMuonCoordsS3, &b_nMuonCoordsS3);
   fChain->SetBranchAddress("nMuonCoordsS4", &nMuonCoordsS4, &b_nMuonCoordsS4);
   fChain->SetBranchAddress("nMuonTracks", &nMuonTracks, &b_nMuonTracks);
   Notify();
}

Bool_t BDecayTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BDecayTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BDecayTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BDecayTree_cxx
