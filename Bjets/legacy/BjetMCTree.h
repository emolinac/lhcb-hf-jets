//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 14 19:20:22 2022 by ROOT version 6.24/06
// from TChain mcjettuple/MCJetTree/
//////////////////////////////////////////////////////////

#ifndef BjetMCTree_h
#define BjetMCTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BjetMCTree
{
public:
  TTree *fChain;  //! pointer to the analyzed TTree or TChain
  Int_t fCurrent; //! current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Int_t NumBHads;
  Bool_t hasb;
  Bool_t hasbbar;
  Bool_t hasc;
  Bool_t hascbar;
  Double_t MCJet_truth_mup_PX;
  Double_t MCJet_truth_mup_PY;
  Double_t MCJet_truth_mup_PZ;
  Double_t MCJet_truth_mup_PE;
  Double_t MCJet_truth_mup_PT;
  Double_t MCJet_truth_mup_PHI;
  Double_t MCJet_truth_mup_ETA;
  Double_t MCJet_truth_mup_M;
  Double_t MCJet_truth_mup_XVTX;
  Double_t MCJet_truth_mup_YVTX;
  Double_t MCJet_truth_mup_ZVTX;
  Double_t MCJet_truth_mum_PX;
  Double_t MCJet_truth_mum_PY;
  Double_t MCJet_truth_mum_PZ;
  Double_t MCJet_truth_mum_PE;
  Double_t MCJet_truth_mum_PT;
  Double_t MCJet_truth_mum_PHI;
  Double_t MCJet_truth_mum_ETA;
  Double_t MCJet_truth_mum_M;
  Double_t MCJet_truth_mum_XVTX;
  Double_t MCJet_truth_mum_YVTX;
  Double_t MCJet_truth_mum_ZVTX;
  Double_t MCJet_truth_K_PX;
  Double_t MCJet_truth_K_PY;
  Double_t MCJet_truth_K_PZ;
  Double_t MCJet_truth_K_PE;
  Double_t MCJet_truth_K_PT;
  Double_t MCJet_truth_K_PHI;
  Double_t MCJet_truth_K_ETA;
  Double_t MCJet_truth_K_M;
  Double_t MCJet_truth_K_XVTX;
  Double_t MCJet_truth_K_YVTX;
  Double_t MCJet_truth_K_ZVTX;
  Double_t MCJet_truth_Jpsi_PX;
  Double_t MCJet_truth_Jpsi_PY;
  Double_t MCJet_truth_Jpsi_PZ;
  Double_t MCJet_truth_Jpsi_PE;
  Double_t MCJet_truth_Jpsi_PT;
  Double_t MCJet_truth_Jpsi_PHI;
  Double_t MCJet_truth_Jpsi_ETA;
  Double_t MCJet_truth_Jpsi_M;
  Double_t MCJet_truth_B_PX;
  Double_t MCJet_truth_B_PY;
  Double_t MCJet_truth_B_PZ;
  Double_t MCJet_truth_B_PE;
  Double_t MCJet_truth_B_PT;
  Double_t MCJet_truth_B_PHI;
  Double_t MCJet_truth_B_ETA;
  Double_t MCJet_truth_B_M;

  Double_t MCJet_recojet_deltaR;
  Double_t MCJet_recojet_PT;
  Int_t MCJet_recojet_nrecodtrs;
  Double_t MCJet_recojet_PX;
  Double_t MCJet_recojet_PY;
  Double_t MCJet_recojet_PZ;
  Double_t MCJet_recojet_PE;
  Double_t MCJet_recojet_ETA;
  Double_t MCJet_recojet_PHI;
  Double_t MCJet_recojet_M;
  Double_t MCJet_recojet_XVTX;
  Double_t MCJet_recojet_YVTX;
  Double_t MCJet_recojet_CPF;
  Double_t MCJet_recojet_MTF;
  Double_t MCJet_recojet_NHasPVInfo;
  Double_t MCJet_recojet_MPT;
  Int_t MCJet_NDtr;
  Float_t MCJet_recojet_Dtr_ID[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ThreeCharge[1024];      //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_IsMeson[1024];          //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_IsBaryon[1024];         //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_IsStrange[1024];        //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PIDe[1024];             //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PIDmu[1024];            //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PIDk[1024];             //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PIDp[1024];             //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PIDpi[1024];            //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ProbNNe[1024];          //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ProbNNmu[1024];         //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ProbNNk[1024];          //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ProbNNp[1024];          //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ProbNNpi[1024];         //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ProbNNghost[1024];      //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PX[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PY[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PZ[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_E[1024];                //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_P[1024];                //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PT[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_M[1024];                //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PHI[1024];              //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_ETA[1024];              //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_VX[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_VY[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_VZ[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_Z[1024];                //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_JT[1024];               //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_R[1024];                //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_TrackType[1024];        //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_TrackKey[1024];         //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_TrackChi2[1024];        //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_TrackNDF[1024];         //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_TrackPCHI2[1024];       //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PFType[1024];           //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PFDaughterType[1024];   //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PFclustype[1024];       //[MCJet_NDtr]
  Float_t MCJet_recojet_Dtr_PFclusE[1024];          //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_PX[1024];   //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_PY[1024];   //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_PZ[1024];   //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_E[1024];    //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_P[1024];    //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_PT[1024];   //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_M[1024];    //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_PHI[1024];  //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_ETA[1024];  //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_XVTX[1024]; //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_YVTX[1024]; //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_ZVTX[1024]; //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_PID[1024];  //[MCJet_NDtr]
  Float_t MCJet_recojet_matchedtruthdtr_KEY[1024];  //[MCJet_NDtr]
  Double_t MCJet_PT;
  Double_t MCJet_P;
  Double_t MCJet_PX;
  Double_t MCJet_PY;
  Double_t MCJet_PZ;
  Double_t MCJet_PE;
  Double_t MCJet_ETA;
  Double_t MCJet_PHI;
  Double_t MCJet_M;
  Int_t MCJet_FLAV;
  Int_t MCJet_Dtr_nmcdtrs;
  Int_t MCJet_Dtr_nmcdtr_recomatches;
  Float_t MCJet_Dtr_mcrecomatch_p[1024];            //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_pt[1024];           //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_pide[1024];         //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_pidmu[1024];        //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_pidk[1024];         //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_pidp[1024];         //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_pidpi[1024];        //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_probnn_e[1024];     //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_probnn_mu[1024];    //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_probnn_k[1024];     //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_probnn_p[1024];     //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_probnn_pi[1024];    //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_probnn_ghost[1024]; //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_track_key[1024];    //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_track_ndf[1024];    //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_track_type[1024];   //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_track_pchi2[1024];  //[MCJet_NDtr]
  Float_t MCJet_Dtr_mcrecomatch_track_chi2[1024];   //[MCJet_NDtr]
  Float_t MCJet_Dtr_ID[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_Key[1024];                      //[MCJet_NDtr]
  Float_t MCJet_Dtr_ThreeCharge[1024];              //[MCJet_NDtr]
  Float_t MCJet_Dtr_IsMeson[1024];                  //[MCJet_NDtr]
  Float_t MCJet_Dtr_IsBaryon[1024];                 //[MCJet_NDtr]
  Float_t MCJet_Dtr_IsStrange[1024];                //[MCJet_NDtr]
  Float_t MCJet_Dtr_PX[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_PY[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_PZ[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_E[1024];                        //[MCJet_NDtr]
  Float_t MCJet_Dtr_P[1024];                        //[MCJet_NDtr]
  Float_t MCJet_Dtr_R[1024];                        //[MCJet_NDtr]
  Float_t MCJet_Dtr_JT[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_Z[1024];                        //[MCJet_NDtr]
  Float_t MCJet_Dtr_PT[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_M[1024];                        //[MCJet_NDtr]
  Float_t MCJet_Dtr_PHI[1024];                      //[MCJet_NDtr]
  Float_t MCJet_Dtr_ETA[1024];                      //[MCJet_NDtr]
  Float_t MCJet_Dtr_VX[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_VY[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_VZ[1024];                       //[MCJet_NDtr]
  Float_t MCJet_Dtr_TopMotherID[1024];              //[MCJet_NDtr]
  Double_t CELLjet_TRUEP_E;
  Double_t CELLjet_TRUEP_X;
  Double_t CELLjet_TRUEP_Y;
  Double_t CELLjet_TRUEP_Z;
  Double_t CELLjet_TRUEPT;
  Double_t CELLjet_TRUEORIGINVERTEX_X;
  Double_t CELLjet_TRUEORIGINVERTEX_Y;
  Double_t CELLjet_TRUEORIGINVERTEX_Z;
  Double_t CELLjet_TRUEENDVERTEX_X;
  Double_t CELLjet_TRUEENDVERTEX_Y;
  Double_t CELLjet_TRUEENDVERTEX_Z;
  Bool_t CELLjet_TRUEISSTABLE;
  Double_t CELLjet_TRUETAU;
  UInt_t nCandidate;
  ULong64_t totCandidates;
  ULong64_t EventInSequence;
  UInt_t EVT_Int__I;
  Double_t EVT_Int__Mean;
  Double_t EVT_Int__Prob;
  UInt_t runNumber;
  ULong64_t eventNumber;
  UInt_t BCID;
  Int_t BCType;
  UInt_t OdinTCK;
  UInt_t L0DUTCK;
  UInt_t HLT1TCK;
  UInt_t HLT2TCK;
  ULong64_t GpsTime;
  Short_t Polarity;
  Int_t nPVs;
  Int_t nTracks;
  Int_t nLongTracks;
  Int_t nDownstreamTracks;
  Int_t nUpstreamTracks;
  Int_t nVeloTracks;
  Int_t nTTracks;
  Int_t nBackTracks;
  Int_t nRich1Hits;
  Int_t nRich2Hits;
  Int_t nVeloClusters;
  Int_t nITClusters;
  Int_t nTTClusters;
  Int_t nOTClusters;
  Int_t nSPDHits;
  Int_t nMuonCoordsS0;
  Int_t nMuonCoordsS1;
  Int_t nMuonCoordsS2;
  Int_t nMuonCoordsS3;
  Int_t nMuonCoordsS4;
  Int_t nMuonTracks;

  // List of branches
  TBranch *b_NumBHads;
  TBranch *b_hasb;
  TBranch *b_hasbbar;
  TBranch *b_hasc;
  TBranch *b_hascbar;
  TBranch *b_MCJet_truth_mup_PX;   //!
  TBranch *b_MCJet_truth_mup_PY;   //!
  TBranch *b_MCJet_truth_mup_PZ;   //!
  TBranch *b_MCJet_truth_mup_PE;   //!
  TBranch *b_MCJet_truth_mup_PT;   //!
  TBranch *b_MCJet_truth_mup_PHI;  //!
  TBranch *b_MCJet_truth_mup_ETA;  //!
  TBranch *b_MCJet_truth_mup_M;    //!
  TBranch *b_MCJet_truth_mup_XVTX; //!
  TBranch *b_MCJet_truth_mup_YVTX; //!
  TBranch *b_MCJet_truth_mup_ZVTX; //!
  TBranch *b_MCJet_truth_mum_PX;   //!
  TBranch *b_MCJet_truth_mum_PY;   //!
  TBranch *b_MCJet_truth_mum_PZ;   //!
  TBranch *b_MCJet_truth_mum_PE;   //!
  TBranch *b_MCJet_truth_mum_PT;   //!
  TBranch *b_MCJet_truth_mum_PHI;  //!
  TBranch *b_MCJet_truth_mum_ETA;  //!
  TBranch *b_MCJet_truth_mum_M;    //!
  TBranch *b_MCJet_truth_mum_XVTX; //!
  TBranch *b_MCJet_truth_mum_YVTX; //!
  TBranch *b_MCJet_truth_mum_ZVTX; //!
  TBranch *b_MCJet_truth_K_PX;     //!
  TBranch *b_MCJet_truth_K_PY;     //!
  TBranch *b_MCJet_truth_K_PZ;     //!
  TBranch *b_MCJet_truth_K_PE;     //!
  TBranch *b_MCJet_truth_K_PT;     //!
  TBranch *b_MCJet_truth_K_PHI;    //!
  TBranch *b_MCJet_truth_K_ETA;    //!
  TBranch *b_MCJet_truth_K_M;      //!
  TBranch *b_MCJet_truth_K_XVTX;   //!
  TBranch *b_MCJet_truth_K_YVTX;   //!
  TBranch *b_MCJet_truth_K_ZVTX;
  TBranch *b_MCJet_truth_Jpsi_PX;  //!
  TBranch *b_MCJet_truth_Jpsi_PY;  //!
  TBranch *b_MCJet_truth_Jpsi_PZ;  //!
  TBranch *b_MCJet_truth_Jpsi_PE;  //!
  TBranch *b_MCJet_truth_Jpsi_PT;  //!
  TBranch *b_MCJet_truth_Jpsi_PHI; //!
  TBranch *b_MCJet_truth_Jpsi_ETA; //!
  TBranch *b_MCJet_truth_Jpsi_M;
  TBranch *b_MCJet_truth_B_PX;                   //!
  TBranch *b_MCJet_truth_B_PY;                   //!
  TBranch *b_MCJet_truth_B_PZ;                   //!
  TBranch *b_MCJet_truth_B_PE;                   //!
  TBranch *b_MCJet_truth_B_PT;                   //!
  TBranch *b_MCJet_truth_B_PHI;                  //!
  TBranch *b_MCJet_truth_B_ETA;                  //!
  TBranch *b_MCJet_truth_B_M;                    //!
  TBranch *b_MCJet_recojet_deltaR;               //!
  TBranch *b_MCJet_recojet_PT;                   //!
  TBranch *b_MCJet_recojet_nrecodtrs;            //!
  TBranch *b_MCJet_recojet_PX;                   //!
  TBranch *b_MCJet_recojet_PY;                   //!
  TBranch *b_MCJet_recojet_PZ;                   //!
  TBranch *b_MCJet_recojet_PE;                   //!
  TBranch *b_MCJet_recojet_ETA;                  //!
  TBranch *b_MCJet_recojet_PHI;                  //!
  TBranch *b_MCJet_recojet_M;                    //!
  TBranch *b_MCJet_recojet_XVTX;                 //!
  TBranch *b_MCJet_recojet_YVTX;                 //!
  TBranch *b_MCJet_recojet_CPF;                  //!
  TBranch *b_MCJet_recojet_MTF;                  //!
  TBranch *b_MCJet_recojet_NHasPVInfo;           //!
  TBranch *b_MCJet_recojet_MPT;                  //!
  TBranch *b_MCJet_NDtr;                         //!
  TBranch *b_MCJet_recojet_Dtr_ID;               //!
  TBranch *b_MCJet_recojet_Dtr_ThreeCharge;      //!
  TBranch *b_MCJet_recojet_Dtr_IsMeson;          //!
  TBranch *b_MCJet_recojet_Dtr_IsBaryon;         //!
  TBranch *b_MCJet_recojet_Dtr_IsStrange;        //!
  TBranch *b_MCJet_recojet_Dtr_PIDe;             //!
  TBranch *b_MCJet_recojet_Dtr_PIDmu;            //!
  TBranch *b_MCJet_recojet_Dtr_PIDk;             //!
  TBranch *b_MCJet_recojet_Dtr_PIDp;             //!
  TBranch *b_MCJet_recojet_Dtr_PIDpi;            //!
  TBranch *b_MCJet_recojet_Dtr_ProbNNe;          //!
  TBranch *b_MCJet_recojet_Dtr_ProbNNmu;         //!
  TBranch *b_MCJet_recojet_Dtr_ProbNNk;          //!
  TBranch *b_MCJet_recojet_Dtr_ProbNNp;          //!
  TBranch *b_MCJet_recojet_Dtr_ProbNNpi;         //!
  TBranch *b_MCJet_recojet_Dtr_ProbNNghost;      //!
  TBranch *b_MCJet_recojet_Dtr_PX;               //!
  TBranch *b_MCJet_recojet_Dtr_PY;               //!
  TBranch *b_MCJet_recojet_Dtr_PZ;               //!
  TBranch *b_MCJet_recojet_Dtr_E;                //!
  TBranch *b_MCJet_recojet_Dtr_P;                //!
  TBranch *b_MCJet_recojet_Dtr_PT;               //!
  TBranch *b_MCJet_recojet_Dtr_M;                //!
  TBranch *b_MCJet_recojet_Dtr_PHI;              //!
  TBranch *b_MCJet_recojet_Dtr_ETA;              //!
  TBranch *b_MCJet_recojet_Dtr_VX;               //!
  TBranch *b_MCJet_recojet_Dtr_VY;               //!
  TBranch *b_MCJet_recojet_Dtr_VZ;               //!
  TBranch *b_MCJet_recojet_Dtr_Z;                //!
  TBranch *b_MCJet_recojet_Dtr_JT;               //!
  TBranch *b_MCJet_recojet_Dtr_R;                //!
  TBranch *b_MCJet_recojet_Dtr_TrackType;        //!
  TBranch *b_MCJet_recojet_Dtr_TrackKey;         //!
  TBranch *b_MCJet_recojet_Dtr_TrackChi2;        //!
  TBranch *b_MCJet_recojet_Dtr_TrackNDF;         //!
  TBranch *b_MCJet_recojet_Dtr_TrackPCHI2;       //!
  TBranch *b_MCJet_recojet_Dtr_PFType;           //!
  TBranch *b_MCJet_recojet_Dtr_PFDaughterType;   //!
  TBranch *b_MCJet_recojet_Dtr_PFclustype;       //!
  TBranch *b_MCJet_recojet_Dtr_PFclusE;          //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_PX;   //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_PY;   //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_PZ;   //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_E;    //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_P;    //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_PT;   //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_M;    //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_PHI;  //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_ETA;  //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_XVTX; //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_YVTX; //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_ZVTX; //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_PID;  //!
  TBranch *b_MCJet_recojet_matchedtruthdtr_KEY;  //!
  TBranch *b_MCJet_PT;                           //!
  TBranch *b_MCJet_P;                            //!
  TBranch *b_MCJet_PX;                           //!
  TBranch *b_MCJet_PY;                           //!
  TBranch *b_MCJet_PZ;                           //!
  TBranch *b_MCJet_PE;                           //!
  TBranch *b_MCJet_ETA;                          //!
  TBranch *b_MCJet_PHI;                          //!
  TBranch *b_MCJet_M;                            //!
  TBranch *b_MCJet_FLAV;                         //!
  TBranch *b_MCJet_Dtr_nmcdtrs;                  //!
  TBranch *b_MCJet_Dtr_nmcdtr_recomatches;       //!
  TBranch *b_MCJet_Dtr_mcrecomatch_p;            //!
  TBranch *b_MCJet_Dtr_mcrecomatch_pt;           //!
  TBranch *b_MCJet_Dtr_mcrecomatch_pide;         //!
  TBranch *b_MCJet_Dtr_mcrecomatch_pidmu;        //!
  TBranch *b_MCJet_Dtr_mcrecomatch_pidk;         //!
  TBranch *b_MCJet_Dtr_mcrecomatch_pidp;         //!
  TBranch *b_MCJet_Dtr_mcrecomatch_pidpi;        //!
  TBranch *b_MCJet_Dtr_mcrecomatch_probnn_e;     //!
  TBranch *b_MCJet_Dtr_mcrecomatch_probnn_mu;    //!
  TBranch *b_MCJet_Dtr_mcrecomatch_probnn_k;     //!
  TBranch *b_MCJet_Dtr_mcrecomatch_probnn_p;     //!
  TBranch *b_MCJet_Dtr_mcrecomatch_probnn_pi;    //!
  TBranch *b_MCJet_Dtr_mcrecomatch_probnn_ghost; //!
  TBranch *b_MCJet_Dtr_mcrecomatch_track_key;    //!
  TBranch *b_MCJet_Dtr_mcrecomatch_track_ndf;    //!
  TBranch *b_MCJet_Dtr_mcrecomatch_track_type;   //!
  TBranch *b_MCJet_Dtr_mcrecomatch_track_pchi2;  //!
  TBranch *b_MCJet_Dtr_mcrecomatch_track_chi2;   //!
  TBranch *b_MCJet_Dtr_ID;                       //!
  TBranch *b_MCJet_Dtr_Key;                      //!
  TBranch *b_MCJet_Dtr_ThreeCharge;              //!
  TBranch *b_MCJet_Dtr_IsMeson;                  //!
  TBranch *b_MCJet_Dtr_IsBaryon;                 //!
  TBranch *b_MCJet_Dtr_IsStrange;                //!
  TBranch *b_MCJet_Dtr_PX;                       //!
  TBranch *b_MCJet_Dtr_PY;                       //!
  TBranch *b_MCJet_Dtr_PZ;                       //!
  TBranch *b_MCJet_Dtr_E;                        //!
  TBranch *b_MCJet_Dtr_P;                        //!
  TBranch *b_MCJet_Dtr_R;                        //!
  TBranch *b_MCJet_Dtr_JT;                       //!
  TBranch *b_MCJet_Dtr_Z;                        //!
  TBranch *b_MCJet_Dtr_PT;                       //!
  TBranch *b_MCJet_Dtr_M;                        //!
  TBranch *b_MCJet_Dtr_PHI;                      //!
  TBranch *b_MCJet_Dtr_ETA;                      //!
  TBranch *b_MCJet_Dtr_VX;                       //!
  TBranch *b_MCJet_Dtr_VY;                       //!
  TBranch *b_MCJet_Dtr_VZ;                       //!
  TBranch *b_MCJet_Dtr_TopMotherID;              //!
  TBranch *b_CELLjet_TRUEP_E;                    //!
  TBranch *b_CELLjet_TRUEP_X;                    //!
  TBranch *b_CELLjet_TRUEP_Y;                    //!
  TBranch *b_CELLjet_TRUEP_Z;                    //!
  TBranch *b_CELLjet_TRUEPT;                     //!
  TBranch *b_CELLjet_TRUEORIGINVERTEX_X;         //!
  TBranch *b_CELLjet_TRUEORIGINVERTEX_Y;         //!
  TBranch *b_CELLjet_TRUEORIGINVERTEX_Z;         //!
  TBranch *b_CELLjet_TRUEENDVERTEX_X;            //!
  TBranch *b_CELLjet_TRUEENDVERTEX_Y;            //!
  TBranch *b_CELLjet_TRUEENDVERTEX_Z;            //!
  TBranch *b_CELLjet_TRUEISSTABLE;               //!
  TBranch *b_CELLjet_TRUETAU;                    //!
  TBranch *b_nCandidate;                         //!
  TBranch *b_totCandidates;                      //!
  TBranch *b_EventInSequence;                    //!
  TBranch *b_EVT_Int__I;                         //!
  TBranch *b_EVT_Int__Mean;                      //!
  TBranch *b_EVT_Int__Prob;                      //!
  TBranch *b_runNumber;                          //!
  TBranch *b_eventNumber;                        //!
  TBranch *b_BCID;                               //!
  TBranch *b_BCType;                             //!
  TBranch *b_OdinTCK;                            //!
  TBranch *b_L0DUTCK;                            //!
  TBranch *b_HLT1TCK;                            //!
  TBranch *b_HLT2TCK;                            //!
  TBranch *b_GpsTime;                            //!
  TBranch *b_Polarity;                           //!
  TBranch *b_nPVs;                               //!
  TBranch *b_nTracks;                            //!
  TBranch *b_nLongTracks;                        //!
  TBranch *b_nDownstreamTracks;                  //!
  TBranch *b_nUpstreamTracks;                    //!
  TBranch *b_nVeloTracks;                        //!
  TBranch *b_nTTracks;                           //!
  TBranch *b_nBackTracks;                        //!
  TBranch *b_nRich1Hits;                         //!
  TBranch *b_nRich2Hits;                         //!
  TBranch *b_nVeloClusters;                      //!
  TBranch *b_nITClusters;                        //!
  TBranch *b_nTTClusters;                        //!
  TBranch *b_nOTClusters;                        //!
  TBranch *b_nSPDHits;                           //!
  TBranch *b_nMuonCoordsS0;                      //!
  TBranch *b_nMuonCoordsS1;                      //!
  TBranch *b_nMuonCoordsS2;                      //!
  TBranch *b_nMuonCoordsS3;                      //!
  TBranch *b_nMuonCoordsS4;                      //!
  TBranch *b_nMuonTracks;                        //!

  BjetMCTree(TTree *tree = 0, int dataset = 93539, bool isData = false, bool onlysim9 = false);
  virtual ~BjetMCTree();
  virtual Int_t Cut(Long64_t entry);
  virtual Int_t GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void Init(TTree *tree);
  virtual void Loop();
  virtual Bool_t Notify();
  virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef BjetMCTree_cxx
BjetMCTree::BjetMCTree(TTree *tree = 0, int dataset = 93539, bool isData = false, bool onlysim9 = false) : fChain(0)
{
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0)
  {

#ifdef SINGLE_TREE
    // The following code should be used if you want this class to access
    // a single tree instead of a chain
    TFile *f = (TFile *)gROOT->GetListOfFiles()->FindObject("Memory Directory");
    if (!f || !f->IsOpen())
    {
      f = new TFile("Memory Directory");
    }
    f->GetObject("mcjettuple/MCJetTree", tree);

#else  // SINGLE_TREE

    // The following code should be used if you want this class to access a chain
    // of trees.
    TChain *chain = new TChain("mcjettuple/MCJetTree", "");
    TString data_dir = "../../data/";
    if (dataset == 61590)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");  
        }      
    }
    else if (dataset == 61591)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
     
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root/MCJets/MCJetTree"); 
        }       
    }
    else if (dataset == 61599)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)        
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        }
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {        
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");   
        }     
    }
    else if (dataset == 71590)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        }        
    }
    else if (dataset == 71591)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");
        }        
    }
    else if (dataset == 71599)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        }
            
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");
        }        
    }
    else if (dataset == 81590)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        }        
    }
    else if (dataset == 81591)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");        
        }
    }
    else if (dataset == 81599)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        }
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        if (!onlysim9)
        {
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");
        }
    }
    else if (dataset == 91599)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");

        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/MCJets/MCJetTree"); 

        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");

        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root/MCJets/MCJetTree"); 

        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");                  
    }
    else if (dataset == 91590)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");

        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");

        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_05152024_full.root/MCJets/MCJetTree");
        
       
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");

        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");

        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MD_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MD_02092025_full.root/MCJets/MCJetTree");        
    }
    else if (dataset == 91591)
    {
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");         

        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_05152024_full.root/MCJets/MCJetTree");
        //chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_05152024_full.root/MCJets/MCJetTree");
        
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2016_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");
        
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2017_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");         

        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09h_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09i_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09k_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim09l_MU_02092025_full.root/MCJets/MCJetTree");
        chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_HighPT_2018_Sim10a_MU_02092025_full.root/MCJets/MCJetTree");        
    }    

    tree = chain;
#endif // SINGLE_TREE
  }
  Init(tree);
}

BjetMCTree::~BjetMCTree()
{
  if (!fChain)
    return;
  delete fChain->GetCurrentFile();
}

Int_t BjetMCTree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain)
    return 0;
  return fChain->GetEntry(entry);
}
Long64_t BjetMCTree::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain)
    return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0)
    return centry;
  if (fChain->GetTreeNumber() != fCurrent)
  {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void BjetMCTree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set branch addresses and branch pointers
  if (!tree)
    return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("NumBHads", &NumBHads, &b_NumBHads);
  fChain->SetBranchAddress("hasb", &hasb, &b_hasb);
  fChain->SetBranchAddress("hasbbar", &hasbbar, &b_hasbbar);
  fChain->SetBranchAddress("hasc", &hasc, &b_hasc);
  fChain->SetBranchAddress("hascbar", &hascbar, &b_hascbar);

  fChain->SetBranchAddress("MCJet_truth_mup_PX", &MCJet_truth_mup_PX, &b_MCJet_truth_mup_PX);
  fChain->SetBranchAddress("MCJet_truth_mup_PY", &MCJet_truth_mup_PY, &b_MCJet_truth_mup_PY);
  fChain->SetBranchAddress("MCJet_truth_mup_PZ", &MCJet_truth_mup_PZ, &b_MCJet_truth_mup_PZ);
  fChain->SetBranchAddress("MCJet_truth_mup_PE", &MCJet_truth_mup_PE, &b_MCJet_truth_mup_PE);
  fChain->SetBranchAddress("MCJet_truth_mup_PT", &MCJet_truth_mup_PT, &b_MCJet_truth_mup_PT);
  fChain->SetBranchAddress("MCJet_truth_mup_PHI", &MCJet_truth_mup_PHI, &b_MCJet_truth_mup_PHI);
  fChain->SetBranchAddress("MCJet_truth_mup_ETA", &MCJet_truth_mup_ETA, &b_MCJet_truth_mup_ETA);
  fChain->SetBranchAddress("MCJet_truth_mup_M", &MCJet_truth_mup_M, &b_MCJet_truth_mup_M);
  fChain->SetBranchAddress("MCJet_truth_mup_XVTX", &MCJet_truth_mup_XVTX, &b_MCJet_truth_mup_XVTX);
  fChain->SetBranchAddress("MCJet_truth_mup_YVTX", &MCJet_truth_mup_YVTX, &b_MCJet_truth_mup_YVTX);
  fChain->SetBranchAddress("MCJet_truth_mup_ZVTX", &MCJet_truth_mup_ZVTX, &b_MCJet_truth_mup_ZVTX);
  fChain->SetBranchAddress("MCJet_truth_mum_PX", &MCJet_truth_mum_PX, &b_MCJet_truth_mum_PX);
  fChain->SetBranchAddress("MCJet_truth_mum_PY", &MCJet_truth_mum_PY, &b_MCJet_truth_mum_PY);
  fChain->SetBranchAddress("MCJet_truth_mum_PZ", &MCJet_truth_mum_PZ, &b_MCJet_truth_mum_PZ);
  fChain->SetBranchAddress("MCJet_truth_mum_PE", &MCJet_truth_mum_PE, &b_MCJet_truth_mum_PE);
  fChain->SetBranchAddress("MCJet_truth_mum_PT", &MCJet_truth_mum_PT, &b_MCJet_truth_mum_PT);
  fChain->SetBranchAddress("MCJet_truth_mum_PHI", &MCJet_truth_mum_PHI, &b_MCJet_truth_mum_PHI);
  fChain->SetBranchAddress("MCJet_truth_mum_ETA", &MCJet_truth_mum_ETA, &b_MCJet_truth_mum_ETA);
  fChain->SetBranchAddress("MCJet_truth_mum_M", &MCJet_truth_mum_M, &b_MCJet_truth_mum_M);
  fChain->SetBranchAddress("MCJet_truth_mum_XVTX", &MCJet_truth_mum_XVTX, &b_MCJet_truth_mum_XVTX);
  fChain->SetBranchAddress("MCJet_truth_mum_YVTX", &MCJet_truth_mum_YVTX, &b_MCJet_truth_mum_YVTX);
  fChain->SetBranchAddress("MCJet_truth_mum_ZVTX", &MCJet_truth_mum_ZVTX, &b_MCJet_truth_mum_ZVTX);
  fChain->SetBranchAddress("MCJet_truth_K_PX", &MCJet_truth_K_PX, &b_MCJet_truth_K_PX);
  fChain->SetBranchAddress("MCJet_truth_K_PY", &MCJet_truth_K_PY, &b_MCJet_truth_K_PY);
  fChain->SetBranchAddress("MCJet_truth_K_PZ", &MCJet_truth_K_PZ, &b_MCJet_truth_K_PZ);
  fChain->SetBranchAddress("MCJet_truth_K_PE", &MCJet_truth_K_PE, &b_MCJet_truth_K_PE);
  fChain->SetBranchAddress("MCJet_truth_K_PT", &MCJet_truth_K_PT, &b_MCJet_truth_K_PT);
  fChain->SetBranchAddress("MCJet_truth_K_PHI", &MCJet_truth_K_PHI, &b_MCJet_truth_K_PHI);
  fChain->SetBranchAddress("MCJet_truth_K_ETA", &MCJet_truth_K_ETA, &b_MCJet_truth_K_ETA);
  fChain->SetBranchAddress("MCJet_truth_K_M", &MCJet_truth_K_M, &b_MCJet_truth_K_M);
  fChain->SetBranchAddress("MCJet_truth_K_XVTX", &MCJet_truth_K_XVTX, &b_MCJet_truth_K_XVTX);
  fChain->SetBranchAddress("MCJet_truth_K_YVTX", &MCJet_truth_K_YVTX, &b_MCJet_truth_K_YVTX);
  fChain->SetBranchAddress("MCJet_truth_K_ZVTX", &MCJet_truth_K_ZVTX, &b_MCJet_truth_K_ZVTX);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_PX", &MCJet_truth_Jpsi_PX, &b_MCJet_truth_Jpsi_PX);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_PY", &MCJet_truth_Jpsi_PY, &b_MCJet_truth_Jpsi_PY);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_PZ", &MCJet_truth_Jpsi_PZ, &b_MCJet_truth_Jpsi_PZ);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_PE", &MCJet_truth_Jpsi_PE, &b_MCJet_truth_Jpsi_PE);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_PT", &MCJet_truth_Jpsi_PT, &b_MCJet_truth_Jpsi_PT);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_PHI", &MCJet_truth_Jpsi_PHI, &b_MCJet_truth_Jpsi_PHI);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_ETA", &MCJet_truth_Jpsi_ETA, &b_MCJet_truth_Jpsi_ETA);
  fChain->SetBranchAddress("MCJet_truth_Jpsi_M", &MCJet_truth_Jpsi_M, &b_MCJet_truth_Jpsi_M);
  fChain->SetBranchAddress("MCJet_truth_B_PX", &MCJet_truth_B_PX, &b_MCJet_truth_B_PX);
  fChain->SetBranchAddress("MCJet_truth_B_PY", &MCJet_truth_B_PY, &b_MCJet_truth_B_PY);
  fChain->SetBranchAddress("MCJet_truth_B_PZ", &MCJet_truth_B_PZ, &b_MCJet_truth_B_PZ);
  fChain->SetBranchAddress("MCJet_truth_B_PE", &MCJet_truth_B_PE, &b_MCJet_truth_B_PE);
  fChain->SetBranchAddress("MCJet_truth_B_PT", &MCJet_truth_B_PT, &b_MCJet_truth_B_PT);
  fChain->SetBranchAddress("MCJet_truth_B_PHI", &MCJet_truth_B_PHI, &b_MCJet_truth_B_PHI);
  fChain->SetBranchAddress("MCJet_truth_B_ETA", &MCJet_truth_B_ETA, &b_MCJet_truth_B_ETA);
  fChain->SetBranchAddress("MCJet_truth_B_M", &MCJet_truth_B_M, &b_MCJet_truth_B_M);
  fChain->SetBranchAddress("MCJet_recojet_deltaR", &MCJet_recojet_deltaR, &b_MCJet_recojet_deltaR);
  fChain->SetBranchAddress("MCJet_recojet_PT", &MCJet_recojet_PT, &b_MCJet_recojet_PT);
  fChain->SetBranchAddress("MCJet_recojet_nrecodtrs", &MCJet_recojet_nrecodtrs, &b_MCJet_recojet_nrecodtrs);
  fChain->SetBranchAddress("MCJet_recojet_PX", &MCJet_recojet_PX, &b_MCJet_recojet_PX);
  fChain->SetBranchAddress("MCJet_recojet_PY", &MCJet_recojet_PY, &b_MCJet_recojet_PY);
  fChain->SetBranchAddress("MCJet_recojet_PZ", &MCJet_recojet_PZ, &b_MCJet_recojet_PZ);
  fChain->SetBranchAddress("MCJet_recojet_PE", &MCJet_recojet_PE, &b_MCJet_recojet_PE);
  fChain->SetBranchAddress("MCJet_recojet_ETA", &MCJet_recojet_ETA, &b_MCJet_recojet_ETA);
  fChain->SetBranchAddress("MCJet_recojet_PHI", &MCJet_recojet_PHI, &b_MCJet_recojet_PHI);
  fChain->SetBranchAddress("MCJet_recojet_M", &MCJet_recojet_M, &b_MCJet_recojet_M);
  fChain->SetBranchAddress("MCJet_recojet_XVTX", &MCJet_recojet_XVTX, &b_MCJet_recojet_XVTX);
  fChain->SetBranchAddress("MCJet_recojet_YVTX", &MCJet_recojet_YVTX, &b_MCJet_recojet_YVTX);
  fChain->SetBranchAddress("MCJet_recojet_CPF", &MCJet_recojet_CPF, &b_MCJet_recojet_CPF);
  fChain->SetBranchAddress("MCJet_recojet_MTF", &MCJet_recojet_MTF, &b_MCJet_recojet_MTF);
  fChain->SetBranchAddress("MCJet_recojet_NHasPVInfo", &MCJet_recojet_NHasPVInfo, &b_MCJet_recojet_NHasPVInfo);
  fChain->SetBranchAddress("MCJet_recojet_MPT", &MCJet_recojet_MPT, &b_MCJet_recojet_MPT);
  fChain->SetBranchAddress("MCJet_NDtr", &MCJet_NDtr, &b_MCJet_NDtr);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ID", MCJet_recojet_Dtr_ID, &b_MCJet_recojet_Dtr_ID);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ThreeCharge", MCJet_recojet_Dtr_ThreeCharge, &b_MCJet_recojet_Dtr_ThreeCharge);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ThreeCharge", MCJet_recojet_Dtr_ThreeCharge, &b_MCJet_recojet_Dtr_ThreeCharge);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_IsMeson", MCJet_recojet_Dtr_IsMeson, &b_MCJet_recojet_Dtr_IsMeson);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_IsBaryon", MCJet_recojet_Dtr_IsBaryon, &b_MCJet_recojet_Dtr_IsBaryon);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_IsStrange", MCJet_recojet_Dtr_IsStrange, &b_MCJet_recojet_Dtr_IsStrange);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PIDe", MCJet_recojet_Dtr_PIDe, &b_MCJet_recojet_Dtr_PIDe);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PIDmu", MCJet_recojet_Dtr_PIDmu, &b_MCJet_recojet_Dtr_PIDmu);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PIDk", MCJet_recojet_Dtr_PIDk, &b_MCJet_recojet_Dtr_PIDk);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PIDp", MCJet_recojet_Dtr_PIDp, &b_MCJet_recojet_Dtr_PIDp);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PIDpi", MCJet_recojet_Dtr_PIDpi, &b_MCJet_recojet_Dtr_PIDpi);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ProbNNe", MCJet_recojet_Dtr_ProbNNe, &b_MCJet_recojet_Dtr_ProbNNe);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ProbNNmu", MCJet_recojet_Dtr_ProbNNmu, &b_MCJet_recojet_Dtr_ProbNNmu);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ProbNNk", MCJet_recojet_Dtr_ProbNNk, &b_MCJet_recojet_Dtr_ProbNNk);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ProbNNp", MCJet_recojet_Dtr_ProbNNp, &b_MCJet_recojet_Dtr_ProbNNp);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ProbNNpi", MCJet_recojet_Dtr_ProbNNpi, &b_MCJet_recojet_Dtr_ProbNNpi);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ProbNNghost", MCJet_recojet_Dtr_ProbNNghost, &b_MCJet_recojet_Dtr_ProbNNghost);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PX", MCJet_recojet_Dtr_PX, &b_MCJet_recojet_Dtr_PX);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PY", MCJet_recojet_Dtr_PY, &b_MCJet_recojet_Dtr_PY);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PZ", MCJet_recojet_Dtr_PZ, &b_MCJet_recojet_Dtr_PZ);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_E", MCJet_recojet_Dtr_E, &b_MCJet_recojet_Dtr_E);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_P", MCJet_recojet_Dtr_P, &b_MCJet_recojet_Dtr_P);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PT", MCJet_recojet_Dtr_PT, &b_MCJet_recojet_Dtr_PT);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_M", MCJet_recojet_Dtr_M, &b_MCJet_recojet_Dtr_M);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PHI", MCJet_recojet_Dtr_PHI, &b_MCJet_recojet_Dtr_PHI);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_ETA", MCJet_recojet_Dtr_ETA, &b_MCJet_recojet_Dtr_ETA);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_VX", MCJet_recojet_Dtr_VX, &b_MCJet_recojet_Dtr_VX);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_VY", MCJet_recojet_Dtr_VY, &b_MCJet_recojet_Dtr_VY);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_VZ", MCJet_recojet_Dtr_VZ, &b_MCJet_recojet_Dtr_VZ);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_Z", MCJet_recojet_Dtr_Z, &b_MCJet_recojet_Dtr_Z);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_JT", MCJet_recojet_Dtr_JT, &b_MCJet_recojet_Dtr_JT);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_R", MCJet_recojet_Dtr_R, &b_MCJet_recojet_Dtr_R);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_TrackType", MCJet_recojet_Dtr_TrackType, &b_MCJet_recojet_Dtr_TrackType);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_TrackKey", MCJet_recojet_Dtr_TrackKey, &b_MCJet_recojet_Dtr_TrackKey);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_TrackChi2", MCJet_recojet_Dtr_TrackChi2, &b_MCJet_recojet_Dtr_TrackChi2);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_TrackNDF", MCJet_recojet_Dtr_TrackNDF, &b_MCJet_recojet_Dtr_TrackNDF);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_TrackPCHI2", MCJet_recojet_Dtr_TrackPCHI2, &b_MCJet_recojet_Dtr_TrackPCHI2);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PFType", MCJet_recojet_Dtr_PFType, &b_MCJet_recojet_Dtr_PFType);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PFDaughterType", MCJet_recojet_Dtr_PFDaughterType, &b_MCJet_recojet_Dtr_PFDaughterType);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PFclustype", MCJet_recojet_Dtr_PFclustype, &b_MCJet_recojet_Dtr_PFclustype);
  fChain->SetBranchAddress("MCJet_recojet_Dtr_PFclusE", MCJet_recojet_Dtr_PFclusE, &b_MCJet_recojet_Dtr_PFclusE);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_PX", MCJet_recojet_matchedtruthdtr_PX, &b_MCJet_recojet_matchedtruthdtr_PX);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_PY", MCJet_recojet_matchedtruthdtr_PY, &b_MCJet_recojet_matchedtruthdtr_PY);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_PZ", MCJet_recojet_matchedtruthdtr_PZ, &b_MCJet_recojet_matchedtruthdtr_PZ);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_E", MCJet_recojet_matchedtruthdtr_E, &b_MCJet_recojet_matchedtruthdtr_E);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_P", MCJet_recojet_matchedtruthdtr_P, &b_MCJet_recojet_matchedtruthdtr_P);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_PT", MCJet_recojet_matchedtruthdtr_PT, &b_MCJet_recojet_matchedtruthdtr_PT);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_M", MCJet_recojet_matchedtruthdtr_M, &b_MCJet_recojet_matchedtruthdtr_M);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_PHI", MCJet_recojet_matchedtruthdtr_PHI, &b_MCJet_recojet_matchedtruthdtr_PHI);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_ETA", MCJet_recojet_matchedtruthdtr_ETA, &b_MCJet_recojet_matchedtruthdtr_ETA);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_XVTX", MCJet_recojet_matchedtruthdtr_XVTX, &b_MCJet_recojet_matchedtruthdtr_XVTX);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_YVTX", MCJet_recojet_matchedtruthdtr_YVTX, &b_MCJet_recojet_matchedtruthdtr_YVTX);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_ZVTX", MCJet_recojet_matchedtruthdtr_ZVTX, &b_MCJet_recojet_matchedtruthdtr_ZVTX);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_PID", MCJet_recojet_matchedtruthdtr_PID, &b_MCJet_recojet_matchedtruthdtr_PID);
  fChain->SetBranchAddress("MCJet_recojet_matchedtruthdtr_KEY", MCJet_recojet_matchedtruthdtr_KEY, &b_MCJet_recojet_matchedtruthdtr_KEY);

  fChain->SetBranchAddress("MCJet_PT", &MCJet_PT, &b_MCJet_PT);
  fChain->SetBranchAddress("MCJet_P", &MCJet_P, &b_MCJet_P);
  fChain->SetBranchAddress("MCJet_PX", &MCJet_PX, &b_MCJet_PX);
  fChain->SetBranchAddress("MCJet_PY", &MCJet_PY, &b_MCJet_PY);
  fChain->SetBranchAddress("MCJet_PZ", &MCJet_PZ, &b_MCJet_PZ);
  fChain->SetBranchAddress("MCJet_PE", &MCJet_PE, &b_MCJet_PE);
  fChain->SetBranchAddress("MCJet_ETA", &MCJet_ETA, &b_MCJet_ETA);
  fChain->SetBranchAddress("MCJet_PHI", &MCJet_PHI, &b_MCJet_PHI);
  fChain->SetBranchAddress("MCJet_M", &MCJet_M, &b_MCJet_M);
  fChain->SetBranchAddress("MCJet_FLAV", &MCJet_FLAV, &b_MCJet_FLAV);
  fChain->SetBranchAddress("MCJet_Dtr_nmcdtrs", &MCJet_Dtr_nmcdtrs, &b_MCJet_Dtr_nmcdtrs);
  fChain->SetBranchAddress("MCJet_Dtr_nmcdtr_recomatches", &MCJet_Dtr_nmcdtr_recomatches, &b_MCJet_Dtr_nmcdtr_recomatches);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_p", MCJet_Dtr_mcrecomatch_p, &b_MCJet_Dtr_mcrecomatch_p);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_pt", MCJet_Dtr_mcrecomatch_pt, &b_MCJet_Dtr_mcrecomatch_pt);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_pide", MCJet_Dtr_mcrecomatch_pide, &b_MCJet_Dtr_mcrecomatch_pide);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_pidmu", MCJet_Dtr_mcrecomatch_pidmu, &b_MCJet_Dtr_mcrecomatch_pidmu);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_pidk", MCJet_Dtr_mcrecomatch_pidk, &b_MCJet_Dtr_mcrecomatch_pidk);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_pidp", MCJet_Dtr_mcrecomatch_pidp, &b_MCJet_Dtr_mcrecomatch_pidp);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_pidpi", MCJet_Dtr_mcrecomatch_pidpi, &b_MCJet_Dtr_mcrecomatch_pidpi);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_probnn_e", MCJet_Dtr_mcrecomatch_probnn_e, &b_MCJet_Dtr_mcrecomatch_probnn_e);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_probnn_mu", MCJet_Dtr_mcrecomatch_probnn_mu, &b_MCJet_Dtr_mcrecomatch_probnn_mu);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_probnn_k", MCJet_Dtr_mcrecomatch_probnn_k, &b_MCJet_Dtr_mcrecomatch_probnn_k);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_probnn_p", MCJet_Dtr_mcrecomatch_probnn_p, &b_MCJet_Dtr_mcrecomatch_probnn_p);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_probnn_pi", MCJet_Dtr_mcrecomatch_probnn_pi, &b_MCJet_Dtr_mcrecomatch_probnn_pi);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_probnn_ghost", MCJet_Dtr_mcrecomatch_probnn_ghost, &b_MCJet_Dtr_mcrecomatch_probnn_ghost);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_track_key", MCJet_Dtr_mcrecomatch_track_key, &b_MCJet_Dtr_mcrecomatch_track_key);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_track_ndf", MCJet_Dtr_mcrecomatch_track_ndf, &b_MCJet_Dtr_mcrecomatch_track_ndf);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_track_type", MCJet_Dtr_mcrecomatch_track_type, &b_MCJet_Dtr_mcrecomatch_track_type);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_track_pchi2", MCJet_Dtr_mcrecomatch_track_pchi2, &b_MCJet_Dtr_mcrecomatch_track_pchi2);
  fChain->SetBranchAddress("MCJet_Dtr_mcrecomatch_track_chi2", MCJet_Dtr_mcrecomatch_track_chi2, &b_MCJet_Dtr_mcrecomatch_track_chi2);
  fChain->SetBranchAddress("MCJet_Dtr_ID", MCJet_Dtr_ID, &b_MCJet_Dtr_ID);
  fChain->SetBranchAddress("MCJet_Dtr_Key", MCJet_Dtr_Key, &b_MCJet_Dtr_Key);
  fChain->SetBranchAddress("MCJet_Dtr_ThreeCharge", MCJet_Dtr_ThreeCharge, &b_MCJet_Dtr_ThreeCharge);
  fChain->SetBranchAddress("MCJet_Dtr_IsMeson", MCJet_Dtr_IsMeson, &b_MCJet_Dtr_IsMeson);
  fChain->SetBranchAddress("MCJet_Dtr_IsBaryon", MCJet_Dtr_IsBaryon, &b_MCJet_Dtr_IsBaryon);
  fChain->SetBranchAddress("MCJet_Dtr_IsStrange", MCJet_Dtr_IsStrange, &b_MCJet_Dtr_IsStrange);
  fChain->SetBranchAddress("MCJet_Dtr_PX", MCJet_Dtr_PX, &b_MCJet_Dtr_PX);
  fChain->SetBranchAddress("MCJet_Dtr_PY", MCJet_Dtr_PY, &b_MCJet_Dtr_PY);
  fChain->SetBranchAddress("MCJet_Dtr_PZ", MCJet_Dtr_PZ, &b_MCJet_Dtr_PZ);
  fChain->SetBranchAddress("MCJet_Dtr_E", MCJet_Dtr_E, &b_MCJet_Dtr_E);
  fChain->SetBranchAddress("MCJet_Dtr_P", MCJet_Dtr_P, &b_MCJet_Dtr_P);
  fChain->SetBranchAddress("MCJet_Dtr_R", MCJet_Dtr_R, &b_MCJet_Dtr_R);
  fChain->SetBranchAddress("MCJet_Dtr_JT", MCJet_Dtr_JT, &b_MCJet_Dtr_JT);
  fChain->SetBranchAddress("MCJet_Dtr_Z", MCJet_Dtr_Z, &b_MCJet_Dtr_Z);
  fChain->SetBranchAddress("MCJet_Dtr_PT", MCJet_Dtr_PT, &b_MCJet_Dtr_PT);
  fChain->SetBranchAddress("MCJet_Dtr_M", MCJet_Dtr_M, &b_MCJet_Dtr_M);
  fChain->SetBranchAddress("MCJet_Dtr_PHI", MCJet_Dtr_PHI, &b_MCJet_Dtr_PHI);
  fChain->SetBranchAddress("MCJet_Dtr_ETA", MCJet_Dtr_ETA, &b_MCJet_Dtr_ETA);
  fChain->SetBranchAddress("MCJet_Dtr_VX", MCJet_Dtr_VX, &b_MCJet_Dtr_VX);
  fChain->SetBranchAddress("MCJet_Dtr_VY", MCJet_Dtr_VY, &b_MCJet_Dtr_VY);
  fChain->SetBranchAddress("MCJet_Dtr_VZ", MCJet_Dtr_VZ, &b_MCJet_Dtr_VZ);
  fChain->SetBranchAddress("MCJet_Dtr_TopMotherID", MCJet_Dtr_TopMotherID, &b_MCJet_Dtr_TopMotherID);
  fChain->SetBranchAddress("CELLjet_TRUEP_E", &CELLjet_TRUEP_E, &b_CELLjet_TRUEP_E);
  fChain->SetBranchAddress("CELLjet_TRUEP_X", &CELLjet_TRUEP_X, &b_CELLjet_TRUEP_X);
  fChain->SetBranchAddress("CELLjet_TRUEP_Y", &CELLjet_TRUEP_Y, &b_CELLjet_TRUEP_Y);
  fChain->SetBranchAddress("CELLjet_TRUEP_Z", &CELLjet_TRUEP_Z, &b_CELLjet_TRUEP_Z);
  fChain->SetBranchAddress("CELLjet_TRUEPT", &CELLjet_TRUEPT, &b_CELLjet_TRUEPT);
  fChain->SetBranchAddress("CELLjet_TRUEORIGINVERTEX_X", &CELLjet_TRUEORIGINVERTEX_X, &b_CELLjet_TRUEORIGINVERTEX_X);
  fChain->SetBranchAddress("CELLjet_TRUEORIGINVERTEX_Y", &CELLjet_TRUEORIGINVERTEX_Y, &b_CELLjet_TRUEORIGINVERTEX_Y);
  fChain->SetBranchAddress("CELLjet_TRUEORIGINVERTEX_Z", &CELLjet_TRUEORIGINVERTEX_Z, &b_CELLjet_TRUEORIGINVERTEX_Z);
  fChain->SetBranchAddress("CELLjet_TRUEENDVERTEX_X", &CELLjet_TRUEENDVERTEX_X, &b_CELLjet_TRUEENDVERTEX_X);
  fChain->SetBranchAddress("CELLjet_TRUEENDVERTEX_Y", &CELLjet_TRUEENDVERTEX_Y, &b_CELLjet_TRUEENDVERTEX_Y);
  fChain->SetBranchAddress("CELLjet_TRUEENDVERTEX_Z", &CELLjet_TRUEENDVERTEX_Z, &b_CELLjet_TRUEENDVERTEX_Z);
  fChain->SetBranchAddress("CELLjet_TRUEISSTABLE", &CELLjet_TRUEISSTABLE, &b_CELLjet_TRUEISSTABLE);
  fChain->SetBranchAddress("CELLjet_TRUETAU", &CELLjet_TRUETAU, &b_CELLjet_TRUETAU);
  fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
  fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
  fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
  fChain->SetBranchAddress("EVT_Int__I", &EVT_Int__I, &b_EVT_Int__I);
  fChain->SetBranchAddress("EVT_Int__Mean", &EVT_Int__Mean, &b_EVT_Int__Mean);
  fChain->SetBranchAddress("EVT_Int__Prob", &EVT_Int__Prob, &b_EVT_Int__Prob);
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

Bool_t BjetMCTree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void BjetMCTree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain)
    return;
  fChain->Show(entry);
}
Int_t BjetMCTree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
#endif // #ifdef BjetMCTree_cxx
