//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb  8 21:09:48 2024 by ROOT version 6.24/08
// from TChain Jpsi2MuMu/DecayTree/
//////////////////////////////////////////////////////////

#ifndef JpsiTree_h
#define JpsiTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class JpsiTree
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxJ_psi_1S_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxJ_psi_1S_OWNPV_COV = 1;
   static constexpr Int_t kMaxJ_psi_1S0_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxJ_psi_1S0_OWNPV_COV = 1;

   // Declaration of leaf types
   Double_t J_psi_1S_ENDVERTEX_X;
   Double_t J_psi_1S_ENDVERTEX_Y;
   Double_t J_psi_1S_ENDVERTEX_Z;
   Double_t J_psi_1S_ENDVERTEX_XERR;
   Double_t J_psi_1S_ENDVERTEX_YERR;
   Double_t J_psi_1S_ENDVERTEX_ZERR;
   Double_t J_psi_1S_ENDVERTEX_CHI2;
   Int_t J_psi_1S_ENDVERTEX_NDOF;
   Float_t J_psi_1S_ENDVERTEX_COV_[3][3];
   Double_t J_psi_1S_OWNPV_X;
   Double_t J_psi_1S_OWNPV_Y;
   Double_t J_psi_1S_OWNPV_Z;
   Double_t J_psi_1S_OWNPV_XERR;
   Double_t J_psi_1S_OWNPV_YERR;
   Double_t J_psi_1S_OWNPV_ZERR;
   Double_t J_psi_1S_OWNPV_CHI2;
   Int_t J_psi_1S_OWNPV_NDOF;
   Float_t J_psi_1S_OWNPV_COV_[3][3];
   Double_t J_psi_1S_IP_OWNPV;
   Double_t J_psi_1S_IPCHI2_OWNPV;
   Double_t J_psi_1S_FD_OWNPV;
   Double_t J_psi_1S_FDCHI2_OWNPV;
   Double_t J_psi_1S_DIRA_OWNPV;
   Double_t J_psi_1S_P;
   Double_t J_psi_1S_PT;
   Double_t J_psi_1S_PE;
   Double_t J_psi_1S_PX;
   Double_t J_psi_1S_PY;
   Double_t J_psi_1S_PZ;
   Double_t J_psi_1S_MM;
   Double_t J_psi_1S_MMERR;
   Double_t J_psi_1S_M;
   Int_t J_psi_1S_ID;
   Bool_t J_psi_1S_L0Global_Dec;
   Bool_t J_psi_1S_L0Global_TIS;
   Bool_t J_psi_1S_L0Global_TOS;
   Bool_t J_psi_1S_Hlt1Global_Dec;
   Bool_t J_psi_1S_Hlt1Global_TIS;
   Bool_t J_psi_1S_Hlt1Global_TOS;
   Bool_t J_psi_1S_Hlt1Phys_Dec;
   Bool_t J_psi_1S_Hlt1Phys_TIS;
   Bool_t J_psi_1S_Hlt1Phys_TOS;
   Bool_t J_psi_1S_Hlt2Global_Dec;
   Bool_t J_psi_1S_Hlt2Global_TIS;
   Bool_t J_psi_1S_Hlt2Global_TOS;
   Bool_t J_psi_1S_Hlt2Phys_Dec;
   Bool_t J_psi_1S_Hlt2Phys_TIS;
   Bool_t J_psi_1S_Hlt2Phys_TOS;
   Bool_t J_psi_1S_L0MuonDecision_Dec;
   Bool_t J_psi_1S_L0MuonDecision_TIS;
   Bool_t J_psi_1S_L0MuonDecision_TOS;
   Bool_t J_psi_1S_L0PhotonDecision_Dec;
   Bool_t J_psi_1S_L0PhotonDecision_TIS;
   Bool_t J_psi_1S_L0PhotonDecision_TOS;
   Bool_t J_psi_1S_L0HadronDecision_Dec;
   Bool_t J_psi_1S_L0HadronDecision_TIS;
   Bool_t J_psi_1S_L0HadronDecision_TOS;
   Bool_t J_psi_1S_L0ElectronDecision_Dec;
   Bool_t J_psi_1S_L0ElectronDecision_TIS;
   Bool_t J_psi_1S_L0ElectronDecision_TOS;
   Bool_t J_psi_1S_L0DiMuonDecision_Dec;
   Bool_t J_psi_1S_L0DiMuonDecision_TIS;
   Bool_t J_psi_1S_L0DiMuonDecision_TOS;
   Bool_t J_psi_1S_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t J_psi_1S_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t J_psi_1S_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t J_psi_1S_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t J_psi_1S_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t J_psi_1S_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS;
   Double_t J_psi_1S0_ENDVERTEX_X;
   Double_t J_psi_1S0_ENDVERTEX_Y;
   Double_t J_psi_1S0_ENDVERTEX_Z;
   Double_t J_psi_1S0_ENDVERTEX_XERR;
   Double_t J_psi_1S0_ENDVERTEX_YERR;
   Double_t J_psi_1S0_ENDVERTEX_ZERR;
   Double_t J_psi_1S0_ENDVERTEX_CHI2;
   Int_t J_psi_1S0_ENDVERTEX_NDOF;
   Float_t J_psi_1S0_ENDVERTEX_COV_[3][3];
   Double_t J_psi_1S0_OWNPV_X;
   Double_t J_psi_1S0_OWNPV_Y;
   Double_t J_psi_1S0_OWNPV_Z;
   Double_t J_psi_1S0_OWNPV_XERR;
   Double_t J_psi_1S0_OWNPV_YERR;
   Double_t J_psi_1S0_OWNPV_ZERR;
   Double_t J_psi_1S0_OWNPV_CHI2;
   Int_t J_psi_1S0_OWNPV_NDOF;
   Float_t J_psi_1S0_OWNPV_COV_[3][3];
   Double_t J_psi_1S0_IP_OWNPV;
   Double_t J_psi_1S0_IPCHI2_OWNPV;
   Double_t J_psi_1S0_FD_OWNPV;
   Double_t J_psi_1S0_FDCHI2_OWNPV;
   Double_t J_psi_1S0_DIRA_OWNPV;
   Double_t J_psi_1S0_P;
   Double_t J_psi_1S0_PT;
   Double_t J_psi_1S0_PE;
   Double_t J_psi_1S0_PX;
   Double_t J_psi_1S0_PY;
   Double_t J_psi_1S0_PZ;
   Double_t J_psi_1S0_MM;
   Double_t J_psi_1S0_MMERR;
   Double_t J_psi_1S0_M;
   Int_t J_psi_1S0_ID;
   Bool_t J_psi_1S0_L0Global_Dec;
   Bool_t J_psi_1S0_L0Global_TIS;
   Bool_t J_psi_1S0_L0Global_TOS;
   Bool_t J_psi_1S0_Hlt1Global_Dec;
   Bool_t J_psi_1S0_Hlt1Global_TIS;
   Bool_t J_psi_1S0_Hlt1Global_TOS;
   Bool_t J_psi_1S0_Hlt1Phys_Dec;
   Bool_t J_psi_1S0_Hlt1Phys_TIS;
   Bool_t J_psi_1S0_Hlt1Phys_TOS;
   Bool_t J_psi_1S0_Hlt2Global_Dec;
   Bool_t J_psi_1S0_Hlt2Global_TIS;
   Bool_t J_psi_1S0_Hlt2Global_TOS;
   Bool_t J_psi_1S0_Hlt2Phys_Dec;
   Bool_t J_psi_1S0_Hlt2Phys_TIS;
   Bool_t J_psi_1S0_Hlt2Phys_TOS;
   Bool_t J_psi_1S0_L0MuonDecision_Dec;
   Bool_t J_psi_1S0_L0MuonDecision_TIS;
   Bool_t J_psi_1S0_L0MuonDecision_TOS;
   Bool_t J_psi_1S0_L0DiMuonDecision_Dec;
   Bool_t J_psi_1S0_L0DiMuonDecision_TIS;
   Bool_t J_psi_1S0_L0DiMuonDecision_TOS;
   Bool_t J_psi_1S0_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t J_psi_1S0_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t J_psi_1S0_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedDecision_Dec;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedDecision_TIS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedDecision_TOS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_Dec;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TIS;
   Bool_t J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TOS;
   Bool_t J_psi_1S0_Hlt2DiMuonJPsiDecision_Dec;
   Bool_t J_psi_1S0_Hlt2DiMuonJPsiDecision_TIS;
   Bool_t J_psi_1S0_Hlt2DiMuonJPsiDecision_TOS;
   Bool_t J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_Dec;
   Bool_t J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TIS;
   Bool_t J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TOS;
   UInt_t nCandidate;
   ULong64_t totCandidates;
   ULong64_t EventInSequence;
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
   UInt_t StrippingBs2MuMuLinesBu2JPsiKFullLineDecision;

   // List of branches
   TBranch *b_J_psi_1S_ENDVERTEX_X;                          //!
   TBranch *b_J_psi_1S_ENDVERTEX_Y;                          //!
   TBranch *b_J_psi_1S_ENDVERTEX_Z;                          //!
   TBranch *b_J_psi_1S_ENDVERTEX_XERR;                       //!
   TBranch *b_J_psi_1S_ENDVERTEX_YERR;                       //!
   TBranch *b_J_psi_1S_ENDVERTEX_ZERR;                       //!
   TBranch *b_J_psi_1S_ENDVERTEX_CHI2;                       //!
   TBranch *b_J_psi_1S_ENDVERTEX_NDOF;                       //!
   TBranch *b_J_psi_1S_ENDVERTEX_COV_;                       //!
   TBranch *b_J_psi_1S_OWNPV_X;                              //!
   TBranch *b_J_psi_1S_OWNPV_Y;                              //!
   TBranch *b_J_psi_1S_OWNPV_Z;                              //!
   TBranch *b_J_psi_1S_OWNPV_XERR;                           //!
   TBranch *b_J_psi_1S_OWNPV_YERR;                           //!
   TBranch *b_J_psi_1S_OWNPV_ZERR;                           //!
   TBranch *b_J_psi_1S_OWNPV_CHI2;                           //!
   TBranch *b_J_psi_1S_OWNPV_NDOF;                           //!
   TBranch *b_J_psi_1S_OWNPV_COV_;                           //!
   TBranch *b_J_psi_1S_IP_OWNPV;                             //!
   TBranch *b_J_psi_1S_IPCHI2_OWNPV;                         //!
   TBranch *b_J_psi_1S_FD_OWNPV;                             //!
   TBranch *b_J_psi_1S_FDCHI2_OWNPV;                         //!
   TBranch *b_J_psi_1S_DIRA_OWNPV;                           //!
   TBranch *b_J_psi_1S_P;                                    //!
   TBranch *b_J_psi_1S_PT;                                   //!
   TBranch *b_J_psi_1S_PE;                                   //!
   TBranch *b_J_psi_1S_PX;                                   //!
   TBranch *b_J_psi_1S_PY;                                   //!
   TBranch *b_J_psi_1S_PZ;                                   //!
   TBranch *b_J_psi_1S_MM;                                   //!
   TBranch *b_J_psi_1S_MMERR;                                //!
   TBranch *b_J_psi_1S_M;                                    //!
   TBranch *b_J_psi_1S_ID;                                   //!
   TBranch *b_J_psi_1S_L0Global_Dec;                         //!
   TBranch *b_J_psi_1S_L0Global_TIS;                         //!
   TBranch *b_J_psi_1S_L0Global_TOS;                         //!
   TBranch *b_J_psi_1S_Hlt1Global_Dec;                       //!
   TBranch *b_J_psi_1S_Hlt1Global_TIS;                       //!
   TBranch *b_J_psi_1S_Hlt1Global_TOS;                       //!
   TBranch *b_J_psi_1S_Hlt1Phys_Dec;                         //!
   TBranch *b_J_psi_1S_Hlt1Phys_TIS;                         //!
   TBranch *b_J_psi_1S_Hlt1Phys_TOS;                         //!
   TBranch *b_J_psi_1S_Hlt2Global_Dec;                       //!
   TBranch *b_J_psi_1S_Hlt2Global_TIS;                       //!
   TBranch *b_J_psi_1S_Hlt2Global_TOS;                       //!
   TBranch *b_J_psi_1S_Hlt2Phys_Dec;                         //!
   TBranch *b_J_psi_1S_Hlt2Phys_TIS;                         //!
   TBranch *b_J_psi_1S_Hlt2Phys_TOS;                         //!
   TBranch *b_J_psi_1S_L0MuonDecision_Dec;                   //!
   TBranch *b_J_psi_1S_L0MuonDecision_TIS;                   //!
   TBranch *b_J_psi_1S_L0MuonDecision_TOS;                   //!
   TBranch *b_J_psi_1S_L0PhotonDecision_Dec;                 //!
   TBranch *b_J_psi_1S_L0PhotonDecision_TIS;                 //!
   TBranch *b_J_psi_1S_L0PhotonDecision_TOS;                 //!
   TBranch *b_J_psi_1S_L0HadronDecision_Dec;                 //!
   TBranch *b_J_psi_1S_L0HadronDecision_TIS;                 //!
   TBranch *b_J_psi_1S_L0HadronDecision_TOS;                 //!
   TBranch *b_J_psi_1S_L0ElectronDecision_Dec;               //!
   TBranch *b_J_psi_1S_L0ElectronDecision_TIS;               //!
   TBranch *b_J_psi_1S_L0ElectronDecision_TOS;               //!
   TBranch *b_J_psi_1S_L0DiMuonDecision_Dec;                 //!
   TBranch *b_J_psi_1S_L0DiMuonDecision_TIS;                 //!
   TBranch *b_J_psi_1S_L0DiMuonDecision_TOS;                 //!
   TBranch *b_J_psi_1S_Hlt1DiMuonHighMassDecision_Dec;       //!
   TBranch *b_J_psi_1S_Hlt1DiMuonHighMassDecision_TIS;       //!
   TBranch *b_J_psi_1S_Hlt1DiMuonHighMassDecision_TOS;       //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedDecision_Dec;       //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedDecision_TIS;       //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedDecision_TOS;       //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_Dec;  //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TIS;  //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TOS;  //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_Dec;  //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TIS;  //!
   TBranch *b_J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TOS;  //!
   TBranch *b_J_psi_1S_Hlt2DiMuonJPsiDecision_Dec;           //!
   TBranch *b_J_psi_1S_Hlt2DiMuonJPsiDecision_TIS;           //!
   TBranch *b_J_psi_1S_Hlt2DiMuonJPsiDecision_TOS;           //!
   TBranch *b_J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_Dec;     //!
   TBranch *b_J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TIS;     //!
   TBranch *b_J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS;     //!
   TBranch *b_J_psi_1S0_ENDVERTEX_X;                         //!
   TBranch *b_J_psi_1S0_ENDVERTEX_Y;                         //!
   TBranch *b_J_psi_1S0_ENDVERTEX_Z;                         //!
   TBranch *b_J_psi_1S0_ENDVERTEX_XERR;                      //!
   TBranch *b_J_psi_1S0_ENDVERTEX_YERR;                      //!
   TBranch *b_J_psi_1S0_ENDVERTEX_ZERR;                      //!
   TBranch *b_J_psi_1S0_ENDVERTEX_CHI2;                      //!
   TBranch *b_J_psi_1S0_ENDVERTEX_NDOF;                      //!
   TBranch *b_J_psi_1S0_ENDVERTEX_COV_;                      //!
   TBranch *b_J_psi_1S0_OWNPV_X;                             //!
   TBranch *b_J_psi_1S0_OWNPV_Y;                             //!
   TBranch *b_J_psi_1S0_OWNPV_Z;                             //!
   TBranch *b_J_psi_1S0_OWNPV_XERR;                          //!
   TBranch *b_J_psi_1S0_OWNPV_YERR;                          //!
   TBranch *b_J_psi_1S0_OWNPV_ZERR;                          //!
   TBranch *b_J_psi_1S0_OWNPV_CHI2;                          //!
   TBranch *b_J_psi_1S0_OWNPV_NDOF;                          //!
   TBranch *b_J_psi_1S0_OWNPV_COV_;                          //!
   TBranch *b_J_psi_1S0_IP_OWNPV;                            //!
   TBranch *b_J_psi_1S0_IPCHI2_OWNPV;                        //!
   TBranch *b_J_psi_1S0_FD_OWNPV;                            //!
   TBranch *b_J_psi_1S0_FDCHI2_OWNPV;                        //!
   TBranch *b_J_psi_1S0_DIRA_OWNPV;                          //!
   TBranch *b_J_psi_1S0_P;                                   //!
   TBranch *b_J_psi_1S0_PT;                                  //!
   TBranch *b_J_psi_1S0_PE;                                  //!
   TBranch *b_J_psi_1S0_PX;                                  //!
   TBranch *b_J_psi_1S0_PY;                                  //!
   TBranch *b_J_psi_1S0_PZ;                                  //!
   TBranch *b_J_psi_1S0_MM;                                  //!
   TBranch *b_J_psi_1S0_MMERR;                               //!
   TBranch *b_J_psi_1S0_M;                                   //!
   TBranch *b_J_psi_1S0_ID;                                  //!
   TBranch *b_J_psi_1S0_L0Global_Dec;                        //!
   TBranch *b_J_psi_1S0_L0Global_TIS;                        //!
   TBranch *b_J_psi_1S0_L0Global_TOS;                        //!
   TBranch *b_J_psi_1S0_Hlt1Global_Dec;                      //!
   TBranch *b_J_psi_1S0_Hlt1Global_TIS;                      //!
   TBranch *b_J_psi_1S0_Hlt1Global_TOS;                      //!
   TBranch *b_J_psi_1S0_Hlt1Phys_Dec;                        //!
   TBranch *b_J_psi_1S0_Hlt1Phys_TIS;                        //!
   TBranch *b_J_psi_1S0_Hlt1Phys_TOS;                        //!
   TBranch *b_J_psi_1S0_Hlt2Global_Dec;                      //!
   TBranch *b_J_psi_1S0_Hlt2Global_TIS;                      //!
   TBranch *b_J_psi_1S0_Hlt2Global_TOS;                      //!
   TBranch *b_J_psi_1S0_Hlt2Phys_Dec;                        //!
   TBranch *b_J_psi_1S0_Hlt2Phys_TIS;                        //!
   TBranch *b_J_psi_1S0_Hlt2Phys_TOS;                        //!
   TBranch *b_J_psi_1S0_L0MuonDecision_Dec;                  //!
   TBranch *b_J_psi_1S0_L0MuonDecision_TIS;                  //!
   TBranch *b_J_psi_1S0_L0MuonDecision_TOS;                  //!
   TBranch *b_J_psi_1S0_L0DiMuonDecision_Dec;                //!
   TBranch *b_J_psi_1S0_L0DiMuonDecision_TIS;                //!
   TBranch *b_J_psi_1S0_L0DiMuonDecision_TOS;                //!
   TBranch *b_J_psi_1S0_Hlt1DiMuonHighMassDecision_Dec;      //!
   TBranch *b_J_psi_1S0_Hlt1DiMuonHighMassDecision_TIS;      //!
   TBranch *b_J_psi_1S0_Hlt1DiMuonHighMassDecision_TOS;      //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedDecision_Dec;      //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedDecision_TIS;      //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedDecision_TOS;      //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_Dec; //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TIS; //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TOS; //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_Dec;  //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TIS;  //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TOS;  //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_Dec; //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TIS; //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TOS; //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonJPsiDecision_Dec;          //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonJPsiDecision_TIS;          //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonJPsiDecision_TOS;          //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_Dec;    //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TIS;    //!
   TBranch *b_J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TOS;    //!
   TBranch *b_nCandidate;                                    //!
   TBranch *b_totCandidates;                                 //!
   TBranch *b_EventInSequence;                               //!
   TBranch *b_runNumber;                                     //!
   TBranch *b_eventNumber;                                   //!
   TBranch *b_BCID;                                          //!
   TBranch *b_BCType;                                        //!
   TBranch *b_OdinTCK;                                       //!
   TBranch *b_L0DUTCK;                                       //!
   TBranch *b_HLT1TCK;                                       //!
   TBranch *b_HLT2TCK;                                       //!
   TBranch *b_GpsTime;                                       //!
   TBranch *b_Polarity;                                      //!
   TBranch *b_nPVs;                                          //!
   TBranch *b_nTracks;                                       //!
   TBranch *b_nLongTracks;                                   //!
   TBranch *b_nDownstreamTracks;                             //!
   TBranch *b_nUpstreamTracks;                               //!
   TBranch *b_nVeloTracks;                                   //!
   TBranch *b_nTTracks;                                      //!
   TBranch *b_nBackTracks;                                   //!
   TBranch *b_nRich1Hits;                                    //!
   TBranch *b_nRich2Hits;                                    //!
   TBranch *b_nVeloClusters;                                 //!
   TBranch *b_nITClusters;                                   //!
   TBranch *b_nTTClusters;                                   //!
   TBranch *b_nOTClusters;                                   //!
   TBranch *b_nSPDHits;                                      //!
   TBranch *b_nMuonCoordsS0;                                 //!
   TBranch *b_nMuonCoordsS1;                                 //!
   TBranch *b_nMuonCoordsS2;                                 //!
   TBranch *b_nMuonCoordsS3;                                 //!
   TBranch *b_nMuonCoordsS4;                                 //!
   TBranch *b_nMuonTracks;                                   //!
   TBranch *b_StrippingBs2MuMuLinesBu2JPsiKFullLineDecision; //!

   JpsiTree(TTree *tree = 0, int dataset = 91599, bool isData = true);
   virtual ~JpsiTree();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef JpsiTree_cxx
JpsiTree::JpsiTree(TTree *tree, int dataset = 93539, bool isData = true) : fChain(0)
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
      f->GetObject("Jpsi2MuMu/DecayTree", tree);

#else  // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain *chain = new TChain("Jpsi2MuMu/DecayTree", "");
      TString data_dir = "../../data/";

      if (dataset == 61590)
      {
         if (isData)
         {

            chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_Random_2016_MD_02122024.root/Jpsi2MuMu/DecayTree");

            // chain->Add(data_dir + "Bjet_betas_2016_MD_02122024_url1245_0_800.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_betas_2016_MD_02122024_url1245_800_1000.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2016_MD_12192023_url1159_0_800.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2016_MD_12192023_url1159_800_1023.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            chain->Add(data_dir + "Jpsi2MuMu_MC_2016_Sim09k_MD_02142024_full.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 61591)
      {
         if (isData)
         {

            chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_Random_2016_MU_02122024.root/Jpsi2MuMu/DecayTree");

            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2016_MU_12192023_url1160_0_800.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2016_MU_12192023_url1160_800_973.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            // chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim09k_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim09l_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim10a_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 61599)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_betas_2016_MU_02122024.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_betas_2016_MD_02122024_url1245_0_800.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_betas_2016_MD_02122024_url1245_800_1000.root/Jpsi2MuMu/DecayTree");
         }
         else
         {

            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim09k_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim09l_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim10a_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim09k_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim09l_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2016_Sim10a_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 71590)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_Random_2017_MD_02122024.root/Jpsi2MuMu/DecayTree");

            // chain->Add(data_dir + "Bjet_betas_2017_MD_02122024.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2017_MD_12192023_url1161_0_721.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            chain->Add(data_dir + "Jpsi2MuMu_2017_Sim09g_MD_02152024.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 71591)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_Random_2017_MU_02122024.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2017_MU_12192023_url1162_0_701.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            chain->Add(data_dir + "Jpsi2MuMu_2018_Sim09g_MD_02162024.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 71599)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_Dimuon_2017_MU_Data.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_Dimuon_2017_MD_Data.root/Jpsi2MuMu/DecayTree");
         }
         else
         {

            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2017_Sim09k_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2017_Sim09l_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2017_Sim10a_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2017_Sim09k_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2017_Sim09l_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2017_Sim10a_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 81590)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_Random_2018_MD_02122024.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2018_MD_12192023_url1163_0_800.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2018_MD_12192023_url1163_800_826.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            chain->Add(data_dir + "Jpsi2MuMu_2018_Sim09g_MD_02162024.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 81591)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_Jpsi2MuMu_Data_Random_2018_MU_02122024.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2018_MD_12192023_url1164_0_800.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_Dimuon_nojetid_DTF_2018_MD_12192023_url1164_800_876.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            // chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim09k_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim09l_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            // chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim10a_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
         }
      }
      else if (dataset == 81599)
      {
         if (isData)
         {
            chain->Add(data_dir + "Bjet_Dimuon_2018_MU_Data.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_Dimuon_2018_MD_Data.root/Jpsi2MuMu/DecayTree");
         }
         else
         {

            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim09k_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim09l_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim10a_MU_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim09k_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim09l_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Bjet_MC_Jpsi2MuMu_2018_Sim10a_MD_02102024_full.root/Jpsi2MuMu/DecayTree");
         }
      }

      if (dataset == 91599)
      {
         // data_dir = "/afs/cern.ch/work/i/ichahrou/DeadCone/DaVinci/Bu2JpsiK/";
         if (isData)
         {
            chain->Add(data_dir + "Jpsi2MuMu_Data_HighPT_2016_MD_01202025.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_Data_HighPT_2016_MU_01212025.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_Data_HighPT_2017_MD_01212025.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_Data_HighPT_2017_MU_01212025.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_Data_HighPT_2018_MD_01212025.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_Data_HighPT_2018_MU_01212025.root/Jpsi2MuMu/DecayTree");
         }
         else
         {
            chain->Add(data_dir + "Jpsi2MuMu_MC_2016_Sim09l_MD_01242025_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_MC_2016_Sim09l_MU_01242025_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_MC_2017_Sim09l_MD_01242025_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_MC_2017_Sim09l_MU_01242025_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_MC_2018_Sim09l_MD_01242025_full.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_MC_2018_Sim09l_MU_01242025_full.root/Jpsi2MuMu/DecayTree");

            chain->Add(data_dir + "Jpsi2MuMu_2017_Sim09g_MD_01252025.root/Jpsi2MuMu/DecayTree");
            chain->Add(data_dir + "Jpsi2MuMu_2017_Sim09g_MU_01252025.root/Jpsi2MuMu/DecayTree");

            // chain->Add(data_dir + "Jpsi2MuMu_MC_2016_Sim09k_MD_02142024_full.root/Jpsi2MuMu/DecayTree");
         }
      }
      tree = chain;
#endif // SINGLE_TREE
   }
   Init(tree);
}

JpsiTree::~JpsiTree()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t JpsiTree::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t JpsiTree::LoadTree(Long64_t entry)
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

void JpsiTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_X", &J_psi_1S_ENDVERTEX_X, &b_J_psi_1S_ENDVERTEX_X);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_Y", &J_psi_1S_ENDVERTEX_Y, &b_J_psi_1S_ENDVERTEX_Y);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_Z", &J_psi_1S_ENDVERTEX_Z, &b_J_psi_1S_ENDVERTEX_Z);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_XERR", &J_psi_1S_ENDVERTEX_XERR, &b_J_psi_1S_ENDVERTEX_XERR);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_YERR", &J_psi_1S_ENDVERTEX_YERR, &b_J_psi_1S_ENDVERTEX_YERR);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_ZERR", &J_psi_1S_ENDVERTEX_ZERR, &b_J_psi_1S_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_CHI2", &J_psi_1S_ENDVERTEX_CHI2, &b_J_psi_1S_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_NDOF", &J_psi_1S_ENDVERTEX_NDOF, &b_J_psi_1S_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("J_psi_1S_ENDVERTEX_COV_", J_psi_1S_ENDVERTEX_COV_, &b_J_psi_1S_ENDVERTEX_COV_);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_X", &J_psi_1S_OWNPV_X, &b_J_psi_1S_OWNPV_X);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_Y", &J_psi_1S_OWNPV_Y, &b_J_psi_1S_OWNPV_Y);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_Z", &J_psi_1S_OWNPV_Z, &b_J_psi_1S_OWNPV_Z);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_XERR", &J_psi_1S_OWNPV_XERR, &b_J_psi_1S_OWNPV_XERR);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_YERR", &J_psi_1S_OWNPV_YERR, &b_J_psi_1S_OWNPV_YERR);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_ZERR", &J_psi_1S_OWNPV_ZERR, &b_J_psi_1S_OWNPV_ZERR);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_CHI2", &J_psi_1S_OWNPV_CHI2, &b_J_psi_1S_OWNPV_CHI2);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_NDOF", &J_psi_1S_OWNPV_NDOF, &b_J_psi_1S_OWNPV_NDOF);
   fChain->SetBranchAddress("J_psi_1S_OWNPV_COV_", J_psi_1S_OWNPV_COV_, &b_J_psi_1S_OWNPV_COV_);
   fChain->SetBranchAddress("J_psi_1S_IP_OWNPV", &J_psi_1S_IP_OWNPV, &b_J_psi_1S_IP_OWNPV);
   fChain->SetBranchAddress("J_psi_1S_IPCHI2_OWNPV", &J_psi_1S_IPCHI2_OWNPV, &b_J_psi_1S_IPCHI2_OWNPV);
   fChain->SetBranchAddress("J_psi_1S_FD_OWNPV", &J_psi_1S_FD_OWNPV, &b_J_psi_1S_FD_OWNPV);
   fChain->SetBranchAddress("J_psi_1S_FDCHI2_OWNPV", &J_psi_1S_FDCHI2_OWNPV, &b_J_psi_1S_FDCHI2_OWNPV);
   fChain->SetBranchAddress("J_psi_1S_DIRA_OWNPV", &J_psi_1S_DIRA_OWNPV, &b_J_psi_1S_DIRA_OWNPV);
   fChain->SetBranchAddress("J_psi_1S_P", &J_psi_1S_P, &b_J_psi_1S_P);
   fChain->SetBranchAddress("J_psi_1S_PT", &J_psi_1S_PT, &b_J_psi_1S_PT);
   fChain->SetBranchAddress("J_psi_1S_PE", &J_psi_1S_PE, &b_J_psi_1S_PE);
   fChain->SetBranchAddress("J_psi_1S_PX", &J_psi_1S_PX, &b_J_psi_1S_PX);
   fChain->SetBranchAddress("J_psi_1S_PY", &J_psi_1S_PY, &b_J_psi_1S_PY);
   fChain->SetBranchAddress("J_psi_1S_PZ", &J_psi_1S_PZ, &b_J_psi_1S_PZ);
   fChain->SetBranchAddress("J_psi_1S_MM", &J_psi_1S_MM, &b_J_psi_1S_MM);
   fChain->SetBranchAddress("J_psi_1S_MMERR", &J_psi_1S_MMERR, &b_J_psi_1S_MMERR);
   fChain->SetBranchAddress("J_psi_1S_M", &J_psi_1S_M, &b_J_psi_1S_M);
   fChain->SetBranchAddress("J_psi_1S_ID", &J_psi_1S_ID, &b_J_psi_1S_ID);
   fChain->SetBranchAddress("J_psi_1S_L0Global_Dec", &J_psi_1S_L0Global_Dec, &b_J_psi_1S_L0Global_Dec);
   fChain->SetBranchAddress("J_psi_1S_L0Global_TIS", &J_psi_1S_L0Global_TIS, &b_J_psi_1S_L0Global_TIS);
   fChain->SetBranchAddress("J_psi_1S_L0Global_TOS", &J_psi_1S_L0Global_TOS, &b_J_psi_1S_L0Global_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt1Global_Dec", &J_psi_1S_Hlt1Global_Dec, &b_J_psi_1S_Hlt1Global_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt1Global_TIS", &J_psi_1S_Hlt1Global_TIS, &b_J_psi_1S_Hlt1Global_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt1Global_TOS", &J_psi_1S_Hlt1Global_TOS, &b_J_psi_1S_Hlt1Global_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt1Phys_Dec", &J_psi_1S_Hlt1Phys_Dec, &b_J_psi_1S_Hlt1Phys_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt1Phys_TIS", &J_psi_1S_Hlt1Phys_TIS, &b_J_psi_1S_Hlt1Phys_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt1Phys_TOS", &J_psi_1S_Hlt1Phys_TOS, &b_J_psi_1S_Hlt1Phys_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2Global_Dec", &J_psi_1S_Hlt2Global_Dec, &b_J_psi_1S_Hlt2Global_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2Global_TIS", &J_psi_1S_Hlt2Global_TIS, &b_J_psi_1S_Hlt2Global_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2Global_TOS", &J_psi_1S_Hlt2Global_TOS, &b_J_psi_1S_Hlt2Global_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2Phys_Dec", &J_psi_1S_Hlt2Phys_Dec, &b_J_psi_1S_Hlt2Phys_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2Phys_TIS", &J_psi_1S_Hlt2Phys_TIS, &b_J_psi_1S_Hlt2Phys_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2Phys_TOS", &J_psi_1S_Hlt2Phys_TOS, &b_J_psi_1S_Hlt2Phys_TOS);
   fChain->SetBranchAddress("J_psi_1S_L0MuonDecision_Dec", &J_psi_1S_L0MuonDecision_Dec, &b_J_psi_1S_L0MuonDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_L0MuonDecision_TIS", &J_psi_1S_L0MuonDecision_TIS, &b_J_psi_1S_L0MuonDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_L0MuonDecision_TOS", &J_psi_1S_L0MuonDecision_TOS, &b_J_psi_1S_L0MuonDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_L0PhotonDecision_Dec", &J_psi_1S_L0PhotonDecision_Dec, &b_J_psi_1S_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_L0PhotonDecision_TIS", &J_psi_1S_L0PhotonDecision_TIS, &b_J_psi_1S_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_L0PhotonDecision_TOS", &J_psi_1S_L0PhotonDecision_TOS, &b_J_psi_1S_L0PhotonDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_L0ElectronDecision_Dec", &J_psi_1S_L0ElectronDecision_Dec, &b_J_psi_1S_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_L0ElectronDecision_TIS", &J_psi_1S_L0ElectronDecision_TIS, &b_J_psi_1S_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_L0ElectronDecision_TOS", &J_psi_1S_L0ElectronDecision_TOS, &b_J_psi_1S_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_L0HadronDecision_Dec", &J_psi_1S_L0HadronDecision_Dec, &b_J_psi_1S_L0HadronDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_L0HadronDecision_TIS", &J_psi_1S_L0HadronDecision_TIS, &b_J_psi_1S_L0HadronDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_L0HadronDecision_TOS", &J_psi_1S_L0HadronDecision_TOS, &b_J_psi_1S_L0HadronDecision_TOS);

   fChain->SetBranchAddress("J_psi_1S_L0DiMuonDecision_Dec", &J_psi_1S_L0DiMuonDecision_Dec, &b_J_psi_1S_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_L0DiMuonDecision_TIS", &J_psi_1S_L0DiMuonDecision_TIS, &b_J_psi_1S_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_L0DiMuonDecision_TOS", &J_psi_1S_L0DiMuonDecision_TOS, &b_J_psi_1S_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt1DiMuonHighMassDecision_Dec", &J_psi_1S_Hlt1DiMuonHighMassDecision_Dec, &b_J_psi_1S_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt1DiMuonHighMassDecision_TIS", &J_psi_1S_Hlt1DiMuonHighMassDecision_TIS, &b_J_psi_1S_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt1DiMuonHighMassDecision_TOS", &J_psi_1S_Hlt1DiMuonHighMassDecision_TOS, &b_J_psi_1S_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedDecision_Dec", &J_psi_1S_Hlt2DiMuonDetachedDecision_Dec, &b_J_psi_1S_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedDecision_TIS", &J_psi_1S_Hlt2DiMuonDetachedDecision_TIS, &b_J_psi_1S_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedDecision_TOS", &J_psi_1S_Hlt2DiMuonDetachedDecision_TOS, &b_J_psi_1S_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_Dec", &J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TIS", &J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TOS", &J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_J_psi_1S_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_Dec", &J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TIS", &J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TOS", &J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_J_psi_1S_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_Dec", &J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TIS", &J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TOS", &J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_J_psi_1S_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonJPsiDecision_Dec", &J_psi_1S_Hlt2DiMuonJPsiDecision_Dec, &b_J_psi_1S_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonJPsiDecision_TIS", &J_psi_1S_Hlt2DiMuonJPsiDecision_TIS, &b_J_psi_1S_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonJPsiDecision_TOS", &J_psi_1S_Hlt2DiMuonJPsiDecision_TOS, &b_J_psi_1S_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_Dec", &J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TIS", &J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS", &J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_J_psi_1S_Hlt2DiMuonJPsiHighPTDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_X", &J_psi_1S0_ENDVERTEX_X, &b_J_psi_1S0_ENDVERTEX_X);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_Y", &J_psi_1S0_ENDVERTEX_Y, &b_J_psi_1S0_ENDVERTEX_Y);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_Z", &J_psi_1S0_ENDVERTEX_Z, &b_J_psi_1S0_ENDVERTEX_Z);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_XERR", &J_psi_1S0_ENDVERTEX_XERR, &b_J_psi_1S0_ENDVERTEX_XERR);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_YERR", &J_psi_1S0_ENDVERTEX_YERR, &b_J_psi_1S0_ENDVERTEX_YERR);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_ZERR", &J_psi_1S0_ENDVERTEX_ZERR, &b_J_psi_1S0_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_CHI2", &J_psi_1S0_ENDVERTEX_CHI2, &b_J_psi_1S0_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_NDOF", &J_psi_1S0_ENDVERTEX_NDOF, &b_J_psi_1S0_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("J_psi_1S0_ENDVERTEX_COV_", J_psi_1S0_ENDVERTEX_COV_, &b_J_psi_1S0_ENDVERTEX_COV_);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_X", &J_psi_1S0_OWNPV_X, &b_J_psi_1S0_OWNPV_X);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_Y", &J_psi_1S0_OWNPV_Y, &b_J_psi_1S0_OWNPV_Y);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_Z", &J_psi_1S0_OWNPV_Z, &b_J_psi_1S0_OWNPV_Z);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_XERR", &J_psi_1S0_OWNPV_XERR, &b_J_psi_1S0_OWNPV_XERR);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_YERR", &J_psi_1S0_OWNPV_YERR, &b_J_psi_1S0_OWNPV_YERR);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_ZERR", &J_psi_1S0_OWNPV_ZERR, &b_J_psi_1S0_OWNPV_ZERR);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_CHI2", &J_psi_1S0_OWNPV_CHI2, &b_J_psi_1S0_OWNPV_CHI2);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_NDOF", &J_psi_1S0_OWNPV_NDOF, &b_J_psi_1S0_OWNPV_NDOF);
   fChain->SetBranchAddress("J_psi_1S0_OWNPV_COV_", J_psi_1S0_OWNPV_COV_, &b_J_psi_1S0_OWNPV_COV_);
   fChain->SetBranchAddress("J_psi_1S0_IP_OWNPV", &J_psi_1S0_IP_OWNPV, &b_J_psi_1S0_IP_OWNPV);
   fChain->SetBranchAddress("J_psi_1S0_IPCHI2_OWNPV", &J_psi_1S0_IPCHI2_OWNPV, &b_J_psi_1S0_IPCHI2_OWNPV);
   fChain->SetBranchAddress("J_psi_1S0_FD_OWNPV", &J_psi_1S0_FD_OWNPV, &b_J_psi_1S0_FD_OWNPV);
   fChain->SetBranchAddress("J_psi_1S0_FDCHI2_OWNPV", &J_psi_1S0_FDCHI2_OWNPV, &b_J_psi_1S0_FDCHI2_OWNPV);
   fChain->SetBranchAddress("J_psi_1S0_DIRA_OWNPV", &J_psi_1S0_DIRA_OWNPV, &b_J_psi_1S0_DIRA_OWNPV);
   fChain->SetBranchAddress("J_psi_1S0_P", &J_psi_1S0_P, &b_J_psi_1S0_P);
   fChain->SetBranchAddress("J_psi_1S0_PT", &J_psi_1S0_PT, &b_J_psi_1S0_PT);
   fChain->SetBranchAddress("J_psi_1S0_PE", &J_psi_1S0_PE, &b_J_psi_1S0_PE);
   fChain->SetBranchAddress("J_psi_1S0_PX", &J_psi_1S0_PX, &b_J_psi_1S0_PX);
   fChain->SetBranchAddress("J_psi_1S0_PY", &J_psi_1S0_PY, &b_J_psi_1S0_PY);
   fChain->SetBranchAddress("J_psi_1S0_PZ", &J_psi_1S0_PZ, &b_J_psi_1S0_PZ);
   fChain->SetBranchAddress("J_psi_1S0_MM", &J_psi_1S0_MM, &b_J_psi_1S0_MM);
   fChain->SetBranchAddress("J_psi_1S0_MMERR", &J_psi_1S0_MMERR, &b_J_psi_1S0_MMERR);
   fChain->SetBranchAddress("J_psi_1S0_M", &J_psi_1S0_M, &b_J_psi_1S0_M);
   fChain->SetBranchAddress("J_psi_1S0_ID", &J_psi_1S0_ID, &b_J_psi_1S0_ID);
   fChain->SetBranchAddress("J_psi_1S0_L0Global_Dec", &J_psi_1S0_L0Global_Dec, &b_J_psi_1S0_L0Global_Dec);
   fChain->SetBranchAddress("J_psi_1S0_L0Global_TIS", &J_psi_1S0_L0Global_TIS, &b_J_psi_1S0_L0Global_TIS);
   fChain->SetBranchAddress("J_psi_1S0_L0Global_TOS", &J_psi_1S0_L0Global_TOS, &b_J_psi_1S0_L0Global_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1Global_Dec", &J_psi_1S0_Hlt1Global_Dec, &b_J_psi_1S0_Hlt1Global_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1Global_TIS", &J_psi_1S0_Hlt1Global_TIS, &b_J_psi_1S0_Hlt1Global_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1Global_TOS", &J_psi_1S0_Hlt1Global_TOS, &b_J_psi_1S0_Hlt1Global_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1Phys_Dec", &J_psi_1S0_Hlt1Phys_Dec, &b_J_psi_1S0_Hlt1Phys_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1Phys_TIS", &J_psi_1S0_Hlt1Phys_TIS, &b_J_psi_1S0_Hlt1Phys_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1Phys_TOS", &J_psi_1S0_Hlt1Phys_TOS, &b_J_psi_1S0_Hlt1Phys_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2Global_Dec", &J_psi_1S0_Hlt2Global_Dec, &b_J_psi_1S0_Hlt2Global_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2Global_TIS", &J_psi_1S0_Hlt2Global_TIS, &b_J_psi_1S0_Hlt2Global_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2Global_TOS", &J_psi_1S0_Hlt2Global_TOS, &b_J_psi_1S0_Hlt2Global_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2Phys_Dec", &J_psi_1S0_Hlt2Phys_Dec, &b_J_psi_1S0_Hlt2Phys_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2Phys_TIS", &J_psi_1S0_Hlt2Phys_TIS, &b_J_psi_1S0_Hlt2Phys_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2Phys_TOS", &J_psi_1S0_Hlt2Phys_TOS, &b_J_psi_1S0_Hlt2Phys_TOS);
   fChain->SetBranchAddress("J_psi_1S0_L0MuonDecision_Dec", &J_psi_1S0_L0MuonDecision_Dec, &b_J_psi_1S0_L0MuonDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_L0MuonDecision_TIS", &J_psi_1S0_L0MuonDecision_TIS, &b_J_psi_1S0_L0MuonDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_L0MuonDecision_TOS", &J_psi_1S0_L0MuonDecision_TOS, &b_J_psi_1S0_L0MuonDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_L0DiMuonDecision_Dec", &J_psi_1S0_L0DiMuonDecision_Dec, &b_J_psi_1S0_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_L0DiMuonDecision_TIS", &J_psi_1S0_L0DiMuonDecision_TIS, &b_J_psi_1S0_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_L0DiMuonDecision_TOS", &J_psi_1S0_L0DiMuonDecision_TOS, &b_J_psi_1S0_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1DiMuonHighMassDecision_Dec", &J_psi_1S0_Hlt1DiMuonHighMassDecision_Dec, &b_J_psi_1S0_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1DiMuonHighMassDecision_TIS", &J_psi_1S0_Hlt1DiMuonHighMassDecision_TIS, &b_J_psi_1S0_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt1DiMuonHighMassDecision_TOS", &J_psi_1S0_Hlt1DiMuonHighMassDecision_TOS, &b_J_psi_1S0_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedDecision_Dec", &J_psi_1S0_Hlt2DiMuonDetachedDecision_Dec, &b_J_psi_1S0_Hlt2DiMuonDetachedDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedDecision_TIS", &J_psi_1S0_Hlt2DiMuonDetachedDecision_TIS, &b_J_psi_1S0_Hlt2DiMuonDetachedDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedDecision_TOS", &J_psi_1S0_Hlt2DiMuonDetachedDecision_TOS, &b_J_psi_1S0_Hlt2DiMuonDetachedDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_Dec", &J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TIS", &J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TOS", &J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_J_psi_1S0_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_Dec", &J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TIS", &J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TOS", &J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_J_psi_1S0_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_Dec", &J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_Dec, &b_J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TIS", &J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TIS, &b_J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TOS", &J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TOS, &b_J_psi_1S0_Hlt2DiMuonDetachedPsi2SDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonJPsiDecision_Dec", &J_psi_1S0_Hlt2DiMuonJPsiDecision_Dec, &b_J_psi_1S0_Hlt2DiMuonJPsiDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonJPsiDecision_TIS", &J_psi_1S0_Hlt2DiMuonJPsiDecision_TIS, &b_J_psi_1S0_Hlt2DiMuonJPsiDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonJPsiDecision_TOS", &J_psi_1S0_Hlt2DiMuonJPsiDecision_TOS, &b_J_psi_1S0_Hlt2DiMuonJPsiDecision_TOS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_Dec", &J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_Dec, &b_J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_Dec);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TIS", &J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TIS, &b_J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TIS);
   fChain->SetBranchAddress("J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TOS", &J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TOS, &b_J_psi_1S0_Hlt2DiMuonJPsiHighPTDecision_TOS);
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
   fChain->SetBranchAddress("StrippingBs2MuMuLinesBu2JPsiKFullLineDecision", &StrippingBs2MuMuLinesBu2JPsiKFullLineDecision, &b_StrippingBs2MuMuLinesBu2JPsiKFullLineDecision);
   Notify();
}

Bool_t JpsiTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void JpsiTree::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t JpsiTree::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef JpsiTree_cxx

