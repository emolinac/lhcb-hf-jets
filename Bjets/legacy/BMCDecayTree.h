//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug 31 19:39:21 2023 by ROOT version 6.24/08
// from TChain MCDecay/MCDecayTree/
//////////////////////////////////////////////////////////

#ifndef BMCDecayTree_h
#define BMCDecayTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class BMCDecayTree
{
public:
   TTree *fChain;  //! pointer to the analyzed TTree or TChain
   Int_t fCurrent; //! current Tree number in a TChain

   // Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t B_ETA;
   Double_t B_ID;
   Double_t B_P;
   Double_t B_PHI;
   Double_t B_ThreeCharge;
   Double_t B_TRUEP_E;
   Double_t B_TRUEP_X;
   Double_t B_TRUEP_Y;
   Double_t B_TRUEP_Z;
   Double_t B_TRUEPT;
   Double_t B_TRUEORIGINVERTEX_X;
   Double_t B_TRUEORIGINVERTEX_Y;
   Double_t B_TRUEORIGINVERTEX_Z;
   Double_t B_TRUEENDVERTEX_X;
   Double_t B_TRUEENDVERTEX_Y;
   Double_t B_TRUEENDVERTEX_Z;
   Bool_t B_TRUEISSTABLE;
   Double_t B_TRUETAU;
   Int_t B_Reconstructible;
   Int_t B_Reconstructed;
   Int_t B_ProtoParticles;
   Float_t B_PP_PX[20];       //[B_ProtoParticles]
   Float_t B_PP_PY[20];       //[B_ProtoParticles]
   Float_t B_PP_PZ[20];       //[B_ProtoParticles]
   Float_t B_PP_Weight[20];   //[B_ProtoParticles]
   Float_t B_PP_tr_pchi2[20]; //[B_ProtoParticles]
   Float_t B_PP_DLLe[20];     //[B_ProtoParticles]
   Float_t B_PP_DLLk[20];     //[B_ProtoParticles]
   Float_t B_PP_DLLp[20];     //[B_ProtoParticles]
   Float_t B_PP_DLLmu[20];    //[B_ProtoParticles]
   Double_t Bplus0_ETA;
   Double_t Bplus0_ID;
   Double_t Bplus0_P;
   Double_t Bplus0_PHI;
   Double_t Bplus0_ThreeCharge;
   Double_t Bplus0_TRUEP_E;
   Double_t Bplus0_TRUEP_X;
   Double_t Bplus0_TRUEP_Y;
   Double_t Bplus0_TRUEP_Z;
   Double_t Bplus0_TRUEPT;
   Double_t Bplus0_TRUEORIGINVERTEX_X;
   Double_t Bplus0_TRUEORIGINVERTEX_Y;
   Double_t Bplus0_TRUEORIGINVERTEX_Z;
   Double_t Bplus0_TRUEENDVERTEX_X;
   Double_t Bplus0_TRUEENDVERTEX_Y;
   Double_t Bplus0_TRUEENDVERTEX_Z;
   Bool_t Bplus0_TRUEISSTABLE;
   Double_t Bplus0_TRUETAU;
   Int_t Bplus0_Reconstructible;
   Int_t Bplus0_Reconstructed;
   Int_t Bplus0_ProtoParticles;
   Float_t Bplus0_PP_PX[20];       //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_PY[20];       //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_PZ[20];       //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_Weight[20];   //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_tr_pchi2[20]; //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_DLLe[20];     //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_DLLk[20];     //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_DLLp[20];     //[Bplus0_ProtoParticles]
   Float_t Bplus0_PP_DLLmu[20];    //[Bplus0_ProtoParticles]
   Double_t Jpsi_ETA;
   Double_t Jpsi_ID;
   Double_t Jpsi_P;
   Double_t Jpsi_PHI;
   Double_t Jpsi_ThreeCharge;
   Double_t Jpsi_TRUEP_E;
   Double_t Jpsi_TRUEP_X;
   Double_t Jpsi_TRUEP_Y;
   Double_t Jpsi_TRUEP_Z;
   Double_t Jpsi_TRUEPT;
   Double_t Jpsi_TRUEORIGINVERTEX_X;
   Double_t Jpsi_TRUEORIGINVERTEX_Y;
   Double_t Jpsi_TRUEORIGINVERTEX_Z;
   Double_t Jpsi_TRUEENDVERTEX_X;
   Double_t Jpsi_TRUEENDVERTEX_Y;
   Double_t Jpsi_TRUEENDVERTEX_Z;
   Bool_t Jpsi_TRUEISSTABLE;
   Double_t Jpsi_TRUETAU;
   Int_t Jpsi_Reconstructible;
   Int_t Jpsi_Reconstructed;
   Int_t Jpsi_ProtoParticles;
   Float_t Jpsi_PP_PX[20];       //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_PY[20];       //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_PZ[20];       //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_Weight[20];   //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_tr_pchi2[20]; //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_DLLe[20];     //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_DLLk[20];     //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_DLLp[20];     //[Jpsi_ProtoParticles]
   Float_t Jpsi_PP_DLLmu[20];    //[Jpsi_ProtoParticles]
   Double_t mup_ETA;
   Double_t mup_ID;
   Double_t mup_P;
   Double_t mup_PHI;
   Double_t mup_ThreeCharge;
   Double_t mup_TRUEP_E;
   Double_t mup_TRUEP_X;
   Double_t mup_TRUEP_Y;
   Double_t mup_TRUEP_Z;
   Double_t mup_TRUEPT;
   Double_t mup_TRUEORIGINVERTEX_X;
   Double_t mup_TRUEORIGINVERTEX_Y;
   Double_t mup_TRUEORIGINVERTEX_Z;
   Double_t mup_TRUEENDVERTEX_X;
   Double_t mup_TRUEENDVERTEX_Y;
   Double_t mup_TRUEENDVERTEX_Z;
   Bool_t mup_TRUEISSTABLE;
   Double_t mup_TRUETAU;
   Int_t mup_Reconstructible;
   Int_t mup_Reconstructed;
   Int_t mup_ProtoParticles;
   Float_t mup_PP_PX[20];       //[mup_ProtoParticles]
   Float_t mup_PP_PY[20];       //[mup_ProtoParticles]
   Float_t mup_PP_PZ[20];       //[mup_ProtoParticles]
   Float_t mup_PP_Weight[20];   //[mup_ProtoParticles]
   Float_t mup_PP_tr_pchi2[20]; //[mup_ProtoParticles]
   Float_t mup_PP_DLLe[20];     //[mup_ProtoParticles]
   Float_t mup_PP_DLLk[20];     //[mup_ProtoParticles]
   Float_t mup_PP_DLLp[20];     //[mup_ProtoParticles]
   Float_t mup_PP_DLLmu[20];    //[mup_ProtoParticles]
   Double_t mum_ETA;
   Double_t mum_ID;
   Double_t mum_P;
   Double_t mum_PHI;
   Double_t mum_ThreeCharge;
   Double_t mum_TRUEP_E;
   Double_t mum_TRUEP_X;
   Double_t mum_TRUEP_Y;
   Double_t mum_TRUEP_Z;
   Double_t mum_TRUEPT;
   Double_t mum_TRUEORIGINVERTEX_X;
   Double_t mum_TRUEORIGINVERTEX_Y;
   Double_t mum_TRUEORIGINVERTEX_Z;
   Double_t mum_TRUEENDVERTEX_X;
   Double_t mum_TRUEENDVERTEX_Y;
   Double_t mum_TRUEENDVERTEX_Z;
   Bool_t mum_TRUEISSTABLE;
   Double_t mum_TRUETAU;
   Int_t mum_Reconstructible;
   Int_t mum_Reconstructed;
   Int_t mum_ProtoParticles;
   Float_t mum_PP_PX[20];       //[mum_ProtoParticles]
   Float_t mum_PP_PY[20];       //[mum_ProtoParticles]
   Float_t mum_PP_PZ[20];       //[mum_ProtoParticles]
   Float_t mum_PP_Weight[20];   //[mum_ProtoParticles]
   Float_t mum_PP_tr_pchi2[20]; //[mum_ProtoParticles]
   Float_t mum_PP_DLLe[20];     //[mum_ProtoParticles]
   Float_t mum_PP_DLLk[20];     //[mum_ProtoParticles]
   Float_t mum_PP_DLLp[20];     //[mum_ProtoParticles]
   Float_t mum_PP_DLLmu[20];    //[mum_ProtoParticles]
   Double_t K_ETA;
   Double_t K_ID;
   Double_t K_P;
   Double_t K_PHI;
   Double_t K_ThreeCharge;
   Double_t K_TRUEP_E;
   Double_t K_TRUEP_X;
   Double_t K_TRUEP_Y;
   Double_t K_TRUEP_Z;
   Double_t K_TRUEPT;
   Double_t K_TRUEORIGINVERTEX_X;
   Double_t K_TRUEORIGINVERTEX_Y;
   Double_t K_TRUEORIGINVERTEX_Z;
   Double_t K_TRUEENDVERTEX_X;
   Double_t K_TRUEENDVERTEX_Y;
   Double_t K_TRUEENDVERTEX_Z;
   Bool_t K_TRUEISSTABLE;
   Double_t K_TRUETAU;
   Int_t K_Reconstructible;
   Int_t K_Reconstructed;
   Int_t K_ProtoParticles;
   Float_t K_PP_PX[20];       //[K_ProtoParticles]
   Float_t K_PP_PY[20];       //[K_ProtoParticles]
   Float_t K_PP_PZ[20];       //[K_ProtoParticles]
   Float_t K_PP_Weight[20];   //[K_ProtoParticles]
   Float_t K_PP_tr_pchi2[20]; //[K_ProtoParticles]
   Float_t K_PP_DLLe[20];     //[K_ProtoParticles]
   Float_t K_PP_DLLk[20];     //[K_ProtoParticles]
   Float_t K_PP_DLLp[20];     //[K_ProtoParticles]
   Float_t K_PP_DLLmu[20];    //[K_ProtoParticles]
   UInt_t nCandidate;
   ULong64_t totCandidates;
   ULong64_t EventInSequence;

   // List of branches
   TBranch *b_B_ETA;                     //!
   TBranch *b_B_ID;                      //!
   TBranch *b_B_P;                       //!
   TBranch *b_B_PHI;                     //!
   TBranch *b_B_ThreeCharge;             //!
   TBranch *b_B_TRUEP_E;                 //!
   TBranch *b_B_TRUEP_X;                 //!
   TBranch *b_B_TRUEP_Y;                 //!
   TBranch *b_B_TRUEP_Z;                 //!
   TBranch *b_B_TRUEPT;                  //!
   TBranch *b_B_TRUEORIGINVERTEX_X;      //!
   TBranch *b_B_TRUEORIGINVERTEX_Y;      //!
   TBranch *b_B_TRUEORIGINVERTEX_Z;      //!
   TBranch *b_B_TRUEENDVERTEX_X;         //!
   TBranch *b_B_TRUEENDVERTEX_Y;         //!
   TBranch *b_B_TRUEENDVERTEX_Z;         //!
   TBranch *b_B_TRUEISSTABLE;            //!
   TBranch *b_B_TRUETAU;                 //!
   TBranch *b_B_Reconstructible;         //!
   TBranch *b_B_Reconstructed;           //!
   TBranch *b_B_ProtoParticles;          //!
   TBranch *b_B_PP_PX;                   //!
   TBranch *b_B_PP_PY;                   //!
   TBranch *b_B_PP_PZ;                   //!
   TBranch *b_B_PP_Weight;               //!
   TBranch *b_B_PP_tr_pchi2;             //!
   TBranch *b_B_PP_DLLe;                 //!
   TBranch *b_B_PP_DLLk;                 //!
   TBranch *b_B_PP_DLLp;                 //!
   TBranch *b_B_PP_DLLmu;                //!
   TBranch *b_Bplus0_ETA;                //!
   TBranch *b_Bplus0_ID;                 //!
   TBranch *b_Bplus0_P;                  //!
   TBranch *b_Bplus0_PHI;                //!
   TBranch *b_Bplus0_ThreeCharge;        //!
   TBranch *b_Bplus0_TRUEP_E;            //!
   TBranch *b_Bplus0_TRUEP_X;            //!
   TBranch *b_Bplus0_TRUEP_Y;            //!
   TBranch *b_Bplus0_TRUEP_Z;            //!
   TBranch *b_Bplus0_TRUEPT;             //!
   TBranch *b_Bplus0_TRUEORIGINVERTEX_X; //!
   TBranch *b_Bplus0_TRUEORIGINVERTEX_Y; //!
   TBranch *b_Bplus0_TRUEORIGINVERTEX_Z; //!
   TBranch *b_Bplus0_TRUEENDVERTEX_X;    //!
   TBranch *b_Bplus0_TRUEENDVERTEX_Y;    //!
   TBranch *b_Bplus0_TRUEENDVERTEX_Z;    //!
   TBranch *b_Bplus0_TRUEISSTABLE;       //!
   TBranch *b_Bplus0_TRUETAU;            //!
   TBranch *b_Bplus0_Reconstructible;    //!
   TBranch *b_Bplus0_Reconstructed;      //!
   TBranch *b_Bplus0_ProtoParticles;     //!
   TBranch *b_Bplus0_PP_PX;              //!
   TBranch *b_Bplus0_PP_PY;              //!
   TBranch *b_Bplus0_PP_PZ;              //!
   TBranch *b_Bplus0_PP_Weight;          //!
   TBranch *b_Bplus0_PP_tr_pchi2;        //!
   TBranch *b_Bplus0_PP_DLLe;            //!
   TBranch *b_Bplus0_PP_DLLk;            //!
   TBranch *b_Bplus0_PP_DLLp;            //!
   TBranch *b_Bplus0_PP_DLLmu;           //!
   TBranch *b_Jpsi_ETA;                  //!
   TBranch *b_Jpsi_ID;                   //!
   TBranch *b_Jpsi_P;                    //!
   TBranch *b_Jpsi_PHI;                  //!
   TBranch *b_Jpsi_ThreeCharge;          //!
   TBranch *b_Jpsi_TRUEP_E;              //!
   TBranch *b_Jpsi_TRUEP_X;              //!
   TBranch *b_Jpsi_TRUEP_Y;              //!
   TBranch *b_Jpsi_TRUEP_Z;              //!
   TBranch *b_Jpsi_TRUEPT;               //!
   TBranch *b_Jpsi_TRUEORIGINVERTEX_X;   //!
   TBranch *b_Jpsi_TRUEORIGINVERTEX_Y;   //!
   TBranch *b_Jpsi_TRUEORIGINVERTEX_Z;   //!
   TBranch *b_Jpsi_TRUEENDVERTEX_X;      //!
   TBranch *b_Jpsi_TRUEENDVERTEX_Y;      //!
   TBranch *b_Jpsi_TRUEENDVERTEX_Z;      //!
   TBranch *b_Jpsi_TRUEISSTABLE;         //!
   TBranch *b_Jpsi_TRUETAU;              //!
   TBranch *b_Jpsi_Reconstructible;      //!
   TBranch *b_Jpsi_Reconstructed;        //!
   TBranch *b_Jpsi_ProtoParticles;       //!
   TBranch *b_Jpsi_PP_PX;                //!
   TBranch *b_Jpsi_PP_PY;                //!
   TBranch *b_Jpsi_PP_PZ;                //!
   TBranch *b_Jpsi_PP_Weight;            //!
   TBranch *b_Jpsi_PP_tr_pchi2;          //!
   TBranch *b_Jpsi_PP_DLLe;              //!
   TBranch *b_Jpsi_PP_DLLk;              //!
   TBranch *b_Jpsi_PP_DLLp;              //!
   TBranch *b_Jpsi_PP_DLLmu;             //!
   TBranch *b_mup_ETA;                   //!
   TBranch *b_mup_ID;                    //!
   TBranch *b_mup_P;                     //!
   TBranch *b_mup_PHI;                   //!
   TBranch *b_mup_ThreeCharge;           //!
   TBranch *b_mup_TRUEP_E;               //!
   TBranch *b_mup_TRUEP_X;               //!
   TBranch *b_mup_TRUEP_Y;               //!
   TBranch *b_mup_TRUEP_Z;               //!
   TBranch *b_mup_TRUEPT;                //!
   TBranch *b_mup_TRUEORIGINVERTEX_X;    //!
   TBranch *b_mup_TRUEORIGINVERTEX_Y;    //!
   TBranch *b_mup_TRUEORIGINVERTEX_Z;    //!
   TBranch *b_mup_TRUEENDVERTEX_X;       //!
   TBranch *b_mup_TRUEENDVERTEX_Y;       //!
   TBranch *b_mup_TRUEENDVERTEX_Z;       //!
   TBranch *b_mup_TRUEISSTABLE;          //!
   TBranch *b_mup_TRUETAU;               //!
   TBranch *b_mup_Reconstructible;       //!
   TBranch *b_mup_Reconstructed;         //!
   TBranch *b_mup_ProtoParticles;        //!
   TBranch *b_mup_PP_PX;                 //!
   TBranch *b_mup_PP_PY;                 //!
   TBranch *b_mup_PP_PZ;                 //!
   TBranch *b_mup_PP_Weight;             //!
   TBranch *b_mup_PP_tr_pchi2;           //!
   TBranch *b_mup_PP_DLLe;               //!
   TBranch *b_mup_PP_DLLk;               //!
   TBranch *b_mup_PP_DLLp;               //!
   TBranch *b_mup_PP_DLLmu;              //!
   TBranch *b_mum_ETA;                   //!
   TBranch *b_mum_ID;                    //!
   TBranch *b_mum_P;                     //!
   TBranch *b_mum_PHI;                   //!
   TBranch *b_mum_ThreeCharge;           //!
   TBranch *b_mum_TRUEP_E;               //!
   TBranch *b_mum_TRUEP_X;               //!
   TBranch *b_mum_TRUEP_Y;               //!
   TBranch *b_mum_TRUEP_Z;               //!
   TBranch *b_mum_TRUEPT;                //!
   TBranch *b_mum_TRUEORIGINVERTEX_X;    //!
   TBranch *b_mum_TRUEORIGINVERTEX_Y;    //!
   TBranch *b_mum_TRUEORIGINVERTEX_Z;    //!
   TBranch *b_mum_TRUEENDVERTEX_X;       //!
   TBranch *b_mum_TRUEENDVERTEX_Y;       //!
   TBranch *b_mum_TRUEENDVERTEX_Z;       //!
   TBranch *b_mum_TRUEISSTABLE;          //!
   TBranch *b_mum_TRUETAU;               //!
   TBranch *b_mum_Reconstructible;       //!
   TBranch *b_mum_Reconstructed;         //!
   TBranch *b_mum_ProtoParticles;        //!
   TBranch *b_mum_PP_PX;                 //!
   TBranch *b_mum_PP_PY;                 //!
   TBranch *b_mum_PP_PZ;                 //!
   TBranch *b_mum_PP_Weight;             //!
   TBranch *b_mum_PP_tr_pchi2;           //!
   TBranch *b_mum_PP_DLLe;               //!
   TBranch *b_mum_PP_DLLk;               //!
   TBranch *b_mum_PP_DLLp;               //!
   TBranch *b_mum_PP_DLLmu;              //!
   TBranch *b_K_ETA;                     //!
   TBranch *b_K_ID;                      //!
   TBranch *b_K_P;                       //!
   TBranch *b_K_PHI;                     //!
   TBranch *b_K_ThreeCharge;             //!
   TBranch *b_K_TRUEP_E;                 //!
   TBranch *b_K_TRUEP_X;                 //!
   TBranch *b_K_TRUEP_Y;                 //!
   TBranch *b_K_TRUEP_Z;                 //!
   TBranch *b_K_TRUEPT;                  //!
   TBranch *b_K_TRUEORIGINVERTEX_X;      //!
   TBranch *b_K_TRUEORIGINVERTEX_Y;      //!
   TBranch *b_K_TRUEORIGINVERTEX_Z;      //!
   TBranch *b_K_TRUEENDVERTEX_X;         //!
   TBranch *b_K_TRUEENDVERTEX_Y;         //!
   TBranch *b_K_TRUEENDVERTEX_Z;         //!
   TBranch *b_K_TRUEISSTABLE;            //!
   TBranch *b_K_TRUETAU;                 //!
   TBranch *b_K_Reconstructible;         //!
   TBranch *b_K_Reconstructed;           //!
   TBranch *b_K_ProtoParticles;          //!
   TBranch *b_K_PP_PX;                   //!
   TBranch *b_K_PP_PY;                   //!
   TBranch *b_K_PP_PZ;                   //!
   TBranch *b_K_PP_Weight;               //!
   TBranch *b_K_PP_tr_pchi2;             //!
   TBranch *b_K_PP_DLLe;                 //!
   TBranch *b_K_PP_DLLk;                 //!
   TBranch *b_K_PP_DLLp;                 //!
   TBranch *b_K_PP_DLLmu;                //!
   TBranch *b_nCandidate;                //!
   TBranch *b_totCandidates;             //!
   TBranch *b_EventInSequence;           //!

   BMCDecayTree(TTree *tree = 0, int dataset = 93539, bool isData = true);
   virtual ~BMCDecayTree();
   virtual Int_t Cut(Long64_t entry);
   virtual Int_t GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void Init(TTree *tree);
   virtual void Loop();
   virtual Bool_t Notify();
   virtual void Show(Long64_t entry = -1);
};

#endif

#ifdef BMCDecayTree_cxx
BMCDecayTree::BMCDecayTree(TTree *tree = 0, int dataset = 93539, bool isData = true) : fChain(0)
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
      f->GetObject("MCDecay/MCDecayTree", tree);

#else  // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain *chain = new TChain("MCDecay/MCDecayTree", "");
      TString data_dir = "/Users/josearias18/Desktop/QCDc2/pythia/Data/";
      
      if (dataset == 93599)
      {
         cout << "Got all data" << endl;
         if (isData)
         {
            chain->Add(data_dir + "Bjet_b_2016_MD_09082022_url333_0_800.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2016_MD_09082022_url333_800_1429.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2016_MU_09082022_url334_0_800.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2016_MU_09082022_url334_800_1299.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2017_MD_09082022_url335_0_800.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2017_MD_09082022_url335_800_882.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2017_MU_09082022_url336_0_800.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2017_MU_09082022_url336_800_846.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2018_MD_09082022_url337_0_800.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2018_MD_09082022_url337_800_988.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2018_MU_09082022_url338_0_800.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_b_2018_MU_09082022_url338_800_1042.root/MCDecay/MCDecayTree");
         }
         else
         {
            // chain->Add("/afs/cern.ch/work/i/ichahrou/DeadCone/DaVinci/Run2_BuJpsiK/BuPsiK_in_jet_13TeV.root/MCDecay/MCDecayTree");
            // chain->Add(data_dir + "Bjet_MC_Sim09_2016_MD_09292022_full.root/MCDecay/MCDecayTree");
            // chain->Add(data_dir + "Bjet_MC_Sim09_2016_MU_09292022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim09b_MD_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim09b_MU_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim09g_MU_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim10a_MD_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim09g_MD_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim09k_MD_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim10a_MU_10242022_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "Bjet_MC_Sim09k_MU_10242022_full.root/MCDecay/MCDecayTree");
         }
      }
      else if (dataset == 91599)
      {
         if (isData)
         {
         }
         else
         {
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09l_MD_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09l_MU_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim10a_MD_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim10a_MU_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09l_MD_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09l_MU_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim10a_MD_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim10a_MU_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09l_MD_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09l_MU_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim10a_MD_09012023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim10a_MU_09012023_full.root/MCDecay/MCDecayTree");
         }
      }
      else if (dataset == 91590)
      {
         if (isData)
         {
         }
         else
         {
            // data_dir = "/afs/cern.ch/work/i/ichahrou/DeadCone/DaVinci/Run2_BuJpsiK/";
            // chain->Add(data_dir + "truth_BuPsiK_in_jet_13TeV_nojetid.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            // chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");

            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            // chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");

            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            // chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");
         }
      }
      else if (dataset == 91591)
      {
         if (isData)
         {
         }
         else
         {
            // data_dir = "/afs/cern.ch/work/i/ichahrou/DeadCone/DaVinci/Run2_BuJpsiK/";
            // chain->Add(data_dir + "truth_BuPsiK_in_jet_13TeV_nojetid.root/MCDecay/MCDecayTree");
            // chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2016_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");

            // chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2017_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");

            // chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09k_MD_08312023_full.root/MCDecay/MCDecayTree");
            chain->Add(data_dir + "BDecay_MC_DFT_2018_Sim09k_MU_08312023_full.root/MCDecay/MCDecayTree");
         }
      }
      else if (dataset == 99599)
      {
         if (isData)
         {
         }
         else
         {
            data_dir = "/afs/cern.ch/work/i/ichahrou/DeadCone/DaVinci/Run2_BuJpsiK/";

            chain->Add(data_dir + "BuPsiK_in_jet_13TeV_nojetid_dtf.root/MCDecay/MCDecayTree");

            // chain->Add(data_dir + "Bjet_MC_nojetid_DFT_2016_Sim09k_MD_11082023_full.root/MCDecay/MCDecayTree");
         }
      }

      tree = chain;
#endif // SINGLE_TREE
   }
   Init(tree);
}

BMCDecayTree::~BMCDecayTree()
{
   if (!fChain)
      return;
   delete fChain->GetCurrentFile();
}

Int_t BMCDecayTree::GetEntry(Long64_t entry)
{
   // Read contents of entry.
   if (!fChain)
      return 0;
   return fChain->GetEntry(entry);
}
Long64_t BMCDecayTree::LoadTree(Long64_t entry)
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

void BMCDecayTree::Init(TTree *tree)
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

   fChain->SetBranchAddress("B_ETA", &B_ETA, &b_B_ETA);
   fChain->SetBranchAddress("B_ID", &B_ID, &b_B_ID);
   fChain->SetBranchAddress("B_P", &B_P, &b_B_P);
   fChain->SetBranchAddress("B_PHI", &B_PHI, &b_B_PHI);
   fChain->SetBranchAddress("B_ThreeCharge", &B_ThreeCharge, &b_B_ThreeCharge);
   fChain->SetBranchAddress("B_TRUEP_E", &B_TRUEP_E, &b_B_TRUEP_E);
   fChain->SetBranchAddress("B_TRUEP_X", &B_TRUEP_X, &b_B_TRUEP_X);
   fChain->SetBranchAddress("B_TRUEP_Y", &B_TRUEP_Y, &b_B_TRUEP_Y);
   fChain->SetBranchAddress("B_TRUEP_Z", &B_TRUEP_Z, &b_B_TRUEP_Z);
   fChain->SetBranchAddress("B_TRUEPT", &B_TRUEPT, &b_B_TRUEPT);
   fChain->SetBranchAddress("B_TRUEORIGINVERTEX_X", &B_TRUEORIGINVERTEX_X, &b_B_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("B_TRUEORIGINVERTEX_Y", &B_TRUEORIGINVERTEX_Y, &b_B_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("B_TRUEORIGINVERTEX_Z", &B_TRUEORIGINVERTEX_Z, &b_B_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("B_TRUEENDVERTEX_X", &B_TRUEENDVERTEX_X, &b_B_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("B_TRUEENDVERTEX_Y", &B_TRUEENDVERTEX_Y, &b_B_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("B_TRUEENDVERTEX_Z", &B_TRUEENDVERTEX_Z, &b_B_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("B_TRUEISSTABLE", &B_TRUEISSTABLE, &b_B_TRUEISSTABLE);
   fChain->SetBranchAddress("B_TRUETAU", &B_TRUETAU, &b_B_TRUETAU);
   fChain->SetBranchAddress("B_Reconstructible", &B_Reconstructible, &b_B_Reconstructible);
   fChain->SetBranchAddress("B_Reconstructed", &B_Reconstructed, &b_B_Reconstructed);
   fChain->SetBranchAddress("B_ProtoParticles", &B_ProtoParticles, &b_B_ProtoParticles);
   fChain->SetBranchAddress("B_PP_PX", B_PP_PX, &b_B_PP_PX);
   fChain->SetBranchAddress("B_PP_PY", B_PP_PY, &b_B_PP_PY);
   fChain->SetBranchAddress("B_PP_PZ", B_PP_PZ, &b_B_PP_PZ);
   fChain->SetBranchAddress("B_PP_Weight", B_PP_Weight, &b_B_PP_Weight);
   fChain->SetBranchAddress("B_PP_tr_pchi2", B_PP_tr_pchi2, &b_B_PP_tr_pchi2);
   fChain->SetBranchAddress("B_PP_DLLe", B_PP_DLLe, &b_B_PP_DLLe);
   fChain->SetBranchAddress("B_PP_DLLk", B_PP_DLLk, &b_B_PP_DLLk);
   fChain->SetBranchAddress("B_PP_DLLp", B_PP_DLLp, &b_B_PP_DLLp);
   fChain->SetBranchAddress("B_PP_DLLmu", B_PP_DLLmu, &b_B_PP_DLLmu);
   fChain->SetBranchAddress("Bplus0_ETA", &Bplus0_ETA, &b_Bplus0_ETA);
   fChain->SetBranchAddress("Bplus0_ID", &Bplus0_ID, &b_Bplus0_ID);
   fChain->SetBranchAddress("Bplus0_P", &Bplus0_P, &b_Bplus0_P);
   fChain->SetBranchAddress("Bplus0_PHI", &Bplus0_PHI, &b_Bplus0_PHI);
   fChain->SetBranchAddress("Bplus0_ThreeCharge", &Bplus0_ThreeCharge, &b_Bplus0_ThreeCharge);
   fChain->SetBranchAddress("Bplus0_TRUEP_E", &Bplus0_TRUEP_E, &b_Bplus0_TRUEP_E);
   fChain->SetBranchAddress("Bplus0_TRUEP_X", &Bplus0_TRUEP_X, &b_Bplus0_TRUEP_X);
   fChain->SetBranchAddress("Bplus0_TRUEP_Y", &Bplus0_TRUEP_Y, &b_Bplus0_TRUEP_Y);
   fChain->SetBranchAddress("Bplus0_TRUEP_Z", &Bplus0_TRUEP_Z, &b_Bplus0_TRUEP_Z);
   fChain->SetBranchAddress("Bplus0_TRUEPT", &Bplus0_TRUEPT, &b_Bplus0_TRUEPT);
   fChain->SetBranchAddress("Bplus0_TRUEORIGINVERTEX_X", &Bplus0_TRUEORIGINVERTEX_X, &b_Bplus0_TRUEORIGINVERTEX_X);
   fChain->SetBranchAddress("Bplus0_TRUEORIGINVERTEX_Y", &Bplus0_TRUEORIGINVERTEX_Y, &b_Bplus0_TRUEORIGINVERTEX_Y);
   fChain->SetBranchAddress("Bplus0_TRUEORIGINVERTEX_Z", &Bplus0_TRUEORIGINVERTEX_Z, &b_Bplus0_TRUEORIGINVERTEX_Z);
   fChain->SetBranchAddress("Bplus0_TRUEENDVERTEX_X", &Bplus0_TRUEENDVERTEX_X, &b_Bplus0_TRUEENDVERTEX_X);
   fChain->SetBranchAddress("Bplus0_TRUEENDVERTEX_Y", &Bplus0_TRUEENDVERTEX_Y, &b_Bplus0_TRUEENDVERTEX_Y);
   fChain->SetBranchAddress("Bplus0_TRUEENDVERTEX_Z", &Bplus0_TRUEENDVERTEX_Z, &b_Bplus0_TRUEENDVERTEX_Z);
   fChain->SetBranchAddress("Bplus0_TRUEISSTABLE", &Bplus0_TRUEISSTABLE, &b_Bplus0_TRUEISSTABLE);
   fChain->SetBranchAddress("Bplus0_TRUETAU", &Bplus0_TRUETAU, &b_Bplus0_TRUETAU);
   fChain->SetBranchAddress("Bplus0_Reconstructible", &Bplus0_Reconstructible, &b_Bplus0_Reconstructible);
   fChain->SetBranchAddress("Bplus0_Reconstructed", &Bplus0_Reconstructed, &b_Bplus0_Reconstructed);
   fChain->SetBranchAddress("Bplus0_ProtoParticles", &Bplus0_ProtoParticles, &b_Bplus0_ProtoParticles);
   fChain->SetBranchAddress("Bplus0_PP_PX", Bplus0_PP_PX, &b_Bplus0_PP_PX);
   fChain->SetBranchAddress("Bplus0_PP_PY", Bplus0_PP_PY, &b_Bplus0_PP_PY);
   fChain->SetBranchAddress("Bplus0_PP_PZ", Bplus0_PP_PZ, &b_Bplus0_PP_PZ);
   fChain->SetBranchAddress("Bplus0_PP_Weight", Bplus0_PP_Weight, &b_Bplus0_PP_Weight);
   fChain->SetBranchAddress("Bplus0_PP_tr_pchi2", Bplus0_PP_tr_pchi2, &b_Bplus0_PP_tr_pchi2);
   fChain->SetBranchAddress("Bplus0_PP_DLLe", Bplus0_PP_DLLe, &b_Bplus0_PP_DLLe);
   fChain->SetBranchAddress("Bplus0_PP_DLLk", Bplus0_PP_DLLk, &b_Bplus0_PP_DLLk);
   fChain->SetBranchAddress("Bplus0_PP_DLLp", Bplus0_PP_DLLp, &b_Bplus0_PP_DLLp);
   fChain->SetBranchAddress("Bplus0_PP_DLLmu", Bplus0_PP_DLLmu, &b_Bplus0_PP_DLLmu);
   fChain->SetBranchAddress("Jpsi_ETA", &Jpsi_ETA, &b_Jpsi_ETA);
   fChain->SetBranchAddress("Jpsi_ID", &Jpsi_ID, &b_Jpsi_ID);
   fChain->SetBranchAddress("Jpsi_P", &Jpsi_P, &b_Jpsi_P);
   fChain->SetBranchAddress("Jpsi_PHI", &Jpsi_PHI, &b_Jpsi_PHI);
   fChain->SetBranchAddress("Jpsi_ThreeCharge", &Jpsi_ThreeCharge, &b_Jpsi_ThreeCharge);
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
   fChain->SetBranchAddress("Jpsi_Reconstructible", &Jpsi_Reconstructible, &b_Jpsi_Reconstructible);
   fChain->SetBranchAddress("Jpsi_Reconstructed", &Jpsi_Reconstructed, &b_Jpsi_Reconstructed);
   fChain->SetBranchAddress("Jpsi_ProtoParticles", &Jpsi_ProtoParticles, &b_Jpsi_ProtoParticles);
   fChain->SetBranchAddress("Jpsi_PP_PX", Jpsi_PP_PX, &b_Jpsi_PP_PX);
   fChain->SetBranchAddress("Jpsi_PP_PY", Jpsi_PP_PY, &b_Jpsi_PP_PY);
   fChain->SetBranchAddress("Jpsi_PP_PZ", Jpsi_PP_PZ, &b_Jpsi_PP_PZ);
   fChain->SetBranchAddress("Jpsi_PP_Weight", Jpsi_PP_Weight, &b_Jpsi_PP_Weight);
   fChain->SetBranchAddress("Jpsi_PP_tr_pchi2", Jpsi_PP_tr_pchi2, &b_Jpsi_PP_tr_pchi2);
   fChain->SetBranchAddress("Jpsi_PP_DLLe", Jpsi_PP_DLLe, &b_Jpsi_PP_DLLe);
   fChain->SetBranchAddress("Jpsi_PP_DLLk", Jpsi_PP_DLLk, &b_Jpsi_PP_DLLk);
   fChain->SetBranchAddress("Jpsi_PP_DLLp", Jpsi_PP_DLLp, &b_Jpsi_PP_DLLp);
   fChain->SetBranchAddress("Jpsi_PP_DLLmu", Jpsi_PP_DLLmu, &b_Jpsi_PP_DLLmu);
   fChain->SetBranchAddress("mup_ETA", &mup_ETA, &b_mup_ETA);
   fChain->SetBranchAddress("mup_ID", &mup_ID, &b_mup_ID);
   fChain->SetBranchAddress("mup_P", &mup_P, &b_mup_P);
   fChain->SetBranchAddress("mup_PHI", &mup_PHI, &b_mup_PHI);
   fChain->SetBranchAddress("mup_ThreeCharge", &mup_ThreeCharge, &b_mup_ThreeCharge);
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
   fChain->SetBranchAddress("mup_Reconstructible", &mup_Reconstructible, &b_mup_Reconstructible);
   fChain->SetBranchAddress("mup_Reconstructed", &mup_Reconstructed, &b_mup_Reconstructed);
   fChain->SetBranchAddress("mup_ProtoParticles", &mup_ProtoParticles, &b_mup_ProtoParticles);
   fChain->SetBranchAddress("mup_PP_PX", mup_PP_PX, &b_mup_PP_PX);
   fChain->SetBranchAddress("mup_PP_PY", mup_PP_PY, &b_mup_PP_PY);
   fChain->SetBranchAddress("mup_PP_PZ", mup_PP_PZ, &b_mup_PP_PZ);
   fChain->SetBranchAddress("mup_PP_Weight", mup_PP_Weight, &b_mup_PP_Weight);
   fChain->SetBranchAddress("mup_PP_tr_pchi2", mup_PP_tr_pchi2, &b_mup_PP_tr_pchi2);
   fChain->SetBranchAddress("mup_PP_DLLe", mup_PP_DLLe, &b_mup_PP_DLLe);
   fChain->SetBranchAddress("mup_PP_DLLk", mup_PP_DLLk, &b_mup_PP_DLLk);
   fChain->SetBranchAddress("mup_PP_DLLp", mup_PP_DLLp, &b_mup_PP_DLLp);
   fChain->SetBranchAddress("mup_PP_DLLmu", mup_PP_DLLmu, &b_mup_PP_DLLmu);
   fChain->SetBranchAddress("mum_ETA", &mum_ETA, &b_mum_ETA);
   fChain->SetBranchAddress("mum_ID", &mum_ID, &b_mum_ID);
   fChain->SetBranchAddress("mum_P", &mum_P, &b_mum_P);
   fChain->SetBranchAddress("mum_PHI", &mum_PHI, &b_mum_PHI);
   fChain->SetBranchAddress("mum_ThreeCharge", &mum_ThreeCharge, &b_mum_ThreeCharge);
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
   fChain->SetBranchAddress("mum_Reconstructible", &mum_Reconstructible, &b_mum_Reconstructible);
   fChain->SetBranchAddress("mum_Reconstructed", &mum_Reconstructed, &b_mum_Reconstructed);
   fChain->SetBranchAddress("mum_ProtoParticles", &mum_ProtoParticles, &b_mum_ProtoParticles);
   fChain->SetBranchAddress("mum_PP_PX", mum_PP_PX, &b_mum_PP_PX);
   fChain->SetBranchAddress("mum_PP_PY", mum_PP_PY, &b_mum_PP_PY);
   fChain->SetBranchAddress("mum_PP_PZ", mum_PP_PZ, &b_mum_PP_PZ);
   fChain->SetBranchAddress("mum_PP_Weight", mum_PP_Weight, &b_mum_PP_Weight);
   fChain->SetBranchAddress("mum_PP_tr_pchi2", mum_PP_tr_pchi2, &b_mum_PP_tr_pchi2);
   fChain->SetBranchAddress("mum_PP_DLLe", mum_PP_DLLe, &b_mum_PP_DLLe);
   fChain->SetBranchAddress("mum_PP_DLLk", mum_PP_DLLk, &b_mum_PP_DLLk);
   fChain->SetBranchAddress("mum_PP_DLLp", mum_PP_DLLp, &b_mum_PP_DLLp);
   fChain->SetBranchAddress("mum_PP_DLLmu", mum_PP_DLLmu, &b_mum_PP_DLLmu);
   fChain->SetBranchAddress("K_ETA", &K_ETA, &b_K_ETA);
   fChain->SetBranchAddress("K_ID", &K_ID, &b_K_ID);
   fChain->SetBranchAddress("K_P", &K_P, &b_K_P);
   fChain->SetBranchAddress("K_PHI", &K_PHI, &b_K_PHI);
   fChain->SetBranchAddress("K_ThreeCharge", &K_ThreeCharge, &b_K_ThreeCharge);
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
   fChain->SetBranchAddress("K_Reconstructible", &K_Reconstructible, &b_K_Reconstructible);
   fChain->SetBranchAddress("K_Reconstructed", &K_Reconstructed, &b_K_Reconstructed);
   fChain->SetBranchAddress("K_ProtoParticles", &K_ProtoParticles, &b_K_ProtoParticles);
   fChain->SetBranchAddress("K_PP_PX", K_PP_PX, &b_K_PP_PX);
   fChain->SetBranchAddress("K_PP_PY", K_PP_PY, &b_K_PP_PY);
   fChain->SetBranchAddress("K_PP_PZ", K_PP_PZ, &b_K_PP_PZ);
   fChain->SetBranchAddress("K_PP_Weight", K_PP_Weight, &b_K_PP_Weight);
   fChain->SetBranchAddress("K_PP_tr_pchi2", K_PP_tr_pchi2, &b_K_PP_tr_pchi2);
   fChain->SetBranchAddress("K_PP_DLLe", K_PP_DLLe, &b_K_PP_DLLe);
   fChain->SetBranchAddress("K_PP_DLLk", K_PP_DLLk, &b_K_PP_DLLk);
   fChain->SetBranchAddress("K_PP_DLLp", K_PP_DLLp, &b_K_PP_DLLp);
   fChain->SetBranchAddress("K_PP_DLLmu", K_PP_DLLmu, &b_K_PP_DLLmu);
   fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
   fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
   fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
   Notify();
}

Bool_t BMCDecayTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BMCDecayTree::Show(Long64_t entry)
{
   // Print contents of entry.
   // If entry is not specified, print current entry
   if (!fChain)
      return;
   fChain->Show(entry);
}
Int_t BMCDecayTree::Cut(Long64_t entry)
{
   // This function may be called from Loop.
   // returns  1 if entry is accepted.
   // returns -1 otherwise.
   return 1;
}
#endif // #ifdef BMCDecayTree_cxx
