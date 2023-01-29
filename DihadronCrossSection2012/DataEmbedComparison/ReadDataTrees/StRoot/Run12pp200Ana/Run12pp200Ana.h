

/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 22 17:04:53 2010 by ROOT version 5.22/00
// from TTree ftree/LRC
// found on file: 00E5DE07942C320D8D35F5F8AD40C261_7total.root
//////////////////////////////////////////////////////////

#ifndef Run12pp200Ana_h
#define Run12pp200Ana_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class Run12pp200Ana
{
public:
        TChain *fChain; //! pointer to the analyzed TTree or TChain
        Int_t fCurrent; //! current Tree number in a TChain

        // Declaration of leaf types
        Int_t frefmult;
        Int_t fmaxpar;
        Int_t ffillNum;
        Int_t frunNum;
        Int_t fspinconfig;
        vector<unsigned int> *ftrigger;
        Double_t fVZ;
        float fvpdVz; // Included by Babu 11/22/2019
        unsigned int fevTime;
        Double_t fverRank;              // included by Babu  11/22/2019
        Double_t fpT[2422];             //[fmaxpar]
        Double_t fp[2422];              //[fmaxpar]
        Double_t fpX[2422];             //[fmaxpar]
        Double_t fpY[2422];             //[fmaxpar]
        Double_t fpZ[2422];             //[fmaxpar]
        Double_t feta[2422];            //[fmaxpar]
        Double_t fphi[2422];            //[fmaxpar]
        Short_t fcharge[2422];          //[fmaxpar]
        Double_t fnSigmaPion[2422];     //[fmaxpar]
        Double_t fnSigmaKaon[2422];     //[fmaxpar]
        Double_t fnSigmaProton[2422];   //[fmaxpar]
        Double_t fnSigmaElectron[2422]; //[fmaxpar]
        Double_t fdEdx[2422];           //[fmaxpar]
        Double_t fdca[2422];            //[fmaxpar]
        UShort_t ffitPts[2422];         //[fmaxpar]
        UShort_t ffitPtsPoss[2422];     //[fmaxpar]
        UShort_t fhitsdedx[2422];       //[fmaxpar]
        Double_t fBetaToF[2422];        //[fmaxpar]

        // List of branches
        TBranch *b_frefmult;    //!
        TBranch *b_fmaxpar;     //!
        TBranch *b_ffillNum;    //!
        TBranch *b_frunNum;     //!
        TBranch *b_fspinconfig; //!
        TBranch *b_ftrigger;    //!
        TBranch *b_fVZ;
        TBranch *b_fvpdVz;
        TBranch *b_fevTime;
        TBranch *b_fverRank;        //!
        TBranch *b_fpT;             //!
        TBranch *b_fp;              //!
        TBranch *b_fpX;             //!
        TBranch *b_fpY;             //!
        TBranch *b_fpZ;             //!
        TBranch *b_feta;            //!
        TBranch *b_fphi;            //!
        TBranch *b_fcharge;         //!
        TBranch *b_fnSigmaPion;     //!
        TBranch *b_fnSigmaKaon;     //!
        TBranch *b_fnSigmaProton;   //!
        TBranch *b_fnSigmaElectron; //!
        TBranch *b_fdEdx;           //!
        TBranch *b_fdca;            //!
        TBranch *b_ffitPts;         //!
        TBranch *b_ffitPtsPoss;     //!
        TBranch *b_fhitsdedx;       //!
        TBranch *b_fBetaToF;        //!

        Run12pp200Ana(char *ifile);
        virtual ~Run12pp200Ana();
        virtual Int_t Cut(Long64_t entry);
        virtual Int_t GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void Init();
        virtual void Loop();
        virtual Bool_t Notify();
        virtual void Finish(char *ofile);
        virtual void Show(Long64_t entry = -1);
        virtual Bool_t dcaCut(Double_t, Double_t);
        // virtual void     loadPol();
        inline void setHome(TDirectory *pHome)
        {
                cout << "setting home.." << endl;
                home = pHome;
        }

        static const int n = 23;
        int trigger;

        TH1D *hVZ[3];

        TH1D *hxsecMJP0;
        TH1D *hxsecMJP1;
        TH1D *hxsecMJP2;

        TH1D *htrackpT[3];
        TH1D *htrackEta[3]; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *htrackPhi[3];
        TH1D *htracknSigmaPion[3];

        TH1D *hpionpT[3];
        TH1D *hpionEta[3]; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *hpionPhi[3];

        // pair quantities
        TH1D *hCone[3];
        TH1D *hMinv[3];
        TH1D *hpTPair[3];
        TH1D *hEtaPair[3];
        TH1D *hPhiRY[3];
        TH1D *hPhiRB[3];
        // for combined trigger
        TH1D *hVZJP;

        TH1D *htrackpTJP;
        TH1D *htrackEtaJP; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *htrackPhiJP;
        TH1D *htracknSigmaPionJP;

        TH1D *hpionpTJP;
        TH1D *hpionEtaJP; // 3 = JP0, JP1, JP2, 2 = +, -
        TH1D *hpionPhiJP;

        // pair quantities
        TH1D *hConeJP;
        TH1D *hMinvJP;
        TH1D *hpTPairJP;
        TH1D *hEtaPairJP;
        TH1D *hPhiRYJP;
        TH1D *hPhiRBJP;

protected:
        TDirectory *home;

private:
        ClassDef(Run12pp200Ana, 1);
};

#endif
