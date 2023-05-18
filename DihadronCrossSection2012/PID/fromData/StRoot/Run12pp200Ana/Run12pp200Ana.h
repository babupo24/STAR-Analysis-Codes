

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
#include "TH2.h"
#include "TNtuple.h"
#include <iostream>
#include <fstream>

using namespace std;

class Run12pp200Ana
{
public:
        TChain *fChain;
        // pointer to the analyzed TTree or TChain
        Int_t fCurrent;
        // current Tree number in a TChain

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
        Double_t fverRank; // included by Babu  11/22/2019
        Double_t fpT[2422];
        Double_t fp[2422];
        Double_t fpX[2422];
        Double_t fpY[2422];
        Double_t fpZ[2422];
        Double_t feta[2422];
        Double_t fphi[2422];
        Short_t fcharge[2422];
        Double_t fnSigmaPion[2422];
        Double_t fnSigmaKaon[2422];
        Double_t fnSigmaProton[2422];
        Double_t fnSigmaElectron[2422];

        Double_t fnSigmaPionTof[2422];
        Double_t fnSigmaKaonTof[2422];
        Double_t fnSigmaProtonTof[2422];
        Double_t fnSigmaElectronTof[2422];

        Double_t fdEdx[2422];
        Double_t fdca[2422];
        UShort_t ffitPts[2422];
        UShort_t ffitPtsPoss[2422];
        UShort_t fhitsdedx[2422];
        Double_t fBetaToF[2422];

        // List of branches
        TBranch *b_frefmult;
        TBranch *b_fmaxpar;
        TBranch *b_ffillNum;
        TBranch *b_frunNum;
        TBranch *b_fspinconfig;
        TBranch *b_ftrigger;
        TBranch *b_fVZ;
        TBranch *b_fvpdVz;
        TBranch *b_fevTime;
        TBranch *b_fverRank;
        TBranch *b_fpT;
        TBranch *b_fp;
        TBranch *b_fpX;
        TBranch *b_fpY;
        TBranch *b_fpZ;
        TBranch *b_feta;
        TBranch *b_fphi;
        TBranch *b_fcharge;
        TBranch *b_fnSigmaPion;
        TBranch *b_fnSigmaKaon;
        TBranch *b_fnSigmaProton;
        TBranch *b_fnSigmaElectron;
        TBranch *b_fnSigmaPionTof;
        TBranch *b_fnSigmaKaonTof;
        TBranch *b_fnSigmaProtonTof;
        TBranch *b_fnSigmaElectronTof;

        TBranch *b_fdEdx;
        TBranch *b_fdca;
        TBranch *b_ffitPts;
        TBranch *b_ffitPtsPoss;
        TBranch *b_fhitsdedx;
        TBranch *b_fBetaToF;

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
        // virtual void     loadPol();
        inline void setHome(TDirectory *pHome)
        {
                cout << "setting home.." << endl;
                home = pHome;
        }

        static const int n = 23;
        int trigger;

protected:
        TDirectory *home;

        TH2D *hSigmaPivsK[2][17];
        TH2D *hSigmaPivsP[2][17];
        TH2D *hSigmaPivsE[2][17];

        TH1D *hSigmaPion[2][17];
        TH2D *hPidEdxVsP[2][17];
        TH2D *hdEdxVsP[2][17];

        TH1D *hSigmaPionTOF[2][17]; // for tof tracks
        TH1D *hSigmaPionTPC[2][17]; // for tof tracks

        TH2D *hSigmaPionTPCvsTOF[2][17]; // for tof tracks

        TH1D *hp[2][17];

        TH2D *hSigmaPivsK_M0[5];
        TH2D *hSigmaPivsP_M0[5];
        TH2D *hSigmaPivsE_M0[5];
        TH1D *hSigmaPion_M0[5];
       
        TH2D *hSigmaPivsKtof_M0[5];
        TH2D *hSigmaPivsPtof_M0[5];
        TH2D *hSigmaPivsEtof_M0[5];
        TH1D *hSigmaPiontof_M0[5];
       
        TH2D *hSigmaPivsKboth_M0[5];
        TH2D *hSigmaPivsPboth_M0[5];
        TH2D *hSigmaPivsEboth_M0[5];
        TH1D *hSigmaPionboth_M0[5];

private:
        ClassDef(Run12pp200Ana, 1);
};

#endif
