

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
        TChain *fChain; // pointer to the analyzed TTree or TChain
        Int_t fCurrent; // current Tree number in a TChain

        // Declaration of leaf types
        Int_t frefmult;
        Int_t fmaxpar;
        Int_t fmaxpar1;
        Int_t ffillNum;
        Int_t frunNum;
        Int_t fspinconfig;
        vector<int> *ftrigger;
        Int_t trigJP0;
        Int_t trigJP1;
        Int_t trigJP2;
        Double_t fVZ;
        Double_t fvpdVz;
        unsigned int fevTime;
        Double_t fverRank;
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
        Double_t fdEdx[2422];
        Double_t fdca[2422];
        UShort_t ffitPts[2422];
        UShort_t ffitPtsPoss[2422];
        UShort_t fhitsdedx[2422];
        Double_t fBetaToF[2422];
        Int_t idTruth[2422];

        // Part level info
        Double_t partonicPtBin;
        Int_t bPartId;
        Double_t bPartX1;
        Double_t bPartE;
        Double_t bPartEta;
        Double_t bPartPt;
        Double_t bPartPhi;

        Int_t yPartId;
        Double_t yPartX2;
        Double_t yPartE;
        Double_t yPartPhi;
        Double_t yPartEta;
        Double_t yPartPt;

        Int_t fPart1Id;
        Double_t fPart1Pt;
        Double_t fPart1Eta;
        Double_t fPart1Phi;
        Double_t fPart1E;
        Double_t fPart1X;

        Int_t fPart2Id;
        Double_t fPart2Pt;
        Double_t fPart2Eta;
        Double_t fPart2Phi;
        Double_t fPart2E;
        Double_t fPart2X;
        Double_t ParticPtBin;
        Double_t fVZ_pyth;
        Double_t fpT_pyth[2422];
        Double_t fp_pyth[2422];
        Double_t feta_pyth[2422];
        Double_t fphi_pyth[2422];
        Double_t fpId_pyth[2422];
        Int_t statusCode_pyth[2422];
        Bool_t isPrim_pyth[2422];

        Double_t fVZ_mc;
        Double_t fpT_mc[2422];
        Double_t fp_mc[2422];
        Double_t feta_mc[2422];
        Double_t fphi_mc[2422];
        Int_t fpId_mc[2422];
        Int_t fId_mc[2422];

        // List of branches
        TBranch *b_frefmult;
        TBranch *b_fmaxpar;
        TBranch *b_fmaxpar1;
        TBranch *b_ffillNum;
        TBranch *b_frunNum;
        TBranch *b_fspinconfig;
        TBranch *b_ftrigger;
        TBranch *b_trigJP0;
        TBranch *b_trigJP1;
        TBranch *b_trigJP2;
        TBranch *b_fVZ;
        TBranch *b_fvpdVz;
        TBranch *b_fevTime;
        TBranch *b_fverRank;
        TBranch *b_fpT;
        TBranch *b_fp;
        TBranch *b_idTruth;

        TBranch *b_feta;
        TBranch *b_fphi;
        TBranch *b_fcharge;
        TBranch *b_fnSigmaPion;
        TBranch *b_fnSigmaKaon;
        TBranch *b_fnSigmaProton;
        TBranch *b_fnSigmaElectron;
        TBranch *b_fdEdx;
        TBranch *b_fdca;
        TBranch *b_ffitPts;
        TBranch *b_ffitPtsPoss;
        TBranch *b_fhitsdedx;
        TBranch *b_fBetaToF;

        TBranch *b_fVZ_pyth;
        TBranch *b_bPartId;
        TBranch *b_bPartEta;
        TBranch *b_bPartPt;
        TBranch *b_bPartPhi;
        TBranch *b_bPartX1;
        TBranch *b_bPartE;

        TBranch *b_yPartId;
        TBranch *b_yPartX2;
        TBranch *b_yPartPhi;
        TBranch *b_yPartEta;

        TBranch *b_yPartE;
        TBranch *b_yPartPt;

        TBranch *b_fPart1Id;
        TBranch *b_fPart1Eta;
        TBranch *b_fPart1Phi;
        TBranch *b_fPart1E;
        TBranch *b_fPart1X;
        TBranch *b_fPart1Pt;

        TBranch *b_fPart2X;
        TBranch *b_fPart2Id;
        TBranch *b_fPart2Eta;
        TBranch *b_fPart2Phi;
        TBranch *b_fPart2E;
        TBranch *b_fPart2Pt;
        TBranch *b_partonicPtBin;

        TBranch *b_isPrim_pyth;
        TBranch *b_statusCode_pyth;
        TBranch *b_fpT_pyth;
        TBranch *b_fp_pyth;
        TBranch *b_feta_pyth;
        TBranch *b_fphi_pyth;
        TBranch *b_fpId_pyth;

        TBranch *b_fVZ_mc;
        TBranch *b_fpT_mc;
        TBranch *b_fp_mc;
        TBranch *b_feta_mc;
        TBranch *b_fphi_mc;
        TBranch *b_fId_mc;
        TBranch *b_fpId_mc;

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
        inline void setHome(TDirectory *pHome)
        {
                cout << "setting home.." << endl;
                home = pHome;
        }
        virtual void matchFinder(double, double, double, double, int, double[], double[], double[], double[], int[], double &, double &, double &, double &, double &, double &, double &, double &, double &);

        static const int n = 23;
        int trigger;
        TH1D *hVZJP0;
        TH1D *hVZJP1;
        TH1D *hVZJP2;

        TH1D *hpartonicPt;
        // histograms for resolution
        TH1D *hMDiffInt;
        TH1D *hMGenRec[12];
        TH1D *hMGenRecU[12];

        TH1D *hxsecMinvGen;
        TH1D *hxsecMinvRecJP0;
        TH1D *hxsecMinvRecJP1;
        TH1D *hxsecMinvRecJP2;
        // histograms for unfolding
        TH2D *hMGenVsRecJP0;
        TH2D *hMGenVsRecJP1;
        TH2D *hMGenVsRecJP2;

        TH2D *hMGenVsRecJP0NoUF;
        TH2D *hMGenVsRecJP1NoUF;
        TH2D *hMGenVsRecJP2NoUF;
        // migration matrix in pseudorapidity bins
        TH2D *hMGenVsRecJP0Gt;
        TH2D *hMGenVsRecJP1Gt;
        TH2D *hMGenVsRecJP2Gt;

        TH2D *hMGenVsRecJP0Lt;
        TH2D *hMGenVsRecJP1Lt;
        TH2D *hMGenVsRecJP2Lt;

        TH1D *hMGenJP0;
        TH1D *hMGenJP1;
        TH1D *hMGenJP2;

        TH1D *hMGenJP0f;
        TH1D *hMGenJP1f;
        TH1D *hMGenJP2f;

        TH1D *hMRecJP2;
        TH1D *hMRecJP1;
        TH1D *hMRecJP0;

        TH1D *hMRecJP2All;
        TH1D *hMRecJP1All;
        TH1D *hMRecJP0All;

        TH1D *hMRecJP2AllLt;
        TH1D *hMRecJP1AllLt;
        TH1D *hMRecJP0AllLt;

        TH1D *hMRecJP2AllGt;
        TH1D *hMRecJP1AllGt;
        TH1D *hMRecJP0AllGt;

        TH1D *hMRecJP2x;
        TH1D *hMRecJP1x;
        TH1D *hMRecJP0x;

        TH1D *hMRecJP2AllXbin;
        TH1D *hMRecJP1AllXbin;
        TH1D *hMRecJP0AllXbin;

        TH1D *hMRecJP2bk1;
        TH1D *hMRecJP1bk1;
        TH1D *hMRecJP0bk1;

        TH1D *hMRecJP2bk2;
        TH1D *hMRecJP1bk2;
        TH1D *hMRecJP0bk2;

        // total background in pseudorapidity bins
        TH1D *hMRecJP2bkGt;
        TH1D *hMRecJP1bkGt;
        TH1D *hMRecJP0bkGt;

        TH1D *hMRecJP2bkLt;
        TH1D *hMRecJP1bkLt;
        TH1D *hMRecJP0bkLt;

        TH1D *hMRecJP2bkx; // background in x-section bins
        TH1D *hMRecJP1bkx;
        TH1D *hMRecJP0bkx;
        // histograms for unfolding test
        TH2D *hMGenVsRecTestM;
        TH1D *hMRecTestInput;
        TH1D *hMGenTestInput;
        TH1D *hMRecTestControlPlot1st;
        TH1D *hMGenTestControlPlot1st;
        TH1D *hMRecTestControlPlot2nd;
        TH1D *hMGenTestControlPlot2nd;

protected:
        TDirectory *home;

private:
        ClassDef(Run12pp200Ana, 1);
};

#endif
