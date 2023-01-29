

/////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 22 17:04:53 2010 by ROOT version 5.22/00
// from TTree ftree/LRC 
// found on file: 00E5DE07942C320D8D35F5F8AD40C261_7total.root
//////////////////////////////////////////////////////////

#ifndef Iff2015_h
#define Iff2015_h

#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TNtuple.h"
#include <iostream>
#include <fstream>

using namespace std;


class Iff2015 {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           frefmult;
   Int_t           fmaxpar;
   Int_t           ffillNum;
   Int_t           frunNum;
   Int_t           fspinconfig;
   vector<double>   *ftrigger;
   Double_t        fVZ;
   float           fvpdVz;// Included by Babu 11/22/2019
   unsigned int    fevTime;
   Double_t        fverRank;// included by Babu  11/22/2019
   Double_t        fpT[2422];   //[fmaxpar]
   Double_t        fp[2422];   //[fmaxpar]
   Double_t        feta[2422];   //[fmaxpar]
   Double_t        fphi[2422];   //[fmaxpar]
   Short_t         fcharge[2422];   //[fmaxpar]
   Double_t        fnSigmaPion[2422];   //[fmaxpar]
   Double_t        fnSigmaKaon[2422];   //[fmaxpar]
   Double_t        fnSigmaProton[2422];   //[fmaxpar]
   Double_t        fnSigmaElectron[2422];   //[fmaxpar]
   Double_t        fnSigmaPionTof[2422];   //[fmaxpar]
   Double_t        fnSigmaKaonTof[2422];   //[fmaxpar]
   Double_t        fnSigmaProtonTof[2422];   //[fmaxpar]
   Double_t        fnSigmaElectronTof[2422];   //[fmaxpar]
   Double_t        fdEdx[2422];   //[fmaxpar]
   Double_t        fdca[2422];   //[fmaxpar]
   UShort_t        ffitPts[2422];   //[fmaxpar]
   UShort_t        ffitPtsPoss[2422];   //[fmaxpar]
   UShort_t        fhitsdedx[2422];   //[fmaxpar]
   Double_t        fBetaToF[2422];//[fmaxpar]

   // List of branches
   TBranch        *b_frefmult;   //!
   TBranch        *b_fmaxpar;   //!
   TBranch        *b_ffillNum;  //!
   TBranch        *b_frunNum;  //!
   TBranch        *b_fspinconfig;  //!
   TBranch        *b_ftrigger;  //!
   TBranch        *b_fVZ;
   TBranch        *b_fvpdVz;
   TBranch        *b_fevTime;
   TBranch        *b_fverRank;   //!
   TBranch        *b_fpT;   //!
   TBranch        *b_fp;   //!
   TBranch        *b_feta;   //!
   TBranch        *b_fphi;   //!
   TBranch        *b_fcharge;   //!
   TBranch        *b_fnSigmaPion;   //!
   TBranch        *b_fnSigmaKaon;   //!
   TBranch        *b_fnSigmaProton;   //!
   TBranch        *b_fnSigmaElectron;   //!
   TBranch        *b_fnSigmaPionTof;   //!
   TBranch        *b_fnSigmaKaonTof;   //!
   TBranch        *b_fnSigmaProtonTof;   //!
   TBranch        *b_fnSigmaElectronTof;   //!
   TBranch        *b_fdEdx;   //!
   TBranch        *b_fdca;   //!
   TBranch        *b_ffitPts;   //!
   TBranch        *b_ffitPtsPoss;   //!
   TBranch        *b_fhitsdedx;   //!
   TBranch        *b_fBetaToF;  //!

   Iff2015(char* ifile);
   virtual ~Iff2015();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Finish(char* ofile);
   virtual void     Show(Long64_t entry = -1);
   //virtual void     loadPol();
   inline   void setHome(TDirectory* pHome)
     {
       cout <<"setting home.." <<endl;
       home=pHome;
     }

static const int n=23;
	int trigger;
 	TH2D* hTpcVsTofM[2][5];
        TH2D* hTpcVsTofpT[2][5];
        TH2D* hTpcVsTofMBkg[2][5];
        TH2D* hTpcVsTofpTBkg[2][5];



/*int trigid[23];
        trigid[0] = 480003; trigid[1] = 480004;  trigid[2] = 480005;  
        trigid[3] = 480904; trigid[4] = 480007;  trigid[5] = 480201;  
        trigid[6] = 480202; trigid[7] = 480203;  trigid[8] = 480204;  
        trigid[9] = 480205; trigid[10]= 480206; trigid[11]= 480301;  
        trigid[12]= 480401; trigid[13]= 480402;  trigid[14]= 480403;  
        trigid[15]= 480404; trigid[16]= 480405;  trigid[17]= 480406;  
        trigid[18]= 480714; trigid[19]= 480411;  trigid[20]= 480414;  
        trigid[21]= 480415; trigid[22]= 480501; 
char trigname[23];
        trigname[0] = "BBCMB";	        trigname[1] = "VPDMB-novtx";         	trigname[2] = "ZDCMB-trgonly";
        trigname[3] = "VPDMB30";     	trigname[4] = "VPDMB-5-trgonly";      	trigname[5] = "BHT0*VPDMB-5";
        trigname[6] = "BHT1_VPDMB30";	trigname[7] = "BHT0*BBCMB";	 	trigname[8] = "BHT1*BBCMB";
        trigname[9] = "BHT2_BBCMB";   	trigname[10]= "BHT1*VPDMB-30-nobsmd";	trigname[11]= "EHT0";
        trigname[12]= "JP2";	        trigname[13]= "JP2-bsmd";           	trigname[14]= "AJP";
        trigname[15]= "JP1";          	trigname[16]= "JP2*L2JetHigh";       	trigname[17]= "BHT2*BJP1*L2Bgamma";
        trigname[18]= "RP_CPEI";	trigname[19]= "JP2";                 	trigname[20]= "JP1";
        trigname[21]= "JP2_L2JetHigh";	trigname[22]= "EHT0*EJP1*L2Egamma";

*/

 protected:
	//ntuples for TPC pid Cut
   TNtuple *ntuple1tpc;   
   TNtuple *ntuple2tpc;
   TNtuple *ntuple3tpc;
   TNtuple *ntuple4tpc;
   TNtuple *ntuple5tpc;
	//ntuples for TPC and TOF Cut(IF TOF avalilable)
   TNtuple *ntuple1f;   
   TNtuple *ntuple2f;
   TNtuple *ntuple3f;
   TNtuple *ntuple4f;
   TNtuple *ntuple5f;
   	//ntuples for TPC and TOF cut (only tracks with TOF)
   TNtuple *ntuple1tof;   
   TNtuple *ntuple2tof;
   TNtuple *ntuple3tof;
   TNtuple *ntuple4tof;
   TNtuple *ntuple5tof;

   	//ntuples for background  (only tracks with TOF)
   TNtuple *ntuple1tofbkg;   
   TNtuple *ntuple2tofbkg;
   TNtuple *ntuple3tofbkg;
   TNtuple *ntuple4tofbkg;
   TNtuple *ntuple5tofbkg;
   	//ntuples for background after TPC and TOF cut (only tracks with TOF)
   TNtuple *ntuple1bkg;   
   TNtuple *ntuple2bkg;
   TNtuple *ntuple3bkg;
   TNtuple *ntuple4bkg;
   TNtuple *ntuple5bkg;


   TNtuple *ntuple6;

  
    TDirectory *home;


 private:
   ClassDef(Iff2015,1);


};

#endif





