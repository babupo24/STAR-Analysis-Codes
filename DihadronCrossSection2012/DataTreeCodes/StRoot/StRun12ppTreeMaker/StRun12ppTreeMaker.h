/*
 * StRun12ppTreeMaker.h
 * Babu Pokhrel, Temple University
 * May 1, 2022
 *
 */

#ifndef STAR_StRun12ppTreeMaker
#define STAR_StRun12ppTreeMaker

#ifndef StMaker_H
#include "StMaker.h"
#endif
#include <vector>
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"

#include "StRoot/StTriggerUtilities/StTriggerSimuMaker.h"

class TTree;
class TFile;
class TH2F;
class TChain;
class StEvent;
class StPicoDst;
class StPicoEvent;
class StPicoTrack;
class StPicoDstMaker;
class StPicoBTofPidTraits;
class StTriggerSimuMaker;

class StEventInfo;
class StEmcCollection;
class StBemcTables;
class EEmcGeomSimple;
class StEEmcDbMaker;
class StEEmcDb;
class EEmcDbItem;
class StPythiaEvent;
class StTriggerSimuMaker;

class StRun12ppTreeMaker : public StMaker
{
private:
  StPicoDst *picoDst;
  StPicoDstMaker *mPicoDstMaker;
  StPicoEvent *picoEvent;
  StPicoTrack *mTrack;
  StPicoBTofPidTraits *trait;

  StBemcTables *mBemcTables;
  StEEmcDb *mEemcDb;
  EEmcGeomSimple *mEEmcGeom;

  StTriggerSimuMaker *mTrigSimu;

  TFile *mFile;

protected:
  Double_t pi;
  TString rootFilename;

  TTree *ftree;
  static const int fMaxHit = 3000;
  Int_t ffillNum;
  Int_t fevtNum;
  unsigned int fevTime;
  Double_t fvpdVz;
  Int_t frunNum;
  Int_t frefmult;
  vector<unsigned int> ftrigger;
  vector<int> ftrigger_simu;
  Int_t fmaxpar;
  Int_t nPrimTracks;
  Double_t fVZ;
  Double_t fVY;
  Double_t fVX;
  Double_t fverRank;

  Int_t fTrackCounter;
  Double_t fpT[fMaxHit];
  Double_t fp[fMaxHit];

  Double_t fpX[fMaxHit];
  Double_t fpY[fMaxHit];
  Double_t fpZ[fMaxHit];

  Double_t feta[fMaxHit];
  Double_t fphi[fMaxHit];
  Double_t fnSigmaPion[fMaxHit];
  Double_t fnSigmaKaon[fMaxHit];
  Double_t fnSigmaProton[fMaxHit];
  Double_t fnSigmaElectron[fMaxHit];
  Short_t fcharge[fMaxHit];
  Double_t fdEdx[fMaxHit];
  Double_t fdca[fMaxHit];
  UShort_t ffitPts[fMaxHit];
  UShort_t ffitPtsPoss[fMaxHit];
  UShort_t fhitsdedx[fMaxHit];
  Double_t fBetaToF[fMaxHit];

  Double_t fnSigmaPionTof[fMaxHit];
  Double_t fnSigmaProtonTof[fMaxHit];
  Double_t fnSigmaKaonTof[fMaxHit];
  Double_t fnSigmaElectronTof[fMaxHit];

public:
  StRun12ppTreeMaker(const char *, const char *, StPicoDstMaker *); // constructor
  virtual ~StRun12ppTreeMaker();                                    // destructor
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  Bool_t useTrack(StPicoTrack *);
  Bool_t dcaCut(Double_t, Double_t);
  TDirectory *home;
  inline void setHome(TDirectory *pHome)
  {
    cout << "setting home.." << endl;
    home = pHome;
  }

  ClassDef(StRun12ppTreeMaker, 0) // StAF chain virtual base class for Makers
};

#endif
