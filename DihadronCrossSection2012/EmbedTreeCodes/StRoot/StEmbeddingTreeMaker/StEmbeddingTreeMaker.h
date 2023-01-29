/*
 * StEmbeddingTreeMaker.h
 * J. Kevin Adkins, University of Kentucky
 * January 25, 2016
 *
 */

#ifndef STAR_StEmbeddingTreeMaker
#define STAR_StEmbeddingTreeMaker

#ifndef StMaker_H
#include "StMaker.h"
#endif
#include <vector>

class TTree;
class TFile;
class TH2F;
class TChain;
class StEvent;
class StMuDstMaker;
class StMuDst;
class StMuEvent;
class StRunInfo;
class StMuTrack;
class StMuMcTrack;
class StMuMcVertex;
class StEventInfo;
class StEmcCollection;
class StBemcTables;
class EEmcGeomSimple;
class StEEmcDbMaker;
class StEEmcDb;
class EEmcDbItem;
class StMuPrimaryVertex;
class StPythiaEvent;
class StTriggerSimuMaker;

class StEmbeddingTreeMaker : public StMaker
{
private:
  StMuDstMaker *mMuDstMaker;
  StMuDst *muDst;
  StMuEvent *muEvent;
  StMuPrimaryVertex *mVertex;
  StMuTrack *mTrack;
  StMuMcTrack *mcTrack;
  StMuMcVertex *mcVertex;
  StBemcTables *mBemcTables;
  StEEmcDb *mEemcDb;
  EEmcGeomSimple *mEEmcGeom;
  StPythiaEvent *mPythiaEvent;
  TFile *mFile;
  TChain *mPythiaChain;
  StTriggerSimuMaker *mTriggerSimu;

protected:
  Double_t yPartX2, bPartX1, yPartE, bPartE, yPartPhi, bPartPhi, bPartEta, yPartEta, bPartPt, yPartPt, bPartP, yPartP;
  Int_t bPartId, yPartId;
  Double_t fPart2X, fPart1X, fPart1E, fPart2E, fPart1Phi, fPart2Phi, fPart1Eta, fPart2Eta, fPart1Pt, fPart2Pt, fPart1P, fPart2P;
  Int_t fPart1Id, fPart2Id;
  Int_t nPrimVerts, nPrimTracks;
  Int_t mTrackMulti;
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
  vector<int> ftrigger;

  Int_t trigJP0;
  Int_t trigJP1;
  Int_t trigJP2;

  Int_t fmaxpar1;
  Int_t fmaxpar;

  Double_t fVZ;
  Double_t fVY;
  Double_t fVX;
  Double_t fVZ_mc;
  Double_t fVZ_pyth;
  Double_t fVY_pyth;
  Double_t fVX_pyth;
  Double_t fverRank;

  Int_t fTrackCounter1;
  Int_t fTrackCounter;
  Int_t bPartonId;
  Double_t partonicPtBin;

  Int_t idTruth[fMaxHit];
  Double_t fpT[fMaxHit];
  Double_t fp[fMaxHit];
  Double_t feta[fMaxHit];
  Double_t fphi[fMaxHit];

  Int_t statusCode_pyth[fMaxHit];
  Double_t fpT_pyth[fMaxHit];
  Double_t fp_pyth[fMaxHit];
  Double_t feta_pyth[fMaxHit];
  Double_t fphi_pyth[fMaxHit];
  Double_t fpId_pyth[fMaxHit];
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
  Double_t fPathLengthToF[fMaxHit];
  Double_t fBetaToF[fMaxHit];
  Int_t trueID[fMaxHit];

  Double_t fpT_mc[fMaxHit];
  Double_t fp_mc[fMaxHit];
  Double_t feta_mc[fMaxHit];
  Double_t fphi_mc[fMaxHit];
  Int_t fpId_mc[fMaxHit]; // mc particle identification
  Int_t fId_mc[fMaxHit];  // mc id to match with idTruth

public:
  StEmbeddingTreeMaker(const char *, const char *, StMuDstMaker *, TChain *); // constructor
  virtual ~StEmbeddingTreeMaker();                                            // destructor
  virtual Int_t Init();
  virtual Int_t Make();
  virtual Int_t Finish();
  Bool_t useVertex(StMuPrimaryVertex *);
  Bool_t useTrack(StMuTrack *);

  // Displayed on session exit, leave it as-is please ...
  virtual const char *GetCVS() const
  {
    static const char cvs[] = "Tag $Name:  $ $Id: StEmbeddingTreeMaker.h,v 1.15 2003/09/10 19:47:43 perev Exp $ built " __DATE__ " " __TIME__;
    return cvs;
  }

  ClassDef(StEmbeddingTreeMaker, 0) // StAF chain virtual base class for Makers
};

#endif
