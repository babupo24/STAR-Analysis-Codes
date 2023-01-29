
/*
 * StRun12ppTreeMaker.cxx
 * Babu Pokhrel, Temple University (04/15/2022)
 * Code to create trees from MuDst's for Dihadron Cross-Secction
 * analysis.
 * Code structure extracted from J. Kevin Adkins, University of Kentucky
 *
 */

#include "StRun12ppTreeMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"

#include "TVector3.h"
#include "TSystem.h"
// std
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <math.h>

// root
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TChain.h"

#include "StChain.h"
#include "StEvent/StEvent.h"

ClassImp(StRun12ppTreeMaker)

    StRun12ppTreeMaker::StRun12ppTreeMaker(const char *name, const char *outfilename, StPicoDstMaker *picoMaker) : StMaker(name)
{
  mPicoDstMaker = picoMaker;
  assert(mPicoDstMaker);
  rootFilename = TString(outfilename);
  mFile = 0;
  picoDst = 0;
  picoEvent = 0;
  mTrack = 0;
}

StRun12ppTreeMaker::~StRun12ppTreeMaker()
{
}

Int_t StRun12ppTreeMaker::Init()
{

  pi = TMath::Pi();

  mFile = new TFile(rootFilename, "RECREATE");
  ftree = new TTree("ftree", "dihadron ");

  // tree branches
  ftree->Branch("fevtNum", &fevtNum, "fevtNum/I");
  ftree->Branch("ffillNum", &ffillNum, "ffillNum/I");
  ftree->Branch("frunNum", &frunNum, "frunNum/I");
  ftree->Branch("frefmult", &frefmult, "frefmult/I");
  ftree->Branch("ftrigger", &ftrigger);
  ftree->Branch("trigJP0", &trigJP0, "trigJP0/I");
  ftree->Branch("trigJP1", &trigJP1, "trigJP1/I");
  ftree->Branch("trigJP2", &trigJP2, "trigJP2/I");

  ftree->Branch("fmaxpar", &fmaxpar, "fmaxpar/I");

  ftree->Branch("fVZ", &fVZ, "fVZ/D");
  ftree->Branch("fVX", &fVX, "fVX/D");
  ftree->Branch("fVY", &fVY, "fVY/D");

  ftree->Branch("fverRank", &fverRank, "fverRank/D");

  ftree->Branch("fpT", fpT, "fpT[fmaxpar]/D");
  ftree->Branch("fp", fp, "fp[fmaxpar]/D");
  ftree->Branch("fpX", fpX, "fpX[fmaxpar]/D");
  ftree->Branch("fpY", fpY, "fpY[fmaxpar]/D");
  ftree->Branch("fpZ", fpZ, "fpZ[fmaxpar]/D");
  // tpc information
  ftree->Branch("feta", feta, "feta[fmaxpar]/D");
  ftree->Branch("fphi", fphi, "fphi[fmaxpar]/D");
  ftree->Branch("fcharge", fcharge, "fcharge[fmaxpar]/S");
  ftree->Branch("fnSigmaPion", fnSigmaPion, "fnSigmaPion[fmaxpar]/D");
  ftree->Branch("fnSigmaKaon", fnSigmaKaon, "fnSigmaKaon[fmaxpar]/D");
  ftree->Branch("fnSigmaProton", fnSigmaProton, "fnSigmaProton[fmaxpar]/D");
  ftree->Branch("fnSigmaElectron", fnSigmaElectron, "fnSigmaElectron[fmaxpar]/D");
  ftree->Branch("fdEdx", fdEdx, "fdEdx[fmaxpar]/D");
  ftree->Branch("fdca", fdca, "fdca[fmaxpar]/D");
  ftree->Branch("ffitPts", ffitPts, "ffitPts[fmaxpar]/s");
  ftree->Branch("ffitPtsPoss", ffitPtsPoss, "ffitPtsPoss[fmaxpar]/s");
  ftree->Branch("fhitsdedx", fhitsdedx, "fhitsdedx[fmaxpar]/s");
  ftree->Branch("fevTime", &fevTime, "fevTime/i");
  ftree->Branch("fvpdVz", &fvpdVz, "fvpdVz/D");
  // tof information
  ftree->Branch("fBetaToF", fBetaToF, "fBetaToF[fmaxpar]/D");
  ftree->Branch("fnSigmaPionTof", fnSigmaPionTof, "fnSigmaPionTof[fmaxpar]/D");
  ftree->Branch("fnSigmaProtonTof", fnSigmaProtonTof, "fnSigmaProtonTof[fmaxpar]/D");
  ftree->Branch("fnSigmaKaonTof", fnSigmaKaonTof, "fnSigmaKaonTof[fmaxpar]/D");
  ftree->Branch("fnSigmaElectronTof", fnSigmaElectronTof, "fnSigmaElectronTof[fmaxpar]/D");

  return StMaker::Init();
}

Int_t StRun12ppTreeMaker::Make()
{
  assert(mPicoDstMaker);

  picoDst = mPicoDstMaker->picoDst();
  assert(picoDst);

  picoEvent = picoDst->event();
  assert(picoEvent);

  ftrigger = picoEvent->triggerIds();
  trigJP0 = 0, trigJP1 = 0, trigJP2 = 0;
  //cout << "Trigger Size: " << ftrigger.size() << endl;
  for (unsigned int i = 0; i < ftrigger.size(); i++)
  {
    if (ftrigger.at(i) == 370601)
      trigJP0 = 1;
    if (ftrigger.at(i) == 370611)
      trigJP1 = 1;
    if (ftrigger.at(i) == 370621)
      trigJP2 = 1;

    cout << "Ids: " << ftrigger.at(i) << endl;
  }

  fverRank = picoEvent->ranking();
  TVector3 vertexPos = picoEvent->primaryVertex();
  fVZ = vertexPos.Z();
  fVX = vertexPos.X();
  fVY = vertexPos.Y();
  fvpdVz = picoEvent->vzVpd();
  fevTime = picoEvent->time();

  ffillNum = (Int_t)(picoEvent->fillId());
  frunNum = picoEvent->runId();
  fevtNum = picoEvent->eventId();
  frefmult = picoEvent->refMult();

  nPrimTracks = 0;
  nPrimTracks = picoDst->numberOfTracks();
  fTrackCounter = 0;
  for (Int_t pTrack = 0; pTrack < nPrimTracks; pTrack++)
  {
    mTrack = picoDst->track(pTrack);
    assert(mTrack);

    if (!useTrack(mTrack))
      continue;
    Double_t gdca = fabs(mTrack->gDCA(fVX, fVY, fVZ));
    if (!dcaCut(mTrack->pPt(), gdca))
      continue;

    fpT[fTrackCounter] = mTrack->pPt();
    fp[fTrackCounter] = mTrack->pMom().Mag();
    fpX[fTrackCounter] = mTrack->pMom().X();
    fpY[fTrackCounter] = mTrack->pMom().Y();
    fpZ[fTrackCounter] = mTrack->pMom().Z();
    feta[fTrackCounter] = mTrack->pMom().Eta();
    fphi[fTrackCounter] = mTrack->pMom().Phi();
    fnSigmaPion[fTrackCounter] = mTrack->nSigmaPion();
    fnSigmaKaon[fTrackCounter] = mTrack->nSigmaKaon();
    fnSigmaProton[fTrackCounter] = mTrack->nSigmaProton();
    fnSigmaElectron[fTrackCounter] = mTrack->nSigmaElectron();
    fcharge[fTrackCounter] = mTrack->charge();
    fdEdx[fTrackCounter] = mTrack->dEdx();
    fdca[fTrackCounter] = gdca;
    ffitPts[fTrackCounter] = mTrack->nHitsFit();
    fhitsdedx[fTrackCounter] = mTrack->nHitsDedx();
    ffitPtsPoss[fTrackCounter] = mTrack->nHitsPoss();
    // get TOF information for TOF tracks
    Int_t tofIndex = mTrack->bTofPidTraitsIndex(); // returns -1 if track doesn't have TOF hits
    trait = picoDst->btofPidTraits(tofIndex);
    if (tofIndex >= 0)
    { // selects only track that has TOF hits
      fBetaToF[fTrackCounter] = trait->btofBeta();
      fnSigmaPionTof[fTrackCounter] = trait->nSigmaPion();
      fnSigmaProtonTof[fTrackCounter] = trait->nSigmaProton();
      fnSigmaKaonTof[fTrackCounter] = trait->nSigmaKaon();
      fnSigmaElectronTof[fTrackCounter] = trait->nSigmaElectron();
    }
    else if (tofIndex < 0)
    {
      fBetaToF[fTrackCounter] = -999;
      fnSigmaPionTof[fTrackCounter] = -999;
      fnSigmaProtonTof[fTrackCounter] = -999;
      fnSigmaKaonTof[fTrackCounter] = -999;
      fnSigmaElectronTof[fTrackCounter] = -999;
    }
    // cout << "Tof Index: " << tofIndex << ",  SigmaPion TOF: " << fnSigmaPionTof[fTrackCounter] << endl;

    fTrackCounter++;
  }
  // Fill track multiplicity
  fmaxpar = fTrackCounter;
  cout << "Max track count: " << fmaxpar << endl;
  if (fmaxpar > 1)
    ftree->Fill();

  return kStOK;
} // maker

Bool_t StRun12ppTreeMaker::useTrack(StPicoTrack *track)
{
  if (track->isPrimary())
  {
    if (track->nHitsFit() > 15)
    {
      if (track->pPt() > 0.2 && fabs(track->pMom().Eta()) < 1.0)
      {
        return true;
      }
    }
  }
  return false;
}

Bool_t StRun12ppTreeMaker::dcaCut(Double_t pT, Double_t dca)
{
  if ((pT < 0.5 && dca < 2.0) || (pT > 0.5 && pT < 1.5 && dca < ((-1) * pT + 2.5)) || (pT > 1.5 && dca < 1.0))
  {
    return true;
  }
  return false;
}

Int_t StRun12ppTreeMaker::Finish()
{
  ftree->SetDirectory(mFile);
  ftree->Write();
  mFile->Write();
  mFile->Close();
  delete mFile;

  return kStOk;
}
