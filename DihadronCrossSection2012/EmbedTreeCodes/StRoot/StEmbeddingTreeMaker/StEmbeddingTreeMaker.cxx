/*
 * StEmbeddingTreeMaker.cxx
 * Babu Pokhrel, Temple University (04/15/2022)
 * Code to create trees from embedding MuDst and pythia files for Dihadron Cross-Secction
 * analysis.
 * Code structure extracted from J. Kevin Adkins, University of Kentucky
 *
 */

#include "StEmbeddingTreeMaker.h"
#include "TDataSetIter.h"
#include "StDAQMaker/StDAQReader.h"

#include "StEventTypes.h"
#include "TVector3.h"

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

// STAR
#include "StMessMgr.h"

// StMuDstMaker
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuEmcCollection.h"
#include "StMuDSTMaker/COMMON/StMuDst.h"
#include "StMuDSTMaker/COMMON/StMuEmcUtil.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcTrack.h"
#include "StMuDSTMaker/COMMON/StMuMcVertex.h"

// Vertex & Trigger Info
#include "StMuDSTMaker/COMMON/StMuPrimaryVertex.h"
#include "StMuDSTMaker/COMMON/StMuTriggerIdCollection.h"

// StEmc
#include "StEmcClusterCollection.h"
#include "StEmcPoint.h"
#include "StEmcUtil/geometry/StEmcGeom.h"
#include "StEmcUtil/others/emcDetectorName.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/StEmcRawMaker.h"
#include "StEmcRawMaker/defines.h"
#include "StEmcUtil/database/StBemcTables.h"
#include "StEvent/StEmcCollection.h"

#include "StChain.h"
#include "StMcEventMaker/StMcEventMaker.h"
#include "StMcEventTypes.hh"
#include "StMcEvent.hh"
#include "StEvent/StEvent.h"

#include "TDataSetIter.h"

#include "StEEmcUtil/database/EEmcDbItem.h"
#include "StEEmcUtil/database/StEEmcDb.h"
#include "StEEmcUtil/EEmcGeom/EEmcGeomDefs.h"
#include "StEEmcUtil/EEmcGeom/EEmcGeomSimple.h"

#include "St_DataSet.h"
#include "St_DataSetIter.h"
#include "StTriggerUtilities/StTriggerSimuMaker.h"
// include pythia event
#include "StSpinPool/StJetSkimEvent/StPythiaEvent.h"

ClassImp(StEmbeddingTreeMaker)

    StEmbeddingTreeMaker::StEmbeddingTreeMaker(const char *name, const char *outfilename, StMuDstMaker *muDstMakerName, TChain *pythiaChain) : StMaker(name), mMuDstMaker(muDstMakerName), mBemcTables(0)
{
  assert(mMuDstMaker);
  rootFilename = TString(outfilename);
  mPythiaChain = pythiaChain;
  mEemcDb = 0;
  mFile = 0;
  muDst = 0;
  muEvent = 0;
  mVertex = 0;
  mTrack = 0;
  mcTrack = 0;
  mcVertex = 0;
  mEEmcGeom = new EEmcGeomSimple();
  if (StEmcADCtoEMaker *adc2e = (StEmcADCtoEMaker *)StMaker::GetChain()->GetMakerInheritsFrom("StEmcADCtoEMaker"))
    mBemcTables = adc2e->getBemcData()->getTables();
}

StEmbeddingTreeMaker::~StEmbeddingTreeMaker()
{
  if (mMuDstMaker)
    delete mMuDstMaker;
  if (muEvent)
    delete muEvent;
  if (muDst)
    delete muDst;
  if (mEemcDb)
    delete mEemcDb;
  if (mEEmcGeom)
    delete mEEmcGeom;
  if (mBemcTables)
    delete mBemcTables;
  if (mVertex)
    delete mVertex;
  if (mTrack)
    delete mTrack;
  if (mcVertex)
    delete mcVertex;
  if (mcTrack)
    delete mcTrack;
  if (mPythiaEvent)
    delete mPythiaEvent;
  if (mPythiaChain)
    delete mPythiaChain;
  if (mFile)
    delete mFile;
}

Int_t StEmbeddingTreeMaker::Init()
{
  mEemcDb = (StEEmcDb *)this->GetDataSet("StEEmcDb");
  assert(mEemcDb);
  mTriggerSimu = dynamic_cast<StTriggerSimuMaker *>(GetMakerInheritsFrom("StTriggerSimuMaker"));
  pi = TMath::Pi();

  mFile = new TFile(rootFilename, "RECREATE");
  ftree = new TTree("ftree", "dihadron ");

  mPythiaChain->SetBranchAddress("PythiaBranch", &mPythiaEvent);

  // tree branches
  ftree->Branch("fevtNum", &fevtNum, "fevtNum/I");
  ftree->Branch("ffillNum", &ffillNum, "ffillNum/I");
  ftree->Branch("frunNum", &frunNum, "frunNum/I");
  ftree->Branch("frefmult", &frefmult, "frefmult/I");
  // ftree->Branch("ftrigger", &ftrigger, "ftrigger/I");
  ftree->Branch("ftrigger", &ftrigger);
  ftree->Branch("trigJP0", &trigJP0, "trigJP0/I");
  ftree->Branch("trigJP1", &trigJP1, "trigJP1/I");
  ftree->Branch("trigJP2", &trigJP2, "trigJP2/I");

  ftree->Branch("fmaxpar1", &fmaxpar1, "fmaxpar1/I");
  ftree->Branch("fmaxpar", &fmaxpar, "fmaxpar/I");

  ftree->Branch("fVZ", &fVZ, "fVZ/D");
  ftree->Branch("fVX", &fVX, "fVX/D");
  ftree->Branch("fVY", &fVY, "fVY/D");
  ftree->Branch("fVZ_mc", &fVZ_mc, "fVZ_mc/D");
  ftree->Branch("fVZ_pyth", &fVZ_pyth, "fVZ_pyth/D");
  ftree->Branch("fVY_pyth", &fVY_pyth, "fVY_pyth/D");
  ftree->Branch("fVX_pyth", &fVX_pyth, "fVX_pyth/D");

  ftree->Branch("fverRank", &fverRank, "fverRank/D");

  ftree->Branch("partonicPtBin", &partonicPtBin, "partonicPtBin/D");
  // partons before collision
  ftree->Branch("bPartX1", &bPartX1, "bPartX1/D");
  ftree->Branch("bPartId", &bPartId, "bPartId/I");
  ftree->Branch("bPartE", &bPartE, "bPartE/D");
  ftree->Branch("bPartP", &bPartP, "bPartP/D");
  ftree->Branch("bPartPhi", &bPartPhi, "bPartPhi/D");
  ftree->Branch("bPartEta", &bPartEta, "bPartEta/D");
  ftree->Branch("bPartPt", &bPartPt, "bPartPt/D");

  ftree->Branch("yPartX2", &yPartX2, "yPartX2/D");
  ftree->Branch("yPartId", &yPartId, "yPartId/I");
  ftree->Branch("yPartE", &yPartE, "yPartE/D");
  ftree->Branch("yPartP", &yPartP, "yPartP/D");
  ftree->Branch("yPartPhi", &yPartPhi, "yPartPhi/D");
  ftree->Branch("yPartEta", &yPartEta, "yPartEta/D");
  ftree->Branch("yPartPt", &yPartPt, "yPartPt/D");
  // parton after collision
  ftree->Branch("fPart1Id", &fPart1Id, "fPart1Id/I");
  ftree->Branch("fPart1Eta", &fPart1Eta, "fPart1Eta/D");
  ftree->Branch("fPart1Phi", &fPart1Phi, "fPart1Phi/D");
  ftree->Branch("fPart1Pt", &fPart1Pt, "fPart1Pt/D");
  ftree->Branch("fPart1X", &fPart1X, "fPart1X/D");
  ftree->Branch("fPart1E", &fPart1E, "fPart1E/D");

  ftree->Branch("fPart2Id", &fPart2Id, "fPart2Id/I");
  ftree->Branch("fPart2Eta", &fPart2Eta, "fPart2Eta/D");
  ftree->Branch("fPart2Phi", &fPart2Phi, "fPart2Phi/D");
  ftree->Branch("fPart2Pt", &fPart2Pt, "fPart2Pt/D");
  ftree->Branch("fPart2X", &fPart2X, "fPart2X/D");
  ftree->Branch("fPart2E", &fPart2E, "fPart2E/D");

  ftree->Branch("idTruth", idTruth, "idTruth[fmaxpar]/I");
  ftree->Branch("fpT", fpT, "fpT[fmaxpar]/D");
  ftree->Branch("fp", fp, "fp[fmaxpar]/D");
  ftree->Branch("feta", feta, "feta[fmaxpar]/D");
  ftree->Branch("fphi", fphi, "fphi[fmaxpar]/D");
  ftree->Branch("fcharge", fcharge, "fcharge[fmaxpar]/S");

  ftree->Branch("fpT_pyth", fpT_pyth, "fpT_pyth[fmaxpar1]/D");
  ftree->Branch("fp_pyth", fp_pyth, "fp_pyth[fmaxpar1]/D");
  ftree->Branch("feta_pyth", feta_pyth, "feta_pyth[fmaxpar1]/D");
  ftree->Branch("fphi_pyth", fphi_pyth, "fphi_pyth[fmaxpar1]/D");
  ftree->Branch("fpId_pyth", fpId_pyth, "fpId_pyth[fmaxpar1]/D");
  ftree->Branch("fnSigmaPion", fnSigmaPion, "fnSigmaPion[fmaxpar]/D");
  ftree->Branch("fnSigmaKaon", fnSigmaKaon, "fnSigmaKaon[fmaxpar]/D");
  ftree->Branch("fnSigmaProton", fnSigmaProton, "fnSigmaProton[fmaxpar]/D");
  ftree->Branch("fnSigmaElectron", fnSigmaElectron, "fnSigmaElectron[fmaxpar]/D");
  ftree->Branch("fdEdx", fdEdx, "fdEdx[fmaxpar]/D");
  ftree->Branch("fdca", fdca, "fdca[fmaxpar]/D");
  ftree->Branch("ffitPts", ffitPts, "ffitPts[fmaxpar]/s");
  ftree->Branch("ffitPtsPoss", ffitPtsPoss, "ffitPtsPoss[fmaxpar]/s");
  ftree->Branch("fhitsdedx", fhitsdedx, "fhitsdedx[fmaxpar]/s");
  ftree->Branch("fPathLengthToF", fPathLengthToF, "fPathLengthToF[fmaxpar]/D");
  ftree->Branch("fBetaToF", fBetaToF, "fBetaToF[fmaxpar]/D");
  ftree->Branch("fevTime", &fevTime, "fevTime/i");
  ftree->Branch("fvpdVz", &fvpdVz, "fvpdVz/D");

  // matched mc tracks with the reconstructed tracks. unmatched flag == -999
  ftree->Branch("fpT_mc", fpT_mc, "fpT_mc[fmaxpar]/D");
  ftree->Branch("fp_mc", fp_mc, "fp_mc[fmaxpar]/D");
  ftree->Branch("feta_mc", feta_mc, "feta_mc[fmaxpar]/D");
  ftree->Branch("fphi_mc", fphi_mc, "fphi_mc[fmaxpar]/D");
  ftree->Branch("fpId_mc", fpId_mc, "fpId_mc[fmaxpar]/I");
  ftree->Branch("fId_mc", fId_mc, "fId_mc[fmaxpar]/I");

  return StMaker::Init();
}

Int_t StEmbeddingTreeMaker::Make()
{
  muDst = mMuDstMaker->muDst();
  assert(muDst);
  muEvent = muDst->event();
  assert(muEvent);
  assert(mPythiaEvent);

  if ((muEvent->runId() != mPythiaEvent->runId()) || (muEvent->eventNumber() != mPythiaEvent->eventId()))
  {
    cout << muEvent->runId() << " " << mPythiaEvent->runId() << " " << muEvent->eventNumber() << " " << mPythiaEvent->eventId() << "<--c BAD ONE HERE!!!!!!!!" << endl;
    return kStWarn;
  }

  assert((muEvent->runId() == mPythiaEvent->runId()) && (muEvent->eventNumber() == mPythiaEvent->eventId()));

  // mc vertex and tracks selection
  TClonesArray *MuMcVertices = muDst->mcArray(0);
  TClonesArray *MuMcTracks = muDst->mcArray(1);
  Int_t NoMuMcVertices = MuMcVertices->GetEntriesFast();
  Int_t NoMuMcTracks = MuMcTracks->GetEntriesFast();

  if (0)
  {
    cout << "\t" << StMuArrays::mcArrayTypes[0] << " " << NoMuMcVertices << endl;
    cout << "\t" << StMuArrays::mcArrayTypes[1] << " " << NoMuMcTracks << endl;
  }
  if (!NoMuMcVertices || !NoMuMcTracks)
  {
    cout << "This event has no MC information ==> skip it" << endl;
    return kStOK;
  }

  //  trigger info
  trigJP0 = 0, trigJP1 = 0, trigJP2 = 0;

  ftrigger = mTriggerSimu->triggerIds();

  if (mTriggerSimu->isTrigger(370601))
    trigJP0 = 1;
  if (mTriggerSimu->isTrigger(370611))
    trigJP1 = 1;
  if (mTriggerSimu->isTrigger(370621))
    trigJP2 = 1;

  // MC information
  partonicPtBin = mPythiaEvent->pt();

  const TParticle *line5particle = mPythiaEvent->particle(4); // parton from blue beam
  const TParticle *line6particle = mPythiaEvent->particle(5); // parton from Yellow beam
  const TParticle *line7particle = mPythiaEvent->particle(6); // parton from line7
  const TParticle *line8particle = mPythiaEvent->particle(7); // parton from line8

  bPartX1 = mPythiaEvent->x1();
  bPartId = line5particle->GetPdgCode();
  bPartE = line5particle->Energy();
  bPartP = line5particle->P();
  bPartPhi = line5particle->Phi();
  if (bPartPhi > pi)
    bPartPhi -= 2 * pi;
  bPartEta = line5particle->Eta();
  bPartPt = line5particle->Pt();

  yPartX2 = mPythiaEvent->x2();
  yPartId = line6particle->GetPdgCode();
  yPartE = line6particle->Energy();
  yPartP = line6particle->P();
  yPartPhi = line6particle->Phi();
  if (yPartPhi > pi)
    yPartPhi -= 2 * pi;
  yPartEta = line6particle->Eta();
  yPartPt = line6particle->Pt();

  fPart1X = mPythiaEvent->x1();
  fPart1Id = line7particle->GetPdgCode();
  fPart1E = line7particle->Energy();
  fPart1P = line7particle->P();
  fPart1Phi = line7particle->Phi();
  if (fPart1Phi > pi)
    fPart1Phi -= 2 * pi;
  fPart1Eta = line7particle->Eta();
  fPart1Pt = line7particle->Pt();

  fPart2X = mPythiaEvent->x2();
  fPart2Id = line8particle->GetPdgCode();
  fPart2E = line8particle->Energy();
  fPart2P = line8particle->P();
  fPart2Phi = line8particle->Phi();
  if (fPart2Phi > pi)
    fPart2Phi -= 2 * pi;
  fPart2Eta = line8particle->Eta();
  fPart2Pt = line8particle->Pt();

  const TVector3 pythVertex = mPythiaEvent->vertex();
  fVZ_pyth = pythVertex.Z();
  fVY_pyth = pythVertex.Y();
  fVX_pyth = pythVertex.X();

  //------------pythia level particle loop----------------------
  int ntrk_pythia = mPythiaEvent->numberOfParticles();
  fTrackCounter1 = 0;
  for (int j = 0; j < ntrk_pythia; j++)
  {
    if (fabs(mPythiaEvent->particle(j)->Eta()) < 4)
    {
      statusCode_pyth[fTrackCounter1] = mPythiaEvent->particle(j)->GetStatusCode();
      fpT_pyth[fTrackCounter1] = mPythiaEvent->particle(j)->Pt();   // track->p().perp();
      fp_pyth[fTrackCounter1] = mPythiaEvent->particle(j)->P();     // track->p().mag();
      feta_pyth[fTrackCounter1] = mPythiaEvent->particle(j)->Eta(); // track->eta();
      fphi_pyth[fTrackCounter1] = mPythiaEvent->particle(j)->Phi(); // track->phi();
      if (fphi_pyth[fTrackCounter1] > pi)
        fphi_pyth[fTrackCounter1] -= 2 * pi;
      fpId_pyth[fTrackCounter1] = mPythiaEvent->particle(j)->GetPdgCode(); // Pdg code
      fTrackCounter1++;
    }
  } // pythia track loop ended;
  fmaxpar1 = fTrackCounter1;
  //---------------------
  // Fill and Run info
  frefmult = muEvent->refMult();
  StRunInfo &runInfo = muEvent->runInfo();
  // int blue = 1; // for fill info
  // int yellow = 0;

  ffillNum = (Int_t)(runInfo.beamFillNumber(blue));
  frunNum = muEvent->runNumber();
  fevtNum = muEvent->eventNumber();
  // cout << "Fill Number: " << ffillNum << endl;
  StEventInfo &evInfo = muEvent->eventInfo();
  fevTime = evInfo.time();
  /*************** TPC Track Information Begin ***************/
  nPrimVerts = muDst->numberOfPrimaryVertices();
  for (Int_t vert = 0; vert < nPrimVerts; ++vert)
  { // Primary Vertex Loop
    nPrimTracks = 0;

    mVertex = muDst->primaryVertex(vert);
    assert(mVertex);
    muDst->setVertexIndex(vert);

    if (vert > 0)
      continue; // Zeroth vertex is highest rank, or "best" for analysis
    if (!useVertex(mVertex))
      continue;

    fverRank = mVertex->ranking();
    fVZ = mVertex->position().z();
    fVY = mVertex->position().y();
    fVX = mVertex->position().x();
    fvpdVz = muDst->btofHeader()->vpdVz();

    nPrimTracks = muDst->numberOfPrimaryTracks();
    fTrackCounter = 0;
    for (Int_t pTrack = 0; pTrack < nPrimTracks; pTrack++)
    {
      mTrack = muDst->primaryTracks(pTrack);
      assert(mTrack);

      if (!useTrack(mTrack))
        continue;

      idTruth[fTrackCounter] = mTrack->idTruth();
      fpT[fTrackCounter] = mTrack->p().perp();
      fp[fTrackCounter] = mTrack->p().mag();
      feta[fTrackCounter] = mTrack->eta();
      fphi[fTrackCounter] = mTrack->phi();
      if (fphi[fTrackCounter] > pi)
        fphi[fTrackCounter] -= 2 * pi;
      fnSigmaPion[fTrackCounter] = mTrack->nSigmaPion();
      fnSigmaKaon[fTrackCounter] = mTrack->nSigmaKaon();
      fnSigmaProton[fTrackCounter] = mTrack->nSigmaProton();
      fnSigmaElectron[fTrackCounter] = mTrack->nSigmaElectron();
      fcharge[fTrackCounter] = mTrack->charge();
      fdEdx[fTrackCounter] = (mTrack->dEdx() * pow(10, 6));
      fdca[fTrackCounter] = mTrack->dcaGlobal().mag();
      ffitPts[fTrackCounter] = mTrack->nHitsFit();
      fhitsdedx[fTrackCounter] = mTrack->nHitsDedx();
      ffitPtsPoss[fTrackCounter] = mTrack->nHitsPoss();
      fPathLengthToF[fTrackCounter] = mTrack->btofPidTraits().pathLength();
      fBetaToF[fTrackCounter] = mTrack->btofPidTraits().beta();

      // Initialize flag for the mc tracks when the match is not found
      fId_mc[fTrackCounter] = -999;
      fpT_mc[fTrackCounter] = -999;
      fp_mc[fTrackCounter] = -999;
      feta_mc[fTrackCounter] = -999;
      fphi_mc[fTrackCounter] = -999;
      fpId_mc[fTrackCounter] = -999;

      // pythia - reconstructed tracks matching
      if (mTrack->idTruth() > 0)
      {
        mcTrack = (StMuMcTrack *)MuMcTracks->UncheckedAt(mTrack->idTruth() - 1);
        if (mcTrack && (mcTrack->Id() == mTrack->idTruth()))
        {
          fId_mc[fTrackCounter] = mcTrack->Id();
          fpT_mc[fTrackCounter] = mcTrack->Pxyz().perp();
          fp_mc[fTrackCounter] = mcTrack->Pxyz().mag();
          feta_mc[fTrackCounter] = mcTrack->Pxyz().pseudoRapidity();
          fphi_mc[fTrackCounter] = mcTrack->Pxyz().phi();
          if (fphi_mc[fTrackCounter] > pi)
            fphi_mc[fTrackCounter] -= 2 * pi;
          fpId_mc[fTrackCounter] = mcTrack->GePid(); // Pdg code
        }
      }

      fTrackCounter++;
    }
    // Fill track multiplicity
    fmaxpar = fTrackCounter;
  }
  if (fmaxpar > 1 && fmaxpar1 > 1)
    ftree->Fill();

  /**************** TPC Track Information End *****************/

  return kStOK;
} // maker

Bool_t StEmbeddingTreeMaker::useVertex(StMuPrimaryVertex *vert)
{
  if (vert)
  {
    if (vert->ranking() > 1e6)
    {
      if (fabs(vert->position().z()) < 100.)
      {
        return true;
      }
    }
  }
  return false;
}

Bool_t StEmbeddingTreeMaker::useTrack(StMuTrack *track)
{
  if (track->flag() >= 0)
  {
    if (track->nHitsFit() > 15)
    {
      if (track->pt() > 0.2 && fabs(track->eta()) < 1.0 && track->dcaGlobal().mag() < 2.0)
      {
        return true;
      }
    }
  }
  return false;
}

Int_t StEmbeddingTreeMaker::Finish()
{
  ftree->Write();
  mFile->Write();
  mFile->Close();
  delete mFile;

  return kStOk;
}
