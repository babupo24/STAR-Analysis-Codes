#include <TSystem>
class StMaker;
class StChain;
class StPicoDstMaker;
class StTriggerFilterMaker;

void Run12ppTreeMaker(Int_t nEvents = 100,
                      const char *picoFile = "/star/data27/reco/pp200_production_2012/ReversedFullField/P12id.SL21d/2012/047/13047042/st_physics_13047042_raw_2010001.picoDst.root",
                      const char *outFile = "testTree.root")
{
  TDirectory *home = gDirectory;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StDaqLib");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StEmcSimulatorMaker");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEEmcDbMaker");
  gSystem->Load("StSpinDbMaker");
  gSystem->Load("StEmcTriggerMaker");
  gSystem->Load("StTriggerUtilities");
  gSystem->Load("StMCAsymMaker");
  gSystem->Load("StRandomSelector");
  gSystem->Load("StJetSkimEvent");

  gSystem->Load("StRun12ppTreeMaker");

  StChain *chain = new StChain("StChain");

  StPicoDstMaker *picoMaker = new StPicoDstMaker(2, picoFile, "picoDst"); // set up maker in read mode(Bassam's approach)
  StPicoDstMaker *picoDstDb = StPicoDstMaker::instance();

  St_db_Maker *mDbMk = new St_db_Maker("StarDb", "MySQL:StarDb"); // connect to the STAR database

  StEEmcDbMaker *eemcDb = new StEEmcDbMaker;

  StEmcADCtoEMaker *adc = new StEmcADCtoEMaker;
  adc->saveAllStEvent(true);

  StTriggerSimuMaker *simuTrig = new StTriggerSimuMaker;
  simuTrig->useOnlineDB();
  simuTrig->setMC(false);
  simuTrig->useBemc();
  simuTrig->useEemc();
  simuTrig->bemc->setConfig(StBemcTriggerSimu::kOnline);

  StRun12ppTreeMaker *picoAna = new StRun12ppTreeMaker("StRun12ppTreeMaker", outFile, picoMaker);
  picoAna->setHome(home);
  // Initialize chain
  chain->Init();
  chain->PrintInfo();
  if (nEvents == 0)
    nEvents = 10000000; // get all events

  for (int iEnt = 0; iEnt < nEvents; ++iEnt)
  {
    if (iEnt % 100 == 0)
      cout << "Working on event number: " << iEnt << endl;
    chain->Clear();
    int iRet = chain->Make(iEnt);
    if (iRet)
    {
      cout << " BAD return code!" << endl;
      break;
    }
  } // End event loop
  chain->Finish();
}
