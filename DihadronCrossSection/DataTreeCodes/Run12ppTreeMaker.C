#include <TSystem>
class StMaker;
class StChain;
class StPicoDstMaker;

void Run12ppTreeMaker(Int_t nEvents = 1000,
                      const char *picoFile = "/star/data27/reco/pp200_production_2012/ReversedFullField/P12id.SL21d/2012/041/13041102/st_physics_13041102_raw_4010002.picoDst.root",
                      const char *outFile = "testTree.root")
{
  TDirectory *home = gDirectory;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");
  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");

  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");

  gSystem->Load("StRun12ppTreeMaker");

  StChain *chain = new StChain("StChain");

  StPicoDstMaker *picoMaker = new StPicoDstMaker(2, picoFile, "picoDst"); // set up maker in read mode(Bassam's approach)
  StPicoDstMaker *picoDstDb = StPicoDstMaker::instance();

  St_db_Maker *mDbMk = new St_db_Maker("StarDb", "MySQL:StarDb"); // connect to the STAR database

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
