void RunEmbeddingTreeMaker(Int_t nEvents =100000,
                           const char *mudstfile = "/star/embed/embedding/pp200_production_2012/v2/Pythia6_pt11_15_100_20212001/P12id.SL12d/2012/047/13047003/st_zerobias_adc_13047003_raw_0570001_r0.MuDst.root",
                           const char *pythiafile = "/star/embed/embedding/pp200_production_2012/v2/Pythia6_pt11_15_100_20212001/P12id.SL12d/2012/047/13047003/st_zerobias_adc_13047003_raw_0570001_r0.pythia.root",
                           const char *outFile = "embedTestTree.root")
// void RunEmbeddingTreeMaker(Int_t nEvents = 2e2,
//  const char *mudstfile = "/star/embed/embedding/pp200_production_2012/pt2_3_100_20153801/P12id.SL12d/2012/047/13047003/pt2_3_13047003_1.MuDst.root",
//  const char *pythiafile = "/star/embed/embedding/pp200_production_2012/pt2_3_100_20153801/P12id.SL12d/2012/047/13047003/pt2_3_13047003_1.pythia.root",
//  const char *outFile = "embedTestTree.root")
{
  // Load libraries
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");
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
  gSystem->Load("StEmbeddingTreeMaker");
  // cout << "Loading libraries done....!" << endl;

  StChain *chain = new StChain("StChain");
  // cout << "Stchain loaded..!" << endl;
  TChain *pythiaChain = new TChain("PythiaTree");
  pythiaChain->Add(pythiafile);
  // cout << "Pythia chain loaded..!" << endl;

  StMuDstMaker *muDstMaker = new StMuDstMaker(0, 0, "", mudstfile, "", 1e6, "MuDst");
  // muDstMaker->SetStatus("*", 0);               // Turn off all branches
  // muDstMaker->SetStatus("MuEvent", 1);         // Turn on the Event data (esp. Event number)
  // muDstMaker->SetStatus("PrimaryVertices", 1); // Turn on the primary track data
  // muDstMaker->SetStatus("PrimaryTracks", 1);   // Turn on the primary track data
  // muDstMaker->SetStatus("MCAll", 1);           // Turn on the McVertex/McTrack data
  // muDstMaker->SetStatus("BTofHeader", 1);      // Turn on the btof data
  // muDstMaker->SetDebug(0);
  //  StMuDbReader *muDstDb = StMuDbReader::instance();

  St_db_Maker *mDbMk = new St_db_Maker("StarDb", "MySQL:StarDb"); // connect to the STAR database
  // mDbMk->SetMaxEntryTime(20160423,000000);
  // cout << "St_db_Maker initialized" << endl;

  StEEmcDbMaker *eemcDb = new StEEmcDbMaker;
  // cout << "St_EEmcdb_Maker initialized" << endl;

  StEmcADCtoEMaker *adc = new StEmcADCtoEMaker;
  adc->saveAllStEvent(true);
  // cout << "St_EEmcADCtoE_Maker initialized" << endl;

  StTriggerSimuMaker *simuTrig = new StTriggerSimuMaker;
  simuTrig->useOfflineDB();
  simuTrig->setMC(2);
  simuTrig->useBemc();
  simuTrig->useEemc();
  simuTrig->bemc->setConfig(StBemcTriggerSimu::kOffline);
  // cout << "StTriggerSimu_Maker initialized" << endl;

  StEmbeddingTreeMaker *eventQA = new StEmbeddingTreeMaker("StEmbeddingTreeMaker", outFile, muDstMaker, pythiaChain);
  // cout << "Tree maker called..!" << endl;
  //  Initialize chain
  chain->Init();
  // chain->PrintInfo();

  for (int iEnt = 0; iEnt < nEvents; ++iEnt)
  {
    if (pythiaChain->GetEvent(iEnt) <= 0)
      break;

    if (iEnt % 100 == 0)
      cout << "Working on event number: " << iEnt << endl;
    chain->Clear();
    int iRet = chain->Make(iEnt);
    if (iRet)
    {
      // cout << " BAD return code!" << endl;
      break;
    }
  } // End event loop
  chain->Finish();
}
