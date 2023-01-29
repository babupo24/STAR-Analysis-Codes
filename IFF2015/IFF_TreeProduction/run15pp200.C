#include <iostream>
#include <vector>
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TInterpreter.h"
#include "TSystem.h"
//#include "StRoot/StSpinPool/StSpinDbMaker/StSpinDbMaker.h"

class StMuDstMaker;
class StMuDebug;
class StSpinDbMaker;
class StChain;
class StMaker;

//StMaker *spinDbMaker;
StMuDstMaker* maker;
//StSpinDbMaker* spDbMaker = dynamic_cast<StSpinDbMaker*>(spinDbMaker);


//void run15pp200(const char* list, const char* oFile, const char* jobid)
void run15pp200(const char* list, const char* oFile) //interactive mode
{
	//char *list = "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_pp200trans_2015/ReversedFullField/P16id/2015/087/16087052/st_physics_16087052_raw_1000003.MuDst.root";
	//char *list = "root://xrdstar.rcf.bnl.gov:1095//home/starlib/home/starreco/reco/production_pp200trans_2015/ReversedFullField/P16id/2015/065/16065023/st_physics_16065023_raw_1000003.MuDst.root";
	//gInterpreter->GenerateDictionary("vector<int>","vector");
	gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
	loadSharedLibraries();
	gSystem->Load("StSpinDbMaker");
	gSystem->Load("StDbBroker");
	gSystem->Load("St_db_Maker");

	//Add libraries to impliment Start-less TOF (Bassam presentation)10/21/2020 
	gSystem->Load("StBTofUtil");
	gSystem->Load("StVpdCalibMaker");
	gSystem->Load("StBTofCalibMaker");
	///-----------------------------------------


	// some definitions
	// Float_t->F, Int_t->I,Double_t->D, UShort_t ->s,Short_t ->S
	// Define in /StRoot/StMuDSTMaker/COMMON/StMuTrack.h
	const int fMaxHit = 4000;
	const int ntrig = 23;
	unsigned int frunTime; 
	Int_t ffillNum;
	Int_t frunId;
	unsigned int fevTime;
	float fvpdVz;
	Int_t fevtNum;
	Int_t frunNum;
	Int_t fspinBits;
	Int_t fspinconfig;
	Int_t frefmult;
	vector<double> ftrigger;
	Int_t fmaxpar;
	Double_t fVZ;
	Double_t fverRank;
	Int_t fTrackCounter;
	Double_t fpT[fMaxHit];
	Double_t fp[fMaxHit];
	Double_t feta[fMaxHit];
	Double_t fphi[fMaxHit];
	Double_t fnSigmaPion[fMaxHit];
	Double_t fnSigmaKaon[fMaxHit];
	Double_t fnSigmaProton[fMaxHit];
	Double_t fnSigmaElectron[fMaxHit];  
	Double_t fnSigmaPionFit[fMaxHit];
	Double_t fnSigmaKaonFit[fMaxHit];
	Double_t fnSigmaProtonFit[fMaxHit];
	Double_t fnSigmaElectronFit[fMaxHit];  
	Short_t fcharge[fMaxHit];
	Double_t fdEdx[fMaxHit];
	Double_t  fdca[fMaxHit];
	UShort_t ffitPts[fMaxHit];
	UShort_t ffitPtsPoss[fMaxHit];
	UShort_t fhitsdedx[fMaxHit];
	Double_t fPathLengthToF[fMaxHit];
	Double_t fBetaToF[fMaxHit];

	TFile *fout  = new TFile(oFile,"recreate");  
	TTree *ftree = new TTree("ftree","dihadron ");
	ftree->Branch("fevtNum",&fevtNum,"fevtNum/I");
	ftree->Branch("ffillNum",&ffillNum,"ffillNum/I");
	ftree->Branch("frunId",&frunId,"frunId/I");
	ftree->Branch("frunNum",&frunNum,"frunNum/I");
	ftree->Branch("fspinBits",&fspinBits,"fspinBits/I");
	ftree->Branch("fspinconfig",&fspinconfig,"fspinconfig/I");
	ftree->Branch("frefmult",&frefmult,"frefmult/I");
	ftree->Branch("fmaxpar",&fmaxpar,"fmaxpar/I");
	ftree->Branch("fVZ",&fVZ,"fVZ/D");
	ftree->Branch("fverRank",&fverRank,"fverRank/D");
	ftree->Branch("fverRank",&fverRank,"fverRank/D");
	ftree->Branch("ftrigger",&ftrigger);
	ftree->Branch("fpT",fpT,"fpT[fmaxpar]/D");
	ftree->Branch("fp",fp,"fp[fmaxpar]/D");
	ftree->Branch("feta",feta,"feta[fmaxpar]/D");
	ftree->Branch("fphi",fphi,"fphi[fmaxpar]/D");
	ftree->Branch("fcharge",fcharge,"fcharge[fmaxpar]/S");
	ftree->Branch("fnSigmaPion",fnSigmaPion,"fnSigmaPion[fmaxpar]/D");
	ftree->Branch("fnSigmaKaon",fnSigmaKaon,"fnSigmaKaon[fmaxpar]/D");
	ftree->Branch("fnSigmaProton",fnSigmaProton,"fnSigmaProton[fmaxpar]/D");
	ftree->Branch("fnSigmaElectron",fnSigmaElectron,"fnSigmaElectron[fmaxpar]/D");
	ftree->Branch("fdEdx",fdEdx,"fdEdx[fmaxpar]/D");
	ftree->Branch("fdca",fdca,"fdca[fmaxpar]/D");
	ftree->Branch("ffitPts",ffitPts,"ffitPts[fmaxpar]/s");
	ftree->Branch("ffitPtsPoss",ffitPtsPoss,"ffitPtsPoss[fmaxpar]/s"); 
	ftree->Branch("fhitsdedx",fhitsdedx,"fhitsdedx[fmaxpar]/s"); 
	ftree->Branch("fPathLengthToF",fPathLengthToF,"fPathLengthToF[fmaxpar]/D");
	ftree->Branch("fBetaToF",fBetaToF,"fBetaToF[fmaxpar]/D");
	ftree->Branch("frunTime", &frunTime, "frunTime/i");
	ftree->Branch("fevTime", &fevTime, "fevTime/i");
	ftree->Branch("fvpdVz",&fvpdVz, "fvpdVz/F");

	//trigger ids
	int trigid[23];
	trigid[0] = 480003; /*BBCMB*/        trigid[1] = 480004;  /*VPDMB-novtx*/         trigid[2] = 480005;  /*ZDCMB-trgonly*/
	trigid[3] = 480904; /*VPDMB-30*/     trigid[4] = 480007; /*VPDMB-5-trgonly*/      trigid[5] = 480201;  /*BHT0*VPDMB-5*/
	trigid[6] = 480202; /*BHT1*VPDMB-30*/trigid[7] = 480203;  /*BHT0*BBCMB*/          trigid[8] = 480204;  /*BHT1*BBCMB*/
	trigid[9] = 480205; /*BHT2*BBCMB*/   trigid[10]= 480206;  /*BHT1*VPDMB-30-nobsmd*/trigid[11]= 480301;  /*EHT0*/
	trigid[12]= 480401; /*JP2 */         trigid[13]= 480402;   /*JP2-bsmd*/           trigid[14]= 480403;  /*AJP*/
	trigid[15]= 480404; /*JP1*/          trigid[16]= 480405;  /*JP2*L2JetHigh*/       trigid[17]= 480406;   /*BHT2*BJP1*L2Bgamma*/
	trigid[18]= 480714; /*RP_CPEI*/      trigid[19]= 480411;  /*JP2*/                 trigid[20]= 480414;  /*JP1*/
	trigid[21]= 480415; /*JP2*L2JetHigh*/trigid[22]= 480501;  /*EHT0*EJP1*L2Egamma*/
	// int trigid[4];
	//trigid[0] = 480003; /*BBCMB*/
	//trigid[1] = 480004;  /*VPDMB-novtx*/
	//trigid[2] = 480005;  /*ZDCMB-trgonly*/
	//trigid[3] = 480904; /*VPDMB-30*/
	//trigid[4] = 480007; /*VPDMB-5-trgonly*/
	//trigid[0] = 480401; /*JP2 */
	//trigid[1] = 480404; /*JP1*/
	//trigid[2] = 480411;  /*JP2*/
	//trigid[3] = 480414;  /*JP1*/

	chain = new StChain("StChain");
	chain->SetDebug(0);

	StMuDebug::setLevel(0);  // switch of some debug output
	char theFilter[80];
	sprintf(theFilter,".MuDst.root:MuDst.root");
	//maker = new StMuDstMaker(0,0,"", list, theFilter, 1e6, "MuDst");// set up maker in read mode(Bassam's approach)
	maker = new StMuDstMaker(0,0,"", list, theFilter, 1000, "MuDst");// set up maker in read mode(original approach)

	St_db_Maker *stDb = new St_db_Maker("StarDb", "MySQL:StarDb");
	StSpinDbMaker* spindb = new StSpinDbMaker;
	//Added from Bassam's presentation-------------
	StVpdCalibMaker *vpdCalib = new StVpdCalibMaker();
 	vpdCalib->setMuDstIn();
	StBTofCalibMaker *btofCalib = new StBTofCalibMaker();
 	btofCalib->setMuDstIn();
 	vpdCalib->setUseVpdStart(kFALSE);
         //----------------------------------------------
	//assert(spindb->isValid());  // all 3 DB records exist , enough to call once in InitRun()
	//assert(spindb->isPolDirTrans());  // you do not want mix Long & Trans by accident, -//-


	Int_t nevents=maker->chain()->GetEntries();
	//cout<<nevents<< " events in chain" << endl;

	int blue = 1;
	int yellow = 0;

	//open database connection for spin configuration
	const char* database = "mysql://db04.star.bnl.gov:3414/RunLog?timeout=60";
	const char* user = "pokhrel";
	const char* pass = "";
	TMySQLServer* mysql = TMySQLServer::Connect(database,user,pass);

	chain->Init();
	while(!chain->Make())
	{
		StMuDst *mu   = maker->muDst();
		StMuEvent *ev = maker->muDst()->event();
		if(!ev)continue;
 		Int_t nevents=maker->chain()->GetEntries();
		int nVertices = mu->numberOfPrimaryVertices();
		//include vertex ranking
		for (int ver =0; ver< nVertices; ver++ ){ 
			StMuPrimaryVertex* vertex = mu->primaryVertex(ver);
			assert(vertex);
			mu->setVertexIndex(ver);
			fverRank = vertex->ranking();
			//cout << "ver rank " << fverRank << endl;
			if(fverRank <1e6) continue;
			//if(ver > 0) continue;
			//cout << "ver rank " << fverRank << endl;
			StZdcTriggerDetector zdc = ev->zdcTriggerDetector();
			StThreeVectorF vertexPos = ev->primaryVertexPosition();
			float fvpdVz = mu->btofHeader()->vpdVz();

			//if(ftrigger==0)continue;
	
			//cout << ftrigger << endl;
			int bx48 =  ev->l0Trigger().bunchCrossingId();
			int fspinconfig = spindb->spin8usingBX48(bx48);

			if(fspinconfig!=51 && fspinconfig!=53 && fspinconfig!=83 && fspinconfig!=85)continue;
			//find trigger ids
			for(int trig=0; trig<23; trig++)
			{if ((ev->triggerIdCollection().nominal().isTrigger(trigid[trig]))) ftrigger.push_back(trigid[trig]);
			}
			//for(int trg=0; trg<ftrigger.size();trg++){
			//cout << "Event: "<< ver<< " fired trigger: "<< ftrigger.at(trg) << endl;	
			//}

			int refMult = ev->refMult();
			int frunId = ev->runId(); 

			StRunInfo& runInfo = ev->runInfo();
			unsigned int  frunTime = runInfo.productionTime();
			//cout << "event time stamp  " << fevTime<< endl;

			StEventInfo& evInfo = ev->eventInfo();
			unsigned int fevTime = evInfo.time();
			//cout << "event time " << fevTime <<endl;

		        int fill = (int) (runInfo.beamFillNumber(blue));
			int ffillNum = fill;
			//cout << "fillnum " << ffillNum << endl;
			int frunNum = ev->runNumber();
			//cout << "run number : " << frunNum << endl;
			int fevtNum =  ev->eventNumber();
			int fspinBits = ev->l0Trigger().spinBits(frunId);
			//cout << "SpinConfig    and     SpinBits  " << endl;
			//cout << fspinconfig << "  \t    "<< fspinBits << endl;
			frefmult = refMult;
			fVZ= vertexPos.z();
			int ntrk = mu->primaryTracks()->GetEntries();
			fTrackCounter=0;
			for (int j = 0; j < ntrk; j++)
			{
				if(mu->primaryTracks(j)->flag() < 0) continue;
				if (fabs(mu->primaryTracks(j)->eta())<1. && mu->primaryTracks(j)->nHitsFit()>12 && fabs(mu->primaryTracks(j)->dcaGlobal().mag())< 1.) 
				{
				//cout << mu->primaryTracks(j)->eta() << ","<<mu->primaryTracks(j)->nHitsFit()<<", "<< mu->primaryTracks(j)->dcaGlobal().mag() << endl;
					fpT[fTrackCounter]         =  mu->primaryTracks(j)->p().perp();
					fp[fTrackCounter]          =  mu->primaryTracks(j)->p().mag();
					feta[fTrackCounter]        =  mu->primaryTracks(j)->eta();
					fphi[fTrackCounter]        =  mu->primaryTracks(j)->phi();
					fnSigmaPion[fTrackCounter] =  mu->primaryTracks(j)->nSigmaPion();
					fnSigmaKaon[fTrackCounter] =  mu->primaryTracks(j)->nSigmaKaon();
					fnSigmaProton[fTrackCounter] =  mu->primaryTracks(j)->nSigmaProton();
					fnSigmaElectron[fTrackCounter] =  mu->primaryTracks(j)->nSigmaElectron();
					fnSigmaPionFit[fTrackCounter] =  mu->primaryTracks(j)->nSigmaPionFit();
					fnSigmaKaonFit[fTrackCounter] =  mu->primaryTracks(j)->nSigmaKaonFit();
					fnSigmaProtonFit[fTrackCounter] =  mu->primaryTracks(j)->nSigmaProtonFit();
					fnSigmaElectronFit[fTrackCounter] =  mu->primaryTracks(j)->nSigmaElectronFit();
					fcharge[fTrackCounter]     =  mu->primaryTracks(j)->charge();
					fdEdx[fTrackCounter]       =  mu->primaryTracks(j)->dEdx();
					fdca[fTrackCounter]       =   mu->primaryTracks(j)->dcaGlobal().mag();
					ffitPts[fTrackCounter]     =  mu->primaryTracks(j)->nHitsFit();
					fhitsdedx[fTrackCounter]     =  mu->primaryTracks(j)->nHitsDedx();
					ffitPtsPoss[fTrackCounter]     =  mu->primaryTracks(j)->nHitsPoss();
					fPathLengthToF[fTrackCounter] = mu->primaryTracks(j)->btofPidTraits().pathLength();
					fBetaToF[fTrackCounter]    =  mu->primaryTracks(j)->btofPidTraits().beta();
					fTrackCounter++;
				}
			}// track loop
			fmaxpar = fTrackCounter;
			if(fmaxpar>1)ftree->Fill();
		
		 }// vertex loop
			ftrigger.clear();
	}// event loop

	fout->cd();
	ftree->Write();     
	fout->Close();
	//cout<<" done"<<endl; flush(cout);
	mysql->Close();
}  
