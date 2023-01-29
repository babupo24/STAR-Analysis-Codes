#include "Iff2015.h"
#include "TH2.h"
#include "TH1.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TVector3.h"
#include "StSpinPool/StSpinDbMaker/StSpinDbMaker.h" 
#include <fstream>
#include <iostream>
#include "TFile.h"
#include <vector>
using namespace std;
ofstream Output;
double pionCut(double sigma_tpc, double sigma_tof);//single pid cut wirks for all bins
double pionCutPt(double sigma_tpc, double sigma_tof,double pt_pair);
double pionCutM(double sigma_tpc, double sigma_tof,double minv);
double pionCutEta(double sigma_tpc, double sigma_tof,double etaPair);
//loadPol function is included in order to compensate polarization decay ,     Babu 11/22/2019
void loadPol(const char* Filename, map<int,unsigned int>& t0Map, map<int,float>& pb0Map,map<int,float>& dpdtbMap, map<int,float>& py0Map, map<int,float>& dpdtyMap);
void Iff2015::Loop()
{
	gROOT->Reset();
	ntuple1tpc   = new TNtuple("ntuple1tpc","PID from TPC","cone:fspinconfig:runNum:fillNum");
	ntuple2tpc   = new TNtuple("ntuple2tpc","PID from TPC","R:ps:Minv:pT_pair:pT_min_pair:eta_pair:psx:psy:psz");
	ntuple3tpc   = new TNtuple("ntuple3tpc","PID from TPC","cosPhiRB:sinPhiRB:cosPhiSB:sinPhiSB:cosPhiRY:sinPhiRY:cosPhiSY:sinPhiSY");
	ntuple4tpc   = new TNtuple("ntuple4tpc","PID from TPC","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple5tpc   = new TNtuple("ntuple5tpc","PID from TPC","fitPts_min_pair:msqr");
	//ntuples with TPC Or TPC+TOF cut (if TOF available). Here, "f" referes to final, as these ntuples will be used for final result 
	ntuple1f   = new TNtuple("ntuple1f","PID from TPC or TOF","cone:fspinconfig:runNum:fillNum");
	ntuple2f   = new TNtuple("ntuple2f","PID from TPC or TOF","R:ps:Minv:pT_pair:pT_min_pair:eta_pair:psx:psy:psz");
	ntuple3f   = new TNtuple("ntuple3f","PID from TPC or TOF","cosPhiRB:sinPhiRB:cosPhiSB:sinPhiSB:cosPhiRY:sinPhiRY:cosPhiSY:sinPhiSY");
	ntuple4f   = new TNtuple("ntuple4f","PID from TPC or TOF","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple5f   = new TNtuple("ntuple5f","PID from TPC or TOF","fitPts_min_pair:msqr");
	//ntuples with TPC && TOF pid cut
	ntuple1tof    = new TNtuple("ntuple1tof","PID from TPC && TOF","cone:fspinconfig:runNum:fillNum");
	ntuple2tof    = new TNtuple("ntuple2tof","PID from TPC && TOF","R:ps:Minv:pT_pair:pT_min_pair:eta_pair:psx:psy:psz");
	ntuple3tof    = new TNtuple("ntuple3tof","PID from TPC && TOF","cosPhiRB:sinPhiRB:cosPhiSB:sinPhiSB:cosPhiRY:sinPhiRY:cosPhiSY:sinPhiSY");
	ntuple4tof    = new TNtuple("ntuple4tof","PID from TPC && TOF","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple5tof    = new TNtuple("ntuple5tof","PID from TPC && TOF","fitPts_min_pair:msqr");
	//Background ntuples with TPC && TOF pid cut
	ntuple1bkg    = new TNtuple("ntuple1bkg","Bkg, PID from TPC && TOF","cone:fspinconfig:runNum:fillNum");
	ntuple2bkg    = new TNtuple("ntuple2bkg","Bkg, PID from TPC && TOF","R:ps:Minv:pT_pair:pT_min_pair:eta_pair:psx:psy:psz");
	ntuple3bkg    = new TNtuple("ntuple3bkg","Bkg, PID from TPC && TOF","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple3bkg    = new TNtuple("ntuple3bkg","Bkg PID from TPC && TOF","cosPhiRB:sinPhiRB:cosPhiSB:sinPhiSB:cosPhiRY:sinPhiRY:cosPhiSY:sinPhiSY");
	ntuple4bkg    = new TNtuple("ntuple4bkg","Bkg PID from TPC && TOF","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple5bkg    = new TNtuple("ntuple5bkg","Bkg PID from TPC && TOF","fitPts_min_pair:msqr");
	//Background ntuples with 
	ntuple1tofbkg    = new TNtuple("ntuple1tofbkg","Bkg, PID from TPC && TOF","cone:fspinconfig:runNum:fillNum");
	ntuple2tofbkg    = new TNtuple("ntuple2tofbkg","Bkg, PID from TPC && TOF","R:ps:Minv:pT_pair:pT_min_pair:eta_pair:psx:psy:psz");
	ntuple3tofbkg    = new TNtuple("ntuple3tofbkg","Bkg, PID from TPC && TOF","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple3tofbkg    = new TNtuple("ntuple3tofbkg","Bkg PID from TPC && TOF","cosPhiRB:sinPhiRB:cosPhiSB:sinPhiSB:cosPhiRY:sinPhiRY:cosPhiSY:sinPhiSY");
	ntuple4tofbkg    = new TNtuple("ntuple4tofbkg","Bkg PID from TPC && TOF","PhiR_cosB:PhiR_sinB:PhiRB:PhiS_cosB:PhiS_sinB:PhiSB:PhiRSB:PhiR_cosY:PhiR_sinY:PhiRY:PhiS_cosY:PhiS_sinY:PhiSY:PhiRSY");
	ntuple5tofbkg    = new TNtuple("ntuple5tofbkg","Bkg PID from TPC && TOF","fitPts_min_pair:msqr");
	//corrected polarization
	ntuple6    = new TNtuple("ntuple6","","polB_corr:polY_corr");//polarization 
	//pT pair bin ranges for Aut vs Minv
 	double pT[6]={2.80,   3.71,  4.303 ,   5.084,   6.404,  15.00};
	//Minv bin ranges for Aut vs pT
	double M[6]={0.20, 0.4795, 0.6714, 0.8452, 1.1100, 4.00};
	//eta bin ranges for aut vs eta
	double eta_range[10]={-1., -0.668, -0.469, -0.281, -0.096, 0.089, 0.275, 0.470, 0.675, 1.};
	//pid control histograms 
	int nCh=2;
	int nBins = 5;
	int nBinsEta = 9;
	const char *charge[2]={"Pos", "Neg"};
	for(int ch =0; ch<2; ch++){
		for (int i =0 ; i<5; i++){
			hTpcVsTofM[ch][i]= new TH2D(Form("hSigmaPionTpcVsTof_%s_Mbin%i",charge[ch],i),"M Bins",100, -5, 10, 100, -5, 10);
			hTpcVsTofpT[ch][i]= new TH2D(Form("hSigmaPionTpcVsTof_%s_pTbin%i",charge[ch],i),"pT Bins",100, -5, 10, 100, -5, 10);
			hTpcVsTofMBkg[ch][i]= new TH2D(Form("hSigmaPionTpcVsTofBkg_%s_Mbin%i",charge[ch],i),"M Bins",100, -5, 10, 100, -5, 10);
			hTpcVsTofpTBkg[ch][i]= new TH2D(Form("hSigmaPionTpcVsTofBkg_%s_pTbin%i",charge[ch],i),"pT Bins",100, -5, 10, 100, -5, 10);
		}
	}


 int trigid[23]={480003,  480004,  480005, 480904,  480007,  480201, 480202,  480203,  480204,  480205,  480206,  480301, 480401,  480402,  480403,480404,  480405,  480406, 480714,  480411,  480414, 480415,  480501};
 TString  trigname[23]={"BBCMB","VPDMB-novtx","ZDCMB-trgonly","VPDMB30","VPDMB-5-trgonly","BHT0*VPDMB-5","BHT1_VPDMB30","BHT0*BBCMB","BHT1*BBCMB","BHT2_BBCMB","BHT1*VPDMB-30-nobsmd","EHT0","JP2_1","JP2-bsmd","AJP","JP1_1","JP2*L2JetHigh_1","BHT2*BJP1*L2Bgamma","RP_CPEI","JP2_2","JP1_2","JP2*L2JetHigh_2","EHT0*EJP1*L2Egamma"};
	double picut1, picut2;
        double picut1pt, picut2pt;	
        double picut1m, picut2m;	
        double picut1eta, picut2eta;	

	double p1x,p2x,p1y,p2y,p1z,p2z,psx,psy,psz,ps,cone,R,Rx,Ry,Rz,R1,Rx1,Ry1,Rz1, p_tr, p_tr2, ptr, ptr2;
	double p1,p2,E1,E2,Minv,pT_pair,pT_min_pair,eta_pair,fitPts_min_pair;
	double cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY;
	double PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY,phiSB,phiRB,phiDiffB,phiRSB;
	int spin51=0, spin53=0, spin83=0, spin85=0, blueUp=0, blueDown=0, yellowUp=0, yellowDown=0; 
	double msqr, msqr1, msqr2;// msqr_tr2;
	double pi=3.14159265359;
	double CONE_CUT_MIN = 0.05;
	double m_pion=0.1396; //GeV
	int fillnum;
	unsigned int evTime;
	//Spin configuration info using spin8(fspinconfig)
	//number  blue  yellow
	//   0    empty empty
	//   3    empty   up
	//   5    empty  down
	//  48     up   empty
	//  51     up     up
	//  53     up    down
	//  80    down  empty
	//  83    down    up
	//  85    down   down
	// 153    unpol unpol

	int paircount =0;
		int ntrack1=0, ntrack2=0;
	if (fChain == 0) return;

	Long64_t nentries = fChain->GetEntries();
	cout << "Entries: " << nentries << endl;
	Long64_t nbytes = 0, nb = 0, trig=0, trig2=0;
	//Event Loop
 for(Long64_t jentry=0; jentry<nentries;jentry++)
 //for (Long64_t jentry=0; jentry<100000;jentry++)
 {       
	 Long64_t ientry = LoadTree(jentry);
	 if (ientry < 0) break;
	 nb = fChain->GetEntry(jentry);   nbytes += nb;
	 //cout << "Working on event : "<< jentry << endl;
	 fillnum = ffillNum-18640;
	 if(fillnum==109 || fillnum==124 ) continue;    //Fills 18749 & 18764 have wrongly assigned fill patterns for Yellow and Blue beam respectively. Don't use in transverse analysis. // Babu 12/22/2019
	 if(fverRank<1e6) continue; //Already applied at tree level
	 if(fabs(fVZ)>60.) continue;
	 //if(fvpdVz==-999) continue;
	 //if(fabs(fVZ-fvpdVz)>6.) continue; // added cut , Carl's suggestion in spin pwg.  Babu 11/22/2019
	 if(fspinconfig!=51 && fspinconfig!=53 && fspinconfig!=83 && fspinconfig!=85)continue;  //applied at tree level 
	 //polarization correction//
	 // I have used event time from which beam fill time(t0) is sbtracted to get dt, the time diffrence,with which polarization correction will be calculated  
	 // The polarization values after corection is tested, which are slightly smaller than polarization values during beam fill, which is expected and resonable.
	 evTime = fevTime;
	 map<int,unsigned int> t0Map;
	 map<int,float> pb0Map;
	 map<int,float> dpdtbMap;
	 map<int,float> py0Map;
	 map<int,float> dpdtyMap;
	 loadPol("/star/u/pokhrel/Run15/Run15_Analysis/IFF_Analysis/iff200_code/AsymWithoutTOF/TransPolzFiles/run15pp200transPol.txt",t0Map,pb0Map,dpdtbMap,py0Map,dpdtyMap);

	 unsigned int   t0    = t0Map[fillnum];
	 float    pb0   = pb0Map[fillnum];
	 float    dpdtb = dpdtbMap[fillnum];
	 float    py0   = py0Map[fillnum];
	 float    dpdty = dpdtyMap[fillnum];
	 //File reads in polarization in percentages, I want absolutes               
	 pb0   = pb0/100.0;
	 dpdtb = dpdtb/100.0;
	 py0   = py0/100.0;
	 dpdty = dpdty/100.0;
	 // set polarization taking into accout decay
	 unsigned int unixTime = evTime;
	 // Set time of event from t0
	 float dt = (unixTime-t0)/3600.0;
	 // calc pol during fill, Using linear decay
	 float    pbD = pb0 + dpdtb*dt;     
	 float    pyD = py0 + dpdty*dt; 

	 //Fill ntuple for corrected polarization values 
	 ntuple6->Fill(pbD, pyD);
	 //Polarization corerction ends!//


	 int jptrig=0;
	 int jp1=0;
	 int jp2=0;
	 for( trig=0; trig<ftrigger->size();trig++){
		 trigger = ftrigger->at(trig);
		 if(trigger==480401 || trigger==480411 || trigger==480404 || trigger==480414) jptrig=1;
		 if(trigger==480401 || trigger==480411) jp2=1;	
		 if(trigger==480404 || trigger==480414) jp1=1;	
	 }
	 if(jptrig!=1)continue;//only JP1 and/or JP2 triggers
	 //cout<<"JP trigger: "<<jptrig<<endl; 
	 //Track Loop
	 for(int tr=0;tr < fmaxpar;tr++)
	 {	ntrack1++;
		 double fitPtsRatio_tr = (double)ffitPts[tr]/(double)ffitPtsPoss[tr];
		 //cout << "nFitPts/nFitPtsPoss tr1: "<< fitPtsRatio_tr << endl;
		 if(fpT[tr]>1.5 && fdca[tr]<1. && ffitPts[tr]>15 && feta[tr]<1. && feta[tr]>-1. && fitPtsRatio_tr>.51) 
		 {
			 //cout << "p="<<fp[tr]<<"\tfnSigmaPion="<<fnSigmaPion[tr]<<"\tdca="<<fdca[tr]<< "\tfitPts="<<ffitPts[tr]<<"\thitsdEdx="<<fhitsdedx[tr]<<"\teta="<<feta[tr]<< "\tfitratio="<<fitPtsRatio_tr<< endl;				
			 double  inverse_pion_beta1 = sqrt(pow(m_pion,2)/pow(fp[tr],2)+1);
			 //apply tof cut only to those tracks that has ToF information. Keep the track if it doesn't have ToF info (flag==-999) 
			 //if(fBetaToF[tr]!=-999 && (abs(1/fBetaToF[tr] - inverse_pion_beta1)>0.03)) continue; 
			 //cout<<"after tof cut: beta: "<<fBetaToF[tr]<< ", invBeta: "<<inverse_pion_beta1<< ", cut value: "<< abs(1/fBetaToF[tr] - inverse_pion_beta1)<<endl;



			 if(fBetaToF[tr]!=-999)msqr1=fp[tr]*fp[tr]*(1/(fBetaToF[tr]*fBetaToF[tr])-1);
			 //Track 2 loop 
			 for(int tr2=tr+1;tr2 < fmaxpar;tr2++)
			 {	
				 ntrack2++;
				 double fitPtsRatio_tr2 = (double)ffitPts[tr2]/(double)ffitPtsPoss[tr2];
				 //cout << "nFitPts/nFitPtsPoss tr2: "<< fitPtsRatio_tr2 << endl;
				 if(fpT[tr2]>1.5 && fdca[tr2]<1. && ffitPts[tr2]>15 && feta[tr2]<1. && feta[tr2]>-1. && fitPtsRatio_tr2 >.51)
				 {
					 double  inverse_pion_beta2 = sqrt(pow(m_pion,2)/pow(fp[tr2],2)+1);
					 //if(fBetaToF[tr2]!=-999 && (abs(1/fBetaToF[tr2] - inverse_pion_beta2)>0.03)) continue;
					 //if(fBetaToF[tr2]!=-999)msqr2=fp[tr2]*fp[tr2]*(1/(fBetaToF[tr2]*fBetaToF[tr2])-1);
					 //cout<<" beta2: "<<fBetaToF[tr2]<< ", invBeta2: "<<inverse_pion_beta2<< ", cut value2: "<< abs(1/fBetaToF[tr2] - inverse_pion_beta2)<<endl;
					 //cout << "p="<<fpT[tr2]<<"\tfnSigmaPion="<<fnSigmaPion[tr2]<<"\tdca="<<fdca[tr2]<< "\tfitPts="<<ffitPts[tr2]<<"\thitsdEdx="<<fhitsdedx[tr2]<<"\teta="<<feta[tr2]<< "\tfitratio="<<fitPtsRatio_tr2<< endl;				
					 //cout << fpT[tr] << ", "<< fpT[tr2]<< endl;
					 if(fcharge[tr]!=fcharge[tr2])
					 {
						 //cout<<"Pair found...."<<endl;

						 double  phiDiff=fphi[tr]-fphi[tr2];
						 if(phiDiff>pi)
							 phiDiff-=(2*pi);
						 if(phiDiff<((-1)*pi))
							 phiDiff+=(2*pi);
						 if(phiDiff>pi || phiDiff <-pi) continue;
						 cone = sqrt(pow(feta[tr]-feta[tr2],2)+pow(phiDiff,2));  //cone cut < 0.3

						 if(cone<0.7) 
						 {
							 paircount++;
							 if(fcharge[tr]>0)
							 {
								 //piplus->Fill(feta[tr],fphi[tr],fpT[tr],fp[tr],ffitPts[tr],fnSigmaPion[tr]);
								 //calculate variables
								 p1x = fpT[tr]*cos(fphi[tr]);
								 p2x = fpT[tr2]*cos(fphi[tr2]);
								 p1y = fpT[tr]*sin(fphi[tr]);
								 p2y = fpT[tr2]*sin(fphi[tr2]);
								 p1z = fpT[tr]*sinh(feta[tr]);
								 p2z = fpT[tr2]*sinh(feta[tr2]);
								 p1 = sqrt(p1x*p1x + p1y*p1y + p1z*p1z);
								 p2 = sqrt(p2x*p2x + p2y*p2y + p2z*p2z);

								 p_tr = fpT[tr]*cosh(feta[tr]);
								 p_tr2 = fpT[tr2]*cosh(feta[tr2]);

								 ptr = fp[tr];
								 ptr2 = fp[tr2];	
								 //cout << "p1x: "<< p1x <<"  "<< "charge: "<<fcharge[tr]<< endl;
							 }

							 if(fcharge[tr]<0)
							 {
								 //piminus->Fill(feta[tr],fphi[tr],fpT[tr],fp[tr],ffitPts[tr],fnSigmaPion[tr]);
								 //this shoud be done for charge ordering. I want the first track to be positive and the second
								 //to be negative. 
								 p1x = fpT[tr2]*cos(fphi[tr2]);
								 p2x = fpT[tr]*cos(fphi[tr]);
								 p1y = fpT[tr2]*sin(fphi[tr2]);
								 p2y = fpT[tr]*sin(fphi[tr]);
								 p1z = fpT[tr2]*sinh(feta[tr2]);
								 p2z = fpT[tr]*sinh(feta[tr]);
								 p1 = sqrt(p1x*p1x + p1y*p1y + p1z*p1z);
								 p2 = sqrt(p2x*p2x + p2y*p2y + p2z*p2z);
								 //cross check momentum ....
								 p_tr = fpT[tr2]*cosh(feta[tr2]);
								 p_tr2 = fpT[tr]*cosh(feta[tr]);

								 ptr = fp[tr2];
								 ptr2 = fp[tr];

								 //cout << "p1x: "<< p1x <<"  "<< "charge: "<<fcharge[tr2]<< endl;	
							 }
							 //cout << "p_tr: "<< p_tr <<",  "<< p1<<", "<<ptr<< endl;
							 //cout << "p_tr2: "<< p_tr2 <<",  "<< p2<<", "<<ptr2<< endl;
							 //momentum all good !!
							 if(fpT[tr]>fpT[tr2])pT_min_pair=fpT[tr2];
							 if(fpT[tr]<fpT[tr2])pT_min_pair=fpT[tr];
							 if(ffitPts[tr]>ffitPts[tr2])fitPts_min_pair=ffitPts[tr2];
							 if(ffitPts[tr]<ffitPts[tr2])fitPts_min_pair=ffitPts[tr];
							 //Components of sum vector
							 psx = p1x+p2x;
							 psy = p1y+p2y;
							 psz = p1z+p2z;
							 //Sum vector
							 ps = sqrt(psx*psx+psy*psy+psz*psz);

							 //Relatve momentum of dihadron system
							 //Charge ordering is important. In R = (Ph1-Ph2)*.5 ,I want  Ph1 to be positive and Ph2 to be negative always. 
							 double Rx1, Rz1,Ry1;
							 if(fcharge[tr]>0){
								 Rx1 = (p1x-p2x);
								 Ry1 = (p1y-p2y);
								 Rz1 = (p1z-p2z);
							 } else if(fcharge[tr2]>0){
								 Rx1 = (p2x-p1x);
								 Ry1 = (p2y-p1y);
								 Rz1 = (p2z-p1z);
							 }

							 Rx = (p1x-p2x);
							 Ry = (p1y-p2y);
							 Rz = (p1z-p2z);
							 R = sqrt(Rx*Rx+Ry*Ry+Rz*Rz);//R and R1 are same
							 R1 = sqrt(Rx1*Rx1+Ry1*Ry1+Rz1*Rz1);
							 //cout << "R1: "<< R1 << ", R: "<< R << endl; //Exact same output 
							 //calculate M,pt,eta
							 E1 = sqrt(m_pion*m_pion+p1*p1);
							 E2 = sqrt(m_pion*m_pion+p2*p2);
							 Minv = sqrt(2*m_pion*m_pion + 2*(E1*E2-p1x*p2x-p1y*p2y-p1z*p2z));
							 pT_pair = sqrt((p1x+p2x)*(p1x+p2x) + (p1y+p2y)*(p1y+p2y));
							 eta_pair = TMath::ASinH((p1z+p2z)/pT_pair); 

							 if(Minv<0. || Minv>4.)continue;
							 if(pT_pair<2.5 || pT_pair>=15.)continue;

							 //cout<<"pT_pair: "<<pT_pair<<", Minv: "<<Minv<<endl;

							 //Vector definations for phiRS calculation
							 TVector3 v_blue(0,0,1); //Blue beam
							 TVector3 v_yell(0,0,-1);//yellow beam
							 TVector3 pSum(psx,psy,psz);//pion pair sum vector
							 TVector3 SUp(0,1,0);//spin UP
							 TVector3 SDn=SUp;//(same as 200GeV analysis to fix y-direction)
							 //TVector3 SDn(0,-1,0);//spin down 
							 TVector3 pSumHat(psx/ps,psy/ps,psz/ps); //pion pair's momentum unit vector
							 TVector3 Rc(Rx,Ry,Rz);//pion pair's momemtum difference vector 

							 //blue beam
							 cosPhiRB=(pSumHat.Cross(v_blue)*((double)1/(pSumHat.Cross(v_blue).Mag()))).Dot((pSumHat.Cross(Rc))*((double)1/(pSumHat.Cross(Rc).Mag())));
							 sinPhiRB=(v_blue.Cross(Rc)).Dot(pSumHat)*((double)1/(pSumHat.Cross(v_blue).Mag()))*((double)1/(pSumHat.Cross(Rc).Mag()));
							 PhiR_cosB = acos(cosPhiRB);
							 PhiR_sinB = asin(sinPhiRB);

							 if(PhiR_sinB>0)PhiRB=PhiR_cosB;
							 if(PhiR_sinB<0)PhiRB=(-1)*PhiR_cosB;
							 //cout << " PhiR_cosB= "<<PhiR_cosB<< " PhiR_sinB="<< PhiR_sinB<< " PhiRB="<<PhiRB<<endl; 
							 if(fspinconfig==51 || fspinconfig==53)//blue up
							 {	blueUp++;
								 if(fspinconfig==51)
									 spin51++;
								 if(fspinconfig==53)
									 spin53++;
								 //cout << "Blue Up : "<< fspinconfig << endl;
								 cosPhiSB=((v_blue.Cross(pSum))*((double)1/(v_blue.Cross(pSum).Mag()))).Dot((v_blue.Cross(SUp))*((double)1/(v_blue.Cross(SUp).Mag())));
								 sinPhiSB=(pSum.Cross(SUp)).Dot(v_blue)*((double)1/(v_blue.Cross(pSum).Mag()))*((double)1/(v_blue.Cross(SUp).Mag()));
							 }else if(fspinconfig==83 || fspinconfig==85)//blue down
							 {
								 blueDown++;
								 if(fspinconfig==83)
									 spin83++;
								 if(fspinconfig==85)
									 spin85++;
								 //cout << "Blue  down: "<< fspinconfig << endl;
								 cosPhiSB=((v_blue.Cross(pSum))*((double)1/(v_blue.Cross(pSum).Mag()))).Dot((v_blue.Cross(SDn))*((double)1/(v_blue.Cross(SDn).Mag())));
								 sinPhiSB=(pSum.Cross(SDn)).Dot(v_blue)*((double)1/(v_blue.Cross(pSum).Mag()))*((double)1/(v_blue.Cross(SDn).Mag()));
							 }else {
								 cout <<"Wrong spin configuration: "<< fspinconfig << " . Exiting...."<< endl;
								 continue;}
							 PhiS_cosB = acos(cosPhiSB);
							 PhiS_sinB = asin(sinPhiSB);
							 if(PhiS_sinB>0)PhiSB=PhiS_cosB;
							 if(PhiS_sinB<0)PhiSB=(-1)*PhiS_cosB;
							 //cout << "PhiRS: "<< PhiSB - PhiRB << endl;
							 //PhiRSB
							 if((PhiSB-PhiRB)<-pi)PhiRSB = (PhiSB-PhiRB)+2*pi;
							 else if((PhiSB-PhiRB)>pi)PhiRSB =(PhiSB-PhiRB)-2*pi;
							 else PhiRSB = PhiSB-PhiRB;
							 //cout << " PhiS_cosB= "<<PhiS_cosB<< " PhiS_sinB="<< PhiS_sinB<< " PhiSB="<<PhiSB<<", PhiRSB: "<< PhiRSB<<endl; 

							 //yellow beam
							 cosPhiRY=(pSumHat.Cross(v_yell)*((double)1/(pSumHat.Cross(v_yell).Mag()))).Dot((pSumHat.Cross(Rc))*((double)1/(pSumHat.Cross(Rc).Mag())));
							 sinPhiRY=(v_yell.Cross(Rc)).Dot(pSumHat)*((double)1/(pSumHat.Cross(v_yell).Mag()))*((double)1/(pSumHat.Cross(Rc).Mag()));
							 PhiR_cosY = acos(cosPhiRY);
							 PhiR_sinY = asin(sinPhiRY);
							 if(PhiR_sinY>0)PhiRY=PhiR_cosY;
							 if(PhiR_sinY<0)PhiRY=-1.*PhiR_cosY;
							 if(fspinconfig==51 || fspinconfig==83)//yellow up
							 {
								 //cout << "Yellow  up: "<< fspinconfig << endl;
								 yellowUp++;
								 cosPhiSY=((v_yell.Cross(pSum))*((double)1/(v_yell.Cross(pSum).Mag()))).Dot((v_yell.Cross(SUp))*((double)1/(v_yell.Cross(SUp).Mag())));
								 sinPhiSY=(pSum.Cross(SUp)).Dot(v_yell)*((double)1/(v_yell.Cross(pSum).Mag()))*((double)1/(v_yell.Cross(SUp).Mag()));
							 }
							 if(fspinconfig==53 || fspinconfig==85)//yellow down
							 {
								 //cout << "Yellow  down: "<< fspinconfig << endl;
								 yellowDown++;
								 cosPhiSY=((v_yell.Cross(pSum))*((double)1/(v_yell.Cross(pSum).Mag()))).Dot((v_yell.Cross(SDn))*((double)1/(v_yell.Cross(SDn).Mag())));
								 sinPhiSY=(pSum.Cross(SDn)).Dot(v_yell)*((double)1/(v_yell.Cross(pSum).Mag()))*((double)1/(v_yell.Cross(SDn).Mag()));
							 }
							 PhiS_cosY = acos(cosPhiSY);
							 PhiS_sinY = asin(sinPhiSY);
							 if(PhiS_sinY>0)PhiSY=PhiS_cosY;
							 if(PhiS_sinY<0)PhiSY=(-1)*PhiS_cosY;

							 //PhiRSY
							 if((PhiSY-PhiRY)<(-pi))PhiRSY = (PhiSY-PhiRY)+(2*pi);
							 else if((PhiSY-PhiRY)>pi)PhiRSY = (PhiSY-PhiRY)-(2*pi);
							 else PhiRSY = PhiSY-PhiRY;


							 picut1=pionCut(fnSigmaPion[tr],fnSigmaPionTof[tr]);							
							 picut2=pionCut(fnSigmaPion[tr2],fnSigmaPionTof[tr2]);
							 picut1pt=pionCutPt(fnSigmaPion[tr],fnSigmaPionTof[tr],pT_pair);							
							 picut2pt=pionCutPt(fnSigmaPion[tr2],fnSigmaPionTof[tr2],pT_pair);
							 picut1m=pionCutPt(fnSigmaPion[tr],fnSigmaPionTof[tr],Minv);							
							 picut2m=pionCutPt(fnSigmaPion[tr2],fnSigmaPionTof[tr2],Minv);
							 picut1eta=pionCutPt(fnSigmaPion[tr],fnSigmaPionTof[tr],eta_pair);							
							 picut2eta=pionCutPt(fnSigmaPion[tr2],fnSigmaPionTof[tr2],eta_pair);

							 if((fnSigmaPion[tr]>-1&&fnSigmaPion[tr]<2) && (fnSigmaPion[tr2]>-1&&fnSigmaPion[tr2]<2)){
								 //PID from TPC only
								 ntuple1tpc->Fill(cone,fspinconfig,frunNum,ffillNum);
								 ntuple2tpc->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
								 ntuple3tpc->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
								 ntuple4tpc->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
								 ntuple5tpc->Fill(fitPts_min_pair,msqr);
								 //----------------------
								 //Polynomial cut is based on PID to remove proton and kaon background based on TOF and TPC.
								 //accept events if both tracks doesn't have TOF
								 if(picut1==-5555 && picut2==-5555){
									 ntuple1f->Fill(cone,fspinconfig,frunNum,ffillNum);
									 ntuple2f->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
									 ntuple3f->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
									 ntuple4f->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
									 ntuple5f->Fill(fitPts_min_pair,msqr);
								 }else if(picut1!=-5555 && fnSigmaPion[tr]>picut1 && picut2!=-5555 &&fnSigmaPion[tr2]>picut2){
									 //accept events if both tracks has TOF info and passed pid cut
									 for(int i=0; i<5; i++){
										 if(Minv>=M[i] && Minv<M[i+1]){
											 if(fcharge[tr]>0){
												 hTpcVsTofM[0][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofM[1][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
											 if(fcharge[tr]<0){
												 hTpcVsTofM[1][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofM[0][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
										 }
										 if(pT_pair>=pT[i] && pT_pair<pT[i+1]){
											 if(fcharge[tr]>0){
												 hTpcVsTofpT[0][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofpT[1][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
											 if(fcharge[tr]<0){
												 hTpcVsTofpT[1][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofpT[0][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
										 }
									 }
									 ntuple1f->Fill(cone,fspinconfig,frunNum,ffillNum);
									 ntuple2f->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
									 ntuple3f->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
									 ntuple4f->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
									 ntuple5f->Fill(fitPts_min_pair,msqr);
								 }else if(picut1==-5555 && picut2!=-5555 && fnSigmaPion[tr2]>picut2){
									 //accept events if first track don't have TOF info, but track second has and passed the pid cut
									 ntuple1f->Fill(cone,fspinconfig,frunNum,ffillNum);
									 ntuple2f->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
									 ntuple3f->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
									 ntuple4f->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
									 ntuple5f->Fill(fitPts_min_pair,msqr);
								 }else if(picut1!=-5555 && fnSigmaPion[tr]>picut1 && picut2==-5555 ){
									 //accept events if first track has TOF info and passed pid cut and track second don't have TOF info
									 ntuple1f->Fill(cone,fspinconfig,frunNum,ffillNum);
									 ntuple2f->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
									 ntuple3f->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
									 ntuple4f->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
									 ntuple5f->Fill(fitPts_min_pair,msqr);
								 }else{cout<<"Skip this event.....!"<<endl; }
								 //-------------------------------------------------


								 //accept events if both tracks have TOF
								 if(picut1!=-5555 && fnSigmaPion[tr]>picut1 && picut2!=-5555 && fnSigmaPion[tr2]>picut2){
									 ntuple1tof->Fill(cone,fspinconfig,frunNum,ffillNum);
									 ntuple2tof->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
									 ntuple3tof->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
									 ntuple4tof->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
									 ntuple5tof->Fill(fitPts_min_pair,msqr);
								 }
								 //accept background events if both tracks have TOF
								 if(picut1!=-5555 && fnSigmaPion[tr]<picut1 && picut2!=-5555 && fnSigmaPion[tr2]<picut2){
									 for(int i=0; i<5; i++){
										 if(Minv>=M[i] && Minv<M[i+1]){
											 if(fcharge[tr]>0){
												 hTpcVsTofMBkg[0][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofMBkg[1][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
											 if(fcharge[tr]<0){
												 hTpcVsTofMBkg[1][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofMBkg[0][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
										 }
										 if(pT_pair>=pT[i] && pT_pair<pT[i+1]){
											 if(fcharge[tr]>0){
												 hTpcVsTofpTBkg[0][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofpTBkg[1][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
											 if(fcharge[tr]<0){
												 hTpcVsTofpTBkg[1][i]->Fill(fnSigmaPionTof[tr], fnSigmaPion[tr]);
												 hTpcVsTofpTBkg[0][i]->Fill(fnSigmaPionTof[tr2], fnSigmaPion[tr2]);
											 }
										 }
									 }
									 ntuple1tofbkg->Fill(cone,fspinconfig,frunNum,ffillNum);
									 ntuple2tofbkg->Fill(R,ps,Minv,pT_pair,pT_min_pair,eta_pair,psx,psy,psz);
									 ntuple3tofbkg->Fill(cosPhiRB,sinPhiRB,cosPhiSB,sinPhiSB,cosPhiRY,sinPhiRY,cosPhiSY,sinPhiSY);
									 ntuple4tofbkg->Fill(PhiR_cosB,PhiR_sinB,PhiRB,PhiS_cosB,PhiS_sinB,PhiSB,PhiRSB,PhiR_cosY,PhiR_sinY,PhiRY,PhiS_cosY,PhiS_sinY,PhiSY,PhiRSY);
									 ntuple5tofbkg->Fill(fitPts_min_pair,msqr);

								 }
							 }//tpc sigma pion cut
								//-------------------------------------------------
							 
						}//cone cut here 
					 }//charge comparison
				 }//track 2 PID and selection cuts 
			 }//track 2 loop	
		 }//track 1 PID and selection cuts   
	 }// track loop
 }//event loop

}//Iff::Loop() 


Iff2015::Iff2015(char* ifile )
{
	fChain=new TChain("ftree");
	fChain->Add(ifile);

	// if parameter tree is not specified (or zero), connect the file
	// used to generate this class and read the Tree.
	Init();
}

Iff2015::~Iff2015()
{
	if (!fChain) return;
	delete fChain->GetCurrentFile();
}

Int_t Iff2015::GetEntry(Long64_t entry)
{
	// Read contents of entry.
	if (!fChain) return 0;
	return fChain->GetEntry(entry);
}
Long64_t Iff2015::LoadTree(Long64_t entry)
{
	// Set the environment to read one entry
	if (!fChain) return -5;
	Long64_t centry = fChain->LoadTree(entry);
	if (centry < 0) return centry;
	if (!fChain->InheritsFrom(TChain::Class()))  return centry;
	TChain *chain = (TChain*)fChain;
	if (chain->GetTreeNumber() != fCurrent) {
		fCurrent = chain->GetTreeNumber();
		Notify();
	}
	return centry;
}

void Iff2015::Init()
{
	// The Init() function is called when the selector needs to initialize
	// a new tree or chain. Typically here the branch addresses and branch
	// pointers of the tree will be set.
	// It is normally not necessary to make changes to the generated
	// code, but the routine can be extended by the user if needed.
	// Init() will be called many times when running on PROOF
	// (once per file to be processed).

	// Set branch addresses and branch pointers

	fCurrent = -1;
	fChain->SetMakeClass(1);

	fChain->SetBranchAddress("frefmult", &frefmult, &b_frefmult);
	fChain->SetBranchAddress("fmaxpar", &fmaxpar, &b_fmaxpar);
	fChain->SetBranchAddress("ffillNum", &ffillNum, &b_ffillNum);
	fChain->SetBranchAddress("frunNum", &frunNum, &b_frunNum);
	fChain->SetBranchAddress("fspinconfig", &fspinconfig, &b_fspinconfig);
	fChain->SetBranchAddress("ftrigger", &ftrigger);
	fChain->SetBranchAddress("fVZ", &fVZ, &b_fVZ);
	fChain->SetBranchAddress("fevTime",&fevTime, &b_fevTime);
	fChain->SetBranchAddress("fverRank", &fverRank, &b_fverRank);
	fChain->SetBranchAddress("fpT", fpT, &b_fpT);
	fChain->SetBranchAddress("fp", fp, &b_fp);
	fChain->SetBranchAddress("feta", feta, &b_feta);
	fChain->SetBranchAddress("fphi", fphi, &b_fphi);
	fChain->SetBranchAddress("fcharge", fcharge, &b_fcharge);
	fChain->SetBranchAddress("fnSigmaPion", fnSigmaPion, &b_fnSigmaPion);
	fChain->SetBranchAddress("fnSigmaKaon", fnSigmaKaon, &b_fnSigmaKaon);
	fChain->SetBranchAddress("fnSigmaProton", fnSigmaProton, &b_fnSigmaProton);
	fChain->SetBranchAddress("fnSigmaElectron", fnSigmaElectron, &b_fnSigmaElectron);
	fChain->SetBranchAddress("fnSigmaPionTof", fnSigmaPionTof, &b_fnSigmaPionTof);
	fChain->SetBranchAddress("fnSigmaKaonTof", fnSigmaKaonTof, &b_fnSigmaKaonTof);
	fChain->SetBranchAddress("fnSigmaProtonTof", fnSigmaProtonTof, &b_fnSigmaProtonTof);
	fChain->SetBranchAddress("fnSigmaElectronTof", fnSigmaElectronTof, &b_fnSigmaElectronTof);
	fChain->SetBranchAddress("fdEdx", fdEdx, &b_fdEdx);
	fChain->SetBranchAddress("fdca", fdca, &b_fdca);
	fChain->SetBranchAddress("ffitPts", ffitPts, &b_ffitPts);
	fChain->SetBranchAddress("ffitPtsPoss", ffitPtsPoss, &b_ffitPtsPoss);
	fChain->SetBranchAddress("fhitsdedx", fhitsdedx, &b_fhitsdedx);
	fChain->SetBranchAddress("fBetaToF", fBetaToF, &b_fBetaToF);
	fChain->SetBranchAddress("fvpdVz", &fvpdVz, &b_fvpdVz);
	Notify();
}

Bool_t Iff2015::Notify()
{
	// The Notify() function is called when a new file is opened. This
	// can be either for a new TTree in a TChain or when when a new TTree
	// is started when using PROOF. It is normally not necessary to make changes
	// to the generated code, but the routine can be extended by the
	// user if needed. The return value is currently not used.

	return kTRUE;
}


//load run15 polarization text file 
  
void loadPol(const char* Filename, map<int,unsigned int>& t0Map, map<int,float>& pb0Map,map<int,float>& dpdtbMap, map<int,float>& py0Map, map<int,float>& dpdtyMap)
//void loadPol(const char* Filename,unsigned int t0Map, float pb0Map,float dpdtbMap, float py0Map, float dpdtyMap)
{ 
	FILE* file = fopen(Filename,"r");
	//assert(file);
	//cout << "file exist "<< Filename << endl;
	int fillnumber; 
	unsigned int t0;
	float pb0, pb0err;
	float dpdtb, dpdtberr;
	float py0, py0err;
	float dpdty, dpdtyerr;
	while (fscanf(file,"%d %d %f +- %f %f +- %f %f +- %f %f +- %f\n",&fillnumber,&t0,&pb0,&pb0err,&dpdtb,&dpdtberr,&py0,&py0err,&dpdty,&dpdtyerr) != EOF){

		int fill=fillnumber-18640;
		t0Map[fill] = t0;
		pb0Map[fill] = pb0;
		dpdtbMap[fill] = dpdtb;
		py0Map[fill] = py0;
		dpdtyMap[fill] = dpdty;
		//cout <<"For fill: "<< fill << " "<< t0Map[fill] <<"  " <<  pb0Map[fill]<<"  "<<dpdtbMap[fill]<<"  "<<py0Map[fill]<<"  "<< dpdtyMap[fill]<<endl;
	}
	fclose(file);
	//cout << "Done with polarization correction."<< endl;
}

//PID cut using TOF and TPC with polynomial cut. This removes proton and kaon background
/*double  pionCutPt(double sigma_tpc,double sigma_tof, double pt_pair){
    double picut;
    double ptPair[6]={2.80,   3.71,  4.303 ,   5.084,   6.404,  15.00};
    if(sigma_tof==-999){
       return -5555;//no tof info, use this track
    }
    else if(sigma_tof!=-999){
     if(pt_pair>=ptPair[0]&&pt_pair<ptPair[1]){
     picut= 0.0156*pow(sigma_tof,3)+0.1187*pow(sigma_tof,2)+0.6996*sigma_tof-2.9738;
     }
     else if(pt_pair>=ptPair[1]&&pt_pair<ptPair[2]){
     picut= 0.034*pow(sigma_tof,3)+0.1822*pow(sigma_tof,2)+0.4719*sigma_tof-3.1947;
     }
     else if(pt_pair>=ptPair[2]&&pt_pair<ptPair[3]){
     picut= 0.0328*pow(sigma_tof,3)+0.197*pow(sigma_tof,2)+0.3386*sigma_tof-2.9323;
     }
     else if(pt_pair>=ptPair[3]&&pt_pair<ptPair[4]){
     picut= 0.038*pow(sigma_tof,3)+0.185*pow(sigma_tof,2)+0.1263*sigma_tof-2.9155;
     }
     else if(pt_pair>=ptPair[4]&&pt_pair<ptPair[5]){
     picut= 0.0303*pow(sigma_tof,3)+0.1482*pow(sigma_tof,2)+0.1978*sigma_tof-2.5664;
     }
     return picut; //polonomial cut value to remove proton-kaon background
   }

}
*/
double  pionCut(double sigma_tpc,double sigma_tof){
    double picut;
    if(sigma_tof==-999){
       return -5555;//no tof info, use this track
    }
    else if(sigma_tof!=-999){
     picut= 0.0156*pow(sigma_tof,3)+0.1187*pow(sigma_tof,2)+0.6996*sigma_tof-2.9738;
     }
     return picut; //polonomial cut value to remove proton-kaon background
}
// same polynomial cut in all pT bins
double  pionCutPt(double sigma_tpc,double sigma_tof, double pt_pair){
    double picut;
    double ptPair[6]={2.80,   3.71,  4.303 ,   5.084,   6.404,  15.00};
    if(sigma_tof==-999){
       return -5555;//no tof info, use this track
    }
    else if(sigma_tof!=-999){
     picut= 0.0156*pow(sigma_tof,3)+0.1187*pow(sigma_tof,2)+0.6996*sigma_tof-2.9738;
     }
     return picut; //polonomial cut value to remove proton-kaon background
}
double  pionCutM(double sigma_tpc,double sigma_tof, double minv){
    double picut;
    if(sigma_tof==-999){
       return -5555;//no tof info, use this track
    }
    else if(sigma_tof!=-999){
     picut= 0.0156*pow(sigma_tof,3)+0.1187*pow(sigma_tof,2)+0.6996*sigma_tof-2.9738;
     }
     return picut; //polonomial cut value to remove proton-kaon background
   }
double  pionCutEta(double sigma_tpc,double sigma_tof, double etaPair){
    double picut;
    if(sigma_tof==-999){
       return -5555;//no tof info, use this track
    }
    else if(sigma_tof!=-999){
     picut= 0.0156*pow(sigma_tof,3)+0.1187*pow(sigma_tof,2)+0.6996*sigma_tof-2.9738;
     }
     return picut; //polonomial cut value to remove proton-kaon background
   }






void Iff2015::Finish(char* ofile)
{
	TFile * fout=new TFile(ofile,"recreate");

	ntuple1tpc->SetDirectory(fout);
	ntuple2tpc->SetDirectory(fout);
	ntuple3tpc->SetDirectory(fout);
	ntuple4tpc->SetDirectory(fout);
	ntuple5tpc->SetDirectory(fout);
	ntuple1tpc->Write();
	ntuple2tpc->Write();
	ntuple3tpc->Write();
	ntuple4tpc->Write();
	ntuple5tpc->Write();
	
	ntuple1f->SetDirectory(fout);
	ntuple2f->SetDirectory(fout);
	ntuple3f->SetDirectory(fout);
	ntuple4f->SetDirectory(fout);
	ntuple5f->SetDirectory(fout);
	ntuple1f->Write();
	ntuple2f->Write();
	ntuple3f->Write();
	ntuple4f->Write();
	ntuple5f->Write();
	
	ntuple1tof->SetDirectory(fout);
	ntuple2tof->SetDirectory(fout);
	ntuple3tof->SetDirectory(fout);
	ntuple4tof->SetDirectory(fout);
 	ntuple5tof->SetDirectory(fout);
	ntuple1tof->Write();
	ntuple2tof->Write();
	ntuple3tof->Write();
	ntuple4tof->Write();
	ntuple5tof->Write();

	ntuple1tofbkg->SetDirectory(fout);
	ntuple2tofbkg->SetDirectory(fout);
	ntuple3tofbkg->SetDirectory(fout);
	ntuple4tofbkg->SetDirectory(fout);
 	ntuple5tofbkg->SetDirectory(fout);
	ntuple1tofbkg->Write();
	ntuple2tofbkg->Write();
	ntuple3tofbkg->Write();
	ntuple4tofbkg->Write();
	ntuple5tofbkg->Write();


	
	ntuple1bkg->SetDirectory(fout);
	ntuple2bkg->SetDirectory(fout);
	ntuple3bkg->SetDirectory(fout);
	ntuple4bkg->SetDirectory(fout);
 	ntuple5bkg->SetDirectory(fout);
	ntuple1bkg->Write();
	ntuple2bkg->Write();
	ntuple3bkg->Write();
	ntuple4bkg->Write();
	ntuple5bkg->Write();
	
	ntuple6->SetDirectory(fout);
	ntuple6->Write();

	for(int ch =0; ch<2; ch++){
		for (int i =0 ; i<5; i++){
			hTpcVsTofM[ch][i]->SetDirectory(fout); 
			hTpcVsTofM[ch][i]->Write(); 
			hTpcVsTofpT[ch][i]->SetDirectory(fout);   
			hTpcVsTofpT[ch][i]->Write();  
			hTpcVsTofMBkg[ch][i]->SetDirectory(fout);  
			hTpcVsTofMBkg[ch][i]->Write();
			hTpcVsTofpTBkg[ch][i]->SetDirectory(fout);  
			hTpcVsTofpTBkg[ch][i]->Write();
		}
	}
	

	fout->Close();
}

void Iff2015::Show(Long64_t entry)
{
	//eventTime-> Write();
	//pions->Write();
	// Print contents of entry.
	// If entry is not specified, print current entry
	if (!fChain) return;
	fChain->Show(entry);
}
Int_t Iff2015::Cut(Long64_t entry)
{
	// This function may be called from Loop.
	// returns  1 if entry is accepted.
	// returns -1 otherwise.
	return 1;
}

ClassImp(Iff2015);
