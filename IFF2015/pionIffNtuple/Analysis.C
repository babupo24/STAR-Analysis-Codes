class Iff2015;

void Analysis(const char* ifile="/star/u/pokhrel/GPFS/IFF_TREES/StartLessTOF/Trees/ifftree_16093018.root", const char* ofile="test.root")//interactive mode
{ 

  LoadLibs();

  TDirectory *home = gDirectory;
//  Iff2015* analysis = new Iff2015(filename);
  Iff2015* analysis = new Iff2015(ifile);
  analysis->setHome(home);
  analysis->Loop();
  analysis->Finish(ofile);
}


void LoadLibs()
{
  //============= replace with your own libraries! =========================
  //  gSystem->Load("/star/u/drach09/pp+pA/pASimus/tppmc/tppmc-read-only_myTest/.libs/libtppmc");
  //  gSystem->Load("/star/u/drach09/pp+pA/pASimus2/eic-smear/.libs/libeicsmear");
  //========================================================================

  //  gSystem->Load("libPhysics");
  //  gSystem->Load("St_base");
  gSystem->Load("Iff2015");
}

  
