class Run12pp200Ana;

void Run12pp200Ana(const char *ifile = "/star/u/pokhrel/GPFS/Run12PicoTrees/run12pp200_13051071.root  ", const char *ofile = "testFile.root") // interactive mode
{

  LoadLibs();

  TDirectory *home = gDirectory;
  //  Run12pp200Ana* analysis = new Run12pp200Ana(filename);
  Run12pp200Ana *analysis = new Run12pp200Ana(ifile);
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
  gSystem->Load("Run12pp200Ana");
}
