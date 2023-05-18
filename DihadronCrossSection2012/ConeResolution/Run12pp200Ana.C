class Run12pp200Ana;

void Run12pp200Ana(const char *ifile = "/star/u/pokhrel/GPFS/Run12EmbeddingTrees/v3/tree_13047003.root", const char *ofile = "test.root") // interactive mode
// void Run12pp200Ana(const char *ifile = "/star/u/pokhrel/GPFS/Run12EmbeddingTrees/NewBatch3/*.root", const char *ofile = "hist4ResAndUnfolding.root") // interactive mode
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
  // gSystem->Load("/star/u/pokhrel/DihadronCrossSection/DataEmbedComparison/ReadEmbedTrees/.sl73_gcc485/lib/Run12pp200Ana.so ");
  //========================================================================

  //  gSystem->Load("libPhysics");
  //  gSystem->Load("St_base");
  gSystem->Load("Run12pp200Ana");
}
