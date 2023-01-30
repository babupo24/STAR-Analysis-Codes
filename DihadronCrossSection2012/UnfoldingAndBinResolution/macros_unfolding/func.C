#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;

void func()
{
    TFile *ifile = new TFile("Histograms/hist4Unfoding341.root", "r");
    if (!ifile)
        cout << "Input file not found!" << endl;

    //get input hsitograms
    TH1D *h1 = (TH1D *)ifile->Get("hMGenVsRecJP0NoUF");
    TH1D *h2 = (TH1D *)ifile->Get("hMGenVsRecJP0NoUF");
}