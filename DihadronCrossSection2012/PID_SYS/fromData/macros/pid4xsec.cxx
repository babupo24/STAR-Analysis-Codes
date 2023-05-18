// Synopsis and updates-------------------
/*
This macro reads histogram file containing pid histograms, which contains charge-separated 2D histogrmas, pion vs (p,k,e) and 1D histogram of nSigmaPion for gaus fit. 2D histogram is used to find peak positions of the other particle species in the pion signal region. Since the PID  being doing in sub-momentum region within mass bin, the file has a histogrmas in five sub-divided momentum bins for the lowest mass bin.

*/
// Babu Pokhrel, Temple University, 2022/07/12
//---------------------------------------
#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"

using namespace std;

void pid4xsec()
{
    TFile *fpid = new TFile("pidhist_Run12.root");
    if (!fpid)
    {
        cout << "File not found!" << endl;
    }
}