// get the average invariant mass and pt of dipion in cross section bins

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include <cmath>
#include <fstream>
#include <iomanip>

using namespace std;

void getAvgQt()
{
    TFile *fhist = new TFile("hist_Xsec.root");
    // TFile *fhist = new TFile("test.root");
    if (!fhist->IsOpen())
        break;

    const int nbins = 13;
    Double_t nBinsEdges[nbins + 1] = {0.27, 0.35, 0.45, 0.60, 0.75, 0.95, 1.15, 1.35, 1.60, 1.90, 2.20, 2.60, 3.20, 4.0};
    TH1D *hM[nbins] = {NULL};
    Double_t avgM[nbins] = {0};
    TH1D *hPt[nbins] = {NULL};
    Double_t avgPt[nbins] = {0};
    TH1D *hEta[nbins] = {NULL};
    Double_t avgEta[nbins] = {0};

    ofstream ftab;
    ftab.open("ResultsMinConeCut/table_avgMpT.txt");
    ofstream fout;
    fout.open("ResultsMinConeCut/avgMpT.txt");
    fout << "#Average dipion invariant mass, p_T, and eta in cross section bins" << endl;
    fout << "----------------------------------------------------------------------------------" << endl;
    fout << "Bin Id \t M_{inv} Range \t <M_{inv}>(GeV/c^2) \t <p_{T}>(GeV/c) \t <Eta> " << endl;
    fout << "----------------------------------------------------------------------------------" << endl;

    ftab << "Bin Id & M_{inv} Bin (GeV/c^2) & <M_{inv}>(GeV/c^2) & <p_{T}>(GeV/c) & <Eta> " << endl;
    for (int i = 0; i < nbins; i++)
    {
        hM[i] = (TH1D *)fhist->Get(Form("hMinv_%i", i));
        avgM[i] = hM[i]->GetMean();

        hPt[i] = (TH1D *)fhist->Get(Form("hPt_%i", i));
        avgPt[i] = hPt[i]->GetMean();

        hEta[i] = (TH1D *)fhist->Get(Form("hEta_%i", i));
        avgEta[i] = hEta[i]->GetMean();

        fout << i + 1 << " \t " << nBinsEdges[i] << " - " << nBinsEdges[i + 1] << " \t   " << avgM[i] << " \t \t   " << avgPt[i] << " \t      " << avgEta[i] << endl;
        ftab << i + 1 << " & " << setprecision(3) << nBinsEdges[i] << " - " << nBinsEdges[i + 1] << " &   " << setprecision(3) << avgM[i] << " & " << setprecision(3) << avgPt[i] << "& " << setprecision(3) << avgEta[i] << "\\\\" << endl;
    }
    fout << "----------------------------------------------------------------------------------" << endl;
    fout.close();
}