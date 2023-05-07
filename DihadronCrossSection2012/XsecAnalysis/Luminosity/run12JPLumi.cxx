/*
Run 12 pp200 Gev Luminosity calculation for JP0, JP1, and JP2
-> Luminosity per trigger per run is grabbed from STAR online
    star.bnl.gov/protected/common/common2012/trigger2012/lumipp200GeV/lum_perrun_JP*.txt

-> Only the QAed run luminosities are added for each trigger for the total delivered luminosity per trigger
-> These values will be used for dihadron cross-section analysis
***** Babu Pokhrel, Temple University (04/20/2022)******
*/

#include <iostream>
#include "TFile.h"
#include "TString.h"
#include <map>

using namespace std;

void run12JPLumi()
{
    const Int_t nruns = 601;
    Int_t run12[nruns] = {0};
    // graphs Luminosity vs run index
    TGraph *gLumVsRunJP0;
    TGraph *gLumVsRunJP1;
    TGraph *gLumVsRunJP2;

    ifstream runFile;
    runFile.open("KevinsRunList.txt");
    if (runFile.is_open())
    {
        Int_t run;
        Int_t ct = 0;
        cout << "Run list found!" << endl;
        while (runFile >> run)
        {
            run12[ct] = run;
            // cout << "Counter: " << ct << ", Run: " << run12[ct] << endl;
            ct++;
        }
    }
    else
    {
        cout << "Run list not found! Exiting...." << endl;
        return;
    }

    ofstream fusedRun[3];
    Int_t nfiles = 3;             // triggers JP0, JP1, and JP2
    Double_t lumJP[3] = {0};      // total delivered luminosity for each triggers
    Double_t prescaleJP[3] = {0}; // 0, 1, 2 average prescale
    string trig[3] = {"JP0", "JP1", "JP2"};

    for (Int_t i = 0; i < nfiles; i++)
    {
        fusedRun[i].open(Form("trig%s.txt", trig[i]));
        Int_t runNum, t1, t2, fill;
        Double_t lum, preScale, a1, a2, a3, a4, a5, a6, a7;
        string bbcmb;

        const char *fname = Form("lum_perrun_JP%i.txt", i);
        ifstream datFile;
        datFile.open(fname, ios::in);
        if (datFile.is_open())
        {
            cout << "datFile: " << fname << endl;
            while (datFile >> runNum >> t1 >> t2 >> fill >> lum >> preScale >> a1 >> bbcmb >> a2 >> a3 >> a4)
            {
                // cout << " runNum : " << runNum << ", lum: " << lum << endl;
                for (Int_t j = 0; j < nruns; j++)
                {
                    if (run12[j] != runNum)
                    {
                        continue;
                    }
                    else
                    {
                        // cout << "Run number matched...." << runNum << endl;
                        fusedRun[i] << j << "\t" << runNum << "\t" << lum << endl;
                        lumJP[i] += lum;
                        prescaleJP[i] += preScale;
                    }
                }
            }
        }
        else
        {
            cout << "File not found....." << endl;
            break;
        }
        fusedRun[i] << Form("**** Total Delivered Luminosity for %s = %g  ****", trig[i], lumJP[i]) << endl;
    }
    cout << "******** Total Delivered Luminosities***********" << endl;
    cout << "Luninosity  JPO: " << lumJP[0] << ",  JP1: " << lumJP[1] << ",  JP2: " << lumJP[2] << endl;
    cout << "Prescale JPO: " << prescaleJP[0] / 601 << ",  JP1: " << prescaleJP[1] / 601 << ",  JP2: " << prescaleJP[2] / 601 << endl;
    cout << "************************************************" << endl;
}