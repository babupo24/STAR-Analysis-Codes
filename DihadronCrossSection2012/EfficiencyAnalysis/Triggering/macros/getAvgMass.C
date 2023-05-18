#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TFile.h"

using namespace std;

void getAvgMass(const char *filename = "hist4TrigEff.root")
{
    TFile *fhist = new TFile(filename);
    if (!fhist->IsOpen())
    {
        cout << "File not found!" << endl;
        break;
    }
    // prescale values for combining triggers
    Double_t presJP0 = 0.00707; //(1. / 141.35);
    Double_t presJP1 = 0.3978;  //(1. / 2.514);
    Double_t presJP2 = 1.0;
    TH1D *hjp12[13];
    TH1D *hjp2[13];
    TH1D *hjp1[13];

    Double_t avgM[13] = {0};
    for (int i = 0; i < 13; i++)
    {
        hjp1[i] = (TH1D *)fhist->Get(Form("havgM_JP1_Bin%i", i));
        hjp2[i] = (TH1D *)fhist->Get(Form("havgM_JP2_Bin%i", i));
        hjp12[i] = (TH1D *)hjp2[i]->Clone();
        hjp12[i]->Add(hjp1[i], presJP1);
        avgM[i] = hjp12[i]->GetMean();
    }
    cout << "Double_t avgM[13] = {";
    for (int i = 0; i < 13; i++)
    {
        if (i < 12)
        {
            cout << avgM[i] << ",";
        }
        else
        {
            cout << avgM[i] << "};" << endl;
        }
    }
}