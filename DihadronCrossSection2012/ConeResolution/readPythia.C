#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "StRoot/StJetSkimEvent/StPythiaEvent.h"

using namespace std;

int readPythia()
{
    // Open the ROOT file containing the event records
    TFile *file = new TFile("/star/embed/embedding/pp200_production_2012/v2/Pythia6_pt11_15_100_20212001/P12id.SL12d/2012/047/13047003/st_zerobias_adc_13047003_raw_0570001_r0.pythia.root", "READ");

    TFile *fout = new TFile("ttree.root", "recreate");
    // Get the TTree from the file and set the branch address for StPythiaEvent
    TTree *tree = (TTree *)file->Get("PythiaEvent");
    StPythiaEvent *pythiaEvent = new StPythiaEvent();
    tree->SetBranchAddress("Pythia.", &pythiaEvent);

    // Loop over all events in the TTree
    for (int ievt = 0; ievt < tree->GetEntries(); ievt++)
    {
        tree->GetEntry(ievt);

        // Get the particle array for this event
        TClonesArray *particles = pythiaEvent->Particles();

        // Loop over all particles in the event
        for (int ipart = 0; ipart < particles->GetEntriesFast(); ipart++)
        {
            TParticle *particle = (TParticle *)particles->At(ipart);

            // Access particle information
            int pdgid = particle->GetPdgCode();
            double px = particle->Px();
            double py = particle->Py();
            double pz = particle->Pz();
            double e = particle->Energy();

            // Access mother particle information
            int motherIndex = particle->GetMother(0);
            if (motherIndex >= 0)
            {
                TParticle *motherParticle = (TParticle *)particles->At(motherIndex);
                int motherPdgid = motherParticle->GetPdgCode();
                double motherPx = motherParticle->Px();
                double motherPy = motherParticle->Py();
                double motherPz = motherParticle->Pz();
                double motherE = motherParticle->Energy();
            }
        }
    }
    fout->Write();
    fout->Close();
    // Close the ROOT file
    file->Close();
    delete file;

    return 0;
}
