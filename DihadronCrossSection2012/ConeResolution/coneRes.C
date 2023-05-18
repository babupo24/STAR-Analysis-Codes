#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include <iomanip>

void coneRes()
{
    // Open the ROOT file containing the event records
    // TFile *file = new TFile("/star/u/pokhrel/GPFS/Run12EmbeddingTrees/v3/tree_13047003.root", "READ");
    TFile *file = new TFile("/star/u/pokhrel/DihadronCrossSection/EmbedTreeProduction/TreeProductionTest/embedTestTree.root", "READ");

    // Get the "ftree" TTree from the file
    TTree *ftree = (TTree *)file->Get("ftree");

    // Set branch addresses for relevant variables
    Int_t fmaxpar1;
    Int_t fpId_pyth[10000];
    Int_t fStatusCode_pyth[10000];
    Int_t fFirstMotherId_pyth[10000];
    ftree->SetBranchAddress("fmaxpar1", &fmaxpar1);
    ftree->SetBranchAddress("fpId_pyth", fpId_pyth);
    ftree->SetBranchAddress("fStatusCode_pyth", fStatusCode_pyth);
    ftree->SetBranchAddress("fFirstMotherId_pyth", fFirstMotherId_pyth);
    int nentries = ftree->GetEntries();
    cout << "Entries: " << ftree->GetEntries() << endl;

    // Loop over all events in the "ftree" TTree
    for (int ievt = 0; ievt < 10; ievt++)
    {
        ftree->GetEntry(ievt);
        cout << "Event Id = " << ievt << "\n";
        cout << "id\tmid\tmstatus \tmpdg\tdpdg\tdstatus" << endl;
        // Loop over all particles in the final state of the event
        for (int ipart = 0; ipart < fmaxpar1; ipart++)
        {
            int pdgid = fpId_pyth[ipart];
            int status = fStatusCode_pyth[ipart];

            // Select only charged pions with status code 1 (stable particles)
            if (abs(pdgid) == 211 && status == 1)
            {
                int mother_index = fFirstMotherId_pyth[ipart];
                int mother_pdgid = fpId_pyth[mother_index];
                int mother_status = fStatusCode_pyth[mother_index];
                // cout << ipart << " \t " << mother_index << " \t " << mother_status << " \t " << mother_pdgid << " \t " << pdgid << " \t " << fStatusCode_pyth[ipart] << endl;

                // Trace the production history of the pion back to the hard-scattered mother parton
                while (mother_pdgid != 21 && mother_pdgid != 1 && mother_pdgid != 2 && mother_pdgid != 3 && mother_pdgid != 4)
                {
                    mother_index = fFirstMotherId_pyth[mother_index];
                    mother_pdgid = fpId_pyth[mother_index];
                }
                /*
                    // Loop over all other pions in the event and check if they share the same mother parton
                    for (int jpart = 0; jpart < fmaxpar1; jpart++)
                    {
                        if (jpart == ipart)
                            continue; // Skip the current pion

                        int jpdgid = fpId_pyth[jpart];
                        int jstatus = fStatusCode_pyth[jpart];
                        if (abs(jpdgid) == 211 && jstatus == 1)
                        {
                            int jmother_index = fFirstMotherId_pyth[jpart];
                            int jmother_pdgid = fpId_pyth[jmother_index];

                            // Trace the production history of the other pion back to the hard-scattered mother parton
                            while (jmother_pdgid != 21 && jmother_pdgid != 1 && jmother_pdgid != 2 && jmother_pdgid != 3 && jmother_pdgid != 4)
                            {
                                jmother_index = fFirstMotherId_pyth[jmother_index];
                                jmother_pdgid = fpId_pyth[jmother_index];
                            }

                            // If the two pions share the same mother parton, they are part of the same dipion pair
                            if (jmother_index == mother_index)
                            {
                                std::cout << "Dipion pair found: (" << pdgid << "," << jpdgid << ")" << std::endl;
                            }
                        }
                    }
                    */
            }
        }
    }

    // Clean up
    file->Close();
}
