#!/bin/bash/

#create all the output directories
mkdir -p /star/u/pokhrel/GPFS/IFFResults_AfterPrelim/Dihad
mkdir -p /star/u/pokhrel/GPFS/IFFResults_AfterPrelim/Dihad/Minv
mkdir -p /star/u/pokhrel/GPFS/IFFResults_AfterPrelim/Dihad/pT
mkdir -p /star/u/pokhrel/GPFS/IFFResults_AfterPrelim/Dihad/Eta
mkdir -p /star/u/pokhrel/GPFS/IFFResults_AfterPrelim/Dihad/intMinv
mkdir -p /star/u/pokhrel/SCRATCH/subFiles
	
star-submit runAsymmetryVsEtaDihad.xml
star-submit runAsymmetryVspTDihad.xml
star-submit runAsymmetryVsMinvDihad.xml
star-submit runIntAsymmetryVsMinvDihad.xml  
