#!/bin/bash/
mkdir -p /star/u/pokhrel/IFF_Work/runIFFfinal/Results/TpcOrTof
mkdir -p /star/u/pokhrel/IFF_Work/runIFFfinal/Results/TpcOrTof/Minv
mkdir -p /star/u/pokhrel/IFF_Work/runIFFfinal/Results/TpcOrTof/pT
mkdir -p /star/u/pokhrel/IFF_Work/runIFFfinal/Results/TpcOrTof/Eta
mkdir -p /star/u/pokhrel/IFF_Work/runIFFfinal/Results/TpcOrTof/intMinv
mkdir -p /star/u/pokhrel/SCRATCH/subFiles



	#star-submit runAsymmetryVsEtaTpcOrTof.xml
	#star-submit runAsymmetryVspTTpcOrTof.xml
	#star-submit runAsymmetryVsMinvTpcOrTof.xml
	star-submit runIntAsymmetryVsMinvTpcOrTof.xml  
