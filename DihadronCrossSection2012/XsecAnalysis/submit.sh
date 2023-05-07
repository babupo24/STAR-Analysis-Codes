#!/bin/bash/

outDir="/star/u/pokhrel/SCRATCH/histDatEmbed/DatHist"
mkdir -p $outDir

treeDir="/star/u/pokhrel/GPFS/Run12PicoTrees2"

for run in `sed -n 2,601p KevinsRunList.list`
do
    histDir=$outDir/$run
    mkdir $histDir
    
    logDir=$histDir/logFiles
    mkdir -p $logDi
    
    subDir=$histDir/subFiles
    mkdir -p $subDir
    
    star-submit-template -template ./Run12pp200Ana.xml -entities runnumber=$run,treeDir=$treeDir,logDir=$logDir,subDir=$subDir,histDir=$histDir
done
