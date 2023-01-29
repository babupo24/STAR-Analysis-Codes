#!/bin/bash
starver SL21d #pp200_production_2012 MuDst->PicoDst production library
SCRATCH=/star/u/pokhrel/SCRATCH/Run12PicoTrees
mkdir -p $SCRATCH
for run in `sed -n 401,601p KevinsRunList.list`; do
    # for run in 13055010; do
    echo $run
    outdir1=$SCRATCH/$run
    mkdir -p $outdir1
    logdir=$outdir1/Log
    mkdir -p $logdir
    # outdir2=$GPFS/$run
    # mkdir -p $outdir2
    star-submit-template -template ./Run12ppTreeMaker.xml -entities runnumber=$run,logdir=$logdir,scratch=$outdir1
    
done





