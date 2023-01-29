#!/bin/bash

outdir="/star/u/pokhrel/GPFS/Run12EmbeddingTrees"
indir="/star/u/pokhrel/SCRATCH/Run12EmbeddingPerRun"
for run in `head -n 100 KevinsRunList.list`
do
    
    if [ -f $outdir/tree_"$run".root ]; then
        echo "File exist ---Exitiing"
        continue
    fi
    
    fCount1=$(ls $indir/$run/*.root | wc -l)
    fCount2=$(cat ./ptListPerRun/$run.list | wc -l)
    echo $run  $fCount1  $fCount2
    if [ $fCount1 -eq $fCount2 ]
    then
        echo "-----------Adding trees for run = $run ----------------"
        star-submit-template -template ./add_trees.xml -entities run=$run,odir=$outdir
    fi
    
done


