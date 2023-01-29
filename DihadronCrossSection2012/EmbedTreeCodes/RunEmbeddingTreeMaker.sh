#!/bin/bash/
wrkDir="/star/u/pokhrel/DihadronCrossSection/EmbedTreeProduction/testDir"
ldir="/star/u/pokhrel/SCRATCH/Run12EmbeddingPerRun"
mkdir -p $ldir
# for ptbin in pt2_3 pt3_4 pt4_5 pt5_7 pt7_9 pt9_11 pt11_15 pt15_20 pt20_25 pt25_35 pt35_45 pt45_55 pt55_-1
for runNum in `sed -n 4,100p KevinsRunList.list`
do
    touch a$runNum.list
    touch b$runNum.list
    # tDir="/star/u/pokhrel/SCRATCH/Run12Embedding/$runNum/"
    tDir="$ldir/$runNum/"
    # locDir="/star/u/pokhrel/SCRATCH/Run12Embedding/$runNum/subFiles/"
    locDir="$ldir/$runNum/subFiles/"
    # outDir="/star/u/pokhrel/SCRATCH/Run12Embedding/$runNum/out/"
    outDir="$ldir/$runNum/out/"
    
    # mkdir -p /star/u/pokhrel/SCRATCH/Run12Embedding/$runNum
    mkdir -p $ldir/$runNum
    # mkdir -p /star/u/pokhrel/SCRATCH/Run12Embedding/$runNum/out
    mkdir -p $ldir/$runNum/out
    # mkdir -p /star/u/pokhrel/SCRATCH/Run12Embedding/$runNum/subFiles
    mkdir -p $ldir/$runNum/subFiles
    #mudst file list
    flist="$wrkDir/ptListPerRun/$runNum.list"
    #filename (just for a character length of filename)
    rofile="st_zerobias_adc_13055001_raw_2590001_r1.MuDst.root"
    #length of filename
    rolength=${#rofile}
    
    for fname in `cat $flist`
    do
        #get length of dir+filename
        length=${#fname}
        #get basename of mudst file
        fbasename=$(echo "$fname" | awk -F '/' '{print $12}' | awk -F '.' '{print $1}')
        #echo "$fbasename"
        #get length dir excluding filename
        dirlen=$(expr $length - $rolength)
        #echo "dirlength: $dirlen  rofile length : $rolength  file length: $length "
        #get directory pointing to the mudst file
        fDir=${fname:0:$dirlen}
        echo $fDir >> a$runNum.list
    done
    #sort the dir list. This step list dir for each unique run
    cat a$runNum.list | sort -u > b$runNum.list
    #loop over every run and submit jobs with wild card in xml for mudst and pythia files within run
    for rundir in `cat b$runNum.list`
    do
        ptbin_lo=$( echo $rundir | awk -F '_' '{print $4}')
        ptbin_up=$( echo $rundir | awk -F '_' '{print $5}')
        ptBin=$( echo $ptbin_lo'_'$ptbin_up )
        echo "pT bin ranges $ptbin_lo $ptbin_up  $ptBin"
        echo "---- Submitting job for:   $run -------- "
        
        star-submit-template -template ./RunEmbeddingTreeMaker.xml -entities runNum=$runNum,loDir=$locDir,oDir=$outDir,treeDir=$tDir,wDir=$wrkDir,runDir=$rundir,ptBin=$ptBin
    done
    #remove temp list files
    rm a$runNum.list
    rm b$runNum.list
    
done
