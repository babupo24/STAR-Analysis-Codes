#!/bin/bash
	dir=/star/data05/scratch/pokhrel/IFF_TREES/StTOF
	#echo "Deleting previous folders in scratch....."
	#rm -rf $dir
	#echo "Deleting done......Now submitting jobs....."
	#mkdir $dir
        report=$dir/subFiles 
	mkdir $dir/ERR_OUT
	mkdir $report
	
for runnumber in `sed -n 4,200p ./MyGoldenRuns.list`; do
	echo "Submitting jobs for : $runnumber"
	outdir=$dir/$runnumber
	eodir=$dir/ERR_OUT/$runnumber
	mkdir -p $outdir
	mkdir -p $eodir
	mkdir -p $eodir/out
	mkdir -p $eodir/err
star-submit-template -template ./run15pp200.xml -entities runnumber=$runnumber,odir=$outdir,eodir=$eodir,rdir=$report
echo "-----------Done with $runnumber -----------" 
done
echo "Job Submission Completed!!!!"
#echo  "fills 18749 (16064077 - 16065018)for YELLOW and 18764 (16069024 - 16069042) for BLUE"
