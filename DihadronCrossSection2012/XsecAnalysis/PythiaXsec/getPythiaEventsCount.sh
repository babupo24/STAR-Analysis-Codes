pt_bin=(2_3 3_4 4_5 5_7 7_9 9_11 11_15 15_20 20_25 25_35 35_45 45_55 55_-1)
#path="/star/embed/embedding/pp200_production_2012/v2/Pythia6_pt2_3_100_20212001/P12id.SL12d/2012/*/*/*.fzd.log.gz"
#nFiles=$(ls $path | wc -l)
#echo $nFiles
fall="EventCounts/pythia_gen_tried_allpT.txt"
touch $fall 
for pt in "${pt_bin[@]}"
do
    sum_gen=0
    sum_tried=0
    fname="pt"$pt"_pythia_gen_tried.txt"
    echo $pt
    #path="/star/embed/embedding/pp200_production_2012/v2/Pythia6_pt"$pt"_100_20212001/P12id.SL12d/2012/*/*/*.fzd.log.gz"
    #zgrep "All included subprocesses" $path | awk -F ' ' '{print $8 " " $9}'>$fname
    
    while read -r line
    do
        read n_gen n_tried <<< $line
        sum_gen=$((sum_gen + n_gen));
        sum_tried=$((sum_tried + n_tried));

    done < "EventCounts/$fname"
    echo $pt $sum_gen $sum_tried>>$fall
done

