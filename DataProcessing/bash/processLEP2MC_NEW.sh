#!/bin/bash

mkdir -p logs
dir=/data/cmcginn/StudyMultSamples/ALEPH/MC/20180402
mkdir -p $dir

#subDirs=(newPathsMC/alephMCRecoAfterCutPaths_1994.txt newPathsMC/alephMCRecoAfterCutPaths_1997_KQQorKWW4F_0.txt newPathsMC/alephMCRecoAfterCutPaths_1998_KQQorKWW4F_0.txt newPathsMC/alephMCRecoAfterCutPaths_1999_KQQorKWW4F_2.txt)

#for i in "${subDirs[@]}"

maxJobs=7

for i in $PWD/inputs/newPathsMC/*.txt
do
    tempName=${i%.txt}
    tempName=${tempName#$PWD/inputs/newPathsMC/}
    echo $i $tempName

    $PWD/bin/scan.exe $i 1 1 $dir/$tempName.root >& logs/$tempName.log &

#    exit

    ps | grep -i scan | wc -l > currentNumberOfJobs.txt
    read lines < currentNumberOfJobs.txt    
    echo $lines
    
    while [ $lines -gt $maxJobs ]
    do
        sleep 60
        ps | grep -i scan | wc -l > currentNumberOfJobs.txt
        read lines < currentNumberOfJobs.txt
        echo $lines
    done
done

