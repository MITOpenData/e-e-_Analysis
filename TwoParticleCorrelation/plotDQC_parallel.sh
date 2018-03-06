#!/bin/sh

#  ridge_parallel.sh
#  
#
#  Created by Anthony Badea on 1/18/18.
#  
#how to run ./ridge.sh inFileName outFileName fileEvents
#!/bin/bash

if [ "$#" -ne 3 ]; 
then
    echo "Usage: sh plotDQC_parallel.sh <inDir> <outDir> <energyCut>"
    exit
fi

inFileLoc=$1
outFileName=$2
energyCut=$3

mkdir -p tempLoc_$outFileName
pos=0

maxRun=9

rm currentNumberOfJobs.txt

for i in $inFileLoc/*.root
do
    echo "./bin/plotDQC.exe $i tempLoc_$outFileName/out_$pos.root $energyCut >& tempLoc_$outFileName/out_$pos.log &"
    ./bin/plotDQC.exe $i tempLoc_$outFileName/out_$pos.root $energyCut >& tempLoc_$outFileName/out_$pos.log &
    pos=$((pos + 1))

    sleep 2

    ps | grep -i plotDQC | wc -l > currentNumberOfJobs.txt
    read lines < currentNumberOfJobs.txt

    while [ $lines -gt $maxRun ]
    do
        sleep 15
        ps | grep -i plotDQC | wc -l > currentNumberOfJobs.txt
        read lines < currentNumberOfJobs.txt
        echo $lines
    done    
done