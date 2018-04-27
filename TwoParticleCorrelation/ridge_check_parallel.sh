#!/bin/sh

#  ridge_parallel.sh
#  
#
#  Created by Anthony Badea on 1/18/18.
#  
#how to run ./ridge.sh inFileName outFileName fileEvents
#!/bin/bash

if [ "$#" -ne 2 ]; 
then
    echo "Usage: sh ridge_check_parallel.sh <inDir> <outDir>"
    exit
fi

inFileLoc=$1
outFileName=$2

mkdir -p tempLoc_$outFileName
pos=0

maxRun=9

rm currentNumberOfJobs.txt

for i in $inFileLoc/*.root
do
    echo "./bin/ridge_check.exe $i tempLoc_$outFileName/out_$pos.root >& tempLoc_$outFileName/out_$pos.log &"
    ./bin/ridge_check.exe $i tempLoc_$outFileName/out_$pos.root >& tempLoc_$outFileName/out_$pos.log &
    pos=$((pos + 1))

    sleep 2

    ps | grep -i ridge_check | wc -l > currentNumberOfJobs.txt
    read lines < currentNumberOfJobs.txt

    while [ $lines -gt $maxRun ]
    do
        sleep 15
        ps | grep -i ridge_check | wc -l > currentNumberOfJobs.txt
        read lines < currentNumberOfJobs.txt
        echo $lines
    done    
done
