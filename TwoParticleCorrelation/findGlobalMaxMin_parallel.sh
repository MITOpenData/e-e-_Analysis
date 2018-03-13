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
    echo "Usage: sh findGlobalMaxMin.sh <inDir> <outDir>"
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
    echo "./bin/findGlobalMaxMin.exe $i tempLoc_$outFileName/out_$pos.root >& tempLoc_$outFileName/out_$pos.txt &"
    ./bin/findGlobalMaxMin.exe $i >& tempLoc_$outFileName/out_$pos.txt &
    pos=$((pos + 1))

    sleep 2

    ps | grep -i findGlobalMaxMin | wc -l > currentNumberOfJobs.txt
    read lines < currentNumberOfJobs.txt

    while [ $lines -gt $maxRun ]
    do
        sleep 15
        ps | grep -i findGlobalMaxMin | wc -l > currentNumberOfJobs.txt
        read lines < currentNumberOfJobs.txt
        echo $lines
    done    
done

for i in tempLoc_$outFileName/*.txt
do
    echo $i;
    while read p; 
    do 
        echo $p; 
    done <$i
done



