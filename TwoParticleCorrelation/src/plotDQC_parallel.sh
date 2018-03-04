#!/bin/sh

#  ridge_parallel.sh
#  
#
#  Created by Anthony Badea on 1/18/18.
#  
#how to run ./ridge.sh inFileName outFileName fileEvents
#!/bin/bash

inFileLoc=$1
outFileName=$2

find $inFileLoc -type f -name "*.root" >inFileList.txt
i=0
mkdir tempLoc_$outFileName

while read line
do
    root -b -q plotDQC.cc\(\"$line\",\"tempLoc_$outFileName/out_$i\"\) &
    let i=i+1
done <inFileList.txt
wait

rm -r inFileList.txt

find tempLoc_$outFileName -type f -name "*.root" > outFileList.txt
root -b -q histoHadd_TH1F.C\(\"outFileList.txt\",\"$outFileName\"\)
rm -r tempLoc_$outFileName
rm -r outFileList.txt



