#!/bin/bash

FULL_PATH=/data/flowex/Datasamples/LEP1_MAIN/LEP1/1991/

cat $1 | while read line
do

suffix=`basename $line`
echo $suffix
suffix_parsed=`echo $suffix|sed "s@\.aleph@\.root@g"`
echo $suffix_parsed
output="/data/abaty/EpEmStudies/ALEPH_ThrustReprocessing_Sept262017/"$suffix_parsed
./bin/scan.exe $line $output

done
