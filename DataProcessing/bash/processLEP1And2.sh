#!/bin/bash

topDir=/data/cmcginn/StudyMultSamples/ALEPH/
dirDate=20180329

mkdir -p logs
mkdir -p $topDir/LEP1/$dirDate
mkdir -p $topDir/LEP2/$dirDate

inputPath=$PWD/inputs

./bin/scan.exe $inputPath/in1992.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1992_recons_aftercut-MERGED.root >& logs/proc1992.log &
./bin/scan.exe $inputPath/in1993.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1993_recons_aftercut-MERGED.root >& logs/proc1993.log &
./bin/scan.exe $inputPath/in1994.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1994_recons_aftercut-MERGED.root >& logs/proc1994.log &
./bin/scan.exe $inputPath/in1995_LEP1.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1995_recons_aftercut-MERGED.root >& logs/proc1995.log &
./bin/scan.exe $inputPath/in1995_LEP2.txt 1 1 $topDir/LEP2/$dirDate/LEP2Data1995_recons_aftercut-MERGED.root >& logs/proc1995.log &
./bin/scan.exe $inputPath/in1996.txt 1 1 $topDir/LEP2/$dirDate/LEP2MC1996YDATAEMINI_recons_aftercut-MERGED.root >& logs/proc1996.log &
./bin/scan.exe $inputPath/in1997.txt 1 1 $topDir/LEP2/$dirDate/LEP2MC1997YDATAEMINI_recons_aftercut-MERGED.root >& logs/proc1997.log &
./bin/scan.exe $inputPath/in1998.txt 1 1 $topDir/LEP2/$dirDate/LEP2MC1998YDATAEMINI_recons_aftercut-MERGED.root >& logs/proc1998.log &
./bin/scan.exe $inputPath/in1999.txt 1 1 $topDir/LEP2/$dirDate/LEP2MC1999YDATAEMINI_recons_aftercut-MERGED.root >& logs/proc1999.log &
./bin/scan.exe $inputPath/in2000.txt 1 1 $topDir/LEP2/$dirDate/LEP2MC2000YDATAEMINI_recons_aftercut-MERGED.root >& logs/proc2000.log &
