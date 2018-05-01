#!/bin/bash

topDir=/data/cmcginn/StudyMultSamples/ALEPH/
dirDate=20180430

mkdir -p logs
mkdir -p $topDir/LEP1/$dirDate

inputPath=$PWD/inputs
exePath=$PWD/bin

$exePath/scan.exe $inputPath/in1992_Single.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1992_recons_aftercut-MERGED.root >& logs/proc1992.log &
$exePath/scan.exe $inputPath/in1993.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1993_recons_aftercut-MERGED.root >& logs/proc1993.log &
$exePath/scan.exe $inputPath/in1994_part1.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1994P1_recons_aftercut-MERGED.root >& logs/proc1994P1.log &
$exePath/scan.exe $inputPath/in1994_part2.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1994P2_recons_aftercut-MERGED.root >& logs/proc1994P2.log &
$exePath/scan.exe $inputPath/in1994_part3.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1994P3_recons_aftercut-MERGED.root >& logs/proc1994P3.log &
$exePath/scan.exe $inputPath/in1995_LEP1.txt 1 1 $topDir/LEP1/$dirDate/LEP1Data1995_recons_aftercut-MERGED.root >& logs/proc1995.log &
