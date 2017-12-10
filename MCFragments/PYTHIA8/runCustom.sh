#!/bin/bash

mkdir -p logs

for i in `seq 0 9`
do
    ./mainCustomZee /data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/20171208/outFile_$i.root 2000000 0 1 >& logs/log_$i.log &

    sleep 2
done