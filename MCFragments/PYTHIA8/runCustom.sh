#!/bin/bash

mkdir -p logs

for i in `seq 0 12`
do
    ./mainCustomZee_RopeWalk /data/cmcginn/GeneratorsHEPMC/PYTHIA8/0p0912/pthat1_Zee_RopeWalk/outFile_$i.root >& logs/log_$i.log &

    sleep 2
done