#!/bin/bash

mkdir -p logs

for i in `seq 40 49`
do
    ./bin/mainCustomZpole_ppAndee.exe PYTHIA8 0 $i 100000 0 0 0 >& logs/log_$i.log &

    sleep 2
done