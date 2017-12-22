#!/bin/bash

for i in `seq 2 5`
do
    echo "Processing 199$i..."
    for j in "199$i"/*.gz
    do
	cp $j .
	gunzip *.gz

	sleep 1

	mv ./*.aleph /data/flowex/Datasamples/LEP1_REPROCESS_20171221/LEP1/199$i/

	rm -f *.gz

	sleep 1
    done
done

for i in `seq 6 9`
do
    echo "Processing 199$i..."
    for j in "199$i"/*.gz
    do
	cp $j .
	gunzip *.gz

	sleep 1

	mv *.aleph /data/flowex/Datasamples/LEP2_REPROCESS_20171221/LEP2/199$i/

	rm -f *.gz

	sleep 1
    done
done


for j in 2000/*.gz
do
    cp $j .
    gunzip *.gz
    
    sleep 1
    
    mv *.aleph /data/flowex/Datasamples/LEP2_REPROCESS_20171221/LEP2/2000/

    rm -f *.gz
    
    sleep 1
done
