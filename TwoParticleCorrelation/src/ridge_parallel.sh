#!/bin/sh

#  ridge_parallel.sh
#  
#
#  Created by Anthony Badea on 1/18/18.
#  
#how to run ./ridge.sh inFileName outFileName fileEvents
#!/bin/bash
NUMFILES=2
mkdir -p tempFiles
cd tempFiles
i=1

# Always round up number of events so we don’t miss events
# To do rounding up in truncating arithmetic, simply add (denom-1) to the numerator
# To do round-to-nearest, add (denom/2) to the numerator (halves will round up)
# we do rounding down here because it is ok to miss a few events but techinically wrong to double count events
let NUMEVENTS=($fileEvents+$NUMFILES)/$NUMFILES

START=1
FINISH=$NUMEVENTS
# divide up the file
while [ $i -le $NUMFILES ]; do
rooteventselector —recreate -f START -l FINISH inFileName:t $i.root
let i=i+1
let START=FINISH+1
let FINISH=$NUMEVENTS*$i
done 
# Divide the file into 20 subfiles we can now run ridge check on them

cd ..
# we are now in TwoParticleCorrelation
for f in tempFiles/*.root; do
ridge_check_parallel.exe "$f" “out_$f” &
done
wait

# delete the temp files
rm -r tempFiles

# hadd the files together
cd ../pdfDir
hadd -f $outFileName.root out_{1..20}.root
# delete the smaller files
rm -r out_*.root
