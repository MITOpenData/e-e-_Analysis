#!/bin/bash
#made from command line option via Anthony - to be fixed
for k in */; 
do 

 for d in "$k"*/; 
 do 

  for i in "$d"*.gz; 
  do 
   gunzip "$i"; 
  done; 

  for f in "$d"*.aleph; 
  do 
   sed -i 's/px=//' "$f"; 
   sed -i 's/py=//' "$f"; 
   sed -i 's/pz=//' "$f"; 
   sed -i 's/charge//' "$f";
   sed -i 's/m=//' "$f";
   sed -i 's/pwflag//' "$f";
   sed -i 's/pid//' "$f";
   sed -i 's/ECM=//' "$f";
   sed -i 's/GEV//' "$f";
   sed -i 's/EVENT//' "$f";
   sed -i 's/pname//' "$f";
   sed -i 's/MC_TRUE_BEFORE_CUT RUN =/   -999. -999. -999. -999./' "$f"; 
   sed -i 's/MC_TRUE_AFTER_CUT RUN =/   -999. -999. -999. -999./' "$f";
   sed -i 's/MC_RECO RUN =/   -999. -999. -999./' "$f";
   sed -i 's/ECM = //' "$f";
   awk '!/END_/' "$f" >> testfile.txt && mv testfile.txt "$f";
   awk '!/END_FILE/' "$f" >> testfile.txt && mv testfile.txt "$f";
  done;
 done;
done;

