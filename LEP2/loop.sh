folder="/data/flowex/Datasamples/LEP2_MAIN/LEP2"
filelist="/data/flowex/Datasamples/LEP2_MAIN/Lists/ALEPH_Data.list"
output="cleaned_ALEPH_Data-all.aleph"
suffix="ALEPH_Data"

rm $filelist
rm $output
rm ROOTfiles/$output.root

cd $folder
rm cleaned*
ls $suffix* >> $filelist
cd ..
mv $folder/$filelist $filelist

while read F  ; do
        echo $F
          rm $folder/cleaned_$F
          rm ROOTfiles/cleaned_$F.root
          cp $folder/$F  $folder/cleaned_$F

          sed -i ''  's/px=//' $folder/cleaned_$F
          sed -i ''  's/py=//' $folder/cleaned_$F
          sed -i ''  's/pz=//' $folder/cleaned_$F
          sed -i ''  's/m=//' $folder/cleaned_$F
          sed -i ''  's/charge//' $folder/cleaned_$F
          sed -i ''  's/pwflag//' $folder/cleaned_$F
          awk '!/END_EVENT/' $folder/cleaned_$F >> testfile.txt && mv testfile.txt $folder/cleaned_$F
          awk '!/END_FILE/' $folder/cleaned_$F >> testfile.txt && mv testfile.txt $folder/cleaned_$F
          sed -i ''  's/EVENT//' $folder/cleaned_$F
          sed -i ''  's/ECM =//' $folder/cleaned_$F
          sed -i ''  's/GEV//' $folder/cleaned_$F
          sed -i ''  's/ALEPH_DATA RUN =/   -999. -999. -999./' $folder/cleaned_$F
           
done <$filelist


cat $folder/cleaned_$suffix* >> $output
echo "-999. -999. -999. -999. -999. -999." >> $output

g++ scan.cc $(root-config --cflags --libs) -g -o scan.exe 
./scan.exe $output
rm scan.exe
