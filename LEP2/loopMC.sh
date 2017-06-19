folder="LEP2_MC/GGUD"
filelist="LEP2MCGGUDY1997E183_mctrue_beforecut-001.list"
output="cleaned_LEP2MCGGUDY1997E183_mctrue_beforecut-001.aleph"
suffix="LEP2MCGGUDY1997E183_mctrue_beforecut-001"


rm $filelist
rm $output
rm ROOTfiles/$output.root

cd $folder
rm cleaned*
ls $suffix* >> $filelist
cd ../../
mv $folder/$filelist $filelist

while read F  ; do
        echo $F
          rm $folder/cleaned_$F
          rm ROOTfiles/cleaned_$F.root
          cp $folder/$F  $folder/cleaned_$F
          sed -i '' 's/MC_TRUE_BEFORE_CUT RUN = / -999. -999. -999. -999./' $folder/cleaned_$F
          sed -i '' 's/px=//'  $folder/cleaned_$F
          sed -i '' 's/py=//'  $folder/cleaned_$F
          sed -i '' 's/pz=//'  $folder/cleaned_$F
          sed -i '' 's/charge//'  $folder/cleaned_$F
          sed -i '' 's/m=//'  $folder/cleaned_$F
          sed -i '' 's/pwflag//'  $folder/cleaned_$F
          sed -i '' 's/pid//'  $folder/cleaned_$F
          sed -i '' 's/ECM =//'  $folder/cleaned_$F
          sed -i '' 's/GEV//'  $folder/cleaned_$F
          sed -i '' 's/EVENT//'  $folder/cleaned_$F
          ex -sc '%s/\(\pname\).*/\1/ | x' $folder/cleaned_$F
          sed -i '' 's/pname//'  $folder/cleaned_$F
          sed -i '' 's/MC_TRUE_AFTER_CUT RUN =/   -999. -999. -999. -999./'  $folder/cleaned_$F
          awk '!/END_/' $folder/cleaned_$F >> testfile.txt && mv testfile.txt $folder/cleaned_$F


           
done <$filelist


cat $folder/cleaned_$suffix* >> $output
echo "-999. -999. -999. -999. -999. -999." >> $output

g++ scanMC.cc $(root-config --cflags --libs) -g -o scanMC.exe 
./scanMC.exe $output
rm scanMC.exe

