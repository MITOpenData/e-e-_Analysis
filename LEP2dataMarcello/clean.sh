sed 's/px=//' ALEPH_Data1998_189GeV_V0.list >> temp_1.txt
sed 's/py=//' temp_1.txt >> temp_2.txt
sed 's/pz=//' temp_2.txt >> temp_3.txt
sed 's/m=//' temp_3.txt >> temp_4.txt
sed 's/charge//' temp_4.txt >> temp_5.txt
sed 's/pwflag//' temp_5.txt >> temp_6.txt
awk '!/END_EVENT/' temp_6.txt >>temp_7.txt
sed 's/EVENT//' temp_7.txt >> temp_8.txt
sed 's/ECM =//' temp_8.txt >> temp_9.txt
sed 's/GEV//' temp_9.txt >> temp_10.txt
sed 's/ALEPH_DATA RUN =/   -999. -999. -999./' temp_10.txt >> cleaned_ALEPH_Data1998_189GeV_V0.txt
#rm temp_*

