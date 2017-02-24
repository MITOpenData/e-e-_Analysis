while read F  ; do
        echo $F
          rm LEP2/cleaned_$F
          rm ROOTfiles/cleaned_$F.root
          sed 's/px=//' LEP2/$F.list >> temp_1.txt
          sed 's/py=//' temp_1.txt >> temp_2.txt
          sed 's/pz=//' temp_2.txt >> temp_3.txt
          sed 's/m=//' temp_3.txt >> temp_4.txt
          sed 's/charge//' temp_4.txt >> temp_5.txt
          sed 's/pwflag//' temp_5.txt >> temp_6.txt
          awk '!/END_EVENT/' temp_6.txt >>temp_7.txt
          awk '!/END_FILE/' temp_7.txt >>temp_8.txt
          sed 's/EVENT//' temp_8.txt >> temp_9.txt
          sed 's/ECM =//' temp_9.txt >> temp_10.txt
          sed 's/GEV//' temp_10.txt >> temp_11.txt
          sed 's/ALEPH_DATA RUN =/   -999. -999. -999./' temp_11.txt >> LEP2/cleaned_$F
          echo "-999. -999. -999. -999. -999. -999." >> LEP2/cleaned_$F

          rm temp_*

          g++ scan.cc $(root-config --cflags --libs) -g -o scan.exe 
          ./scan.exe cleaned_$F
          rm scan.exe
           
done <LEP2events_samples.list
