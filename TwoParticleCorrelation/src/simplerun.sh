DOCENTRAL=1
DOSTUDYVSDIJET=1
DOCOPYPLOTS=1

INPUTDATA="/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180322/LEP1Data1992_recons_aftercut-MERGED.root"
OUTPUT="LEP1Data1992"
mix="0"
overwrite=1
thrust=1
wta=0
perp=0
gen=0
VERBOSE=1
ajrej=0
threejet=0
ajrejcut=0
threejetcut=0

#REGULAR ANALYSIS, CENTRAL VALUES

if [ $DOCENTRAL -eq 1 ]; then      

  OUTPUTROOT=rootfiles/${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut}).root
  OUTPUTHISTO=rootfiles/2PC_${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut}).root
  FOLDERPLOTS=plots/plots_${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut})
  OUTPUTPLOTS=$FOLDERPLOTS/${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut})
  echo $OUTPUTROOT
  
  rm $OUTPUTHISTO
  rm -rf $FOLDERPLOTS
  mkdir $FOLDERPLOTS 
  rm $OUTPUTROOT
  
  root -l -q -b "ridge_check.c+(\"$INPUTDATA\",\"$OUTPUTROOT\",\"$mix\","$overwrite","$thrust","$wta","$perp","$gen","$VERBOSE","${ajrej}","${ajrejcut}","${threejet}","${threejetcut}")"
  root -l -q -b "TPCPlots.cc+(\"$OUTPUTROOT\",\"$OUTPUTHISTO\",\"$OUTPUTPLOTS\")" 
  
fi

ajrej=1
threejet=1  
ajrejcut=0.1


#STUDY VS ASYMMETRY
if [ $DOSTUDYVSDIJET -eq 1 ]; then      

  for threejetcut in 0.02 0.03 0.05 0.1 0.2 0.5
  do 

   OUTPUTROOT=rootfiles/${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut}).root
   OUTPUTHISTO=rootfiles/2PC_${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut}).root
   FOLDERPLOTS=plots/plots_${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut})
   OUTPUTPLOTS=$FOLDERPLOTS/${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut})
    
    echo $OUTPUTROOT
    rm $OUTPUTHISTO
    rm -rf $FOLDERPLOTS
    mkdir $FOLDERPLOTS 
    rm $OUTPUTROOT
    root -l -q -b "ridge_check.c+(\"$INPUTDATA\",\"$OUTPUTROOT\",\"$mix\","$overwrite","$thrust","$wta","$perp","$gen","$VERBOSE","${ajrej}","${ajrejcut}","${threejet}","${threejetcut}")"
    root -l -q -b "TPCPlots.cc+(\"$OUTPUTROOT\",\"$OUTPUTHISTO\",\"$OUTPUTPLOTS\")" 
  done 
fi


if [ $DOCOPYPLOTS -eq 1 ]; then      

  cd plots/
  rm -rf summary
  mkdir summary 
  cd summary

  for etaindex in 0 1 2 
  do 
    mkdir 0_20_$etaindex
    cd 0_20_$etaindex
    cp ../../plots*/*ASummary_0_0_20_$etaindex* .
    rm *.pdf
    cd ..
  done
  
  for etaindex in 0 1 2 
  do 
    mkdir 20_30_$etaindex
    cd 20_30_$etaindex
    cp ../../plots*/*ASummary_0_20_30_$etaindex* .
    rm *.pdf
    cd ..
  done

  for etaindex in 0 1 2 
  do 
    mkdir 30_999_$etaindex
    cd 30_999_$etaindex
    cp ../../plots*/*ASummary_0_30_999_$etaindex* .
    rm *.pdf
    cd ..
  done

  for etaindex in 0 1 2 
  do 
    mkdir 35_999_$etaindex
    cd 35_999_$etaindex
    cp ../../plots*/*ASummary_0_35_999_$etaindex* .
    rm *.pdf
    cd ..
  done

  
  cd ..
fi


function float_to_string()
{
    if [[ $# -ne 1 ]]
    then
        echo -e "${ERRCOLOR}error:${NC} invalid argument number - float_to_string()"
        return 1
    fi
    part1=`echo $1 | awk -F "." '{print $1}'`
    part2=`echo $1 | awk -F "." '{print $2}'`
    rt_float_to_string=${part1:-0}p${part2:-0}
    echo $rt_float_to_string
}

function produce_postfix()
{
    if [[ $# -ne 9 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - produce_postfix()"
        return 1
    fi

    echo thrust${1}_mix${2}_wta${3}_perp${4}_gen${5}_ajrej${6}_ajrejcut$(float_to_string ${7})_threejet${8}_threejetcut$(float_to_string ${9})
}

