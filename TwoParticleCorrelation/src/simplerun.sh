DOCENTRAL=0
DOSTUDYVSDIJET=0
DOCOPYPLOTS=1

mix="0"
overwrite=1
thrust=0
wta=1
perp=0
gen=0
VERBOSE=1
ajrej=0
threejet=0
ajrejcut=0
threejetcut=0
etathrustselection=2.0

#REGULAR ANALYSIS, CENTRAL VALUES
AINPUT=( "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180322/LEP1Data1992_recons_aftercut-MERGED.root" "/data/cmcginn/StudyMultSamples/ALEPH/MC/20180323/alephMCRecoAfterCutPaths_1994.root" )
AOUTPUT=( "LEP1Data1992" "LEP1MC1994_20180323" )

sleep .5 

if [ $DOCENTRAL -eq 1 ]; then       
  
  for i in 0 1
  do
    for optionetasel in 0 1 2
    do  
      INPUTDATA=${AINPUT[$i]}
      OUTPUT=${AOUTPUT[$i]}
      
      echo $INPUTDATA
      echo $OUTPUT
      suffix=${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut} ${optionetasel} ${etathrustselection})
    
      sleep .5 
  
      OUTPUTROOT=rootfiles/${suffix}.root
      OUTPUTHISTO=rootfiles/2PC_${suffix}.root
      FOLDERPLOTS=plots/plots_${suffix}.root
      OUTPUTPLOTS=$FOLDERPLOTS/${suffix}.root
      echo $OUTPUTROOT
  
      rm $OUTPUTHISTO
      rm -rf $FOLDERPLOTS
      mkdir $FOLDERPLOTS 
      rm $OUTPUTROOT
  
      root -l -q -b "ridge_check.c+(\"$INPUTDATA\",\"$OUTPUTROOT\",\"$mix\","$overwrite","$thrust","$wta","$perp","$gen","$VERBOSE","${ajrej}","${ajrejcut}","${threejet}","${threejetcut}","${optionetasel}","${etathrustselection}")"
      root -l -q -b "TPCPlots.cc+(\"$OUTPUTROOT\",\"$OUTPUTHISTO\",\"$OUTPUTPLOTS\")" 
    done 
  done
fi



multlow=( 0 20 30 35)
multhigh=( 20 30 999 999)



if [ $DOCOPYPLOTS -eq 1 ]; then      

  cd plots/
  rm -rf summary
  mkdir summary 
  cd summary


  for etaindex in 0 1 2 
  do 
  
    for indexmult in 0 1 2 3
    do
    folder=${multlow[$indexmult]}_${multhigh[$indexmult]}_$etaindex
    mkdir $folder
    cd $folder
    cp ../../plots*/*ASummary_0_${folder}* .
    rm *.png
    mkdir thrust
    cd thrust
    mv ../*thrust1*.pdf .
    cd ..
    mkdir wta
    cd wta
    mv ../*wta1*.pdf .
    cd ..
    cd ..
  done
  done  
fi
cd ..
 

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
    if [[ $# -ne 11 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - produce_postfix()"
        return 1
    fi

    echo thrust${1}_mix${2}_wta${3}_perp${4}_gen${5}_ajrej${6}_ajrejcut$(float_to_string ${7})_threejet${8}_threejetcut$(float_to_string ${9})_optionetasel${10}_etathrustselection$(float_to_string ${11})
}

