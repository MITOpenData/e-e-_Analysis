DOCENTRAL=1
DOSTUDYVSDIJET=0
DOCOPYPLOTS=1

mix="0"
overwrite=1
VERBOSE=1

listsample=(1) #0=data, 1=mc
listetaselection=(4) #from this list of values (0 0.3 0.5 1.0 2.0) 
listgen=(0 1) #0=no gen selection, 1=gen selection
listaxis=(1) #0=beam, 1=thrust, 2=wta, 3=thrust perp, 4 =wta perp 
listthirdjet=(0 1 2 3) #from this list of values (0 0.05 0.1 0.3) 

################################################################
#### dont change anything below this if you dont know what you are doing #### 
################################################################

# list of axis choices
listthrust=(0 1 0 1 0)
listtwta=(0 0 1 0 1)
listtperp=(0 0 0 1 1)

listetacuts=(0 0.3 0.5 1.0 2.0) 
activateetacut=(0 1 1 1 1) 

listthirdjetcut=(0 0.05 0.1 0.3) 
activatethirdjetcut=(0 1 1 1) 
listactivateajrej=(0 1 1 1) 
listajrejcut=(0 0.1 0.1 0.1) 

#REGULAR ANALYSIS, CENTRAL VALUES
AINPUT=( "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180322/LEP1Data1992_recons_aftercut-MERGED.root" "/data/cmcginn/StudyMultSamples/ALEPH/MC/20180323/alephMCRecoAfterCutPaths_1994.root" )
AOUTPUT=( "LEP1Data1992" "LEP1MC1994_20180323" )

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

    echo thrust${1}_mix${2}_wta${3}_perp${4}_gen${5}_ajrej${6}_ajrejcut$(float_to_string ${7})_threejet${8}_threejetcut$(float_to_string ${9})_optionetasel${10}_etacut$(float_to_string ${11})
}




sleep .5 

if [ $DOCENTRAL -eq 1 ]; then       
  for isample in ${listsample[@]}
    do
    for iaxis in ${listaxis[@]}
      do
      thrust=${listthrust[$iaxis]}
      wta=${listtwta[$iaxis]}
      perp=${listtperp[$iaxis]}
      
      for gen in ${listgen[@]}
        do  
        if [ $isample -eq 0 -a $gen -eq 1  ]; then       
        continue
        fi

        for ietacut in ${listetaselection[@]}
          do 
          etacut=${listetacuts[$ietacut]}
          acivateetacut=${activateetacut[$ietacut]}
          
          if [ $acivateetacut -eq 1 ]; then netaselection=(1 2)
          else netaselection=(0)
          fi
    
          for etacutoption in ${netaselection[@]}
            do 
            
            for ithird in ${listthirdjet[@]}
            do 
              threejetcut=${listthirdjetcut[$ithird]}
              threejet=${activatethirdjetcut[$ithird]}
              ajrejcut=${listajrejcut[$ithird]}
              ajrej=${listactivateajrej[$ithird]}
        
              echo "loop"
              INPUTDATA=${AINPUT[$isample]}
              OUTPUT=${AOUTPUT[$isample]}
              suffix=${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut} ${etacutoption} ${etacut})
              sleep .5  
              OUTPUTROOT=rootfiles/${suffix}.root
              OUTPUTHISTO=rootfiles/2PC_${suffix}.root
              FOLDERPLOTS=plots/plots_${suffix}
              OUTPUTPLOTS=$FOLDERPLOTS/${suffix}
              echo $OUTPUTROOT
              rm $OUTPUTHISTO
              rm -rf $FOLDERPLOTS
              mkdir $FOLDERPLOTS 
              rm $OUTPUTROOT
              root -l -q -b "ridge_check.c+(\"$INPUTDATA\",\"$OUTPUTROOT\",\"$mix\","$overwrite","$thrust","$wta","$perp","$gen","$VERBOSE","${ajrej}","${ajrejcut}","${threejet}","${threejetcut}","${etacutoption}","${etacut}")"
              root -l -q -b "TPCPlots.cc+(\"$OUTPUTROOT\",\"$OUTPUTHISTO\",\"$OUTPUTPLOTS\")" 
              done # with three jets
            done # with eta cut option
          done # with eta selection
        done #done with gen 
      done #done with iaxis
    done #done with samples
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
  cd ..
fi 
