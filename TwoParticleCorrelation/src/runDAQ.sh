AINPUT=( "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180423/LEP1Data1993_recons_aftercut-MERGED.root" "/data/cmcginn/StudyMultSamples/ALEPH/MC/20180423/alephMCRecoAfterCutPaths_1994.root" "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180423/LEP1Data1993_recons_aftercut-MERGED_Mix.root" "/data/cmcginn/StudyMultSamples/ALEPH/MC/20180423/alephMCRecoAfterCutPaths_1994_mix.root" )
AOUTPUT=( "DataQualityLEP1Data1992" "DataQualityLEP1MC1994_20180423", "DataQualityLEP1Data1992_Mixed" "DataQualityLEP1MC1994_20180423_Mixed" )
NPUTDATAMIX="0"
overwrite=1
thrust=1
wta=0
perp=0
gen=0
VERBOSE=1
ajrej=0
ajrejcut=0
threejet=0
threejetcut=0
owbarrel=1
anatyperegion=2
etabarrelcut=2.0
typeEnergyBarrelSel=0
etabarrelcutforEselection=2.0
maxrelenergyinsidebarrel=0.2
typemultiplicity=1

mixedsample=0

listsample=(0 2) #0=data, 1=mc

############################### DONT MODIFY BELOW THIS LINE ###############################

listismixed=(0 0 1 1) 


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
    if [[ $# -ne 16 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - produce_postfix()"
        return 1
    fi

    echo thrust${1}_mixedsample${2}_wta${3}_perp${4}_gen${5}_ajrej${6}_ajrejcut$(float_to_string ${7})_threejet${8}_threejetcut$(float_to_string ${9})_owbarrel${10}_anatyperegion${11}_etabarrelcut$(float_to_string ${12})_typeEnergyBarrelSel${13}_etabarrelcutforEselection$(float_to_string ${14})_maxrelenergyinsidebarrel$(float_to_string ${15})_typemultiplicity${16}
}



  for isample in ${listsample[@]}
  do

    INPUTDATA=${AINPUT[$isample]}
    OUTPUT=${AOUTPUT[$isample]} 
    mixedsample=${listismixed[$isample]}

    suffix=${OUTPUT}_$(produce_postfix ${thrust} ${mixedsample} ${wta} ${perp} ${gen} ${ajrej} ${ajrejcut} ${threejet} ${threejetcut} ${owbarrel} ${anatyperegion} ${etabarrelcut} ${typeEnergyBarrelSel} ${etabarrelcutforEselection} ${maxrelenergyinsidebarrel} ${typemultiplicity})
    sleep .5  
    OUTPUTROOT=rootfilesDataQuality/${suffix}.root
    FOLDERPLOTS=plotsDataQuality/plots_${suffix}
    OUTPUTPLOTS=$FOLDERPLOTS/${suffix}
    echo $OUTPUTROOT
    rm -rf $FOLDERPLOTS
    mkdir $FOLDERPLOTS 
    rm $OUTPUTROOT

    root -l -q -b "DataQuality.c+(\"$INPUTDATA\",\"$OUTPUTROOT\",\"$INPUTDATAMIX\","${mixedsample}","$overwrite","$thrust","$wta","$perp","$gen","$VERBOSE","${ajrej}","${ajrejcut}","${threejet}","${threejetcut}","${owbarrel}","${anatyperegion}","${etabarrelcut}","${typeEnergyBarrelSel}","${etabarrelcutforEselection}","${maxrelenergyinsidebarrel}","${typemultiplicity}")"
  done
