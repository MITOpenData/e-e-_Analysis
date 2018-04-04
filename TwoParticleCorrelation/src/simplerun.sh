INPUTDATA="/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180322/LEP1Data1992_recons_aftercut-MERGED.root"
OUTPUT="LEP1Data1992"
mix="0"
overwrite=1
thrust=1
wta=0
perp=0
gen=0
VERBOSE=1
threejetrej=1

for threejetcut in 0.0 0.1 1.0 10.0
do 

OUTPUTROOT=rootfiles/${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${threejetrej} ${threejetcut}).root
OUTPUTHISTO=rootfiles/2PC_${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${threejetrej} ${threejetcut}).root
FOLDERPLOTS=plots/plots_${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${threejetrej} ${threejetcut})
OUTPUTPLOTS=$FOLDERPLOTS/${OUTPUT}_$(produce_postfix ${thrust} ${mix} ${wta} ${perp} ${gen} ${threejetrej} ${threejetcut})
echo $OUTPUTROOT


rm $OUTPUTHISTO
rm -rf $FOLDERPLOTS
mkdir $FOLDERPLOTS 
rm $OUTPUTROOT
./runHisto.sh $INPUTDATA $OUTPUTROOT $mix $overwrite $thrust $wta $perp $gen $VERBOSE $threejetrej $threejetcut
./runplot.sh $OUTPUTROOT $OUTPUTHISTO $OUTPUTPLOTS 

done 

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
    if [[ $# -ne 7 ]]
    then
        echo -e "\033[1;31merror:${NC} invalid argument number - produce_postfix()"
        return 1
    fi

    echo thrust${1}_mix${2}_wta${3}_perp${4}_gen${5}_threejetrej${6}_threejetcut$(float_to_string ${7})
}
