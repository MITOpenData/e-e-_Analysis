INPUTDATA="/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180322/LEP1Data1992_recons_aftercut-MERGED.root"
OUTPUT="LEP1Data1992_thrust.root"
MIXING=""
OVERWRITE=1
OVERWRITETHRUST=1
OVERWRITEWTA=0
OVERWRITEDOGEN=0

root -l -q -b "ridge_check.c+(\"$INPUTDATA\",\"$OUTPUT\",\"$MIXING\",\"$OVERWRITE\",\"$OVERWRITETHRUST\",0,0)"
root -l -q -b "TPCPlots.cc+(\"$OUTPUT\")"

