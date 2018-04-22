INPUTDATA="/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180402YJTest/LEP1Data1993_recons_aftercut-MERGED.root"
OUTPUTROOT="test.root"
INPUTDATAMIX="0"
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
typeEnergyBarrelSel=1
etabarrelcutforEselection=2.0
maxrelenergyinsidebarrel=0.2
typemultiplicity=1

root -l -q -b "DataQuality.c+(\"$INPUTDATA\",\"$OUTPUTROOT\",\"$INPUTDATAMIX\","$overwrite","$thrust","$wta","$perp","$gen","$VERBOSE","${ajrej}","${ajrejcut}","${threejet}","${threejetcut}","${owbarrel}","${anatyperegion}","${etabarrelcut}","${typeEnergyBarrelSel}","${etabarrelcutforEselection}","${maxrelenergyinsidebarrel}","${typemultiplicity}")"

