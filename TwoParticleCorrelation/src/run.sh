#\"$8\",\"$9\" are respectively the directory and the prefix of the plots
#e.g. "myplots" and "plots_LEP1Data1992_thrust1_mix0_wta0_perp0_gen0_threejetrej1_threejetcut0p03"

root -l -q -b "ridge_check.c+(\"$1\",\"$2\",\"$3\",$4,$5,$6,$7)"
root -l -q -b "TPCPlot.cc+(\"$2\",\"$8\",\"$9\")"