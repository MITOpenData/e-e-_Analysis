#./run.sh "/data/cmcginn/StudyMultSamples/ALEPH/LEP1/20180322/LEP1Data1992_recons_aftercut-MERGED.root" "LEP1Data1992" ""  1 1 0 0 "results.root" "plots_LEP1Data1992_thrust1"

root -l -q -b "ridge_check.c+(\"$1\",\"$2\",\"$3\",$4,$5,$6,$7)"
root -l -q -b "TPCPlot.cc+(\"$2\",\"$8\",\"$9\")"