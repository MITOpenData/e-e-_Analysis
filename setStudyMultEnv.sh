#!/bin/bash

echo "SETTING NEW ENVIRONMENT VARIABLE STUDYMULTDIR AS '$PWD'"
export STUDYMULTDIR=$PWD

export FASTJETDIR="/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install"
echo "SETTING NEW ENVIRONMENT VARIABLE FASTJETDIR AS '$FASTJETDIR'"
echo "   !!!IF YOU ARE NOTE ON LXPLUS, PLEASE POINT FASTJETDIR TO YOUR OWN FASTJET INSTALL!!!"
echo ""

echo "All environment variables set"