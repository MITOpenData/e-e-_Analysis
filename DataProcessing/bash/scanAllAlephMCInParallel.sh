#!/bin/bash

if [ $# -ne 1 ]
  then
    echo "scanAllAlephMCInParallel.sh takes 1 argument, dirpath for output. return"
    exit
fi

mkdir -p logs
DATE=`date +%Y%m%d`

echo "./bin/scan.exe paths/alephMCRecoAfterCutPaths_1997.txt $1/alephMCRecoAfterCut_1997_$DATE.root >& logs/alephMCRecoAfterCutPaths_1997.log &"

echo "./bin/scan.exe paths/alephMCRecoAfterCutPaths_1998.txt $1/alephMCRecoAfterCut_1998_$DATE.root >& logs/alephMCRecoAfterCutPaths_1998.log &"

echo "./bin/scan.exe paths/alephMCRecoAfterCutPaths_1999.txt $1/alephMCRecoAfterCut_1999_$DATE.root >& logs/alephMCRecoAfterCutPaths_1999.log &"
