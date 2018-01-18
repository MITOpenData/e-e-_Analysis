{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf200
{\fonttbl\f0\fnil\fcharset0 Menlo-Regular;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red255\green255\blue255;\red0\green0\blue0;
\red0\green116\blue0;\red0\green0\blue0;\red196\green26\blue22;\red28\green0\blue207;}
{\*\expandedcolortbl;;\csgray\c0;\csgray\c100000;\csgenericrgb\c0\c0\c0;
\csgenericrgb\c0\c45600\c0;\cssrgb\c0\c0\c0;\csgenericrgb\c77000\c10200\c8600;\csgenericrgb\c11000\c0\c81000;}
\margl1440\margr1440\vieww23580\viewh12740\viewkind0
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0

\f0\fs22 \cf2 \cb3 \CocoaLigature0 #how to run ./ridge.sh \cf4 \CocoaLigature1 inFileName\cf2 \CocoaLigature0  outFileName fileEvents\
\pard\tx543\pardeftab543\pardirnatural\partightenfactor0
\cf5 \CocoaLigature1 #!/bin/bash\
\pard\tx560\tx1120\tx1680\tx2240\tx2800\tx3360\tx3920\tx4480\tx5040\tx5600\tx6160\tx6720\pardirnatural\partightenfactor0
\cf2 \CocoaLigature0 NUMFILES = 20\
mkdir -p tempFiles\
cd tempFiles\
i = 1\
\
# Always round up number of events so we don\'92t miss events\
# To do rounding up in truncating arithmetic, simply add (denom-1) to the numerator\
# To do round-to-nearest, add (denom/2) to the numerator (halves will round up)\
let NUMEVENTS = ($fileEvents+$NUMFILES-1)/$NUMFILES\
\
START = 1\
FINISH = $NUMEVENTS\
# divide up the file\
while [$i -le $NUMFILES]; do\
		rooteventselector \'97recreate -f START -l FINISH inFileName:t $i.root\
		let i = i + 1\
		let START = FINISH + 1\
		let FINISH = $NUMEVENTS * $i\
		\
# Divide the file into 20 subfiles we can now run ridge check on them \
\
cd ..\
# we are now in TwoParticleCorrelation\
for f in tempFiles/*.root; do\
   src/ridge_check_parallel.exe \cf6 \CocoaLigature1 "\CocoaLigature0 $f\CocoaLigature1 " \'93out_$f\'94 &\cf2 \CocoaLigature0 \
done\
wait\
\
# delete the temp files\
rm -r tempFiles\
\
# hadd the files together\
cd pdfDir \cf7 \CocoaLigature1 \
\pard\tx543\pardeftab543\pardirnatural\partightenfactor0
\cf4 hadd -f $outFileName.root out_\{\cf8 1..20\cf4 \}.root\
# delete the smaller files\
rm -r out_*.root}