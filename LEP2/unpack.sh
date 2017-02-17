cd LEP2_MC/GGUD
mv DIR*/*.aleph* .
rm -rf DIR0*
rm *mctrue_beforecu*
ls LEP2MC* >>list


while read F  ; do
        echo $F
        gunzip $F           
done <list
cd ../../