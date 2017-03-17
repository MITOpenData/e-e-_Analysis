#!/bin/bash
#analysis(1,0,20,40)
#analysis(1,0,40,20)
#analysis(1,0,50,20)

rm Plots/*.pdf*
rm ROOTfiles/*.root

if (true) then 
time root -b > an.log 2>&1 <<EOI
.L ridge.C+
analysis(1,0,20,40)
.L plot.C+
plot(1,20,0)
EOI
fi

if (true) then 
time root -b > an.log 2>&1 <<EOI
.L ridge.C+
analysis(1,0,40,20)
.L plot.C+
plot(1,40,0)
EOI
fi

if (true) then 
time root -b > an.log 2>&1 <<EOI
.L ridge.C+
analysis(1,0,50,20)
.L plot.C+
plot(1,50,0)
EOI
fi
