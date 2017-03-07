#!/bin/bash
#analysis(1,0,20,40)
#analysis(1,0,40,20)
#analysis(1,0,50,20)

rm Plots/*.pdf*

time root -b > an.log 2>&1 <<EOI
.L ridge.C+
analysis(1,0,20,40)
.L plot.C+
plot(1,20,0)
plot(1,20,1)
EOI

time root -b > an.log 2>&1 <<EOI
.L ridge.C+
analysis(1,0,40,20)
.L plot.C+
plot(1,40,0)
plot(1,40,1)
EOI

time root -b > an.log 2>&1 <<EOI
.L ridge.C+
analysis(1,0,50,20)
.L plot.C+
plot(1,50,0)
plot(1,50,1)
EOI

