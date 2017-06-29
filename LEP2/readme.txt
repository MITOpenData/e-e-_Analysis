1. aleph_analysis.c: plots the eta and the phi of different particles

2. aleph_mult.c: plots the multiplicity of different particles

3. clean.sh: deletes a bunch of files (see code)

4. loop.sh: Goes over datafiles and removes strings like "px =", and leaves pure numbers. That file is then used to create a ntuple

5. loopMC.sh: similar task as loop.sh, but on a different file

6. loopMCreo.sh: slight modification of loop.sh

7. MCvalidation.c: 

8. scan.c: takes file processed from loop.sh or any of its' modifications and creates a ntuple

9. scanMC.c: basically the same file as scan.c

10. unpack.sh: 