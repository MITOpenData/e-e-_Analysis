#MAKEFILE For thrust and reproduction of ALEPH results

#-Werror turned off
CXX = g++
CXXFLAGS = -Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11

ROOT=`root-config --cflags --glibs`

#Fastjet dependence valid and kept pretty up to date on lxplus, replace with own build if running elsewhere
FASTJET=`/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

EVENTPARTICLEDATA=$(PWD)/../DataProcessing/

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_PDFDIR=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdfdir bin/multDist.exe bin/plotMultDist.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

#bin/thrust_distribution.exe: src/thrust_distribution.c
#	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -o bin/thrust_distribution.exe src/thrust_distribution.c

bin/multDist.exe: src/multDist.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -o bin/multDist.exe src/multDist.C

bin/plotMultDist.exe: src/plotMultDist.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -o bin/plotMultDist.exe src/plotMultDist.C

clean:
	rm *~ || true
	rm *# || true
	rm include/*~ || true
	rm include/#*# || true
	rm src/*~ || true
	rm src/#*# || true
	rm src/*.d || true
	rm src/*.so || true
	rm src/*.pcm || true
	rm bin/*.exe || true
	rmdir bin || true
