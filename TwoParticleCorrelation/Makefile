#MAKEFILE For two particle correlation

CXX = g++
CXXFLAGS = -Wall -O2 -Wextra   -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
#-Werror

ROOT=`root-config --cflags --glibs`

#Fastjet dependence valid and kept pretty up to date on lxplus, replace with own build if running elsewhere
FASTJET=`/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

EVENTPARTICLEDATA=$(PWD)/../

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_PDFDIR=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdfdir plotDQC ridge_check histoHadd_TH1F

mkdirBin:
	$(MKDIR_BIN)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

plotDQC: src/plotDQC.cc
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -o bin/plotDQC.exe src/plotDQC.cc

ridge_check: src/ridge_check.c
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -o bin/ridge_check.exe src/ridge_check.c

histoHadd_TH1F: src/histoHadd_TH1F.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -o bin/histoHadd_TH1F.exe src/histoHadd_TH1F.C

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
