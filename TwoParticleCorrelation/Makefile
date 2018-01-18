#MAKEFILE For two particle correlation

CXX = g++
CXXFLAGS = -Wall -O2 -Wextra  -Werror -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11

ROOT=`root-config --cflags --glibs`

#Fastjet dependence valid and kept pretty up to date on lxplus, replace with own build if running elsewhere
FASTJET=`/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

EVENTPARTICLEDATA=$(PWD)/../DataProcessing/

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_PDFDIR=mkdir -p $(PWD)/pdfDir

all: mkdirBin mkdirPdfdir eeplots plotDataQualityChecks ridge_check

mkdirBin:
	$(MKDIR_BIN)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

eeplots: src/eeplots.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -o bin/eeplots.exe src/eeplots.C

plotDataQualityChecks: src/plotDataQualityChecks.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -o bin/plotDataQualityChecks.exe src/plotDataQualityChecks.C

ridge_check: src/ridge_check.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -o bin/ridge_check.exe src/ridge_check.C


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
