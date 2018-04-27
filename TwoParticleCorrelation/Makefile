#MAKEFILE For two particle correlation

CXX = g++
CXXFLAGS = -Wall -O2 -Wextra   -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
#-Werror

ROOT=`root-config --cflags --glibs`

#Fastjet dependence valid and kept pretty up to date on lxplus, replace with own build if running elsewhere
FASTJET=`/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

EVENTPARTICLEDATA=$(PWD)/../

DATAPROCESSING=$(PWD)/../DataProcessing

MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_PDFDIR=mkdir -p $(PWD)/pdfDir

#all: mkdirBin mkdirPdfdir plotDQC ridge_check histoHadd_TH1F findGlobalMaxMin

all: mkdirBin mkdirPdfdir bin/ridge_check.exe bin/histoHadd_TH1F.exe bin/findGlobalMaxMin.exe bin/plotDQC.exe

mkdirBin:
	$(MKDIR_BIN)

mkdirPdfdir:
	$(MKDIR_PDFDIR)

bin/plotDQC.exe: src/plotDQC.cc
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -I $(DATAPROCESSING) -o bin/plotDQC.exe src/plotDQC.cc

bin/ridge_check.exe: src/ridge_check.c
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -I $(DATAPROCESSING) -o bin/ridge_check.exe src/ridge_check.c

bin/histoHadd_TH1F.exe: src/histoHadd_TH1F.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -I $(DATAPROCESSING) -o bin/histoHadd_TH1F.exe src/histoHadd_TH1F.C

bin/findGlobalMaxMin.exe: src/findGlobalMaxMin.C
	$(CXX) $(CXXFLAGS) $(ROOT)  -I $(PWD) -I $(EVENTPARTICLEDATA) -I $(DATAPROCESSING) -o bin/findGlobalMaxMin.exe src/findGlobalMaxMin.C

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
