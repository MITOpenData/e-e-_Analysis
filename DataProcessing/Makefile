#path var + dir commands

LOCAL=$(PWD)
ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

#Public install of latest fastjet code (3.3.0) on lxplus; if not on lxplus or using different version, change this path as appropriate
FASTJET=`/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
MKDIR_BIN=mkdir -p $(LOCAL)/bin
MKDIR_PDF=mkdir -p $(LOCAL)/pdfDir
CHECK_STUDYMULTDIR=printf "\nSTUDYMULTDIR is set to: \'$(STUDYMULTDIR)\'\n   !!!IF STUDYMULTDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

#programs to make

all: mkdirBin mkdirPdf checkGlobal bin/scan.exe bin/mixMain.exe bin/countNEvt.exe bin/doComparison.exe
mkdirBin:
	$(MKDIR_BIN)
mkdirPdf:
	$(MKDIR_PDF)
checkGlobal:
	@$(CHECK_STUDYMULTDIR)

bin/scan.exe: src/scan.cc
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I $(LOCAL) -I $(STUDYMULTDIR) -o bin/scan.exe src/scan.cc
bin/mixMain.exe: src/mixMain.cc
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I $(LOCAL) -I $(STUDYMULTDIR) -o bin/mixMain.exe src/mixMain.cc
bin/doComparison.exe: src/doComparison.cc
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(LOCAL) -I $(STUDYMULTDIR) -o bin/doComparison.exe src/doComparison.cc
bin/countNEvt.exe: src/countNEvt.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I $(LOCAL) -I $(STUDYMULTDIR) -o bin/countNEvt.exe src/countNEvt.C
clean:
	rm -f *~
	rm -f \#*.*#
	rm -f $(LOCAL)/include/#*.*#
	rm -f $(LOCAL)/include/*~
	rm -f $(LOCAL)/src/#*.*#
	rm -f $(LOCAL)/src/*~
	rm -f $(LOCAL)/bash/#*.*#
	rm -f $(LOCAL)/bash/*~
	rm -f $(LOCAL)/paths/#*.*#
	rm -f $(LOCAL)/paths/*~
	rm -f $(LOCAL)/bin/*.exe
	rmdir bin
.PHONY: all

