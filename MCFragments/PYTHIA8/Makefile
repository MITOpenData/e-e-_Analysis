#path var + dir commands
ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
#CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs 
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

#Public install of latest fastjet code (3.3.0) on lxplus; if not on lxplus or using different version, change this path as appropriate
FASTJET=`$(FASTJETDIR)/bin/fastjet-config --cxxflags --libs --plugins`
MKDIR_BIN=mkdir -p $(PWD)/bin
#done for lxplus
SETPYT=export PYTHIA8DATA=/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/pythia8/230/share/Pythia8/xmldoc
#replace with your own pythia8 path, done for lxplus
#PYTPATH=/afs/cern.ch/work/c/cmcginn/public/Packages/PYTHIA8/pythia8230
CHECK_STUDYMULTDIR=printf "\nSTUDYMULTDIR is set to: \'$(STUDYMULTDIR)\'\n   !!!IF STUDYMULTDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"
CHECK_FASTJETDIR=printf "\nFASTJETDIR is set to: \'$(FASTJETDIR)\'\n   !!!IF FASTJETDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"
CHECK_STUDYPYTPATH=printf "\nSTUDYPYTPATH is set to: \'$(STUDYPYTPATH)\'\n   !!!IF STUDYPYTPATH IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"


#programs to make

#all: mkdirBin setpyt mainCustomZpole_ppAndee.exe
all: mkdirBin checkStudyMultDir checkFastJetDir checkStudyPytPath bin/mainCustomZpole_ppAndee.exe
mkdirBin:
	$(MKDIR_BIN)
checkStudyMultDir:
	@$(CHECK_STUDYMULTDIR)
checkFastJetDir:
	@$(CHECK_FASTJETDIR)
checkStudyPytPath:
	@$(CHECK_STUDYPYTPATH)
setpyt:
	$(SETPYT)
bin/mainCustomZpole_ppAndee.exe: src/mainCustomZpole_ppAndee.cc $(STUDYPYTPATH)/lib/libpythia8.a
	$(CXX) src/mainCustomZpole_ppAndee.cc -o bin/mainCustomZpole_ppAndee.exe  -I $(STUDYPYTPATH)/include -pedantic -W $(CXXFLAGS) -Wshadow -fPIC -L$(STUDYPYTPATH)/lib -Wl,-rpath,$(STUDYPYTPATH)/lib -lpythia8 -ldl $(ROOT) $(FASTJET) -I $(PWD) -I $(STUDYMULTDIR)

clean:
	rm -f *~
	rm -f \#*.*#
	rm -f $(PWD)/include/#*.*#
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/#*.*#
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/bash/#*.*#
	rm -f $(PWD)/bash/*~
	rm -f $(PWD)/paths/#*.*#
	rm -f $(PWD)/paths/*~
	rm -f $(PWD)/bin/*.exe
	rmdir bin
.PHONY: all

