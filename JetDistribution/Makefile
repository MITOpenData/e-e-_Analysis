#path var + dir commands

ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs  -Werror -Wno-deprecated-declarations -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

#Public install of latest fastjet code (3.3.0) on lxplus; if not on lxplus or using different version, change this path as appropriate
FASTJET=`$(FASTJETDIR)/bin/fastjet-config --cxxflags --libs --plugins`
MKDIR_BIN=mkdir -p $(PWD)/bin
MKDIR_PDF=mkdir -p $(PWD)/pdfDir
MKDIR_OUTPUT=mkdir -p $(PWD)/output
CHECK_STUDYMULTDIR=printf "\nSTUDYMULTDIR is set to: \'$(STUDYMULTDIR)\'\n   !!!IF STUDYMULTDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"
CHECK_FASTJETDIR=printf "\nFASTJETDIR is set to: \'$(FASTJETDIR)\'\n   !!!IF FASTJETDIR IS EMPTY THE BUILD WILL FAIL!!!\n   Run 'source setStudyMultEnv.sh' in top level StudyMult directory to fix.\n\n"

#programs to make

all: mkdirBin mkdirPdf mkdirOutput checkStudyMultDir checkFastJetDir bin/makeGeneralJetTree.exe bin/checkThrustMajorMinor.exe

mkdirBin:
	$(MKDIR_BIN)
mkdirPdf:
	$(MKDIR_PDF)
mkdirOutput:
	$(MKDIR_OUTPUT)
checkStudyMultDir:
	@$(CHECK_STUDYMULTDIR)
checkFastJetDir:
	@$(CHECK_FASTJETDIR)

bin/makeGeneralJetTree.exe: src/makeGeneralJetTree.C
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I $(STUDYMULTDIR) -o bin/makeGeneralJetTree.exe src/makeGeneralJetTree.C

bin/checkThrustMajorMinor.exe: src/checkThrustMajorMinor.C
	$(CXX) $(CXXFLAGS) $(ROOT) -I $(STUDYMULTDIR) -o bin/checkThrustMajorMinor.exe src/checkThrustMajorMinor.C

clean:
	rm -f *~
	rm -f \#*.*#
	rm -f $(PWD)/include/#*.*#
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/#*.*#
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bash/#*.*#
	rm -f $(PWD)/bash/*~
	rm -f $(PWD)/paths/#*.*#
	rm -f $(PWD)/paths/*~
	rm -f $(PWD)/bin/*.exe
	rmdir bin
.PHONY: all

