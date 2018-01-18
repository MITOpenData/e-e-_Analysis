#path var + dir commands
ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs -Wno-deprecated-declarations -std=c++11
#CXXFLAGS=-Wall -O2 -Wextra -Wno-unused-local-typedefs 
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

#Public install of latest fastjet code (3.3.0) on lxplus; if not on lxplus or using different version, change this path as appropriate
FASTJET=`/afs/cern.ch/work/c/cmcginn/public/Fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
MKDIR_BIN=mkdir -p $(PWD)/bin

#replace with your own pythia8 path
PYTPATH=/afs/cern.ch/work/c/cmcginn/private/Generators/PYTHIA8/pythia8230

PROCPATH=../../DataProcessing/

#programs to make

all: mkdirBin mainCustomZee.exe
mkdirBin:
	$(MKDIR_BIN)
mainCustomZee.exe: src/mainCustomZee.cc $(PYTPATH)/lib/libpythia8.a
	$(CXX) src/mainCustomZee.cc -o bin/mainCustomZee.exe  -I $(PYTPATH)/include -pedantic -W $(CXXFLAGS) -Wshadow -fPIC -L$(PYTPATH)/lib -Wl,-rpath,$(PYTPATH)/lib -lpythia8 -ldl $(ROOT) $(FASTJET) -I $(PWD) -I $(PROCPATH)

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

