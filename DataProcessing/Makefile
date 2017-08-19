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

#programs to make

all: mkdirBin scan.exe
mkdirBin:
	$(MKDIR_BIN)
scan.exe: src/scan.cc
	$(CXX) $(CXXFLAGS) $(ROOT) $(FASTJET) -I $(LOCAL) -o bin/scan.exe src/scan.cc
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

