SUBSYSNAME := phasta
MODULENAME := phE2sms
ISEXECUTABLE := 1
OUTNAME := meshModel-Sim
SIM := 1

ifndef ARCHOS
   ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

dirs := .

ifndef MESHSIM
  MESHSIM := /Install/simmetrix/latest  #7.2-111119
endif

ifeq ($(ARCHOS), sgi_6-64)
  LDFLAGS := $(LDFLAGS) -lftn -lfortran -lm 
  DEFS := $(DEFS) -Dsgi 
endif
ifeq ($(ARCHOS), sgi_6)
  LDFLAGS := $(LDFLAGS) -lftn -lfortran -lm 
  DEFS := $(DEFS) -Dsgi 
endif
ifeq ($(ARCHOS),sun4_5)
  LDFLAGS := $(LDFLAGS) -lfsu -lfai2 -lfsumai -lfminlai -lfminvai -lfmaxvai \
              -lfmaxlai -lF77 -lsunmath -lnsl -lsocket -lrt
  CXXFLAGS := $(CXXFLAGS) -w
  DEFS := $(DEFS) -Dsun4_5
  SIMARCH := sun4_5_sc5
endif
ifeq ($(ARCHOS),sun4_5-64)
  LDFLAGS := $(LDFLAGS) -lfsu -lfai2 -lfsumai -lfminlai -lfminvai -lfmaxvai \
              -lfmaxlai -lF77 -lsunmath -lnsl -lsocket -lrt
  CXXFLAGS := $(CXXFLAGS) -w
  DEFS := $(DEFS) -Dsun4_5
  SIMARCH := sun4_5_sc5-64
endif
ifeq ($(ARCHOS), x86_linux)
  DEFS := $(DEFS) -DLINUX
  SIMARCH := x86_rhel3_gcc32
endif
ifeq ($(ARCHOS), x86_linux-icc)
  DEFS := $(DEFS) -DLINUX
  LDFLAGS := $(LDFLAGS) -lifcore -ldl -lc -lm
  SIMARCH := x86_rhel3_gcc32
endif

ifeq ($(ARCHOS), x86_64_linux-icc)
  DEFS := $(DEFS) -DLINUX -Debug 
  LDFLAGS := $(LDFLAGS) -lifcore -ldl -lc -lm -lpthread -lnsl -Wuninitialized
#  SIMARCH := x64_rhel3_gcc32
  SIMARCH := x64_rhel5_gcc41
endif

INCLUDES := $(INCLUDES) -I/usr/include -I$(MESHSIM)/include -I$(DEVROOT)/phasta/phFortranIO/include -I$(DEVROOT)/phasta/phastaIO/include
DEFS := $(DEFS) -DSIM

DEPS := $(DEPS) phasta/phFortranIO phasta/phastaIO

MESHSIMLIBS := -lSimMeshing -lSimModel -lSimMeshing \
               -lSimMeshTools -lSimDiscrete -lSimModel \
	       -lSimMeshing -lSimMeshTools -lSimModel

ETCLIBS := $(ETCLIBS) -lpthread -lm

LDFLAGS := $(LDFLAGS) -L$(MESHSIM)/lib/$(SIMARCH) $(MESHSIMLIBS) $(ETCLIBS)
LIBS := $(LIBS) phFortranIO phastaIO

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
