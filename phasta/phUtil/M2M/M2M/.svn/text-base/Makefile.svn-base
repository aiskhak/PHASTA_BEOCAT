SUBSYSNAME := phasta/phUtil
MODULENAME := M2M
OUTNAME := M2M
ISEXECUTABLE := 1

ifndef ARCHOS
  ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

dirs := src:utils:solution

ifndef PARASOLID
  PARASOLID := /Install/simmetrix/latest/lib/x64_rhel5_gcc41/psKrnl
#  PARASOLID := /usr/local/parasolid/16.0
endif

ifndef MESHSIM
  MESHSIM := /users/SCOREC/public/meshSim/latest
endif

ifeq ($(ARCHOS),sun4_5)
  SIMARCH := sun4_5_sc5
  ETCFLAGS := $(ETCFLAGS) -lm -lnsl -lsocket -lthread -lsunperf -lF77 -lsunmath -lfsu
  ETCFLAGS := $(ETCFLAGS) -lfsumai -lfui -lfai -lfai2 -lfsumai -lfprodai -lfminlai -lfmaxlai -lfminvai -lfmaxvai -lmopt -lcopt -lz
endif

ifeq ($(ARCHOS),sun4_5-64)
  SIMARCH := sun4_5_sc5-64
  PARASOLID := /usr/local/parasolid/16.0-64
  ETCFLAGS := $(ETCFLAGS) -lm -lnsl -lsocket -lthread -lsunperf -lF77 -lsunmath -lfsu
endif

ifeq ($(ARCHOS),x86_linux)
  SIMARCH := x86_rhel3_gcc32
  ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
endif

ifeq ($(ARCHOS),x86_linux-icc)
  SIMARCH := x86_rhel3_gcc32
  ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
endif

ifeq ($(ARCHOS),x86_64_linux)
  SIMARCH := x64_rhel3_gcc32
  ETCFLAGS := $(ETCFLAGS) -lpthread
endif

ifeq ($(ARCHOS),x86_64_linux-icc)
  SIMARCH := x64_rhel5_gcc41
  ETCFLAGS := $(ETCFLAGS) -lpthread
endif

ifeq ($(ARCHOS),x86_64_linux-IB)
  SIMARCH := x64_rhel5_gcc41
  ETCFLAGS := $(ETCFLAGS) -lpthread
endif

ifeq ($(HOSTNAME),master.medusa.scorec.rpi.edu)
  MESHSIM :=/import/users/SCOREC/public/meshSim/latest/
endif

ifndef SIMARCH
  ifeq ($(ARCHOS),x86_linux)
    SIMARCH := x86_rhel3_gcc32
  endif
  ifeq ($(ARCHOS),x86_linux-icc)
    SIMARCH := x86_rhel3_gcc32
  endif
  ifeq ($(ARCHOS),x86_linux-icc)
    SIMARCH := amd64_rhel3_gcc32
  endif
  ifeq ($(ARCHOS),x86_linux-icc)
    SIMARCH := amd64_rhel3_gcc32
  endif
  ifeq ($(ARCHOS),sun4_5)
    SIMARCH := sun4_5_sc5
  endif
  ifeq ($(ARCHOS),sun4_5-64)
    SIMARCH := sun4_5_sc5-64
  endif
endif

DEPS := $(DEPS) \
	phasta/phastaIO \
        phasta/phUtil/LU

LIBS := $(LIBS) phastaIO LU

DEFS := $(DEFS) -DSIM

#  PARASOLID_LIBS := -lfg -lparasolid
PARASOLID_LIBS := -lpskernel
  ifeq ($(ARCHOS),x86_linux-icc)
SIMPARASOLID_LIBS := -lSimParasolid250
  endif
  ifeq ($(ARCHOS),x86_linux-IB)
SIMPARASOLID_LIBS := -lSimParasolid260
  endif

INCLUDES := $(INCLUDES) \
            -I$(PARASOLID) \
            -I$(MESHSIM)/include 

MESHSIMFLAGS := -L$(MESHSIM)/lib/$(SIMARCH) 
MESHSIMFLAGS := $(MESHSIMFLAGS) \
	        -lSimMeshing -lSimModel -lSimMeshTools $(SIMPARASOLID_LIBS) -lSimDiscrete \
                -lSimModel -lSimMeshing -lSimMeshTools -lSimModel

ADDLIB = $(call libtmpl,$(PARASOLID)) $(PARASOLID_LIBS) -lm

LDFLAGS := $(LDFLAGS) $(MESHSIMFLAGS) $(ETCFLAGS)

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
