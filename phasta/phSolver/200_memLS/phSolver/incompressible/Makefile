SUBSYSNAME := phasta/phSolver/incompressible
MODULENAME := incompressible
BUILDV := 1.0
ISEXECUTABLE := 1
OUTNAME := phastaIC.exe
ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

  PRECOMPILE := copymods
  ARPACKPATH := $(DEVROOT)/ARPACK

  LIBARPACK := $(ARPACKPATH)/libarpack_linux.a
  LIBPARPACK := $(ARPACKPATH)/parpack_MPI-linux.a
  
  ifeq ($(LESLIBDIR), )
       LESLIBDIR := /Install/develop/LIBLES
  endif
  ifeq ($(ARCHOS),x86_64_linux-icc)
	  DEFS := $(DEFS) -DLINUX -DACUSIM_LINUX
	  CFLAGS := -w
	  FFLAGS := -W0 -132 -fpp
#          FFLAGS := -ffixed-line-length-132
	  ifeq ($(AMG),1)
	      LDFLAGS := $(LDFLAGS) $(LIBPARPACK) $(LIBARPACK)
	  endif
#          LDFLAGS := $(LDFLAGS) -lifcore  -lifport -lpthread -lmpl
          LDFLAGS := $(LDFLAGS) -lifcore  -lifport -lpthread -lmpi
          SIMARCH := amd64_rhel3_gcc32
  endif

  ifeq ($(ARCHOS),x86_64_linux-IB)
          DEFS := $(DEFS) -DLINUX -DACUSIM_LINUX
          CFLAGS := -w
          FFLAGS := -W0 -132 -fpp
          ifeq ($(AMG),1)
                 LDFLAGS := $(LDFLAGS) $(LIBPARPACK) $(LIBARPACK)
          endif
#          LDFLAGS := $(LDFLAGS) -lifcore  -lifport -lpthread -lmpl
          LDFLAGS := $(LDFLAGS) -lifcore  -lifport -lpthread -lmpi
          SIMARCH := amd64_rhel3_gcc32
  endif
                                                                 

ifeq ($(AMG),1)
    DEFS := $(DEFS) -DAMG
endif


dirs := .
mods  := stats.f
DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver phasta/phastaIO phasta/phMetis\
        phasta/shapeFunction LIBLES/1.5 memLS
LIBS := common phastaIO metis shapeFunction les memLS

ifeq ($(AMG),1)
    DEPS := $(DEPS) phasta/phSolver/$(VERSION)/phSolver
	LIBS := $(LIBS) AMG
endif

INCLUDES := $(INCLUDES) -I../common
ADDLIB = $(ADDSTDLIB)
# need to include common for intel fortran .mod files
include $(DEVROOT)/Util/buildUtil/make.common

.PHONY: copymods
.PHONY: $(OUTFILE)

copymods:
	cp ../common/*.mod .
ifeq ($(AMG),1)	
	cp ../AMG/*.mod .
endif
