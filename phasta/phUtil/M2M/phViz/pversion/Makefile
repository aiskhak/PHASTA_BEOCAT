SUBSYSNAME := phasta
MODULENAME := phViz
OUTNAME    := phProject
ISEXECUTABLE := 1
MODELER      := Parasolid

ifndef ARCHOS
  ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

dirs := .


ifndef MESHSIM
  MESHSIM := /users/SCOREC/public/meshSim/latest
endif


ifeq ($(ARCHOS),x86_cygwin)
    PARASOL_DIR:=C:\program files\Parasolid\Kernel\v16.0\base
    ETCFLAGS := $(ETCFLAGS) -LIBPATH:"c:\program files\intel\compiler\fortran\9.0\ia32\lib" \
                            Advapi32.lib Ws2_32.lib /subsystem:console /release
endif



ifeq ($(ARCHOS),x86_linux)
     SIMARCH := x86_rh8_gcc32
     ETCFLAGS := $(ETCFLAGS)  -lnsl    -lpthread
endif
ifeq ($(ARCHOS),x86_linux-icc)
     SIMARCH := x86_rhel3_gcc32
     ETCFLAGS := $(ETCFLAGS)  -lnsl   -lpthread 
endif
ifeq ($(ARCHOS),sun4_5)
     SIMARCH := sun4_5_sc5
     ETCFLAGS := $(ETCFLAGS) -lF77 -lsunmath -lnsl -lsocket  -lrt
endif

ifeq ($(ARCHOS),x86_cygwin)
     SIMARCH := x86_win_net
     ETCFLAGS := $(ETCFLAGS)  ws2_32.lib
endif


ifeq ($(ARCHOS),sgi_6)
       SIMARCH := sgi_6
       DEFS := $(DEFS) -Dsgi_6 -Dsgi 
#       fortlib := /usr/lib64/libcomplib.sgimath.so
       LDFLAGS := -v -lc -lfortran $(fortlib) -lftn -lffio\
                  -lfpe -lffio -lfpe -lm -lnsl
       ETCFLAGS := $(ETCFLAGS)    -lfg
endif 

DEFS := $(DEFS) -DSIM

ifeq ($(ARCHOS),x86_cygwin)
  IOLIB_DIR     := $(DEVROOT)/phasta/phastaIO/lib/$(ARCHOS)
  WINIOLIB_DIR  := $(shell cygpath -w $(IOLIB_DIR) )
  LIBS          := -LIBPATH:"$(WINIOLIB_DIR)" libphastaIO.lib
  UTIL_LIBDIR   :=$(DEVROOT)/phasta/phUtil/lib/$(ARCHOS)
  WINUTIL_LIBDIR:=$(shell cygpath -w $(UTIL_LIBDIR) )
  LIBS          := $(LIBS) -LIBPATH:"$(WINUTIL_LIBDIR)" libLU.lib libSolution.lib
  SHAPE_DIR     :=$(DEVROOT)/phasta/phShape/lib/$(ARCHOS)
  WINSHAPE_DIR  :=$(shell cygpath -w $(SHAPE_DIR) )
  LIBS          := $(LIBS) -LIBPATH:"$(WINSHAPE_DIR)" libphShape.lib
else
  DEPS := $(DEPS)  phasta/phastaIO \
          phasta/phUtil/LU \
          phasta/phShape\
          phasta/phUtil/Solution/pversion
  LIBS := $(LIBS) phastaIO LU phShape Solution
endif


ifeq ($(ARCHOS),x86_cygwin)
  INCLUDES := $(INCLUDES) -I"$(MESHSIM)/include" 
  WIN_IOINCLUDE := $(shell cygpath -w $(DEVROOT)/phasta/phastaIO/include )
  INCLUDES := $(INCLUDES) -I"$(WIN_IOINCLUDE)"
else
  INCLUDES := $(INCLUDES) -I$(MESHSIM)/include 
endif


ifeq ($(MODELER),Parasolid)
  MESHSIMMODELVERS = -lSimParasolid180
  INCLUDES :=$(INCLUDES)
  PARASOLID := $(MESHSIM)/lib/$(SIMARCH)/psKrnl
  PSKRNL:= -lpskernel
  ifeq ($(ARCHOS), sun4_5)
    MODELERFLAGS :=  -lnsl -lthread
  endif
  ifeq ($(ARCHOS),x86_linux-icc)
     MODELERFLAGS :=  -lnsl 
  endif
  ifeq ($(ARCHOS),x86_64_linux-icc)
     MODELERFLAGS :=  -lnsl -lpthread
  endif
  ifeq ($(ARCHOS),x86_linux)
     MODELERFLAGS := -lnsl 
  endif
  ifeq ($(ARCHOS),x86_cygwin)
     MESHSIMMODELVERS :=-LIBPATH:"$(MESHSIM)/lib/x86_win_net" SimParasolid160.lib	
     MODELERFLAGS := -LIBPATH:"$(PARASOL_DIR)"  PSKERNEL_ARCHIVE.lib
  endif
endif

ifeq ($(MODELER),Discrete)
 ifeq ($(ARCHOS),sgi_6)
   MESHSIMMODELVERS := -lSimDiscrete
 endif
 ifeq ($(ARCHOS),sgi_6-64)
   MESHSIMMODELVERS := -lSimDiscrete
 endif
 ifeq ($(ARCHOS),sun4_5)
   MESHSIMMODELVERS := -lSimDiscrete -lnsl
   ETCFLAGS := $(ETCFLAGS) -lF77 -lsunmath  -lsocket  -lthread
 endif
 ifeq ($(ARCHOS),x86_linux-icc)
    MESHSIMMODELVERS := -lSimDiscrete
    ETCFLAGS := $(ETCFLAGS)  -lnsl     -lpthread
  endif
 ifeq ($(ARCHOS),x86_cygwin)
   MESHSIMMODELVERS := -LIBPATH:"$(MESHSIM)/lib/x86_win_net" SimDiscrete.lib
 endif
endif



ifeq ($(ARCHOS),x86_cygwin)
  MESHSIMFLAGS :=   $(MESHSIMMODELVERS) \
                    -LIBPATH:"$(MESHSIM)/lib/$(SIMARCH)" SimAdvMeshing.lib \
                     SimModel.lib SimMeshing.lib SimMeshTools.lib \
                     SimModel.lib SimMeshing.lib SimMeshTools.lib 
else
   MESHSIMFLAGS := -L$(MESHSIM)/lib/$(SIMARCH)
   MESHSIMFLAGS := $(MESHSIMFLAGS)  -L$(PARASOLID) \
               -lSimAdvMeshing  -lSimMeshing -lSimMeshTools \
		$(MESHSIMMODELVERS) -lSimModel \
                -lSimMeshTools  $(PSKRNL)   -lm
endif


ifeq ($(ARCHOS),x86_cygwin)
LDFLAGS := $(LDFLAGS) $(ETCFLAGS) $(MESHSIMFLAGS) $(MODELERFLAGS) $(LIBS)
else
LDFLAGS := $(LDFLAGS) $(ETCFLAGS) $(MESHSIMFLAGS) $(MODELERFLAGS)
LIBS := $(LIBS) phShape Solution LU phastaIO
endif


# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
