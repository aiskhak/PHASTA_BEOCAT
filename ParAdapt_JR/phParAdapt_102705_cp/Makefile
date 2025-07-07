SUBSYSNAME := phasta
MODULENAME := phParAdapt
ISEXECUTABLE := 1
SIM := 1
OUTNAME := phParAdapt

ifndef PARASOLID
  PARASOLID := /usr/local/parasolid/latest
endif

ifndef MODELER
  MODELER := Parasolid
endif

OUTNAME:=$(OUTNAME)-$(MODELER)

ifeq ($(HOSTNAME),master.medusa.scorec.rpi.edu)
  DEVROOT := /import/users/$(USER)/develop
endif

ifndef ARCHOS
  ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

ifndef MESHSIM
  MESHSIM := /users/SCOREC/public/meshSim/latest
endif

dirs := .
SIMARCH := $(ARCHOS)

ifeq ($(ARCHOS),sun4_5)
  SIMARCH := sun4_5_sc5
  ETCFLAGS := $(ETCFLAGS) -lnsl -lm -lsocket -lsunperf -lF77 -lsunmath -lthread -lfsu
endif

ifeq ($(ARCHOS),sun4_5-64)
  SIMARCH := sun4_5_sc5-64
  ETCFLAGS := $(ETCFLAGS) -lnsl -lm -lsocket -lsunperf -lF77 -lsunmath -lthread -lfsu
endif

ifeq ($(ARCHOS),x86_linux-icc)
   SIMARCH := x86_rhel3_gcc32
   ETCFLAGS := $(ETCFLAGS) -lnsl -lpthread
endif

ifeq ($(HOSTNAME),master.medusa.scorec.rpi.edu)
    SIMARCH := x86_rhel3_gcc32
    MESHSIM := /import/users/SCOREC/public/meshSim/latest/
    ETCFLAGS := $(ETCFLAGS)  -lnsl -lpthread
endif

ifeq ($(PARALLEL),mpich2)
  ifeq ($(MPICOMM),ch_p4)
    SIMPAR := -mpichp4
  endif
  ifeq ($(MPICOMM),ch_shmem)
    SIMPAR := -mpichp4
  endif
endif
ifeq ($(PARALLEL),sgimpi)
SIMPAR := -mpisgi
endif

ifdef PARALLEL
ETCFLAGS :=  $(ETCFLAGS) -lparmetis$(SIMPAR) -lautopack$(SIMPAR) 
endif

DEPS := $(DEPS) \
	phasta/phShape \
	phasta/phUtil/LU \
	phasta/phastaIO

LIBS := $(LIBS) phShape LU phastaIO


ifeq ($(MODELER),Parasolid)
  ADDMODELLIB = $(call libtmpl,$(PARASOLID)) -lparasolid -lfg
  MESHSIMMODELVERS = -lSimParasolid160 
else
  MODELLIB = -lSim$(MODELER)
endif 


INCLUDES := $(INCLUDES) -I$(MESHSIM)/include
DEFS := $(DEFS) -DSIM -DMODELER_$(shell echo $(MODELER)| tr '[a-z]' '[A-Z]')

MESHSIMFLAGS := -L$(MESHSIM)/lib/$(SIMARCH) 
MESHSIMFLAGS := $(MESHSIMFLAGS)  \
	-lSimMeshing$(SIMPAR) -lSimPMesh$(SIMPAR) -lSimMeshTools$(SIMPAR)\
			   -lSimMeshing$(SIMPAR) -lSimMeshTools$(SIMPAR) \
			   -lSimMeshing$(SIMPAR) -lSimMeshTools$(SIMPAR) \
			   -lSimPMesh$(SIMPAR)  $(MESHSIMMODELVERS) \
				-lSimModel \


ADDLIB = $(ADDMODELLIB)

LDFLAGS := $(LDFLAGS) $(MESHSIMFLAGS) $(ETCFLAGS) -lmetis -lnsl -lm 

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common
