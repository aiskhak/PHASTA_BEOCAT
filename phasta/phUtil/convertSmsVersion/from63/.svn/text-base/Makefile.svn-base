CC := g++
DEFS := -DPARASOLID
MESHSIM := /net/common/meshSim/6.2-080614
SIMARCH := x64_rhel3_gcc32
PARASOLID := $(MESHSIM)/lib/$(SIMARCH)/psKrnl
COPT:=O3

comp:
	$(CC) -c -$(COPT) $(DEFS) -I$(MESHSIM)/include \
           convertSMSVersion.cc -o convertSMSVersion.o
	$(CC) -$(COPT) -o convertSMSVersion convertSMSVersion.o \
           -L$(MESHSIM)/lib/$(SIMARCH) \
           -lSimAdvMeshing -lSimMeshing -lSimMeshTools \
           -lSimParasolid171 -lSimAcis160 -lSimModel \
           -L$(PARASOLID) \
           -lpskernel \
           -lm -lpthread
