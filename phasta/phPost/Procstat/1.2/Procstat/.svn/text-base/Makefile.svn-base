SUBSYSNAME := phasta/phPost/Procstats
MODULENAME := Procstats
ISEXECUTABLE := 1
#OUTNAME :=procstats.exe

dirs := src
#DEPS := $(DEPS) phasta/phasta/phSolver/1.9.5_stats/phSolver/common
DEPS := $(DEPS) phasta/phastaIO
LIBS := phastaIO

#	ifndef ARCHOS
#	  ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
#	endif

	ifeq ($(ARCHOS), )
	    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
	endif

	ifeq ($(ARCHOS),sun4_5)
		DEFS := $(DEFS) -Dsun4_5 -DACUSIM_SUN
	       LDFLAGS := -lfsu -lfai2 -lfsumai -lfminlai -lfminvai \
                  -lfmaxvai -lfmaxlai -lF77 -lsunmath -lthread \
                  -lsocket -lnsl -lrt -lnsl -laio  -lsunperf -lmvec
	endif

include $(DEVROOT)/Util/buildUtil/make.common


