SUBSYSNAME := phasta/phPost/Z-Correlation
MODULENAME := Z-Correlation
ISEXECUTABLE := 1
#OUTNAME :=correlationX.exe

dirs := src
DEPS := $(DEPS) phasta/phastaIO
LIBS := phastaIO

	ifeq ($(ARCHOS), )
	    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
	endif

	ifeq ($(ARCHOS),sun4_5)
		DEFS := $(DEFS) -Dsun4_5 -DACUSIM_SUN
	       LDFLAGS := -lfsu -lfai2 -lfsumai -lfminlai -lfminvai \
                  -lfmaxvai -lfmaxlai -lF77 -lsunmath -lthread \
                  -lsocket -lnsl -lrt -lnsl -laio  -lsunperf -lmvec
	endif


  ifeq ($(ARCHOS),x86_64_linux-icc)
          DEFS := $(DEFS) -DLINUX -DACUSIM_LINUX
          CFLAGS := $(CFLAGS) -w -g
          FFLAGS := $(FFLAGS) -w95 -W0 -cm -80 -fpp -g
          LDFLAGS := $(LDFLAGS) -lifcore -lpthread   -lstdc++
#         LDFLAGS := $(LDFLAGS) -L/usr/local/intel/mkl72/lib/32 -lmkl_ia32 -lmkl_p4 -lmkl_lapack
#         LDFLAGS := $(LDFLAGS) -lintrins -lCEPCF90 -lPEPCF90 -lF90
  endif


include $(DEVROOT)/Util/buildUtil/make.common
REALAR:=ifort -g -o

