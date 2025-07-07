SUBSYSNAME := phasta/phPost
MODULENAME := Reduce
ISEXECUTABLE := 1
OUTNAME := reduce
DEPS := $(DEPS) phasta/phastaIO
LIBS := phastaIO

ifndef ARCHOS
      ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

dirs := src

ifeq ($(ARCHOS), sgi_6-64)
  ADDSTDLIB := -lftn -lfortran -lm 
  DEFS := -Dsgi 
endif
ifeq ($(ARCHOS), sgi_6)
  ADDSTDLIB := -lftn -lfortran -lm 
  DEFS := -Dsgi 
endif
ifeq ($(ARCHOS),sun4_5)
  ADDSTDLIB :=-lfsu -lfai2 -lfsumai -lfminlai -lfminvai -lfmaxvai -lfmaxlai -lF77 -lsunmath
  CXXFLAGS := -w
  DEFS := -Dsun4_5 
  MTFLAG := -mt
endif

ADDLIB = $(ADDSTDLIB)

# Include standard makefile
include $(DEVROOT)/Util/buildUtil/make.common

rclean:
	rm .*.d
