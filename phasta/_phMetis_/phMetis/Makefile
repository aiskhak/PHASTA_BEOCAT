SUBSYSNAME = phasta
MODULENAME = phMetis
OUTNAME = metis
NODEP = 1
NOSHARED = 1
dirs := Lib

ifeq ($(ARCHOS),x86_cygwin)
    DEFS := $(DEFS) -D__VC__ -DINTEL
endif

include $(DEVROOT)/Util/buildUtil/make.common
