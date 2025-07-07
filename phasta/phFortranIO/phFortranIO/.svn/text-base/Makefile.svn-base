SUBSYSNAME := phasta
MODULENAME := phFortranIO
NOSHARED := 1

dirs := include:src
ifeq ($(ARCHOS), x86_linux-icc)
        DEFS := $(DEFS) -DLINUX 
        CFLAGS := -w
        FFLAGS := -w95 -W0 -cm -80
        LDFLAGS := -lifcore -lifport -lpthread -ldl -lc -lm 
endif

include $(DEVROOT)/Util/buildUtil/make.common
