# adjust this to your module’s path
MODDIR     := $(DEVROOT)/phasta/phastaIO/phastaIO
FCMANGLE_H := $(MODDIR)/FCMangle.h
PRECOMPILE := generate_fcmangle
.PHONY: generate_fcmangle
generate_fcmangle: $(FCMANGLE_H)
$(FCMANGLE_H):
	@mkdir -p $(MODDIR)
	@echo "/* Auto‐generated mangling header */"                      >  $@
	@echo "#ifndef FCMANGLE_H"                                        >> $@
	@echo "#define FCMANGLE_H"                                        >> $@
	@echo                                                            >> $@
	@echo "/* lowercase names with trailing underscore */"             >> $@
	@echo "#define FortranCInterface_GLOBAL(name,NAME) name##_"      >> $@
	@echo "#define FortranCInterface_GLOBAL_(name,NAME) name##_"     >> $@
	@echo "#define FortranCInterface_MODULE(name,NAME,mod,MOD) name##_"   >> $@
	@echo "#define FortranCInterface_MODULE_(name,NAME,mod,MOD) name##_"  >> $@
	@echo                                                            >> $@
	@echo "/* per‐routine shorthands */"                                >> $@
	@echo "#define commu_int   FortranCInterface_GLOBAL_(commu_int,COMMU_INT)" >> $@
	@echo "#define gen_ncorp   FortranCInterface_GLOBAL_(gen_ncorp,GEN_NCORP)" >> $@
	@echo                                                            >> $@
	@echo "#endif /* FCMANGLE_H */"                                   >> $@

# previous stuff
SUBSYSNAME := phasta/phUtil
MODULENAME := phastaIO
NODEP := 1
NOSHARED = 1

dirs := .

ifeq ($(ARCHOS), )
    ARCHOS := $(shell $(DEVROOT)/Util/buildUtil/getarch)
endif

include $(DEVROOT)/Util/buildUtil/make.common

