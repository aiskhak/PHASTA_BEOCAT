# Makefile for svLS linear solver (static library)
# Generates:
#   debug:  $(DEVROOT)/phasta/svLS/lib/$(ARCHOS)/libsvLS-$(PARALLEL).a
#   release:$(DEVROOT)/phasta/svLS/lib/$(ARCHOS)/libsvLS-$(PARALLEL)-O.a

# Variables (from environment)
DEVROOT     ?= $(DEVROOT)
ARCHOS      ?= $(ARCHOS)
PARALLEL    ?= $(PARALLEL)
BUILD       ?= $(BUILD)        # "debug" or "release"

# Compiler and tools
FC          ?= $(I_MPI_FC)
AR          ?= ar
ARFLAGS     ?= rcs

# Include directories (set by the top-level script)
INCLUDES    = $(MPI_INCDIR) $(METIS_INCDIR)
FFLAGS      ?= $(FFLAGS)

# Determine library suffix based on build type
ifeq ($(BUILD),debug)
  QUAL =
else
  QUAL = -O
endif

# Directory for output library
LIBDIR      = $(DEVROOT)/phasta/svLS/lib/$(ARCHOS)
LIB         = $(LIBDIR)/libsvLS-$(PARALLEL)$(QUAL).a

# Source and object files
FSRCS       = ADDBCMUL.f \
               BC.f \
               CGRAD.f \
               COMMU.f \
               CPUT.f \
               DOT.f \
               GE.f \
               GMRES.f \
               INCOMMU.f \
               LHS.f \
               LS.f \
               NORM.f \
               NSSOLVER.f \
               PRECOND.f \
               SOLVE.f \
               SPARMUL.f
OBJS        = $(FSRCS:.f=.o)

.PHONY: all clean install

all: $(LIB)

# Create lib directory if needed
$(LIBDIR):
	mkdir -p $(LIBDIR)

# Archive objects into the versioned library
$(LIB): $(OBJS) | $(LIBDIR)
	$(AR) $(ARFLAGS) $@ $^

# Compile Fortran sources
%.o: %.f
	$(FC) $(FFLAGS) $(INCLUDES) -c $< -o $@

# Optional install
install: all
	install -d $(PREFIX)/lib
	install -m 644 $(LIB) $(PREFIX)/lib/
	install -d $(PREFIX)/include
	install -m 644 \
	  svLS_API_c.h \
	  svLS_API.h \
	  svLS_c.h \
	  svLS.h \
	  svLS_STD.h \
	  svLS_STRUCT_c.h \
	  svLS_STRUCT.h \
	  $(PREFIX)/include/

clean:
	rm -f $(OBJS) $(LIB)
