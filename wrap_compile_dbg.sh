# load modules
module purge
module load intel-compilers/2023.2.1
module load impi/2021.10.0-intel-compilers-2023.2.1
module load libtirpc/1.3.3-GCCcore-12.3.0

export VERSION=200_memLS
export PARALLEL=intelmpi
export ARCHOS=x86_64_linux-IB
export DEVROOT="`pwd`"

#––– force Phasta to use your bundled METIS –––
export PHASTA_METIS_INC="${DEVROOT}/phasta/phMetis/phMetis"
export PHASTA_METIS_LIB="${DEVROOT}/phasta/phMetis/lib/${ARCHOS}/libmetis-intelmpi-O.a"
export CFLAGS="-I${PHASTA_METIS_INC} ${CFLAGS}"
export FFLAGS="-I${PHASTA_METIS_INC} ${FFLAGS}"
export CPPFLAGS="-I${PHASTA_METIS_INC} ${CPPFLAGS}"
export INCLUDES="-I${PHASTA_METIS_INC} ${INCLUDES}"
export LDFLAGS="-L${DEVROOT}/phasta/phMetis/lib/${ARCHOS} -lmetis-intelmpi-O ${LDFLAGS}"

export PHASTA_INCLUDE_DIR="${DEVROOT}/phasta/phSolver/200_memLS/phSolver/include"

export I_MPI_CC=icx
export I_MPI_CXX=icpx
export I_MPI_FC=ifx
export CC=icx
export CXX=icpx
export FC=ifx

# Set MPI include/lib paths
export MPICH_INC=${I_MPI_ROOT}/include
export MPICH_LIB=${I_MPI_ROOT}/lib
export FQLIBS=""
ls $MPICH_INC $MPICH_LIB

export SIMMODELER_ROOT="/homes/aiskhak/SimModeler/2025.1-250623dev"
export SIMLIB_DIR="${SIMMODELER_ROOT}/lib/x64_rhel9_gcc11"
export PSDK_LIB="${SIMLIB_DIR}/psKrnl"
export LD_LIBRARY_PATH="${SIMLIB_DIR}:${PSDK_LIB}:${LD_LIBRARY_PATH}"

# Compiler flags
#export CFLAGS="-std=gnu89 -qopenmp -Wno-implicit-function-declaration -Wno-int-conversion -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${CFLAGS}"
#export FFLAGS="-fpp -qopenmp -ffree-form -cpp -w -DLINUX -DNDEBUG -DPARALLEL -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${FFLAGS}"
#export FFLAGS="-g -O0 -DDEBUG -fpp -qopenmp -ffree-form -cpp -w -DLINUX -DPARALLEL -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${FFLAGS}"

#export FFLAGS="-g -O0 -zero -DDEBUG -fpp -qopenmp -ffree-form -cpp -w -heap-arrays -check bounds -check uninit -traceback -DLINUX -DPARALLEL -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${FFLAGS}"

#
export CFLAGS="-g -O0 -DDEBUG -std=gnu89 -qopenmp -Wno-implicit-function-declaration -Wno-int-conversion -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${CFLAGS}"
export CPPFLAGS="-P -traditional-cpp ${CPPFLAGS}"
export FFLAGS="-g -O0 -init=snan -check bounds,uninit,pointer -DDEBUG -fpp -qopenmp -ffree-form -cpp -w -heap-arrays -traceback -DLINUX -DPARALLEL -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${FFLAGS}"
export CXXFLAGS="-g -O0 -DDEBUG ${CXXFLAGS}"

# ASAN
#export CFLAGS="-g -O1 -fsanitize=address -fno-optimize-sibling-calls -std=gnu89 -qopenmp -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include"
#export CPPFLAGS="-P -traditional-cpp ${CPPFLAGS}"
#export FFLAGS="-g -O0 -zero -DDEBUG -fpp -qopenmp -ffree-form -cpp -w -heap-arrays -check bounds -check uninit -traceback -DLINUX -DPARALLEL -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${FFLAGS}"
#export CXXFLAGS="-g -O1 -fsanitize=address -fno-optimize-sibling-calls -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include"
#export LDFLAGS="-fsanitize=address ${LDFLAGS}"

export FFLAGS_FIXED="$FFLAGS"
export FFLAGS_FREE="$FFLAGS"

# Linker and include flags
export INCLUDES="-I${PHASTA_INCLUDE_DIR} -I${MPICH_INC} -I${SIMMODELER_ROOT}/include ${INCLUDES}"
#export LDFLAGS="-L${MPICH_LIB} -qopenmp -Wl,-rpath,${MPICH_LIB} -L${SIMLIB_DIR} -Wl,-rpath,${SIMLIB_DIR} -L${PSDK_LIB} -Wl,-rpath,${PSDK_LIB} -lpskernel -lSimLicense -ltirpc ${LDFLAGS}"
export LDFLAGS="-g -L${MPICH_LIB} -qopenmp -Wl,-rpath,${MPICH_LIB} -L${SIMLIB_DIR} -Wl,-rpath,${SIMLIB_DIR} -L${PSDK_LIB} -Wl,-rpath,${PSDK_LIB} -lpskernel -lSimLicense -ltirpc ${LDFLAGS}"

export SIM_LICENSE_FILE=/homes/aiskhak/SimModeler/SimModeler2025.0-250628/simmodsuite.lic
echo "Using license file: $SIM_LICENSE_FILE"

# uncomment line below to clean build
#isclean="clean"

setup="gmake -j32 VERS=dbg BUILD=debug NODEP=1 setup"
compile="gmake -j32 VERS=dbg BUILD=debug VERSION=200_memLS NODEP=1 NOSHARED=1 $isclean"

dest_path=$DEVROOT/phasta/phastaIO/phastaIO
cd $dest_path
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$setup
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$compile

dest_path=$DEVROOT/phasta/shapeFunction/shapeFunction
cd $dest_path
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$compile

dest_path=$DEVROOT/phasta/phMetis/phMetis
cd $dest_path
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$compile

dest_path=$DEVROOT/phasta/phNSpre/phNSpre
cd $dest_path
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$setup
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$compile

dest_path=$DEVROOT/phasta/phPost/phPost/Reduce
cd $dest_path
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$setup
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$compile

export LDFLAGS="$LDFLAGS \
    /homes/aiskhak/memAlloc_trace.o \
    -Wl,--wrap=malloc \
    -Wl,--wrap=memAlloc \
	-Wl,--wrap=free"
dest_path=$DEVROOT/phasta/phSolver/$VERSION/phSolver/
cd $dest_path
PHASTA_METIS_INC=${PHASTA_METIS_INC} \
PHASTA_METIS_LIB=${PHASTA_METIS_LIB} \
$compile

cd $DEVROOT