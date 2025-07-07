module purge
module load intel-compilers/2023.2.1
module load impi/2021.10.0-intel-compilers-2023.2.1
module load libtirpc/1.3.3-GCCcore-12.3.0

export VERSION=200_memLS
export PARALLEL=intelmpi
export ARCHOS=x86_64_linux-IB
export DEVROOT="`pwd`"

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
export CFLAGS="-std=gnu89 -qopenmp -Wno-implicit-function-declaration -Wno-int-conversion -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${CFLAGS}"
export CPPFLAGS="-P -traditional-cpp ${CPPFLAGS}"
export FFLAGS="-fpp -qopenmp -ffree-form -cpp -diag-disable=5082,5149 -w -DLINUX -DNDEBUG -DPARALLEL -I${PHASTA_INCLUDE_DIR} -I${SIMMODELER_ROOT}/include ${FFLAGS}"

# Linker and include flags
export INCLUDES="-I${PHASTA_INCLUDE_DIR} -I${MPICH_INC} -I${SIMMODELER_ROOT}/include ${INCLUDES}"
export LDFLAGS="-L${MPICH_LIB} -qopenmp -Wl,-rpath,${MPICH_LIB} -L${SIMLIB_DIR} -Wl,-rpath,${SIMLIB_DIR} -L${PSDK_LIB} -Wl,-rpath,${PSDK_LIB} -ltirpc ${LDFLAGS}"

# uncomment line below to clean build
#isclean="clean"

  setup="gmake -j32 VERS=opt NODEP=1 setup"
  compile="gmake -j2 VERS=opt VERSION=200_memLS  NODEP=1 NOSHARED=1 $isclean"

  dest_path=$DEVROOT/phasta/phastaIO/phastaIO
  cd $dest_path
  $setup
  $compile

  dest_path=$DEVROOT/phasta/shapeFunction/shapeFunction
  cd $dest_path
#  $setup
  $compile

  dest_path=$DEVROOT/phasta/phMetis/phMetis
  cd $dest_path
#  $setup
  $compile

  dest_path=$DEVROOT/phasta/phNSpre/phNSpre
  cd $dest_path
  $setup
  $compile

  dest_path=$DEVROOT/phasta/phPost/phPost/Reduce
  cd $dest_path
  $setup
  $compile

  dest_path=$DEVROOT/phasta/phSolver/$VERSION/phSolver/
  cd $dest_path
#  $setup
  $compile

  cd $DEVROOT

