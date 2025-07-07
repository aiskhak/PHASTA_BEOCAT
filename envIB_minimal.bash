ulimit -s unlimited

#source /share/apps/Compiler/11.1/080/bin/ifortvars.sh intel64
#source /share/apps/Compiler/11.1/080/bin/iccvars.sh intel64
#source /opt/intel/impi/4.1.0.024/intel64/bin/mpivars.sh
  export DEVROOT="`pwd`"
#export SIM_LICENSE_FILE=/home/iabolotn/Install/simmetrix/License/NorthCarolinaSU.license
  #export NODEP=1
  #export NOSHARED=1
#export MESHSIM=/home/iabolotn/Install/simmetrix/latest
  #export LESLIBDIR="$DEVROOT/LIBLES"
  #export memLSLIBDIR="$DEVROOT/memLS"
  export PARALLEL=mpich2
  #export MPICH_NO_LOCAL=1
  #export MPICOMM=ch_p4
#export SIMMODELER_HOME=/home/iabolotn/Install/simmetrix/SimModeler-latest
  export ARCHOS=x86_64_linux-IB
  #export MODELER=parasolid
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.1/080/bin/intel64:/opt/intel/Compiler/11.1/080/mkl/lib/em64t    # 11.1 version options
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MESHSIM/lib/x64_rhel5_gcc41:$MESHSIM/lib/x64_rhel5_gcc41/psKrnl:/Install/develop/phasta/phastaIO/lib/$ARCHOS
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/iabolotn/develop/Meshing/BLUtils/lib/$ARCHOS  #:/usr/lib64/mvapich2-psm/lib
#export PATH=/opt/intel/Compiler/11.1/080/bin/intel64:$PATH
