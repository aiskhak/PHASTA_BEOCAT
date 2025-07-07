ulimit -s unlimited

# Came from .bashrc 

LD_LIBRARY_PATH=/share/apps/MPICH2/lib/trace_rlog:/share/apps/MPICH2/lib/:$LD_LIBRARY_PATH;
export LD_LIBRARY_PATH

#Intel version 11.1:  
source /share/apps/Compiler/11.1/080/bin/ifortvars.sh intel64
source /share/apps/Compiler/11.1/080/bin/iccvars.sh intel64

# PHASTA env:
export MPIHOME=/share/apps/MPICH-3.0.2

# ********* end of .bashrc

export DEVROOT="`pwd`"
#export SIMUI_ATT_DEFS=/users/SCOREC/public/meshSim/SpecSim/Config
export SIM_LICENSE_FILE=/home/iabolotn/Install/simmetrix/License/NorthCarolinaSU.license
export NODEP=1
export NOSHARED=1
export MESHSIM=/home/iabolotn/Install/simmetrix/latest
#export MESHSIM=/Install/simmetrix/7.2-120623
#export MESHSIM=/Install/simmetrix/7.2-120829
export LESLIBDIR="$DEVROOT/LIBLES"
export memLSLIBDIR="$DEVROOT/memLS"
#export MPIHOME=/share/apps/MPICH2
export PARALLEL=mpich2
export MPICH_NO_LOCAL=1
export MPICOMM=
# ch_p4 was above
export SIMMODELER_HOME=/home/iabolotn/Install/simmetrix/SimModeler-latest

export ACUSOLVE_LIB=/home/iabolotn/develop/LIBLES/lib/x86_64_linux-icc/libles.a

export ARCHOS=x86_64_linux-icc
#export ARCHOS=x86_64_linux-gnu
export MODELER=discrete
#export MODELER=parasolid
#chmod +x ./Util/buildUtil/getarch
#chmod +x ./Util/buildUtil/static.pl
#gmake setup -C ./Util/buildUtil

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Install/Intel/mkl/lib/intel64:/Install/Intel/lib/intel64   
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.1/080/bin/intel64:/opt/intel/Compiler/11.1/080/mkl/lib/em64t    # 11.1 version options
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MESHSIM/lib/x64_rhel5_gcc41:$MESHSIM/lib/x64_rhel5_gcc41/psKrnl:/Install/develop/phasta/phastaIO/lib/$ARCHOS
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/iabolotn/develop/Meshing/BLUtils/lib/$ARCHOS:/usr/lib64/mvapich2-psm/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/impi/4.1.0.024/intel64/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Install/mpich2-1.2.1p1/src/mpl
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/share/apps/MPICH2/lib


#export PATH=/Install/Intel/bin:/Install/MPICH2/bin:$PATH    # REMOVED for ifort11.1 version !
export PATH=/opt/intel/Compiler/11.1/080/bin/intel64:/share/apps/MPICH2/bin:$PATH
export PATH=/share/apps/MPICH2/bin$PATH


#export PARASOLID=/usr/local/parasolid/latest
#source /usr/local/parasolid/latest/parasolid.env
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/parasolid/latest/shared_object
