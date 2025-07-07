ulimit -s unlimited
export DEVROOT="`pwd`"
#export SIMUI_ATT_DEFS=/users/SCOREC/public/meshSim/SpecSim/Config
export SIM_LICENSE_FILE=/Install/simmetrix/License/NorthCarolinaSU.license
export NODEP=1
export NOSHARED=1
export MESHSIM=/Install/simmetrix/latest
#export MESHSIM=/Install/simmetrix/7.2-120623
#export MESHSIM=/Install/simmetrix/7.2-120829
export LESLIBDIR="$DEVROOT/LIBLES"
export memLSLIBDIR="$DEVROOT/memLS"
export MPIHOME=/Install/MPICH2
export PARALLEL=mpich2
export MPICH_NO_LOCAL=1
export MPICOMM=ch_p4
export SIMMODELER_HOME=/Install/simmetrix/SimModeler-latest

export ARCHOS=x86_64_linux-icc
export MODELER=discrete
#export MODELER=parasolid
#chmod +x ./Util/buildUtil/getarch
#chmod +x ./Util/buildUtil/static.pl
#gmake setup -C ./Util/buildUtil

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Install/Intel/mkl/lib/intel64:/Install/Intel/lib/intel64   
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/Compiler/11.1/080/lib/intel64:/opt/intel/Compiler/11.1/080/mkl/lib/em64t    # 11.1 version options
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MESHSIM/lib/x64_rhel5_gcc41:$MESHSIM/lib/x64_rhel5_gcc41/psKrnl:/Install/develop/phasta/phastaIO/lib/$ARCHOS
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Install/develop/Meshing/BLUtils/lib/$ARCHOS
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Install/MPICH2/lib



#export PATH=/Install/Intel/bin:/Install/MPICH2/bin:$PATH    # REMOVED for ifort11.1 version !
export PATH=/opt/intel/Compiler/11.1/080/bin/intel64:/Install/mpich2-1.2.1p1/bin$PATH
export PATH=/Install/MPICH2/bin$PATH


#export PARASOLID=/usr/local/parasolid/latest
#source /usr/local/parasolid/latest/parasolid.env
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/parasolid/latest/shared_object
