  export VERSION=200_memLS

# uncomment line below to clean build
  isclean="clean"

  setup="gmake NODEP=1 setup"
  compile="gmake VERSION=200_memLS  NODEP=1 NOSHARED=1 $isclean"

  dest_path=$DEVROOT/phasta/phastaIO/phastaIO
  cd $dest_path
  $setup
  $compile

  dest_path=$DEVROOT/phasta/shapeFunction/shapeFunction
  cd $dest_path
#  $setup
  $compile

#  dest_path=$DEVROOT/phasta/phMetis/phMetis
#  cd $dest_path
#  $setup
#  $compile

  dest_path=$DEVROOT/phasta/phNSpre/phNSpre
  cd $dest_path
  $setup
  $compile

  dest_path=$DEVROOT/phasta/phPost/phPost/Reduce
  cd $dest_path
#  $setup
  $compile

  dest_path=$DEVROOT/phasta/phSolver/$VERSION/phSolver/
  cd $dest_path
##  $setup
  $compile

  cd $DEVROOT

