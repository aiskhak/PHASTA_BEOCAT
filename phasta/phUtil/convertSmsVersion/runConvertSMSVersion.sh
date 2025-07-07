#!/bin/sh

SIMARCH=x64_rhel3_gcc32
PARASOLID=/net/common/meshSim/6.2-080614/lib/$SIMARCH/psKrnl

echo
ls $PARASOLID
echo

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PARASOLID

EXE_PATH=/users/behaff/develop/phasta/phUtil/convertSmsVersion

f1=63
v1=61

f2=61
v2=2

#logfile=convertSMSVersion${f1}to${v1}.log
logfile=/dev/stdout

echo
echo "Converting sms version from $f1 to $v1 (make sure to check $logfile)..."
$EXE_PATH/from${f1}/convertSMSVersion geom.xmt_txt geom.sms geom_ver${v1}.sms $v1 > $logfile 2>&1
echo

#logfile=convertSMSVersion${f2}to${v2}.log
logfile=/dev/stdout

echo
echo "Converting sms version from $f2 to $v2 (make sure to check $logfile)..."
$EXE_PATH/from${f2}/convertSMSVersion geom.xmt_txt geom_ver${v1}.sms geom_ver${v2}.sms $v2 > $logfile 2>&1
echo

rm -rf geom_ver${v1}.sms
