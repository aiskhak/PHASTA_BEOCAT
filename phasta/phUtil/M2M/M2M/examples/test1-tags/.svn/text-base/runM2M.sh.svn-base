#!/bin/sh

EXE_PATH=/users/osahni/develop/phasta/phUtil/M2M

logfile=M2M.log
#logfile=/dev/stdout

src_dir=.
dest_dir=.

echo
echo "Running M2M solution transfer (make sure to check $logfile)..."
echo "1 1689 2 1" | $EXE_PATH/bin/x86_64_linux-icc/M2M-O $src_dir/geom.xmt_txt $src_dir/geom.sms $src_dir/restart_src.10000.1 $src_dir/idmap.dat $dest_dir/geom.xmt_txt $dest_dir/geom.sms $dest_dir/restart_dest.10000.1 $dest_dir/idmap.dat 1 0 0 1 > $logfile 2>&1
echo
