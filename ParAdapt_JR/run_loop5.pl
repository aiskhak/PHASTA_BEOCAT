#!/usr/bin/perl
#
# name: run_loop.pl
# created by: Joe Rodriguez, April 2006
####################################################
#
# This Perl script automates the iterative running of
# the PHASTA and phParAdapt codes.  
#
#
#
####################################################
#
# Modules used
use Cwd;
use File::Copy;
#
# Set global variables and paths
#
	$loop_start = 1;	# start number of loop ID
	$nloops = 4;		# number of iterations to perform
	$nprocs = 2;		# number of processors
	$procs_dir = "$nprocs-procs_case";	# name of output directory
#	$rundir = "/home/navy/rodrigu/work/test_ic_loop";
	$rundir = `pwd`;
	$tsleep = 15;
	$phasta_exec = "phastaIC.exe-mpigm-O";
	$phparadapt_exec = "phParAdapt_102705_jmr2-Parasolid-mpigm-O";
#
	$phasta_path = "/home/navy/rodrigu/develop/phasta/phSolver/195_all6/phSolver/bin/x86_linux-icc_8132";
	$mpi_exec = "mpirun.lsf --gm-copy-env";
	$phasta_status_file = "$procs_dir/status.dat";
#
# open output debug file
#
        open(OUTFILE,">loop.out");
#
# Loop over number of loops
#
	for ($iloop=$loop_start; $iloop<=$nloops; $iloop++) {
		print OUTFILE "\n--------------------------------------\n";
		print OUTFILE " LOOP #$iloop\n";
		print  " Loop #$iloop\n";
# run phasta
		run_phasta();
# Find last timestep
		last_timestep();
# Setup for phParAdapt
		setup_adapt($timestep);
# Run phParAdapt
		run_phparadapt($timestep);
# Setup for Phasta
		setup_phasta($timestep);
	}
	close (OUTFILE);


#------------------------------------------------------------------------------
sub run_phasta
{
#
# Filelist to copy to trun directory
#
	$filelist_phasta_1[0] = "dumpnew.dat";
	$filelist_phasta_1[1] = "echo.dat";
	$filelist_phasta_1[2] = "forces.dat";
	$filelist_phasta_1[3] = "fort.76";
	$filelist_phasta_1[4] = "fort.881";
	$filelist_phasta_1[5] = "histor.dat";
	$filelist_phasta_1[6] = "spmass.dat";
	$filelist_phasta_2[0] = "solver.inp";
	$filelist_phasta_2[1] = "myjob.err";
	$filelist_phasta_2[2] = "myjob.out";
	$file_count_1 = @filelist_phasta_1;
	$file_count_2 = @filelist_phasta_2;
#
# remove phasta status file
#	print "status file: $phasta_status_file\n";
	if (-s $phasta_status_file) {
		print OUTFILE "  Deleting $phasta_status_file\n";	
		print " deleting $phasta_status_file\n";
		unlink ($phasta_status_file );
	}
# run phasta
	$strcmd = "$mpi_exec ./$phasta_exec";
	print OUTFILE "\n  Executing: $strcmd\n";
	print "executing: $strcmd\n";
	system($strcmd);
# loop until see return file, phasta_status
	$istop = 0;
	while ($istop == 0) {
		$tmpdir = `pwd`;
		print OUTFILE "    ...waiting for $phasta_status_file to exist\n";
		print "In $tmpdir and waiting for $phasta_status_file to exist\n";
		sleep $tsleep;
		if (-s $phasta_status_file ) {
			print OUTFILE "    ...found $phasta_status_file\n";
			print " found $phasta_status_file\n";
			$istop = 1;
		}
	}
	sleep $tsleep;
#
# gather results to trun directory
# first need to determine the highest trun* directory
#
	$trun_num = 0;
	@trunlist = `ls -dc1 $procs_dir/trun*`;
	$count = @trunlist;
	for ($itrun=1; $itrun<=$count; $itrun++) {
		@trunseg = split(/trun/,$trunlist[$itrun-1]);
# 		print "trunseg = @trunseg\n";
		if ($trunseg[1] > $trun_num) {
			$trun_num = $trunseg[1];
		}
	}
	$trun_num++;
	$trun_dir_name = "$procs_dir/trun$trun_num";
	print OUTFILE "\n  Creating trun directory: $trun_dir_name\n";
#	print "   creating trun directory: $trun_dir_name\n";
#
	system ("mkdir $trun_dir_name");
	print OUTFILE "   ...copying files to $trun_dir_name\n";
	for ($i=1; $i<@filelist_phasta_1; $i++) {
		copy_file($filelist_phasta_1[$i-1],"./$procs_dir",$filelist_phasta_1[$i-1],"./$trun_dir_name");
	}
	for ($i=1; $i<@filelist_phasta_2; $i++) {
		copy_file($filelist_phasta_2[$i-1],".",$filelist_phasta_2[$i-1],"./$trun_dir_name");
	}
	print OUTFILE " done\n";
#
# renaming phasta output directory
#
	$procs_dir_new = "$procs_dir"."_loop_"."$iloop";
	$strcmd = "mv $procs_dir $procs_dir_new";
	system($strcmd);
#
}  # end of subroutine run_phasta


#------------------------------------------------------------------------------
sub last_timestep
{
#
# Determine the last saved restart file
#
	print OUTFILE "\n  Finding last saved time step ..";
	$timestep = 0;
	@filelist = `ls -dc1 $procs_dir_new/restart.[1-9]*.1`;
	$count = @filelist;
#	print "filelist = @filelist\n";
	for ($ifile=1; $ifile<=$count; $ifile++) {
		@fileseg = split(/[.]/,$filelist[$ifile-1]);
#		print "fileseg = $fileseg[0], $fileseg[1], $fileseg[2]\n";
		if ($fileseg[1] > $timestep) {
			$timestep = $fileseg[1];
		}
	}
	print OUTFILE " done: time step = $timestep\n";
}  # end of subroutine last_timestep


#------------------------------------------------------------------------------
sub setup_adapt
{
	print OUTFILE "\n  Setting up files for adaption ... ";
# command line arguments
	my $timestepnum = @_[0];
#
# copy necessary files from PHASTA output 
# restart.$timestepnum.* & restart.$timestepnum.*
#
	for (my $iproc=1; $iproc<=$nprocs; $iproc++) {
		$restartname = "restart.$timestepnum.$iproc";
		$errorname = "error.$timestepnum.$iproc";
		copy_file($restartname,"./$procs_dir_new",$restartname,".");
		copy_file($errorname,"./$procs_dir_new",$errorname,".");
	}
#
# Setup geometry directories
# In general, the existing refined_p.sms directory applies to the previously run
# PHASTA files and should be renamed geom_p.sms for the adaption.  
# Therefore, 
#	1) mv refiend_p.sms geom_p.sms
#	2) cp -r geom_p.sms $procs_dir/  - to save pmesh files
# If this the first adaption, there will not be a refined_p.sms directory.  In
# this case, simply execute step (2) above
	if (-s "refined_p.sms") {
		system (" rm -r geom_p.sms");
		system (" mv refined_p.sms geom_p.sms");
	}
	copy_directory("geom_p.sms",$procs_dir_new);
#
# Modify the timestep number in the adapt.inp file
#
	open (ADAPTINP,"adapt.inp");
	open (ADAPTINPTMP,">adapt.inp_tmp");
	while (my $line = <ADAPTINP>) {
		chomp ($line);
		@lineseg = split(/ +/,$line);
		if ($lineseg[0] eq "timeStepNumber") {
			$line = "timeStepNumber $timestepnum";
		}
		print ADAPTINPTMP "$line\n";
	}
	close (ADAPTINP);
	close (ADAPTINPTMP);
#
	$strcmd = "cp adapt.inp_tmp adapt.inp";
	system ($strcmd);
	unlink ("adapt.inp_tmp");
#
	print OUTFILE " done \n";
}  # end of subroutine setup_adapt


#------------------------------------------------------------------------------
sub run_phparadapt
{
# command line arguments
        my $timestepnum = @_[0];
#
	print OUTFILE "\n  Executing phParAdapt ... \n";
# run phParAdapt
	$strcmd = "$mpi_exec ./$phparadapt_exec";
	print OUTFILE "\n  Executing: $strcmd\n";
	print "executing: $strcmd\n";
	system($strcmd);
# loop until see the geombc.dat.<proc> file created and no longer changing size
	$geomfile = "geombc.dat.$nprocs";
	$phparadapt_status_file = "$timestepnum/$procs_dir/$geomfile";
	$istop = 0;
	$size = 0;
	while ($istop == 0) {
		print OUTFILE "    ...waiting for $phparadapt_status_file to exist\n";
		print "   ...waiting for $phparadapt_status_file to exist\n";
		sleep $tsleep;
		if (-s $phparadapt_status_file) {
			$size_old = $size;
			$size = -s $phparadapt_status_file;
			if ($size == $size_old) {
				print OUTFILE "    ...found $phparadapt_status_file & finished writing\n";
				print " found $phparadapt_status_file & finished writing\n";
				$istop = 1;
			} else {
				print OUTFILE "    ...found $phparadapt_status_file but not finished writing\n";
				print " found $phparadapt_status_file but not finished writing\n";
			}
		}
	}
	sleep $tsleep;
#
# move adaption results to new directory
#
	$refine_dir = "refine_"."$procs_dir"."_loop_$iloop";
	system ("mkdir $refine_dir");
	print OUTFILE "    ...created directory $refine_dir\n";
	print OUTFILE "    ...moving results to $refine_dir ... \n";
# move results to new directory
	for (my $iproc=1; $iproc<=$nprocs; $iproc++) {
		$filename = "phAdapt.$iproc.log";
		move_file($filename,".",$filename,"./$refine_dir");
		$filename = "restart.$timestepnum.$iproc";
		move_file($filename,".",$filename,"./$refine_dir");
		$filename = "error.$timestepnum.$iproc";
		move_file($filename,".",$filename,"./$refine_dir");
	}
	copy_file("adapt.inp",".","adapt.inp","./$refine_dir");
	copy_directory("refined_p.sms",$refine_dir);
	move_directory($timestepnum,$refine_dir);
	print OUTFILE " done \n";
#
}  # end of subroutine run_phparadapt


#------------------------------------------------------------------------------
sub setup_phasta
{
# command line arguments
        my $timestepnum = @_[0];
#
	print OUTFILE "\n  Setting up files for phasta ... \n";
# create run directory
#	system ("mkdir $procs_dir");
# copy contents of new run directory from adapted mesh
	copy_directory("$refine_dir/$timestepnum/$procs_dir",".");
# modify solver.inp for time 
#
	print OUTFILE "\n     ... done \n";
#
}  # end of subroutine setup_phasta


#------------------------------------------------------------------------------
sub copy_file
{
# copies a file and returns control once copy has finished
	my $filename_orig = "@_[1]/@_[0]";
	my $filename_dest = "@_[3]/@_[2]";
	
	print OUTFILE "    copying $filename_orig to $filename_dest";
	
	my $size_orig = -s $filename_orig;
	my $size_dest = 0;
	
	copy($filename_orig,$filename_dest);
	while ($size_dest != $size_orig) {
		sleep (5);
		$size_dest = -s $filename_dest;
	}
	print OUTFILE " ... done\n";
	print OUTFILE "     ... size: orig = $size_orig  dest = $size_dest\n";

}	# end of subroutine copy_file


#------------------------------------------------------------------------------
sub move_file
{
# moves a file and returns control once move has finished
	my $filename_orig = "@_[1]/@_[0]";
	my $filename_dest = "@_[3]/@_[2]";
	
	print OUTFILE "    moving $filename_orig to $filename_dest";
	
	my $size_orig = -s $filename_orig;
	my $size_dest = 0;
	
	move($filename_orig,$filename_dest);
	while ($size_dest != $size_orig) {
		sleep (5);
		$size_dest = -s $filename_dest;
	}
	print OUTFILE " ... done\n";
	print OUTFILE "     ... size: orig = $size_orig  dest = $size_dest\n";

}	# end of subroutine move_file


#------------------------------------------------------------------------------
sub copy_directory
{
# copies a file and returns control once copy has finished
	my $dirname_orig = "@_[0]";
	my $dirname_dest = "@_[1]";
	my $dirname, $filename;
	my $idir, $ifile, $i;
	
	print OUTFILE "   Copying $dirname_orig to $dirname_dest\n";
#
# 
# Get list of all subdirectories in directory
# Using the 'ls -R <dirname> | grep "/"' command
# must remove colon at end of subdirectory names
#
	my @subdir = `ls -R $dirname_orig | grep "/"`;
	for ($idir=1; $idir <= @subdir; $idir++) {
		chomp ($subdir[$idir-1]);
		@subdirseg = split (/[:]/,$subdir[$idir-1]);
		$subdir[$idir-1] = $subdirseg[0];
	}
#
# The @subdir array contains pathnames of the subdirectories w.r.t.
# to the current run directory.  For example:
#	if $dirname_orig = dir1/subdir2/subsubdir3 then the @sublist might
#	look like:
#		dir1/subdir2/subsubdir3
#		dir1/subdir2/subsubdir3/subsubsubdir1
#		dir1/subdir2/subsubdir3/subsubsubdir2
# 	if we are copying to the directory $dirname_dest = dir2 then we want to
#	create the following directories in dir2/:
#		dir2/subsubdir3
#		dir2/subsubdir3/subsubsubdir1
#		dir2/subsubdir3/subsubsubdir2
# Therefore we need a second @sublist2 array that lists the subdirectories w.r.t.
# the bottom directory of $dirname_orig (i.e. sudsubdir3 in our example).
# @sublist2 would look like:
#		subsubdir3
#		subsubdir3/subsubsubdir1
#		subsubdir3/subsubsubdir2
#
# create new @sublist2 array
#
	my @nameseg = split(/\//,$dirname_orig);
	if (@nameseg > 1) {
		$header_name = "$nameseg[0]/";
		for ($i=2; $i<@nameseg; $i++) {
			$header_name = "$header_name$nameseg[$i-1]/";
		}
		print OUTFILE "headername = $header_name\n";
		$dirname_orig2 = truncate_header($dirname_orig,$header_name);
		for ($idir=1; $idir <= @subdir; $idir++) {
			$subdir2[$idir-1] = truncate_header($subdir[$idir-1],$header_name);
		}
	} else {
		$dirname_orig2 = $dirname_orig;
		for ($idir=1; $idir <= @subdir; $idir++) {
			$subdir2[$idir-1] = $subdir[$idir-1];
		}
	}
#
# Create directory structure in new directory
	
	$dirname = "$dirname_dest/$dirname_orig2";
	print OUTFILE "    creating directory $dirname\n";
	system ("mkdir -p $dirname");
	for ($idir=1; $idir <= @subdir2; $idir++) {
		$dirname = "$dirname_dest/$subdir2[$idir-1]";
		print OUTFILE "    creating directory $dirname\n";
		system ("mkdir -p $dirname");
	}

# Copy all files from the main directory

	my @files = glob ("$dirname_orig/*");
#	print OUTFILE "    files in $dirname_orig = @files\n";
	for ($ifile=1; $ifile <= @files; $ifile++) {
		if (-d $files[$ifile-1]) {
			print OUTFILE "    skipping $files[$ifile-1] - is a directory\n";
		} else {
			if ($dirname_orig2 eq $dirname_orig) {
				$filename = $files[$ifile-1];
			} else {
				$filename = truncate_header($files[$ifile-1],$header_name);
			}
			$dirname = "$dirname_dest/$dirname_orig2";
			copy_file($files[$ifile-1],".",$filename,$dirname_dest);
		}
	}
	
# Copy all files from the subdirectories
	for ($idir=1; $idir <= @subdir; $idir++) {
		$strcmd = "$subdir[$idir-1]/*";
		@files = glob ("$strcmd");
#		print OUTFILE "    files in $dirname_orig = @files\n";
		for ($ifile=1; $ifile <= @files; $ifile++) {
			if (-d $files[$ifile-1]) {
				print OUTFILE "    skipping $files[$ifile-1] - is a directory\n";
			} else {
				if ($dirname_orig2 eq $dirname_orig) {
					$filename = $files[$ifile-1];
				} else {
					$filename = truncate_header($files[$ifile-1],$header_name);
				}
				copy_file($files[$ifile-1],".",$filename,$dirname_dest);
			}
		}
	}
}	# end of subroutine copy_directory


#------------------------------------------------------------------------------
sub move_directory
{
# This is identical to the copy_directory subroutine except the directory is
# deleted after the copy is finished.
	my $dirname_orig = "@_[0]";
	my $dirname_dest = "@_[1]";
#
	print OUTFILE "   Moving $dirname_orig to $dirname_dest\n";
# First copy the directory
	print OUTFILE "    ... first copy $dirname_orig to $dirname_dest\n";
	copy_directory($dirname_orig,$dirname_dest);
# now delete
	print OUTFILE "    ... now delete $dirname_orig\n";
	$strcmd = "rm -r $dirname_orig";
	system ("$strcmd");
	
}	# end of subroutine move_directory


#------------------------------------------------------------------------------
sub truncate_header
{
# @_[0] = original string
# @_[1] = string to remove from original string
#
	my $tmpstr = @_[0];
	substr($tmpstr,index(@_[0],@_[1]),length(@_[1]))="";
	print OUTFILE "     ...truncating \"@_[1]\" from \"@_[0]\" = $tmpstr\n";
	return ($tmpstr);
}
# This is identical to the
