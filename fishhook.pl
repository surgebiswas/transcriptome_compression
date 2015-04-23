#!/usr/bin/env perl

# fishhook.pl - v0.1
# by Surge Biswas, Konstantin Kerner

# intended as slave script for srafish.pl. can be used independently however.
# takes an SRR identifier and an output path as input. downloads the corresponding
# .sra file via aspera file transfer by calling ascp or via ftp by calling 
# fastq-dump directly.
# requires ascpera ascp client, and path to ascp executable stored in $PATH.


use strict;
use warnings;
use Data::Dumper;


my $srrid		 = shift;
my $outdir		 = shift;
my $method		 = shift;
my $BANDWIDTH	 = "100m";
my $cmd;

$outdir =~ s/\/$//g;
system("mkdir $outdir/$srrid");

	
if ($method eq "aspera") {
	
	my $OPENSSH_KEY = "/home/konstantin/perl/tradict/asperaweb_id_dsa.openssh";
#	my $OPENSSH_KEY = "/Users/wiggemacbook/Applications/Aspera\\ Connect.app/Contents/Resources/asperaweb_id_dsa.openssh";
	
	my $srr_pre_folder = substr $srrid, 0, 6;
	my $srr_pre_pre_folder = substr $srrid, 0, 3;
	
	# download sra file for specified srr ID with ascp
	$cmd = "ascp -i $OPENSSH_KEY -k1 -Tr -l$BANDWIDTH anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/$srr_pre_pre_folder/$srr_pre_folder/$srrid/$srrid.sra $outdir/$srrid";
	
	print "using aspera ascp client to download run $srrid\n";
	print "EXECUTING: $cmd\n";
	system($cmd);
	
	# unpack sra file with fastq-dump
	$cmd = "fastq-dump -I --split-files --outdir $outdir/$srrid $outdir/$srrid/$srrid.sra";
	print "EXECUTING: $cmd\n";
	system($cmd);
	
	# remove sra file
	system("rm $outdir/$srrid/$srrid.sra");
}

else {
	
	# download directly with fastq-dump (utilizes ftp protocol)
	print  "\nFetching run: $srrid . . .\n";
	$cmd = "fastq-dump -I --split-files --outdir $outdir/$srrid $srrid";
	print "EXECUTING: $cmd\n";
	system($cmd);
	
}
