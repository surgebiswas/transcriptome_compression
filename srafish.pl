#!/usr/bin/perl

# srafish.pl - v0.13
# by Surge Biswas, Konstantin Kerner

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd;
use File::Path;

# Constants
my $EXECUTE = 0;

my $out = getcwd;
my $query_table;
my $do_not_process_list;
my $help;
my $download_protocol;
my $index;
my $ssh_key;
my $nthreads;


GetOptions(
	"out=s"	=> \$out,
	"table=s" => \$query_table,
	"index=s" => \$index,
	"ssh_key=s" => \$ssh_key,
	"nthreads=s" => \$nthreads,
	"do_not_process=s" => \$do_not_process_list,
	"protcl=s" => \$download_protocol,
	"help" => \$help,
);

die <<USAGE

USAGE: srafish.pl 
	-t  Path to a query table (required)
	-i  Sailfish index path (required)
	-s  SSH key for aspera (required)
	-n  Number of threads to use when running Sailfish (default: 6)
	-o  Path to output directory (default: working directory)
	-d  Path to file that contains a list of SRA IDs that are NOT to be processed (optional)
	-p  Download protocol to be used (-d aspera|ftp - default:aspera)
	-h  Print this help message
	
USAGE
if $help;


## 0. Parameter checks and echo.
my $errmessage = "";
$errmessage = $errmessage . "Path to query table required.\n" unless $query_table;
$errmessage = $errmessage . "Sailfish index path required.\n" unless $index;
$errmessage = $errmessage . "SSH key for aspera required.\n" unless $ssh_key or $download_protocol ne "aspera";
die $errmessage if $errmessage;

$nthreads = 6 unless $nthreads;
$download_protocol = "aspera" unless $download_protocol;
$ssh_key = "NA" if $download_protocol ne "aspera";
$out =~ s/\/$//g;

print "Query table:             $query_table\n";
print "Output directory:        $out\n";
print "Sailfish index:          $index\n";
print "SSH Key:                 $ssh_key\n";
print "Do not process list:     $do_not_process_list\n" if $do_not_process_list;
print "Do not process list:     not provided\n" unless $do_not_process_list;
print "Download protocol:       $download_protocol\n";



## BEGIN MAIN CODE BODY ##

## 1. Prepare query table.
my $num_queries = 0;
open QT, "<$query_table" or die $!; <QT>; $num_queries++ while <QT>; close QT;
print  "\nQuery table found at $query_table! $num_queries queries identified\n";
	
## 2. Read in the 'do not process' list.
my $do_not_process = read_do_not_process_list($do_not_process_list);

## 3. Download, fastq-dump, and Sailfish query results.
my $i = 1;
my ($run, $qdetails);
mkdir($out) unless -d $out;

open my $query_results, "$query_table" or die $!; <$query_results>; #skip header
while (<$query_results>) {
	
	$qdetails = query_details($_, $do_not_process, $out, 1);
	$run = $qdetails -> {run_id};
	
	if ($qdetails -> {should_process}){
		        
		run_fishhook($out, $run, $download_protocol, $ssh_key, $EXECUTE);
        
		run_sailfish($out, $run, $index, $nthreads, $EXECUTE);
        		
		if (processed_to_completion($out, $run)) {
			print  "Run $run finished successfully!\n";
		} else {
			print "WARNING: Run $run was unsuccessful. Moving on ...\n";
		}
		
	} else {
		print "Not processing $run.\n";
	}
	
	print "$i of $num_queries records completed\n";
	$i++;
}

## END MAIN CODE BODY ##



################################ SUBROUTINES ##################################
sub read_do_not_process_list {
	my $dnpl = shift;
	my %do_not_process;
	if ($dnpl) {
		open PF, "$dnpl" or die $!;
		
		while (<PF>) {
			
			chomp;
			$do_not_process{$_} = 1;
		}
		close PF;
	}
	return \%do_not_process;
}

sub run_fishhook {
	my $out = shift;
	my $run = shift;
	my $download_protocol = shift;
	my $ssh_key = shift;
	my $EXECUTE = shift;
	my $cmd;
	
	mkdir("$out/$run") unless -d "$out/$run";
	
	# Call fishhook.pl to download and unpack sra files
	$cmd = "fishhook.pl $run $out/$run $download_protocol $ssh_key";
	print "EXECUTING: $cmd\n";
	system($cmd) if $EXECUTE;
}
sub run_sailfish {
	my $out = shift;
	my $run = shift;
	my $index = shift;
	my $nt = shift;
	my $EXECUTE = shift;
	
	my $cmd;
	if (-e "$out/$run/$run\_2.fastq") {
		# Data is paired end.
		$cmd = "sailfish quant -i $index -l 'T=PE:O=><:S=U' -1 $out/$run/$run\_1.fastq -2 $out/$run/$run\_2.fastq -o $out/$run -p $nt";
		print "EXECUTING: $cmd\n";
		system($cmd) if $EXECUTE;
		
		system("rm $out/$run/$run\_1.fastq") if $EXECUTE;
		system("rm $out/$run/$run\_2.fastq") if $EXECUTE;
	}
	else {
		# Data is single end.
		$cmd = "sailfish quant -i $index -l 'T=SE:S=U' -r $out/$run/$run\_1.fastq -o $out/$run -p $nt";
		print "EXECUTING: $cmd\n";
		system($cmd) if $EXECUTE;
		
		system("rm $out/$run/$run\_1.fastq") if $EXECUTE;
	}
	
	# reads.sfc is a reasonably large file produced by sailfish, which we don't need for downstream analysis
	system("rm $out/$run/reads.sfc") if $EXECUTE;
}

sub query_details{
    my $query_table_line = shift;
    my $do_not_process = shift;
    my $out = shift;
    my $print_details = shift;
    my $platformOK = 0;
    my $dnp = 0;
    my $ap = 0;

    
    
    chomp($query_table_line);
    
    # replace dates within quotes with empty string. These interfere with comma parsing.
    $query_table_line =~ s/\".+?\"//g;
    
    my @line = split(/,/, $query_table_line);
    
    my $run = $line[0];
    my $library_name = $line[11];
    my $layout = $line[15];
    my $platform = $line[18];
    
    
    # Make sure its an Illumina run.
    $platformOK = 1 if lc $platform eq lc "ILLUMINA";
        
    # Check if this run is in the "do not process" list.
    $dnp = 1 if exists($do_not_process -> {$run});
    
    # Check if this run already has been processed.
    $ap = processed_to_completion($out, $run);
    
    my %query_details;
    $query_details{run_id} = $run;
    $query_details{platformOK} = $platformOK;
    $query_details{do_not_process} = $dnp;
    $query_details{already_processed} = $ap;
    $query_details{should_process} = $platformOK & !$ap & !$dnp;
    
    if ($print_details) {
	print "\n[$run]\n";
	print "Platform: $platform\n";
	print "Listed as do not process? YES\n" if $dnp;
	print "Listed as do not process? NO\n" unless $dnp;
	print "Already processed? YES\n" if $ap;
	print "Already processed? NO\n" unless $ap;
	print "Should process? YES\n" if $query_details{should_process};
	print "Should process? NO\n" unless $query_details{should_process};
    }
	
    return \%query_details;
}


sub processed_to_completion {
	my $complete = 0;
	my $outdir = shift;
	my $rundir = shift;
	my $filesize;
	
	if (-e "$outdir/$rundir/quant_bias_corrected.sf") {
			
		$filesize =
		$complete = 1 if -s "$outdir/$rundir/quant_bias_corrected.sf" >= 500000;
	}
	return $complete;
}

###############################################################################