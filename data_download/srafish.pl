#!/usr/bin/perl

# srafish.pl - v1.13
# Now runs on kallisto
# by Surge Biswas, Konstantin Kerner

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd;
use File::Path;

# Constants
my $EXECUTE = 1;

my $out = getcwd;
my $query_table;
my $do_not_process_list;
my $help;
my $download_protocol;
my $index;
my $ssh_key;
my $nthreads;
my $min_num_reads;
my $max_num_reads;


GetOptions(
	"out=s"	=> \$out,
	"table=s" => \$query_table,
	"index=s" => \$index,
	"ssh_key=s" => \$ssh_key,
	"nthreads=s" => \$nthreads,
	"do_not_process=s" => \$do_not_process_list,
	"protcl=s" => \$download_protocol,
	"min_num_reads=s" => \$min_num_reads,
	"a_max_num_reads=s" => \$max_num_reads,
	"help" => \$help,
);

die <<USAGE

USAGE: srafish.pl 
	-t  Path to a query table (required)
	-i  Kallisto index path (required)
	-s  SSH key for aspera (required)
	-n  Number of threads to use when running Kallisto (default: 1)
	-o  Path to output directory (default: working directory)
	-d  Path to file that contains a list of SRA IDs that are NOT to be processed (optional)
	-p  Download protocol to be used (-d aspera|ftp - default:aspera)
	-m  Minimum number of reads/pairs a sample should have in order to be processed (default: 4e6)
	-a  Maximum number of reads/pairs a sample should have. Larger ones will be downsampled (default: 40e6)
	-h  Print this help message
	
USAGE
if $help;


## 0. Parameter checks and echo.
my $errmessage = "";
$errmessage = $errmessage . "Path to query table required.\n" unless $query_table;
$errmessage = $errmessage . "Kallisto index path required.\n" unless $index;
$errmessage = $errmessage . "SSH key for aspera required.\n" unless $ssh_key or $download_protocol ne "aspera";
die $errmessage if $errmessage;

$nthreads = 1 unless $nthreads;
$download_protocol = "aspera" unless $download_protocol;
$ssh_key = "NA" if $download_protocol ne "aspera";
$min_num_reads = 4000000 unless $min_num_reads;
$max_num_reads = 40000000 unless $max_num_reads;
$out =~ s/\/$//g;

print "Query table:             $query_table\n";
print "Output directory:        $out\n";
print "Kallisto index:          $index\n";
print "SSH Key:                 $ssh_key\n";
print "Do not process list:     $do_not_process_list\n" if $do_not_process_list;
print "Do not process list:     not provided\n" unless $do_not_process_list;
print "Download protocol:       $download_protocol\n";
print "Min. num. reads:         $min_num_reads\n";
print "Max. num. reads:         $max_num_reads\n";



## BEGIN MAIN CODE BODY ##

## 1. Prepare query table.
my $num_queries = 0;
open QT, "<$query_table" or die $!; <QT>; $num_queries++ while <QT>; close QT;
print  "\nQuery table found at $query_table! $num_queries queries identified\n";
	
## 2. Read in the 'do not process' list.
my $do_not_process = read_do_not_process_list($do_not_process_list);

## 3. Download, fastq-dump, and Kallisto query results.
my $i = 1;
my ($run, $qdetails);
execute_cmd("mkdir $out", $EXECUTE) unless -d $out;

open my $query_results, "$query_table" or die $!; <$query_results>; #skip header
while (<$query_results>) {
	
	$qdetails = query_details($_, $do_not_process, $min_num_reads, $out, 1);
	$run = $qdetails -> {run_id};
	
	if ($qdetails -> {should_process}){
		        
		# Download and unpack the data.
		run_fishhook($out, $run, $download_protocol, $ssh_key, $EXECUTE);
		
		# Subsample reads if the file is too large.
		sub_sample($out, $run, $qdetails, $max_num_reads, $EXECUTE) if $qdetails -> {nreads} > $max_num_reads;
		
		# Run sailfish to quantify transcript abundances.
		#run_sailfish($out, $run, $index, $nthreads, $EXECUTE);
        	run_kallisto($out, $run, $index, $nthreads, $EXECUTE);
		
		if (processed_to_completion($out, $run)) {
            		# Remove the first 6 (header) lines from quant_bias_corrected.sf to
            		# make it easier to parse this file downstream
            		#remove_header_lines_from_file("$out/$run/quant_bias_corrected.sf", 6, $EXECUTE);
            
			print "Run $run finished successfully!\n";
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
sub execute_cmd {
	my $cmd = shift;
	my $exec = shift;
	
	print "EXECUTING: $cmd\n";
	system($cmd) if $exec;
}

sub sub_sample {
	my $out = shift;
	my $run = shift;
	my $qd = shift;
	my $maxn = shift;
	my $EXECUTE = shift;
	
	my $freq = $maxn / $qd->{nreads};

	if (-e "$out/$run/$run\_2.fastq") {
		# Data is paired end
		execute_cmd("subsample_fastq.pl $freq $out/$run/$run\_1.fastq $out/$run/$run\_2.fastq", $EXECUTE);
		
		# Rename the subsampled files back to the original names.
		execute_cmd("mv $out/$run/$run\_1.fastq.sub $out/$run/$run\_1.fastq; mv $out/$run/$run\_2.fastq.sub $out/$run/$run\_2.fastq;", $EXECUTE)
	} else {
		# Data is single end
		execute_cmd("subsample_fastq.pl $freq $out/$run/$run\_1.fastq", $EXECUTE);
		
		# Rename the subsampled file back to the original name.
		execute_cmd("mv $out/$run/$run\_1.fastq.sub $out/$run/$run\_1.fastq", $EXECUTE);
	}
}

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
	
	execute_cmd("mkdir $out/$run", $EXECUTE) unless -d "$out/$run";
	
	# Call fishhook.pl to download and unpack sra files
	execute_cmd("fishhook.pl $run $out/$run $download_protocol $ssh_key", $EXECUTE);
}

sub run_kallisto {
	my $out = shift;
	my $run = shift;
	my $index = shift;
	my $nt = shift;
	my $EXECUTE = shift;
	
	my $cmd;
	if (-e "$out/$run/$run\_2.fastq") {
		# Data is paired end.
		execute_cmd("kallisto quant -i $index -o $out/$run -t $nt $out/$run/$run\_1.fastq $out/$run/$run\_2.fastq", $EXECUTE);
		execute_cmd("rm $out/$run/$run\_1.fastq; rm $out/$run/$run\_2.fastq", $EXECUTE);
	}
	else {
		# Data is single end.
		execute_cmd("kallisto quant -i $index -o $out/$run --single -l 200 -s 40 -t $nt $out/$run/$run\_1.fastq", $EXECUTE);
		execute_cmd("rm $out/$run/$run\_1.fastq", $EXECUTE);
	}
	
	# reads.sfc is a reasonably large file produced by sailfish, which we don't need for downstream analysis
	execute_cmd("rm $out/$run/abundance.h5", $EXECUTE);
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
		execute_cmd("sailfish quant -i $index -l 'T=PE:O=><:S=U' -1 $out/$run/$run\_1.fastq -2 $out/$run/$run\_2.fastq -o $out/$run -p $nt &> $out/$run/sf.out", $EXECUTE);
		execute_cmd("rm $out/$run/$run\_1.fastq; rm $out/$run/$run\_2.fastq; rm $out/$run/sf.out", $EXECUTE);
	}
	else {
		# Data is single end.
		execute_cmd("sailfish quant -i $index -l 'T=SE:S=U' -r $out/$run/$run\_1.fastq -o $out/$run -p $nt &> $out/$run/sf.out", $EXECUTE);
		execute_cmd("rm $out/$run/$run\_1.fastq; rm $out/$run/sf.out", $EXECUTE);
	}
	
	# reads.sfc is a reasonably large file produced by sailfish, which we don't need for downstream analysis
	execute_cmd("rm $out/$run/reads.sfc", $EXECUTE);
}

sub query_details{
    my $query_table_line = shift;
    my $do_not_process = shift;
    my $min_reads = shift;
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
    my $nreads = $line[3];
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
    $query_details{nreads} = $nreads;
    $query_details{should_process} = $platformOK & !$ap & !$dnp & ($nreads >= $min_reads);
    
    if ($print_details) {
	print "\n[$run]\n";
	print "Platform: $platform\n";
	print "Num. reads: $nreads\n";
	print "Listed as do not process? YES\n" if $dnp;
	print "Listed as do not process? NO\n" unless $dnp;
	print "Already processed? YES\n" if $ap;
	print "Already processed? NO\n" unless $ap;
	print "Insufficient reads? NO\n" if $nreads >= $min_reads;
	print "Insufficient reads? YES\n" unless $nreads >= $min_reads;
	print "Should process? YES\n" if $query_details{should_process};
	print "Should process? NO\n" unless $query_details{should_process};
    }
	
    return \%query_details;
}

sub processed_to_completion {
	my $complete = 0;
	my $outdir = shift;
	my $rundir = shift;
	
	if (-e "$outdir/$rundir/quant_bias_corrected.sf") {
		$complete = 1 if -s "$outdir/$rundir/quant_bias_corrected.sf" >= 500000;
	}
	return $complete;
}

sub remove_header_lines_from_file {
	
	my $file = shift; # file you want to change (full path)
	my $remlines = shift; # number of lines you want to remove from the beginning of the file
	my $EXECUTE = shift;

    execute_cmd("tail -n +$remlines $file > $file.tmp", $EXECUTE);
    execute_cmd("rm $file", $EXECUTE);
    execute_cmd("mv $file.tmp $file", $EXECUTE);
}
###############################################################################
