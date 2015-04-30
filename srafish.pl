#!/usr/bin/perl

# srafish.pl - v0.13
# by Surge Biswas, Konstantin Kerner

######### UPDATE HISTORY ######################################################
# (20-04-15) v0.11 - added option to use an already existing query table.
#
#                    added subroutine to identify reference genotypes. if the
#                    keywords for one of the 19 MAGIC parents (e.g. "edi", case
#                    insensitive) is identified in column "library name", the
#                    respective reference matrix will be used. if no genotype
#                    can be identified, col reference will be used.
#
# (21-04-15) v0.12 - Added comments to code.
#
#                    Changed $experiment variable to
#                    $run because 'experiments' are SRX... and can be shared
#                    across samples. 'Runs' however are unique SRR or ERR IDs.
#
#                    Changed chdir calls from $cwd/$out to simply $out. I did
#                    this because if out is a relative path, it will be created
#                    as such (no need to specify the working directory). If
#                    its a full path, then having $cwd/$out doesn't make sense.
#
#                    Kept --split-files on for all fastq-dump calls.
#
#		             Added Sailfish call
#                
#                    added checkpoint for successful sailfish execution based
#                    on the composition of the sailfish logfile.
#
# (22-04-15) v0.13 - changed checkpoint 1: "check if run is already processed"
#                    into a subroutine.
#
#                    added checkpoint 2: "parse sailfish log".
#
#                    re-arranged data transfer. downloading .sra files will
#                    now be managed by the slave script fishhook.pl.
###############################################################################



use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd;
#use IPC::System::Simple qw(system);

my $cwd = getcwd;

my $MAGIC_INDEX_DIR 	= "/nas02/home/s/b/sbiswas/bin/Sailfish-0.6.3-Linux_x86-64/indexes";
my $NTHREADS 			= 6;
my $EXECUTE             = 0;
my $OPEN_SSH_KEY        = "/nas02/home/s/b/sbiswas/tradict/asperaweb_id_dsa.openssh";
my $query;
my $out 				= $cwd;
my $query_table;
my $help;
my $download_protocol;


GetOptions(
	"query=s"			=> \$query,
	"out=s"				=> \$out,
	"table=s"			=> \$query_table,
	"dlprotcl=s"			=> \$download_protocol,
	"help"				=> \$help,
);

die <<USAGE

USAGE: srafish.pl -q 'query' -o /path/to/out(optional)

    -q  specify a query for Entrez search, e.g. '"Homo sapiens"[Organism] AND "strategy rna seq"[Properties]' to build a query table
    OR
    -t  use an already existing query table
	
    -o  path of your output directory
    -d  download protocol to be used (-d aspera|ftp - default:aspera)
    -s  Full path to open ssh key for ascp client
    -h  print this help message
	
USAGE


## 0. Parameter checks.
unless $query or $query_table;
$download_protocol = "aspera" unless $download_protocol;
$out =~ s/\/$//g;



## 1. Prepare query table.
if ($query and $query_table) {
	die "query AND table specified. please choose only one of these two options.";
}


my $queries = 0;
my $query_count;
if ($query_table) {
    # User provided a query table.
    
	open $query_count, "<$query_table" or die $!;
	<$query_count>;
	$queries++ while <$query_count>;
	close $query_count;
	print  "\nquery table found at $query_table! $queries queries identified\n";
}

else {
    # User provided a query that we need to build a table from.
    
	$query_table = "$out/query_results.csv";
	print  "\nBuilding Query table . . .\n\n";	
	system("wget -O $query_table 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$query'") if $EXECUTE;
	
	
	open $query_count, "<$query_table" or die $!;
	<$query_count>;
	$queries++ while <$query_count>;
	close $query_count;
	print  "query table finished! $queries datasets identified\n";
}







## 2. Process (download, fastq-dump, and sailfish) query results.
open my $query_results, "<$query_table" or die $!;

<$query_results>; #skip header

my $i = 1;
my $run;
my $library_name;
my $layout;
my $genotype;
my $cmd;
my $failfish;
my $exists;

mkdir($out) unless -d $out;

while (<$query_results>) {

	chomp;
	
    # replace dates within quotes with empty string. These interfere with comma parsing.
    $_ =~ s/\".+?\"//g;

	my @line = split(/,/, $_);

	$run			= $line[0];
	$library_name	= $line[11];
	$layout			= $line[15];
	
    print "\n[$run]\n";
    
	# checkpoint 1: check if run is already processed
	system("mkdir $out/$run") if $EXECUTE;

	$exists = &already_processed($run);
	
	if ($exists){
		print "run $run already processed! proceeding to next run . . .\n";
	}
	
	else {
		$genotype = &determine_genotype ($library_name);
		
        # Get the data as a FASTQ. Always have --split-files on. In the case the
        # data is SE, it will only create a _1.fastq file.
		# call slave script fishhook.pl to download and unpack sra files
		
        $cmd = "./fishhook.pl $run $out/$run $download_protocol $OPEN_SSH_KEY";
        print "EXECUTING: $cmd\n";
		system($cmd) if $EXECUTE;

        # Run sailfish.
        print  "Running sailfish  . . .\n";
        if (-e "$out/$run/$run\_2.fastq") {
            # Data is paired end.
            $cmd = "sailfish quant -i $MAGIC_INDEX_DIR/$genotype/ -l 'T=PE:O=><:S=U' -1 $out/$run/$run\_1.fastq -2 $out/$run/$run\_2.fastq -o $out/$run -p $NTHREADS";
	    	print "EXECUTING: $cmd\n";
            system($cmd) if $EXECUTE;
	    
	    	system("rm $out/$run/$run\_1.fastq") if $EXECUTE;
	    	system("rm $out/$run/$run\_2.fastq") if $EXECUTE;
        }
        else {
            # Data is single end.
            $cmd = "sailfish quant -i $MAGIC_INDEX_DIR/$genotype/ -l 'T=SE:S=U' -r $out/$run/$run\_1.fastq -o $out/$run -p $NTHREADS";
	    	print "EXECUTING: $cmd\n";
            system($cmd) if $EXECUTE;
	    
	    	system("rm $out/$run/$run\_1.fastq") if $EXECUTE;
        }
        
		# checkpoint 2: parse sailfish log and match for keywords that will always
		# appear in the last line of a logfile from a successful sailfish
		$failfish = &check_logfile($out, $run);
		
		if ($failfish) {
			print "odd looking logfile for run $run. check sailfish logfile at $out/$run/logs.\n";
			print "proceeding anyway. ", ($i / $queries) * 100, " % done.\n";

		}
		else {
			print  "Run $run finished successfully! ", ($i / $queries) * 100, " % done.\n";
		}
		$i++;
	}
}


###############################################################################
###############################################################################

sub already_processed {
	
	my $processed = 0;
	my $rundir = shift @_;
	my @files_in_rundir;
	
	opendir DIR, "$out/$run/" or print "cannot open directory $out/$rundir\n";
	@files_in_rundir = grep { $_ ne '.' && $_ ne '..' } readdir DIR;
	closedir DIR;
	
	for my $file_in_rundir (@files_in_rundir) {
		
		if ($file_in_rundir =~ /quant.sf/) {
			$processed = 1;
			last;
		}
	}
	return $processed;
}


###############################################################################

sub determine_genotype {
	
	my $name = shift @_;

	my %known_genotypes = (
		"bur"	=> "Bur_0",
		"can"	=> "Can_0",
		"ct"	=> "Ct_1",
		"edi"	=> "Edi_0",
		"hi"	=> "Hi_0",
		"kn"	=> "Kn_0",
		"ler"	=> "Ler_0",
		"mt"	=> "Mt_0",
		"no"	=> "No_0",
		"oy"	=> "Oy_0",
		"po"	=> "Po_0",
		"rsch"	=> "Rsch_4",
		"sf"	=> "Sf_2",
		"tsu"	=> "Tsu_0",
		"wil"	=> "Wil_2",
		"ws"	=> "Ws_0",
		"wu"	=> "Wu_0",
		"zu"	=> "Zu_0",
	);
	
	my $final_genotype = "Col_0";
	
    # There are very few non Col-0 datasets. Even for non Col_0 accessions,
    # the mapping rates are mid-70% and mapping these to Col_0 gets you a
    # mapping rate of ~70%. Therefore, for the run on Killdevil we will map
    # to Col_0 exclusively.
    #
    # Commenting out the following:
    
#	for my $known_genotype (keys(%known_genotypes)) {
#		
##		print $known_genotype, "\n";
#		if ($name =~ /$known_genotype[^a-z]/i){
#			$final_genotype = $known_genotypes{$known_genotype};
#			last;
#		}
#	}

	return $final_genotype;
}



###############################################################################

sub check_logfile {
	
	my $logfile_status	 = 1;
	my $outdir		 = shift @_;
	my $logdir		 = shift @_;
	my @files_in_logdir;
	my $logfilename;
	
	if (-e "$outdir/$logdir/$run/logs/") {
			
		opendir LOG, "$outdir/$logdir/logs/" or print "cannot open directory $out/$logdir/logs/\n";
		@files_in_logdir = grep { $_ ne '.' && $_ ne '..' } readdir LOG;
		closedir LOG;
		
		for my $file_in_logdir (@files_in_logdir){
			if ($file_in_logdir =~ /sailfish\.g2log.+\.log/) {
				$logfilename = $file_in_logdir;
				last;
			} 
		}
		
		open my $logfile, "<$outdir/$logdir/logs/$logfilename" or print "cannot find logfile for run $run!\n";
		
		while (<$logfile>){
			$logfile_status = 0 if /g2log\ file\ shutdown/i;
		}
		close $logfile;
	}
	return $logfile_status;
}

###############################################################################
