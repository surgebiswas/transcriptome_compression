#!/usr/bin/perl

# srafish.pl - v0.12
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
###############################################################################



use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd;

my $cwd = getcwd;

my $MAGIC_INDEX_DIR = "/home/surge/applications/Sailfish-0.6.3-Linux_x86-64/indexes/magic";
my $NTHREADS = 6;
my $query;
my $out 				= $cwd;
my $query_table;
my $help;

GetOptions(
	"query=s"			=> \$query,
	"out=s"				=> \$out,
	"table=s"			=> \$query_table,
	"help"				=> \$help,
);

die <<USAGE

USAGE: srafish.pl -q 'query' -o /path/to/out(optional)

	-q		specify a query for Entrez search, e.g. '"Homo sapiens"[Organism] AND "strategy rna seq"[Properties]' to build a query table
	OR
	-t		use an already existing query table
	
	-o		path of your output directory
	-h		print this help message
	
USAGE


## 0. Parameter checks.
unless $query or $query_table;
$out =~ s/\/$//g;



## 1. Prepare query table.
if ($query and $query_table) {
	die "query AND table specified. please choose only one of these two options.";
}


my $queries = 0;
if ($query_table) {
    # User provided a query table.
    
	open my $query_count, "<$query_table" or die $!;
	<$query_count>;
	$queries++ while <$query_count>;
	close $query_count;
	print  "\nquery table found at $query_table! $queries queries identified\n";
}

else {
    # User provided a query that we need to build a table from.
    
	$query_table = "$out/query_results.csv";
	print  "\nBuilding Query table . . .\n\n";	
	system("wget -O $query_table 'http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term=$query'");
	
	
	open my $query_count, "<$query_table" or die $!;
	<$query_count>;
	$queries++ while <$query_count>;
	close $query_count;
	print  "query table finished! $queries datasets identified\n";
}







## 2. Process (fastq-dump and sailfish) query results.
open my $query_results, "<$query_table" or die $!;

<$query_results>; #skip header

my $i = 1;
my $run;
my $library_name;
my $layout;
my $genotype;
my $cmd;
my $check_success;

while (<$query_results>) {

	chomp;
	
    # replace dates within quotes with empty string. These interfere with comma parsing.
    $_ =~ s/\".+?\"//g;

	my @line = split(/,/, $_);

	$run = $line[0];
	$library_name = $line[11];
	$layout = $line[15];
	

 	
	# For now check if quant.sf has been output by sailfish.
	# If so, we'll assume it completed successfully and that this run
	# has already been processed. We might want to implement a more
	# rigorous check of the sailfish log later.
	if (-e "$out/$run/quant.sf") {
		print  "Experiment $run already analysed. skipping . . .\n";
		$i++;
		next;
	}
	
	else {
		$genotype = &determine_genotype ($library_name);
		
        # Get the data as a FASTQ. Always have --split-files on. In the case the
        # data is SE, it will only create a _1.fastq file.
        print  "\nFetching run: $run. Genotype: $genotype . . .\n";
        $cmd = "fastq-dump -I --split-files --outdir $out/$run $run";
        print "EXECUTING: $cmd\n";
        system($cmd);
        
        # Run sailfish.
        print  "Running sailfish  . . .\n";
        if (-e "$out/$run/$run\_2.fastq") {
            # Data is paired end.
            $cmd = "sailfish quant -i $MAGIC_INDEX_DIR/$genotype/ -l T=PE:O=><:S=U -1 $out/$run/$run\_1.fastq -2 $out/$run/$run\_2.fastq -o $out/$run -p $NTHREADS";
	    	print "EXECUTING: $cmd\n";
            system($cmd);
	    
	    	system("rm $out/$run/$run\_1.fastq");
	    	system("rm $out/$run/$run\_2.fastq");
        }
        else {
            # Data is single end.
            $cmd = "sailfish quant -i $MAGIC_INDEX_DIR/$genotype/ -l T=SE:S=U -r $out/$run/$run\_1.fastq -o $out/$run -p $NTHREADS";
	    	print "EXECUTING: $cmd\n";
            system($cmd);
	    
	    system("rm $out/$run/$run\_1.fastq");
        }
        
		$check_success = 0;
		
		if (-e "$out/$run/logs/") {
			
			opendir LOG, "$out/$run/logs/" or die $!;
			my @findlogfile = grep { $_ ne '.' && $_ ne '..' } readdir LOG;
			closedir LOG;
						
			open my $logfile, "<$out/$run/logs/@findlogfile" or print "cannot find logfile for run $run!\n";
			
			while (<$logfile>){
				$check_success = 1 if /g2log\ file\ shutdown/i;
			}
		}
        
		if ($check_success) {
			print  "Run $run finished successfully! ", ($i / $queries) * 100, " % done.\n";
		}
		else {
			print "odd looking logfile for run $run. check sailfish logfile at $out/$run/logs.\n";
			print "proceeding anyway. ", ($i / $queries) * 100, " % done.\n";
		}
		$i++;
	}
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
	
	for my $known_genotype (keys(%known_genotypes)) {
		
#		print $known_genotype, "\n";
		if ($name =~ /$known_genotype[^a-z]/i){
			$final_genotype = $known_genotypes{$known_genotype};
			last;
		}
	}

	return $final_genotype;
}
