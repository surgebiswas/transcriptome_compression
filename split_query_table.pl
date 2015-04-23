#!/usr/bin/env perl

# split_query_table.pl
# by Surge Biswas, Konstantin Kerner

use strict;
use warnings;
use File::Basename;

my $query_table		= shift;
my $subsets			= shift;
my $out				= shift;
my $queries_per_subset;
my $queries			= 0;
my $int_queries;
my $header;
my $query_count;
my $subset_count	= 1;
my $query_table_suffix = basename($query_table);
$query_table_suffix =~ s/\.csv$//;


open my $fh, "<$query_table" or die $!;

$header = <$fh>;

$queries++ while <$fh>;
close  $query_table;

$queries_per_subset = $queries / $subsets;
$int_queries = int($queries_per_subset) +1;

print "$queries queries identified in $query_table. splitting into $subsets subsets of $int_queries queries (last subset might contain less)\n";


open my $filehandle, "<$query_table" or die $!;

<$filehandle>;


while ($subset_count <= $subsets) {
	
	open my $out, ">$out/$query_table_suffix$subset_count.csv" or die $!;
	print $out $header;
	
	$query_count = 1;
	
	while (<$filehandle>) {
		
		print $out $_;
		last if $query_count >= $queries_per_subset;
		$query_count++;
	}
	close $out;
	$subset_count++;
}
close $filehandle;