#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Cwd;
use File::Path;



# Assumes each sequence only takes 1 line in the fastq file.
my $freq = shift; # Proportion to downsample to. A number between 0 and 1. (required)
my $r1 = shift; # Read file, or file for read 1 if paired end.   (required)
my $r2 = shift; # File for read 2 if paired end.
my $seed = shift; # Random number seed.

# Set the random number seed.
srand($seed) if $seed;

my $isPE = $r2;


if ($isPE) {
    
    open R1, $r1;
    open R2, $r2;
    open R1O, ">$r1.sub";
    open R2O, ">$r2.sub";
    
    my (@r1l, @r2l, $i);
    my $itr = 1;
    my $u;
    while (<R1>){
        $r1l[0] = $_;
        $r2l[0] = <R2>;
        
        $r1l[1] = <R1>;
        $r2l[1] = <R2>;
        
        $r1l[2] = <R1>;
        $r2l[2] = <R2>;
        
        $r1l[3] = <R1>;
        $r2l[3] = <R2>;
        
        $u = rand();
        if ($u < $freq) {
            for ($i = 0; $i <= 3; $i++){
                print R1O $r1l[$i];
                print R2O $r2l[$i];
            }
        }
        
        $itr++;
    }
    
    close(R1);
    close(R2);
    close(R1O);
    close(R2O);
    
} else {
    # Single end
    
    open R1, $r1;
    open R1O, ">$r1.sub";
    
    my (@r1l, $i);
    my $itr = 1;
    my $u;
    while (<R1>){
        $r1l[0] = $_;
        
        $r1l[1] = <R1>;
        
        $r1l[2] = <R1>;
        
        $r1l[3] = <R1>;
        
        $u = rand();
        if ($u < $freq) {
            for ($i = 0; $i <= 3; $i++){
                print R1O $r1l[$i];
            }
        }
        
        $itr++;
    }
    
    close(R1);
    close(R1O);    
}




