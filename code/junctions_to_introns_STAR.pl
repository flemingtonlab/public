#!/usr/bin/perl
#junctions_to_introns_STAR.pl by TO'G. Converts SJ.out.tab files from STAR to intron bed files.
#Usage: $ perl <PATH/junctions_to_introns_STAR.pl> </PATH/SJ.out.tab>

use warnings;
use strict;

die "USAGE: 'perl <PATH/junctions_to_introns_STAR.pl> </PATH/SJ.out.tab>'" unless @ARGV == 1;

foreach my $file (@ARGV) {
    open(INF, "<$file" ) or die "couldn't open file";
    open(OUT, ">$file.bed");

    my $line;
    while(<INF>) {
        chomp;
        my @cols = split("\t", $_);

        tr/12/+-/ foreach ($cols[3]);	#change the numeric strand indicators to + or -
        
        my $corr_chrStart = $cols[1] - 1; #change the coordinate to zero-based for bed format
            
        print OUT "$cols[0]\t$corr_chrStart\t$cols[2]\t$cols[4]:$cols[5]:$cols[6]:$cols[7]:$cols[8]\t$cols[6]\t$cols[3]\n"; #print out the line in bed format, with the extra STAR columns as the colon-delimited name

    }
    close(OUT);
    close(INF);
}