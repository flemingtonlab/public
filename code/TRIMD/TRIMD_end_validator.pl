#!/usr/bin/perl

# TRIMD_end_validator.pl
# Copyright (C) 2016 Flemington Lab

#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Accepts a SAM file using Iso-Seq fl data, a SAM file using Illumina data, and a bed file of annotated polyadenylated transcripts. Counts the number of non-clipped Iso-Seq reads with 3' ends at each genomic position and estimates consensus locations of clusters of 3' ends. Extracts Illumina reads containing apparent polyA tails and estimates consensus locations of clusters of polyadenylation sites. Output includes bedgraph files of all 3' ends, bed files of the weighted centers of end clusters, a sam file of reads with polyA tails and a bed file of Iso-Seq 3' ends supported by annotation, Illumina data or a genomic polyadenylation sequence.

#USAGE:
# perl <PATH/TRIMD_end_validator.pl> </PATH/Iso-Seq_sam_file> </PATH/Illumina_sam_file> </PATH/Annotation_bed_file>

use warnings;
use strict;

die "USAGE: 'perl <PATH/TRIMD_end_validator.pl> </PATH/Iso-Seq_sam_file> </PATH/Illumina_sam_file> </PATH/Annotation_bed_file>'" unless @ARGV == 3;

my ($SMRT_file, $ill_file, $ann_file) = @ARGV;

#print "Enter name of viral chromosome (e.g. chrEBV_Akata_inverted): ";
#my $viral_chr = <STDIN>;
#chomp $viral_chr;

my $distance_between_SMRT_peaks;
my $min_As;
my $min_softclip;
my $distance_between_ill_peaks;
my $dist_SMRT_ill_d;
my $dist_SMRT_ill_u;
my $min_SMRT;
my $min_ill;
my $ann_dist;

print "Use default parameters [y/n]? ";
my $answer = <STDIN>;
chomp $answer;

if ($answer eq "y") {
    $distance_between_SMRT_peaks = 8;
    $min_As = 5;
    $min_softclip = 2;
    $distance_between_ill_peaks = 8;
    $dist_SMRT_ill_d = 10;
    $dist_SMRT_ill_u = 4;
    $min_SMRT = 5;
    $min_ill = 1;
    $ann_dist = 10;
    my $gsig_file;
}
else {
    print "Enter desired window for collapsing Iso-Seq 3' ends (e.g. 8): ";
    $distance_between_SMRT_peaks = <STDIN>;
    chomp $distance_between_SMRT_peaks;
    
    print "Enter minimum number of As for Illumina poly(A) tails (e.g. 5): ";
    $min_As = <STDIN>;
    chomp $min_As;
    
    print "Enter minimum number of mismatches for Illumina poly(A) tails (e.g. 2): ";
    $min_softclip = <STDIN>;
    chomp $min_softclip;
    
    print "Enter desired window for collapsing Illumina 3' ends (e.g. 8): ";
    $distance_between_ill_peaks = <STDIN>;
    chomp $distance_between_ill_peaks;
    
    print "Enter number of bases downstream of Iso-Seq ends to look for Illumina support (e.g. 10): ";
    $dist_SMRT_ill_d = <STDIN>;
    chomp $dist_SMRT_ill_d;
    
    print "Enter number of bases upstream of Iso-Seq ends to look for Illumina support (e.g. 4): ";
    $dist_SMRT_ill_u = <STDIN>;
    chomp $dist_SMRT_ill_u;
    
    print "Enter minimum number of Iso_seq reads to report a 3' end (e.g. 5): ";
    $min_SMRT = <STDIN>;
    chomp $min_SMRT;
    
    print "Enter minimum number of Illumina poly(A) tails to support a 3' end (e.g. 1): ";
    $min_ill = <STDIN>;
    chomp $min_ill;
    
    print "Enter maximum distance in bp from an annotated end to be called as 'annotated' (e.g. 10): ";
    $ann_dist = <STDIN>;
    chomp $ann_dist;
    
    print "Consider genomic polyadenylation signals? [n or /PATH/TO/BEDFILE]:\n";
    $gsig_file = <STDIN>;
    chomp $gsig_file;
}

print "================================================\n";

#####----------SMRT FILE CHROMOSOME EXTRACTION-------------######

open(INF, "<$SMRT_file") or die "couldn't open file";

print "Identifying chromosomes...\n";

my %chroms;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if ($cols[0] eq "\@HD" || $cols[0] eq "\@PG" || $cols[0] eq "\@SQ"); #skips SAM file header lines
    next if ($cols[2] eq "*"); #skips unmapped CFLs
    next if (exists $chroms{$cols[2]});
    $chroms{$cols[2]} = 1;
}

close(INF);

#####----------ILLUMINA POLY(A) TAIL DETECTION-------------######

open(INF, "<$ill_file") or die "couldn't open input file";
open(OUT, ">$ill_file.polyA_ends.temp") or die "couldn't open output file";

print "Extracting Illumina reads with at least $min_As As and at least $min_softclip mismatches...\n";

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    next if ($cols[0] eq "\@HD" || $cols[0] eq "\@PG" || $cols[0] eq "\@SQ"); #skips SAM file header lines
    if ($cols[1] == 81 || $cols[1] == 83 || $cols[1] == 89 || $cols[1] == 16) {  #selects reads with FLAG codes indicating they are first in pair on the plus strand
        if (($cols[5] =~ m/\d+S$/) and ($cols[9] =~ m/A{$min_As}$/)) { # selects reads with softclipping and a run of As at the end
            my ($softclips) = $cols[5] =~ m/(\d+)S$/; #pulls out the number of softclipped bases
            if ($softclips > $min_softclip) { #selects reads with at least the specified number of softclipped bases
                print OUT $line, "\n";
            }
        }
    }
    elsif ($cols[1] == 73 || $cols[1] == 97 || $cols[1] == 99 || $cols[1] == 0) {  #selects reads with FLAG codes indicating they are first in pair on the minus strand
        if (($cols[5] =~ m/^\d+S/) and ($cols[9] =~ m/^T{$min_As}/)) { #selects reads with softclipping and a run of Ts at the beginning
            my ($softclips) = $cols[5] =~ m/^(\d+)S/; #pulls out the number of softclipped bases
            if ($softclips > $min_softclip) { #selects reads with at least the specified number of softclipped bases
                print OUT $line, "\n";
            }
        }
    }
}

close(INF);
close(OUT);

system ("sort -k 3,3 -k 4,4n \Q$ill_file\E.polyA_ends.temp > \Q$ill_file\E.polyA_ends.sam");
system ("rm \Q$ill_file\E.polyA_ends.temp");

open(INF, "<$ill_file.polyA_ends.sam") or die "couldn't open file";
open(OUT, ">$ill_file.polyA_sites.temp") or die "couldn't open file";

print "Processing Illumina reads with polyA tails...\n";
#create a file with the coordinates corresponding to the polyA ends of the reads, and sort it by those coordinates

my $cigar_sum;
my $cigar_calc;
my @plus_ends;
my @read_dist;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    if ($cols[1] == 73 || $cols[1] == 97 || $cols[1] == 99 || $cols[1] == 0) { #minus strand
        print OUT "$cols[2]\t$cols[3]\t0\n";
    }
    elsif ($cols[1] == 81 || $cols[1] == 83 || $cols[1] == 89 || $cols[1] == 16) { #plus strand
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
            push (@read_dist, $1);
        }
        $cigar_sum += $_ for @read_dist;
        $cigar_calc = $cols[3] + $cigar_sum - 1; #subtract one to account for start/end inclusion
        $cigar_sum = 0;
        @read_dist = ();
        print OUT "$cols[2]\t$cigar_calc\t1\n";
    }
}
close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$ill_file.polyA_sites.temp\E > \Q$ill_file.polyA_sites.temp\E.sorted");

#create a bedgraph file from the sorted coordinates file

open(INF, "<$ill_file.polyA_sites.temp.sorted") or die "couldn't open file";
open(OUT, ">$ill_file.polyA_sites.temp.bedgraph") or die "couldn't open file";

my $chrom_minus;
my $previous_coordinate_m=0;
my $count_m=0;
my $chrom_plus;
my $previous_coordinate_p=0;
my $count_p=0;

while (my $line = <INF>) {
    
    my @cols = split("\t", $line);
    
    #reads on the plus strand:
    if ($cols[2] == 1) {
        if ($chrom_plus) { #if $chrom_plus has been defined (i.e. there is a previous plus strand read)
            if (($cols[0] eq $chrom_plus) and ($cols[1] == $previous_coordinate_p)) {
                $count_p++;
            }
            else {
                print OUT $chrom_plus, "\t", $previous_coordinate_p-1, "\t", $previous_coordinate_p, "\t", $count_p, "\n"; #prints to output file, converting chrStart to 0-based bedgraph coordinates
                $chrom_plus = $cols[0];
                $previous_coordinate_p = $cols[1];
                $count_p = 1;
            }
        }
        else { #if $chrom_plus has not been defined (i.e. there is no previous plus strand read)
            $chrom_plus = $cols[0];
            $previous_coordinate_p = $cols[1];
            $count_p = 1;
        }
    }
    
    #reads on the minus strand:
    elsif ($cols[2] == 0) {
        if ($chrom_minus) {
            if (($cols[0] eq $chrom_minus) and ($cols[1] == $previous_coordinate_m)) {
                $count_m++;
            }
            else {
                print OUT $chrom_minus, "\t", $previous_coordinate_m-1, "\t", $previous_coordinate_m, "\t-", $count_m, "\n"; #prints to output file, converting chrStart to 0-based bedgraph coordinates
                $chrom_minus = $cols[0];
                $previous_coordinate_m = $cols[1];
                $count_m = 1;
            }
        }
        else {
            $chrom_minus = $cols[0];
            $previous_coordinate_m = $cols[1];
            $count_m = 1;
        }
    }
}
#prints to output file, converting chrStart to 0-based bedgraph coordinates
print OUT $chrom_plus, "\t", $previous_coordinate_p-1, "\t", $previous_coordinate_p, "\t", $count_p, "\n";
print OUT $chrom_minus, "\t", $previous_coordinate_m-1, "\t", $previous_coordinate_m, "\t-", $count_m, "\n";

close(INF);
close(OUT);

system("sort -k 1,1 -k 2,2n \Q$ill_file\E.polyA_sites.temp.bedgraph > \Q$ill_file\E.polyA_sites.bedgraph.noheader");
system("rm \Q$ill_file\E.polyA_sites.temp.bedgraph");
system("rm \Q$ill_file\E.polyA_sites.temp.sorted");
system("rm \Q$ill_file\E.polyA_sites.temp");

#add header to bedgraph file
open(INF, "<$ill_file.polyA_sites.bedgraph.noheader") or die "couldn't open file";
open(OUT, ">$ill_file.polyA_sites.bedgraph") or die "couldn't open file";

print OUT "track type=bedGraph name=\"$ill_file.polyA_sites.bedgraph\" description=\"polyA sites in Illumina reads with at least $min_As As and at least $min_softclip mismatches from end_finder_sam_to_bed.pl\"\n";
while (my $line = <INF>) {
    print OUT $line;
}
close(OUT);
close(INF);

system("rm \Q$ill_file\E.polyA_sites.bedgraph.noheader");

#####----------CHROMOSOME-BY-CHROMOSOME ANALYSIS-------------######

my @total_ends_found;
my @novel_ends_found;
my @SMRT_ends_ann;
my @ann_ends_found;
my @ann_ends;

foreach my $chrom (sort keys %chroms) {
    print "Processing $chrom:\n";
    
    #####----------SMRT FILE PROCESSING-------------######
    system("awk '\$3==\"\Q$chrom\E\"' \Q$SMRT_file\E \| sort -k 4,4n > \Q$SMRT_file\E.sorted.temp");
    system("awk '\$2==0' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.plus.sam.temp");
    system("awk '\$2==16' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.minus.sam.temp");
    system("rm \Q$SMRT_file\E.sorted.temp");
    
    #processing of PLUS sam file
    open(INF, "<$SMRT_file.sorted.plus.sam.temp") or die "couldn't open file";
    open(OUT, ">$SMRT_file.sorted.plus.sam.read_ends.bedgraph.temp") or die "couldn't open file";
    
    my @dist;
    my $sum;
    my %plus_ends;
    print "Processing Iso-Seq plus strand reads...\n";
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        next if ($cols[5] =~ m/\d+S$/); #skips reads soft-clipped at the 3' end
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
            push (@dist, $1);
        }
        $sum += $_ for @dist;
        my $end_coord = $cols[3] + $sum - 1; #subtract 1 to account for start/end inclusion
        my $chr_end_coord = "$cols[2]\:$end_coord"; #combines the chromosome and 3' end coordinate into a key to use for the hash
        $sum = 0;
        @dist = ();
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
        if (exists $plus_ends{$chr_end_coord}) { #if the key is already in the hash, increases the value (count) by 1
            $plus_ends{$chr_end_coord} = $plus_ends{$chr_end_coord} + $split_id[1];
        }
        else {
            $plus_ends{$chr_end_coord} = $split_id[1]; #if the key is not already in the hash, adds it with a value (count) of the read depth
        }
    }
    
    foreach my $chr_end_coord (sort keys %plus_ends) { #prints out a(n inadequately) sorted temporary bedgraph file
        my @split_keys = split("\:", $chr_end_coord);
        print OUT $split_keys[0], "\t", $split_keys[1]-1, "\t", $split_keys[1], "\t", $plus_ends{$chr_end_coord}, "\n"; #prints to output file, converting chrStart to 0-based bedgraph coordinates
    }
    close(INF);
    close(OUT);
    
    system("rm \Q$SMRT_file\E.sorted.plus.sam.temp");
    
    #processing of MINUS sam file
    open(INF, "<$SMRT_file.sorted.minus.sam.temp") or die "couldn't open file";
    open(OUT, ">$SMRT_file.sorted.minus.sam.read_ends.bedgraph.temp") or die "couldn't open file";
    
    my $previous_coordinate=1;
    my $count=0;
    my $previous_chr = "start";
    print "Processing Iso-Seq minus strand reads...\n";
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        next if ($cols[5] =~ m/^\d+S/); #skips reads soft-clipped at the 3' end
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
        if (($cols[2] eq $previous_chr) and ($cols[3] == $previous_coordinate)) {
            $count = $count + $split_id[1]; #increases the count by the read depth for the putative isoform
        }
        else {
            if ($previous_chr eq "start") { #doesn't print out the placeholder first line.
                $previous_chr = $cols[2];	#sets the previous chromosome, previous coordinate and count values
                $previous_coordinate = $cols[3];
                $count = $split_id[1];
            }
            else {
                print OUT $previous_chr, "\t", $previous_coordinate-1, "\t", $previous_coordinate, "\t-", $count, "\n"; #prints to output file, converting chrStart to 0-based bedgraph coordinates
                $previous_chr = $cols[2];
                $previous_coordinate = $cols[3];
                $count = $split_id[1];
            }
        }
    }
    
    print OUT $previous_chr, "\t", $previous_coordinate-1, "\t", $previous_coordinate, "\t-", $count, "\n"; #prints the last start coordinates to output file
    close(INF);
    close(OUT);
    
    system("cat \Q$SMRT_file\E.sorted.plus.sam.read_ends.bedgraph.temp \Q$SMRT_file\E.sorted.minus.sam.read_ends.bedgraph.temp | sort -k2,3n > \Q$SMRT_file\E.\Q$chrom\E.read_ends.bedgraph");
    
    system("rm \Q$SMRT_file\E.sorted.plus.sam.read_ends.bedgraph.temp");
    system("rm \Q$SMRT_file\E.sorted.minus.sam.read_ends.bedgraph.temp");
    system("rm \Q$SMRT_file\E.sorted.minus.sam.temp");
    
    #add header to bedgraph file
    #    open(INF, "<$SMRT_file.$chrom.read_ends.bedgraph.noheader") or die "couldn't open file";
    #    open(OUT, ">$SMRT_file.$chrom.read_ends.bedgraph") or die "couldn't open file";
    #
    #    print OUT "track type=bedGraph name=\"$SMRT_file.$chrom.read_ends.bedgraph\" description=\"3' ends of Iso-Seq reads from end_finder_sam_to_bed.pl\"\n";
    #    while (my $line = <INF>) {
    #        print OUT $line;
    #    }
    #    close(OUT);
    #    close(INF);
    #
    #    system("rm \Q$SMRT_file\E.\Q$chrom\E.read_ends.bedgraph.noheader");
    
    #make a bed file from the SMRT bedgraph file:
    open(INF, "<$SMRT_file.$chrom.read_ends.bedgraph") or die "couldn't open file";
    open(OUT, ">$SMRT_file.ends.temp.bed") or die "couldn't open file";
    
    print "Combining Iso-Seq 3' ends within $distance_between_SMRT_peaks of each other and calculating consensus 3' ends...\n";
    collapse_bedgraph($chrom, $distance_between_SMRT_peaks);
    
    close(INF);
    close(OUT);
    
    system("sort -k 1,1 -k 2,2n \Q$SMRT_file\E.ends.temp.bed > \Q$SMRT_file\E.\Q$chrom\E.SMRT_ends.bed");
    system("rm \Q$SMRT_file.ends.temp.bed\E");
    
    #add header to bed file
    #    open(INF, "<$SMRT_file.ends.bed.noheader") or die "couldn't open file";
    #    open(OUT, ">$SMRT_file.$chrom.SMRT_ends.bed") or die "couldn't open file";
    #
    #    print OUT "track type=bed name=\"$SMRT_file.$chrom.SMRT_ends.bed\" description=\"consensus 3' ends of Iso-Seq reads within $distance_between_SMRT_peaks bp collapsed to weighted center from end_finder_sam_to_bed.pl\"\n";
    #    while (my $line = <INF>) {
    #        print OUT $line;
    #    }
    #    close(OUT);
    #    close(INF);
    #
    #    system("rm \Q$SMRT_file\E.ends.bed.noheader");
    
    #####----------ILLUMINA FILE PROCESSING-------------######
    
    #make a bed file from the Illumina bedgraph file:
    system("awk '\$1==\"\Q$chrom\E\"' \Q$ill_file\E.polyA_sites.bedgraph > \Q$ill_file\E.\Q$chrom\E.polyA_sites.bedgraph.temp");
    open(INF, "<$ill_file.$chrom.polyA_sites.bedgraph.temp") or die "couldn't open file";
    open(OUT, ">$ill_file.$chrom.polyA_sites.temp.bed") or die "couldn't open file";
    
    print "Combining Illumina polyA tails within $distance_between_ill_peaks of each other and calculating consensus 3' ends...\n";
    collapse_bedgraph($chrom, $distance_between_ill_peaks);
    
    close(INF);
    close(OUT);
    
    system("sort -k 1,1 -k 2,2n \Q$ill_file\E.\Q$chrom\E.polyA_sites.temp.bed > \Q$ill_file\E.\Q$chrom\E.polyA_sites.bed");
    system("rm \Q$ill_file\E.\Q$chrom\E.polyA_sites.temp.bed");
    system("rm \Q$ill_file\E.\Q$chrom\E.polyA_sites.bedgraph.temp");
    
    #add header to bed file
    #    open(INF, "<$ill_file.$chrom.polyA_sites.bed.noheader") or die "couldn't open file";
    #    open(OUT, ">$ill_file.$chrom.polyA_sites.bed") or die "couldn't open file";
    #
    #    print OUT "track type=bed name=\"$ill_file.$chrom.polyA_sites.bed\" description=\"consensus polyA sites of Illumina reads with tails of 5 As with 2 mismatches within $distance_between_ill_peaks bp collapsed to weighted centers from end_finder_sam_to_bed.pl\"\n";
    #    while (my $line = <INF>) {
    #        print OUT $line;
    #    }
    #    close(OUT);
    #    close(INF);
    #
    #    system("rm \Q$ill_file\E.\Q$chrom\E.polyA_sites.bed.noheader");
    
    #####----------SEEKING ILLUMINA SUPPORT FOR SMRT ENDS-------------######
    
    open(INF, "<$ill_file.$chrom.polyA_sites.bed" ) or die "couldn't open file";
    
    print "Extracting Iso-Seq 3' ends with Illumina polyA tails within $dist_SMRT_ill_d bases downstream or $dist_SMRT_ill_u upstream...\n";
    
    my %features_ill;
    my $key_combo_ill;
    
    while(my $line = <INF> ) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        if ($cols[5] eq "+") { #for each line in the Illumina polyA reads bed file, creates a key for the hash combining coordinate and strand. Selects chrEnd for ends on the plus strand and chrStart for ends on the minus strand.
            $key_combo_ill = "$cols[2]:$cols[5]";
        }
        if ($cols[5] eq "-") {
            $key_combo_ill = "$cols[1]:$cols[5]";
        }
        $features_ill{$key_combo_ill} = $cols[4]; #enters a count value for the key into the hash
    }
    
    close(INF);
    
    open(INF, "<$SMRT_file.$chrom.SMRT_ends.bed" ) or die "couldn't open file";
    open(OUT, ">$SMRT_file.$chrom.ends.bed.illumina_support.bed.temp");
    
    my $ill_coord;
    my $match_count;
    my $lower_limit;
    my $upper_limit;
    
    while(my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @SMRT_cols = split("\t", $line);
        next if (abs $SMRT_cols[4] < $min_SMRT);
        foreach my $key_combo_ill (keys %features_ill) {
            my @ill_cols = split(":", $key_combo_ill);
            next if (abs $features_ill{$key_combo_ill} < $min_ill);
            
            if ($SMRT_cols[5] eq "+") { #sets boundaries for plus strand support
                $lower_limit = $SMRT_cols[2]-$dist_SMRT_ill_u;
                $upper_limit = $SMRT_cols[2]+$dist_SMRT_ill_d;
            }
            if ($SMRT_cols[5] eq "-") { #sets boundaries for minus strand support
                $lower_limit = $SMRT_cols[1]-$dist_SMRT_ill_d;
                $upper_limit = $SMRT_cols[1]+$dist_SMRT_ill_u;
            }
            if (($SMRT_cols[5] eq $ill_cols[1]) and ($ill_cols[0] >= $lower_limit) and ($ill_cols[0] <= $upper_limit)) {
                
                if ($match_count) { #if more than one Illumina end matches the SMRT end, selects Illumina end with the most reads
                    if ($features_ill{$key_combo_ill} > $match_count){
                        $match_count = $features_ill{$key_combo_ill};
                        $ill_coord = $ill_cols[0];
                    }
                }
                else {
                    $match_count = $features_ill{$key_combo_ill};
                    $ill_coord = $ill_cols[0];
                }
            }
            
        }
        if ($match_count) {
            if ($SMRT_cols[5] eq "+") {
                my $name = "$SMRT_cols[4].IsoSeq_$match_count.Ill";
                my $count = $match_count + $SMRT_cols[4];
                print OUT $SMRT_cols[0], "\t", $ill_coord-1, "\t", $ill_coord, "\t", $name, "\t", $count, "\t", $SMRT_cols[5], "\t", $SMRT_cols[3], "\n"; #prints to output, adjusting chrStart to 0-based
                undef($match_count);
            }
            if ($SMRT_cols[5] eq "-") {
                my $name = "$SMRT_cols[4].IsoSeq_$match_count.Ill";
                my $count = $match_count + $SMRT_cols[4];
                print OUT $SMRT_cols[0], "\t", $ill_coord, "\t", $ill_coord+1, "\t", $name, "\t", $count, "\t", $SMRT_cols[5], "\t", $SMRT_cols[3], "\n"; #prints to output, adujsting chrEnd
                undef($match_count);
            }
        }
        else {
            my @range_cols = split (":", $SMRT_cols[3]); #includes SMRT ends that are not supported by Illumina in this temporary file
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\t$range_cols[2].IsoSeq\t$range_cols[2]\t$SMRT_cols[5]\t$SMRT_cols[3]\n";
        }
    }
    
    close(OUT);
    close(INF);
    
    #####----------COMPARING TO GENOMIC POLYADENYLATION SIGNALS-------------######
    
    if ($gsig_file) {
        
        #get locations of polyadenylation signals:
        open(INF, "<$gsig_file" ) or die "couldn't open file";
        
        print "Extracting Iso-Seq 3' ends within 40 bp downstream of a genomic polyadenylation signal...\n";
        
        my %features_gsig;
        my $key_combo_gsig;
        
        while(my $line = <INF> ) {
            chomp($line);
            next if ($line =~ /^track/); #skips the track definition line
            my @cols = split("\t", $line);
            if ($cols[0] eq $chrom) {
                if ($cols[5] eq "+") { #for each line in the polyadenylation signal bed file, creates a key for the hash combining coordinate and strand. Selects chrEnd for ends on the plus strand and chrStart for ends on the minus strand.
                    $key_combo_gsig = "$cols[2]:$cols[5]";
                }
                if ($cols[5] eq "-") {
                    $key_combo_gsig = "$cols[1]:$cols[5]";
                }
                $features_gsig{$key_combo_gsig} = $cols[4]; #enters a count value for the key into the hash
            }
        }
        
        close(INF);
        
        #then compare the SMRT ends in the temporary file (which already includes information about Illumina support) to the poly(A) signal locations in the hash.
        open(INF, "<$SMRT_file.$chrom.ends.bed.illumina_support.bed.temp" ) or die "couldn't open file";
        open(OUT, ">$SMRT_file.$chrom.ends.bed.illumina_gpolyasig_support.bed.temp");
        
        my $gsig_coord;
        my $gsig_found = 0;
        
        while(my $line = <INF>) {
            chomp($line);
            next if ($line =~ /^track/); #skips any track definition line
            my @SMRT_cols = split("\t", $line);
            
            if ($SMRT_cols[5] eq "+") { #sets boundaries for plus strand support
                $lower_limit = $SMRT_cols[2] - 40;
                $upper_limit = $SMRT_cols[2];
            }
            if ($SMRT_cols[5] eq "-") { #sets boundaries for minus strand support
                $lower_limit = $SMRT_cols[1];
                $upper_limit = $SMRT_cols[1] + 40;
            }
            
            foreach my $key_combo_gsig (keys %features_gsig) {
                my @gsig_cols = split(":", $key_combo_gsig);
                
                if (($SMRT_cols[5] eq $gsig_cols[1]) and ($gsig_cols[0] >= $lower_limit) and ($gsig_cols[0] <= $upper_limit)) {
                    
                    $gsig_found = 1;
                    
                }
                
            }
            
            if ($gsig_found == 1) {
                print OUT $SMRT_cols[0], "\t", $SMRT_cols[1], "\t", $SMRT_cols[2], "\t", $SMRT_cols[3], "_gsig", "\t", $SMRT_cols[4], "\t", $SMRT_cols[5], "\t", $SMRT_cols[6], "\n";
            }
            
            else { #includes SMRT ends that are not supported in this temporary file
                print OUT "$line\n";
            }
            $gsig_found = 0;
            
        }
        
        close(OUT);
        close(INF);
    }
    
    #####----------COMPARING TO ANNOTATED ENDS-------------######
    open(INF, "<$ann_file" ) or die "couldn't open file";
    
    print "Processing annotation file...\n";
    
    #extract 3' ends from the annotation file:
    #annotation file must be sorted by chrStart then chrEnd!
    my @annotated_ends;
    my $plus_prev_coord = 0;
    my $minus_prev_coord = 0;
    
    while(my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @ann_cols = split("\t", $line);
        next if $ann_cols[0] ne $chrom; #skip lines that aren't viral
        if ($ann_cols[5] eq "+") {
            if ($ann_cols[2] != $plus_prev_coord) {
                push (@annotated_ends, "$ann_cols[2]:$ann_cols[5]"); #creates an array with chrEnd and strand
                $plus_prev_coord = $ann_cols[2];
            }
        }
        elsif ($ann_cols[5] eq "-"){
            if ($ann_cols[1] != $minus_prev_coord) {
                push (@annotated_ends, "$ann_cols[1]:$ann_cols[5]"); #creates an array with chrStart and strand
                $minus_prev_coord = $ann_cols[1];
            }
        }
    }
    
    my $annotated = scalar @annotated_ends;
    
    close(INF);
    
    #compare ends in the altered SMRT ends file (that already has info about Illumina ends and genomic polyadenylation signals) with annotated ends
    
    open(OUT, ">$SMRT_file.$chrom.validated_ends.bed");
    
    if ($gsig_file) {
        
        open(INF, "<$SMRT_file.$chrom.ends.bed.illumina_gpolyasig_support.bed.temp" ) or die "couldn't open file";
        
        print "Comparing Iso-Seq ends to annotated ends...\n";
        
        print OUT "track type=bedDetail name=\"$SMRT_file.$chrom.validated_ends.bed\" description=\"validated ends supported by  at least $min_SMRT Iso-Seq read ends within $distance_between_SMRT_peaks bp, with an Illumina polyA site within $dist_SMRT_ill_d bp downstream or $dist_SMRT_ill_u bp upstream, or within 40 bp downstream of a polyadenylation signal, or within $ann_dist bp of an annotated end. Illumina polyA sites have at least $min_ill reads with $min_As As and $min_softclip mismatches, within $distance_between_ill_peaks bp of each other. From end_finder_sam_to_bed.pl\"\n";
    }
    else {
        open(INF, "<$SMRT_file.$chrom.ends.bed.illumina_support.bed.temp" ) or die "couldn't open file";
    
        print "Comparing Iso-Seq ends to annotated ends...\n";
    
        #print OUT "track type=bedDetail name=\"$SMRT_file.$chrom.SMRT_ends.bed.illumina_support.bed\" description=\"validated ends supported by  at least $min_SMRT Iso-Seq read ends within $distance_between_SMRT_peaks bp, with an Illumina polyA site within $dist_SMRT_ill_d bp downstream or $dist_SMRT_ill_u bp upstream, or within $ann_dist bp of an annotated end. Illumina polyA sites have at least $min_ill reads with $min_As As and $min_softclip mismatches, within $distance_between_ill_peaks bp of each other. From end_finder_sam_to_bed.pl\"\n";
    }
    
    my $annotated_found_by_SMRT = 0;
    my $novel_found_by_SMRT_ill = 0;
    my $SMRT_annotated = 0; #this is different than $annotated_found_by_SMRT because depending on input parameters two SMRT ends may correspond to a single annotated end or vice versa.
    
    while(my $line = <INF>) {
        chomp($line);
        my @SMRT_cols = split("\t", $line);
        my $found_flag=0;
        foreach my $ann_end (@annotated_ends) {
            my @ann_cols = split(":", $ann_end);
            my $lower_limit = $ann_cols[0]-$ann_dist;
            my $upper_limit = $ann_cols[0]+$ann_dist;
            if ($ann_cols[1] eq "+") {
                if (($SMRT_cols[5] eq $ann_cols[1]) and ($SMRT_cols[2]>=$lower_limit) and ($SMRT_cols[2]<=$upper_limit)) {
                    if ($found_flag == 0) {
                        print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tann_$SMRT_cols[5]_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\t$SMRT_cols[6]\n";
                        $found_flag = 1;
                        $annotated_found_by_SMRT++;
                        $SMRT_annotated++;
                    }
                    elsif ($found_flag == 1) {
                        $annotated_found_by_SMRT++;
                    }
                }
            }
            if ($ann_cols[1] eq "-") {
                if (($SMRT_cols[5] eq $ann_cols[1]) and ($SMRT_cols[1]>=$lower_limit) and ($SMRT_cols[1]<=$upper_limit)) {
                    if ($found_flag == 0) {
                        print OUT $SMRT_cols[0], "\t", $SMRT_cols[1], "\t", $SMRT_cols[2], "\tann_", $SMRT_cols[5], "_", $SMRT_cols[3], "\t", abs($SMRT_cols[4]), "\t", $SMRT_cols[5], "\t", $SMRT_cols[6], "\n";
                        $found_flag = 1;
                        $annotated_found_by_SMRT++;
                        $SMRT_annotated++;
                    }
                    elsif ($found_flag == 1) {
                        $annotated_found_by_SMRT++;
                    }
                }
            }
        }
        if ($found_flag == 0) {
            if (($SMRT_cols[3] =~ /.+IsoSeq_.+Ill/)||($SMRT_cols[3] =~ /.+_gsig/)) {
                print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tnov_$SMRT_cols[5]_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\t$SMRT_cols[6]\n";
                $novel_found_by_SMRT_ill++;
            }
        }
    }
    
    my $total_found = $SMRT_annotated + $novel_found_by_SMRT_ill;
    
    close(INF);
    close(OUT);
    
    print "------------------------------------------------\n";
    
    push (@total_ends_found, $total_found);
    push (@novel_ends_found, $novel_found_by_SMRT_ill);
    push (@SMRT_ends_ann, $SMRT_annotated);
    push (@ann_ends_found, $annotated_found_by_SMRT);
    push (@ann_ends, $annotated);
    
    #open(OUT, ">${chrom}_validated_ends_stats.txt");
    
    if ($total_found > 0) {
        if ($SMRT_annotated != $annotated_found_by_SMRT) {
            print "$total_found 3' ends found. $novel_found_by_SMRT_ill are novel, $SMRT_annotated are annotated.  $annotated_found_by_SMRT out of $annotated total annotated 3' ends are found.\nNote that two annotated ends may be within $ann_dist bp of a single Iso-Seq end or vice versa.\n";
            #print OUT "$chrom\n\n$total_found 3' ends\n\t$novel_found_by_SMRT_ill novel\n\t$SMRT_annotated annotated\n$annotated 3' ends in annotation\n\t$annotated_found_by_SMRT detected by Iso-Seq\n\ninput files:\n\t$SMRT_file\n\t$ill_file\n\t$ann_file\n";
        }
        else {
            print "$total_found 3' ends found. $novel_found_by_SMRT_ill are novel, $SMRT_annotated are annotated (out of a total of $annotated annotated 3' ends).\n\n";
            #print OUT "$chrom\n\n$total_found 3' ends\n\t$novel_found_by_SMRT_ill novel\n\t$SMRT_annotated annotated\n$annotated 3' ends in annotation\n\ninput files:\n\t$SMRT_file\n\t$ill_file\n\t$ann_file\n";
        }
    }
    else {
        print "No 3' ends validated\n";
        #print OUT "No validated 3' ends found\n\ninput files:\n\t$SMRT_file\n\t$ill_file\n\t$ann_file\n";
    }
    
    #close(OUT);
    
    system("rm \Q$SMRT_file\E.\Q$chrom\E.ends.bed.illumina_support.bed.temp");
    
    print "================================================\n";
}

system("cat \Q$ill_file\E.*.polyA_sites.bed > \Q$ill_file\E.polyA_sites.bed");
#print "consolidating Illumina polyA site bed\n";
#system("cat \Q$ill_file\E.*.polyA_sites.bedgraph > \Q$ill_file\E.polyA_sites.bedgraph");
#print "consolidating Illumina polyA site bedgraph\n";
system("cat \Q$SMRT_file\E.*.SMRT_ends.bed > \Q$SMRT_file\E.SMRT_ends.bed");
#print "consolidating Iso-Seq ends bed\n";
system("cat \Q$SMRT_file\E.*.read_ends.bedgraph > \Q$SMRT_file\E.read_ends.bedgraph");
#print "consolidating Iso-Seq ends bedgraph\n";
system("cat \Q$SMRT_file\E.*.validated_ends.bed > \Q$SMRT_file\E.validated_ends.bed");
#print "consolidating validated ends\n";

system("rm \Q$ill_file\E.*.polyA_sites.bed");
#print "removing Illumina polyA site bed files\n";
#system("rm \Q$ill_file\E.*.polyA_sites.bedgraph");
#print "removing Illumina polyA sites bedgraphs\n";
system("rm \Q$SMRT_file\E.*.SMRT_ends.bed");
#print "removing Iso-Seq ends bed file\n";
system("rm \Q$SMRT_file\E.*.read_ends.bedgraph");
#print "removing Iso_seq ends bedgraph files\n";
system("rm \Q$SMRT_file\E.*.validated_ends.bed");
#print "removing validated ends files\n";

my $sum_total_found = 0;
my $sum_novel_found = 0;
my $sum_SMRT_ann = 0;
my $sum_total_ann = 0;
my $sum_ann_found = 0;

foreach my $total_ends (@total_ends_found) {
    $sum_total_found = $sum_total_found + $total_ends;
}
foreach my $novel_ends (@novel_ends_found) {
    $sum_novel_found = $sum_novel_found + $novel_ends;
}
foreach my $SMRT_ann (@SMRT_ends_ann) {
    $sum_SMRT_ann = $sum_SMRT_ann + $SMRT_ann;
}
foreach my $total_ann (@ann_ends) {
    $sum_total_ann = $sum_total_ann + $total_ann;
}
foreach my $ann_SMRT (@ann_ends_found) {
    $sum_ann_found = $sum_ann_found + $ann_SMRT;
}

open(OUT, ">validated_ends.txt");

print OUT "$sum_total_found 3' ends\n\t$sum_novel_found novel\n\t$sum_SMRT_ann annotated\n$sum_total_ann ends in annotation file\n\t$sum_ann_found detected by Iso-Seq\n\ninput files:\n\t$SMRT_file\n\t$ill_file\n\t$ann_file\n";

close(OUT);

#print %chroms, "\n";

#print @ann_ends, "\n";

#########################
sub collapse_bedgraph {
    my ($chrom, $distance_between_peaks) = @_;
    my $prev_coord_plus = 0;
    my $prev_coord_minus = 0;
    my $count_sum_plus = 0;
    my $count_sum_minus = 0;
    my $weighted_coordinate_sum_plus = 0;
    my $weighted_coordinate_sum_minus = 0;
    my $weighted_average_plus;
    my $weighted_average_minus;
    my $first_plus = 1;
    my $first_minus = 1;
    my @coords_plus;
    my @coords_minus;
    my $chrStart_plus;
    my $chrEnd_plus;
    my $chrStart_minus;
    my $chrEnd_minus;
    
    while (my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        if ($cols[3] > 0) { #if this coordinate has a positive count...
            if ($cols[2] <= $prev_coord_plus + ($distance_between_peaks)) { #if the coordinate is within the specified number of bp of the previous coordinate
                $count_sum_plus = $count_sum_plus + $cols[3]; #adds to the sums to eventually calculate the weighted average
                $weighted_coordinate_sum_plus = $weighted_coordinate_sum_plus + ($cols[2]*$cols[3]);
                push (@coords_plus, $cols[2]);
                $prev_coord_plus = $cols[2]; #sets the current coordinate as the "previous coordinate" before moving on
            }
            else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
                if ($first_plus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                    $count_sum_plus = $cols[3];
                    $weighted_coordinate_sum_plus = $cols[2]*$cols[3];
                    $prev_coord_plus = $cols[2];
                    push (@coords_plus, $cols[2]);
                    $first_plus = 0;
                }
                else {
                    $weighted_average_plus = sprintf("%1.0f", ($weighted_coordinate_sum_plus/$count_sum_plus)); #calculates weighted average
                    $chrStart_plus = $coords_plus[0];
                    $chrEnd_plus = pop(@coords_plus);
                    print OUT $chrom, "\t", $weighted_average_plus-1, "\t", $weighted_average_plus,  "\t", $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, "\t", $count_sum_plus, "\t+\n"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
                    @coords_plus = ($cols[2]);
                    $count_sum_plus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
                    $weighted_coordinate_sum_plus = $cols[2]*$cols[3];
                    $prev_coord_plus = $cols[2];
                }
            }
        }
        elsif ($cols[3] < 0) { #if this coordinate has a negative count...
            if ($cols[1] <= $prev_coord_minus + ($distance_between_peaks)) { #if the coordinate is within the specified number of bp of the previous coordinate
                $count_sum_minus = $count_sum_minus + $cols[3]; #adds to the sums to eventually calculate the weighted average
                $weighted_coordinate_sum_minus = $weighted_coordinate_sum_minus + ($cols[1]*$cols[3]);
                push (@coords_minus, $cols[1]);
                $prev_coord_minus = $cols[1]; #sets the current coordinate as the "previous coordinate" before moving on
            }
            else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
                if ($first_minus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                    $count_sum_minus = $cols[3];
                    $weighted_coordinate_sum_minus = $cols[1]*$cols[3];
                    $prev_coord_minus = $cols[1];
                    push (@coords_minus, $cols[1]);
                    $first_minus = 0;
                }
                else {
                    $weighted_average_minus = sprintf("%1.0f", ($weighted_coordinate_sum_minus/$count_sum_minus)); #calculates weighted average
                    $chrStart_minus = $coords_minus[0];
                    $chrEnd_minus = pop(@coords_minus);
                    print OUT $chrom, "\t", $weighted_average_minus, "\t", $weighted_average_minus+1, "\t", $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, "\t", abs($count_sum_minus), "\t-\n";
                    @coords_minus = ($cols[1]);
                    $count_sum_minus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
                    $weighted_coordinate_sum_minus = $cols[1]*$cols[3];
                    $prev_coord_minus = $cols[1];
                }
            }
        }
    }
    
    if ($count_sum_plus > 0) {#calculates and prints out weighted average for the last feature (plus strand)
        $weighted_average_plus = sprintf("%1.0f", ($weighted_coordinate_sum_plus/$count_sum_plus));
        $chrStart_plus = $coords_plus[0];
        $chrEnd_plus = pop(@coords_plus);
        print OUT $chrom, "\t", $weighted_average_plus-1, "\t", $weighted_average_plus, "\t", $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, "\t", $count_sum_plus, "\t+\n";
    }
    
    if ($count_sum_minus < 0) {#calculates and prints out weighted average for the last feature (minus strand)
        $weighted_average_minus = sprintf("%1.0f", ($weighted_coordinate_sum_minus/$count_sum_minus));
        $chrStart_minus = $coords_minus[0];
        $chrEnd_minus = pop(@coords_minus);
        print OUT $chrom, "\t", $weighted_average_minus, "\t", $weighted_average_minus+1, "\t", $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, "\t", abs($count_sum_minus), "\t-\n";
    }
}