#!/usr/bin/perl

# TRIMD_junction_validator.pl
# Copyright (C) 2016 Flemington Lab

#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Accepts a junctions files from GMAP/Iso-Seq (generated with the -f introns argument), an SJ.out.tab files from STAR/Illumina and an annotation file. Returns 3 bed files: one of SMRT introns, one of Illumina introns and one of introns detected by both methods. Annotation status of validated introns is noted.

#Accepts a junctions files from GMAP/Iso-Seq (generated with the -f introns argument), an SJ.out.tab files from STAR/Illumina and an annotation file. Returns 3 bed files: one of SMRT introns, one of Illumina introns and one of introns detected by both methods. Annotation status of validated introns is noted.

#USAGE:
# perl <PATH/TRIMD_junction_validator.pl> </PATH/Iso-Seq_introns_file> </PATH/Illumina_SJ.out.tab_file> </PATH/transcript_annotation_bed_file> <coordinates_to_ignore_bed_file(optional)>

use warnings;
use strict;

my ($SMRT_jfile, $ill_jfile, $ann_file, $ig_file) = @ARGV;

#print "Enter name of viral chromosome (e.g. chrEBV_Akata_inverted): ";
#my $viral_chr = <STDIN>;
#chomp $viral_chr;

my $min_SMRTj;
my $min_illj;

print "Use default parameters [y/n]? ";
my $answer = <STDIN>;
chomp $answer;

if ($answer eq "y") {
    $min_SMRTj = 1;
    $min_illj = 1;
}
else {
    print "Enter minimum Iso-Seq read depth to report a splice junction (e.g. 1): ";
    $min_SMRTj = <STDIN>;
    chomp $min_SMRTj;
    
    print "Enter minimum Illumina read depth to report a splice junction (e.g. 1): ";
    $min_illj = <STDIN>;
    chomp $min_illj;
}

print "================================================\n";

#####----------SMRT FILE CHROMOSOME EXTRACTION-------------######

open(INF, "<$SMRT_jfile") or die "couldn't open file";

print "Identifying chromosomes...\n";

my %chroms;

while (my $line = <INF>) {
    chomp($line);
    my ($chrom) = $line =~ /\s(.+):/;
    next if (exists $chroms{$chrom});
    $chroms{$chrom} = 1;
}

close(INF);

print "================================================\n";

my @total_juncs_found;
my @novel_juncs_found;
my @ann_juncs_found;
my @ann_juncs;

foreach my $chrom (sort keys %chroms) {
    print "Processing $chrom:\n";
    
    #####----------GMAP/SMRT FILE CONVERSION-------------######
    open(INF, "<$SMRT_jfile");
    open(OUT, ">$SMRT_jfile.temp");
    
    print "Processing Iso-Seq splice junctions...\n";
    
    while(my $line = <INF> ) {
        chomp($line);
        my ($id) = $line =~ /\>(.+)\.i/;
        my ($chr) = $line =~ /\s(.+):/;
        my ($score) = $line =~ /\>.+\/(\d+)\//;
        my ($donor, $acceptor) = $line =~ /:(\d+)\.\.(\d+)/;
        next if $chr ne $chrom;
        if ($acceptor > $donor) {
            print OUT $chr, "\t", $donor, "\t", $acceptor - 1, "\t", $id, "\t", $score, "\t+\n";
        }
        else {
            print OUT $chr, "\t", $acceptor, "\t", $donor - 1, "\t", $id, "\t", $score, "\t-\n";
        }
    }
    close(OUT);
    close(INF);
    
    system("sort -k2,3n \Q$SMRT_jfile\E.temp > \Q$SMRT_jfile\E.sorted.temp"); #sorts so that duplicate introns will be next to each other.
    
    open(INF, "<$SMRT_jfile.sorted.temp" ) or die "couldn't reopen file";
    open(OUT, ">$SMRT_jfile.bed.temp");
    
    my $plus_previous_chr = "start";
    my $plus_count = 0;
    my $plus_previous_start = 0;
    my $plus_previous_end = 0;
    my $minus_previous_chr = "start";
    my $minus_count = 0;
    my $minus_previous_start = 0;
    my $minus_previous_end = 0;
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        if ($cols[5] eq "+") { #plus and minus need to be treated separately in case of introns with the same starts and ends annotated on opposite strands
            if (($cols[0] eq $plus_previous_chr) and ($cols[1] == $plus_previous_start) and ($cols[2] == $plus_previous_end)) { #checks to see if the intron matches the previous intron
                $plus_count = $plus_count + $cols[4];
            }
            else {
                if ($plus_previous_chr eq "start") { #prevents the initial placeholder value from printing out as a line, and sets the values of the first intron
                    $plus_previous_chr = $cols[0];
                    $plus_previous_start = $cols[1];
                    $plus_previous_end = $cols[2];
                    $plus_count = $cols[4];
                }
                else {
                    print OUT "$plus_previous_chr\t$plus_previous_start\t$plus_previous_end\t$plus_count\t$plus_count\t+\n";
                    $plus_previous_chr = $cols[0];
                    $plus_previous_start = $cols[1];
                    $plus_previous_end = $cols[2];
                    $plus_count = $cols[4];
                }
            }
        }
        if ($cols[5] eq "-") {
            if (($cols[0] eq $minus_previous_chr) and ($cols[1] == $minus_previous_start) and ($cols[2] == $minus_previous_end)) {
                $minus_count = $minus_count + $cols[4];
            }
            else {
                if ($minus_previous_chr eq "start") {
                    $minus_previous_chr = $cols[0];
                    $minus_previous_start = $cols[1];
                    $minus_previous_end = $cols[2];
                    $minus_count = $cols[4];
                }
                else {
                    print OUT "$minus_previous_chr\t$minus_previous_start\t$minus_previous_end\t$minus_count\t$minus_count\t-\n"; #prints out in bed format
                    $minus_count = $cols[4];
                    $minus_previous_chr = $cols[0];
                    $minus_previous_start = $cols[1];
                    $minus_previous_end = $cols[2];
                }
            }
        }
    }
    
    print OUT "$plus_previous_chr\t$plus_previous_start\t$plus_previous_end\t$plus_count\t$plus_count\t+\n"; #adds the last plus strand feature
    print OUT "$minus_previous_chr\t$minus_previous_start\t$minus_previous_end\t$minus_count\t$minus_count\t-\n"; #adds the last plus strand feature
    close(OUT);
    close(INF);
    
    system("sort -k2,3n \Q$SMRT_jfile\E.bed.temp > \Q$SMRT_jfile\E.\Q$chrom\E.bed.sorted.temp");
    system("rm \Q$SMRT_jfile\E.temp");
    system("rm \Q$SMRT_jfile.sorted\E.temp");
    system("rm \Q$SMRT_jfile\E.bed.temp");
    
    #if an annotation file of regions to be ignored is supplied, remove the SMRT junctions with a donor or acceptor in those regions:
    if (defined $ig_file) {
        open(INF, "<$ig_file");
        print "Removing Iso-Seq junctions with donor or acceptor in ignored region...\n";
        my @ig_coords;
        while(my $line = <INF>) {
            chomp($line);
            my @cols = split("\t", $line);
            next if $cols[0] ne $chrom;
            my $ig_coord = "$cols[1]:$cols[2]";
            push (@ig_coords, $ig_coord);
        }
        close(INF);
        
        open(INF, "<$SMRT_jfile.$chrom.bed.sorted.temp") or die "couldn't open file";
        open(OUT, ">$SMRT_jfile.$chrom.bed.noheader");
        
        while(my $line = <INF>) {
            chomp($line);
            my @cols = split( "\t", $line );
            my $found_flag=0;
            foreach my $ig_coord (@ig_coords) {
                my ($ig_start, $ig_end) = split (":", $ig_coord);
                if ((($cols[1] >= $ig_start) and ($cols[1] <= $ig_end)) || (($cols[2] >= $ig_start) and ($cols[2] <= $ig_end))) {
                    $found_flag = 1;
                    last;
                }
            }
            if ($found_flag == 0) {
                print OUT $line, "\n";
            }
        }
        
        close(INF);
        close(OUT);
    }
    
    #add header to bed file
    if (defined $ig_file) {
        open(INF, "<$SMRT_jfile.$chrom.bed.noheader") or die "couldn't open file";
    }
    else {
        open(INF, "<$SMRT_jfile.$chrom.bed.sorted.temp") or die "couldn't open file";
    }
    open(OUT, ">$SMRT_jfile.$chrom.bed") or die "couldn't open file";
    
    print OUT "track type=bed name=\"$SMRT_jfile.$chrom.bed\" description=\"Iso-Seq introns from splice_junction_matcher.pl\"\n";
    while (my $line = <INF>) {
        print OUT $line;
    }
    close(OUT);
    close(INF);
    
    if (defined $ig_file) {
        system("rm \Q$SMRT_jfile\E.\Q$chrom\E.bed.noheader");
    }
    system("rm \Q$SMRT_jfile\E.\Q$chrom\E.bed.sorted.temp");
    
    
    #####----------STAR/ILLUMINA FILE CONVERSION-------------######
    
    open(INF, "<$ill_jfile" ) or die "couldn't open file";
    open(OUT, ">$ill_jfile.$chrom.bed");
    
    print "Processing Illumina splice junctions...\n";
    print OUT "track type=bed name=\"$ill_jfile.$chrom.bed\" description=\"Illumina STAR introns from splice_junction_matcher.pl\"\n";
    
    while(my $line = <INF> ) {
        chomp($line);
        my @cols = split( "\t", $line );
        tr/12/+-/ foreach ($cols[3]);	#change the numeric strand indicators to + or -
        next if $cols[0] ne $chrom; #skip lines that aren't viral
        my $chrStart = $cols[1] - 1; #changes the start coordinate to 0-based for bed
        print OUT "$cols[0]\t$chrStart\t$cols[2]\t$cols[4]\t$cols[6]\t$cols[3]\n";
    }
    
    close(OUT);
    close(INF);
    
    #if an annotation file of regions to be ignored is supplied, remove the Illumina junctions with a donor or acceptor in those regions:
    if (defined $ig_file) {
        open(INF, "<$ig_file");
        print "Removing Illumina junctions with donor or acceptor in ignored region...\n";
        my @ig_coords;
        while(my $line = <INF>) {
            chomp($line);
            my @cols = split("\t", $line);
            my $ig_coord = "$cols[1]:$cols[2]";
            push (@ig_coords, $ig_coord);
        }
        close(INF);
        
        open(INF, "<$ill_jfile.$chrom.bed") or die "couldn't open file";
        open(OUT, ">$ill_jfile.$chrom.no_ignored.bed");
        
        print OUT "track type=bed name=\"$ill_jfile.$chrom.no_ignored.bed\" description=\"Illumina STAR introns from splice_junction_matcher.pl\"\n";
        
        while(my $line = <INF>) {
            chomp($line);
            next if ($line =~ /^track/); #skips the track definition line
            my @cols = split( "\t", $line );
            my $found_flag=0;
            foreach my $ig_coord (@ig_coords) {
                my ($ig_start, $ig_end) = split (":", $ig_coord);
                if ((($cols[1] >= $ig_start) and ($cols[1] <= $ig_end)) || (($cols[2] >= $ig_start) and ($cols[2] <= $ig_end))) {
                    $found_flag = 1;
                    last;
                }
            }
            if ($found_flag == 0) {
                print OUT $line, "\n";
            }
        }
        
        close(INF);
        close(OUT);
        system("rm \Q$ill_jfile\E.\Q$chrom\E.bed");
    }
    
    #####----------GMAP/ILLUMINA COMPARISON-------------######
    
    if (defined $ig_file) {
        open(INF, "<$ill_jfile.$chrom.no_ignored.bed" ) or die "couldn't open file";
    }
    else {
        open(INF, "<$ill_jfile.$chrom.bed" ) or die "couldn't open file";
    }
    
    print "Checking for matching splice junctions...\n";
    
    my %ill_junctions;
    
    while(my $line = <INF> ) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        next if ($cols[4] < $min_illj);
        my $ill_key_combo = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in the Illumina file, creates a key for the hash combining chromosome, start coordinate, end coordinate and strand
        $ill_junctions{$ill_key_combo} = $cols[4]; #enters a count value for the key into the hash
    }
    
    close(INF);
    
    open(INF, "<$SMRT_jfile.$chrom.bed" ) or die "couldn't open file";
    open(OUT, ">$SMRT_jfile.$chrom.illumina_support.bed.temp");
    
    while(my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        next if ($cols[4] < $min_SMRTj);
        my $SMRT_key_combo = "$cols[0]$cols[1]$cols[2]$cols[5]"; #for each line in the SMRT file, creates a variable/key combining chromosome, start coordinate, end coordinate and strand
        if (exists $ill_junctions{$SMRT_key_combo}) { #checks to see if the key exists in the Illumina hash: if so, prints it out
            my $junction_depth = $cols[4] + $ill_junctions{$SMRT_key_combo};
            print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$cols[4].IsoSeq_$ill_junctions{$SMRT_key_combo}.Ill\t$junction_depth\t$cols[5]\n";
        }
        else {
            print OUT "$cols[0]\t$cols[1]\t$cols[2]\t$cols[3].IsoSeq\t$cols[4]\t$cols[5]\n";
        }
    }
    close(INF);
    close(OUT);
    
    #####----------ANNOTATION FILE COMPARISON-------------######
    
    #First extract intron coordinates from the annotation file
    open(INF, "<$ann_file");
    
    print "Processing annotation file...\n";
    
    my @intron_start;
    my @intron_end;
    my %ann_intron_coord_pair;
    my $start;
    my $end;
    
    while (my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        next if $cols[0] ne $chrom; #skip lines that aren't viral
        my $intron_number = $cols[9] - 1;
        next if ($intron_number == 0);
        my @block_sizes = split(",", $cols[10]);
        my @block_starts = split(",", $cols[11]);
        for (my $i = 0; $i < $intron_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of intron start sites relative to the genome
            $start = $cols[1] + $block_sizes[$i] + $block_starts[$i];
            push(@intron_start, $start);
        }
        for (my $i2 = 1; $i2 < $cols[9]; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
            $end = $cols[1] + $block_starts[$i2];
            push(@intron_end, $end);
        }
        for (my $i3 = 0; $i3 < $intron_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create a hash of complete intron coordinates relative to the genome
            my $intron_coords = "$cols[0]:$intron_start[$i3]:$intron_end[$i3]:$cols[5]";
            if (exists $ann_intron_coord_pair{$intron_coords}) {
                $ann_intron_coord_pair{$intron_coords} = $ann_intron_coord_pair{$intron_coords} + 1; #if the intron is already in the hash (from another transcript), increase the count
            }
            else {
                $ann_intron_coord_pair{$intron_coords} = 1; #if the intron is not already in the hash, adds it with a value of 1
            }
        }
        @intron_start = ();
        @intron_end = (); #intron starts and ends have been assigned to the %ann_intron_pair hash; empty them for the next transcript
    }
    
    my $ann_count = 0;
    
    if (defined $ig_file) {
        open(INF, "<$ig_file");
        my @ig_coords;
        while(my $line = <INF>) {
            chomp($line);
            my @cols = split("\t", $line);
            my $ig_coord = "$cols[1]:$cols[2]";
            push (@ig_coords, $ig_coord);
        }
        close(INF);
        
        foreach my $ann_intron_coord_pair (keys %ann_intron_coord_pair) {
            my ($ann_chr, $ann_start, $ann_end, $ann_strand) = split (":", $ann_intron_coord_pair);
            my $found_flag = 0;
            foreach my $ig_coord (@ig_coords) {
                my ($ig_start, $ig_end) = split (":", $ig_coord);
                if ((($ann_start >= $ig_start) and ($ann_start <= $ig_end)) || (($ann_end >= $ig_start) and ($ann_end <= $ig_end))) {
                    $found_flag = 1;
                    last;
                }
            }
            if ($found_flag == 0) {
                $ann_count++;
            }
        }
    }
    else {
        $ann_count = scalar (keys %ann_intron_coord_pair);
    }
    
    close(INF);
    
    #Compare introns in the altered (with Illumina data) SMRT file to annotated introns
    
    open(INF, "<$SMRT_jfile.$chrom.illumina_support.bed.temp");
    open(OUT, ">$SMRT_jfile.$chrom.validated_introns.bed");
    
    print "Comparing Iso-Seq junctions to annotation file...\n";
    
    print OUT "track type=bed name=\"$SMRT_jfile.$chrom.validated_introns.bed\" description=\"Introns detected by Iso-Seq with read depth at least $min_SMRTj supported by Illumina-detected junctions with read depth at least $min_illj and/or annotation. From splice_junction_matcher.pl\"\n";
    
    my $val_SMRT_count = 0;
    my $ann_SMRT_count = 0;
    my $nov_SMRT_count = 0;
    
    while (my $line = <INF>) {
        chomp($line);
        my @SMRT_cols = split("\t", $line);
        my $SMRT_intron_coords = "$SMRT_cols[0]:$SMRT_cols[1]:$SMRT_cols[2]:$SMRT_cols[5]"; #creates a key to search the has of annotated introns
        if (exists $ann_intron_coord_pair{$SMRT_intron_coords}) { #if the intron matches an annotated intron, notes that and prints out the line
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tann_$SMRT_cols[5]_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\n";
            $ann_SMRT_count++;
            $val_SMRT_count++;
        }
        else {
            if ($SMRT_cols[3] =~ /.+IsoSeq_.+Ill/) { #if the intron doesn't match an annotated intron but does have Illumina support, notes that and prints out the line
                print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tnov_$SMRT_cols[5]_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\n";
                $nov_SMRT_count++;
                $val_SMRT_count++;
            }
        }
    }
    close(OUT);
    close(INF);
    
    print "------------------------------------------------\n";
    
    push(@total_juncs_found, $val_SMRT_count);
    push(@novel_juncs_found, $nov_SMRT_count);
    push(@ann_juncs_found, $ann_SMRT_count);
    push(@ann_juncs, $ann_count);
    
    #open(OUT, ">${chrom}_validated_introns_stats.txt");
    
    if ($val_SMRT_count > 0) {
        print "$val_SMRT_count validated junctions detected in the Iso-Seq file. $nov_SMRT_count are novel and $ann_SMRT_count are annotated (out of $ann_count annotated junctions).\n";
        if (defined $ig_file) {
            #print OUT "$chrom\n\n$val_SMRT_count validated junctions\n\t$nov_SMRT_count novel\n\t$ann_SMRT_count annotated\n$ann_count junctions in annotation file\n\ninput files:\n\t$SMRT_jfile\n\t$ill_jfile\n\t$ann_file\n\t$ig_file\n";
        }
        else {
            #print OUT "$chrom\n\n$val_SMRT_count validated junctions\n\t$nov_SMRT_count novel\n\t$ann_SMRT_count annotated\n$ann_count junctions in annotation file\n\ninput files:\n\t$SMRT_jfile\n\t$ill_jfile\n\t$ann_file\n";
        }
    }
    else {
        print "No validated junctions found.\n";
        if (defined $ig_file) {
            #print OUT "No validated junctions found.\n\ninput files:\n\t$SMRT_jfile\n\t$ill_jfile\n\t$ann_file\n\t$ig_file\n";
        }
        else{
            #print OUT "No validated junctions found.\n\ninput files:\n\t$SMRT_jfile\n\t$ill_jfile\n\t$ann_file \n";
        }
    }
    
    #close(OUT);
    
    system ("rm \Q$SMRT_jfile\E.\Q$chrom\E.illumina_support.bed.temp");
    
    print "================================================\n";
}

system("cat $SMRT_jfile.*.bed > ${SMRT_jfile}_introns.bed");

if (defined $ig_file) {
    system("cat $ill_jfile*.no_ignored.bed > ${ill_jfile}_introns.bed");
}
else {
    system("cat $ill_jfile.*.bed > ${ill_jfile}_introns.bed");
}


system("cat *.validated_introns.bed > ${SMRT_jfile}_validated_introns.bed");

system("rm $SMRT_jfile.*.bed");
if (defined $ig_file) {
    system("rm $ill_jfile*.no_ignored.bed");
}
system("rm $ill_jfile.*.bed");

my $sum_total_found = 0;
my $sum_novel_found = 0;
my $sum_ann_found = 0;
my $sum_total_ann = 0;


foreach my $total_juncs (@total_juncs_found) {
    $sum_total_found = $sum_total_found + $total_juncs;
}
foreach my $novel_juncs (@novel_juncs_found) {
    $sum_novel_found = $sum_novel_found + $novel_juncs;
}
foreach my $ann_found (@ann_juncs_found) {
    $sum_ann_found = $sum_ann_found + $ann_found;
}
foreach my $total_ann (@ann_juncs) {
    $sum_total_ann = $sum_total_ann + $total_ann;
}

open(OUT, ">validated_introns.txt");

if (defined $ig_file) {
    print OUT "$sum_total_found validated junctions\n\t$sum_novel_found novel\n\t$sum_ann_found annotated\n$sum_total_ann junctions in annotation file\n\ninput files:\n\t$SMRT_jfile\n\t$ill_jfile\n\t$ann_file\n\t$ig_file\n";
}
else {
    print OUT "$sum_total_found validated junctions\n\t$sum_novel_found novel\n\t$sum_ann_found annotated\n$sum_total_ann junctions in annotation file\n\ninput files:\n\t$SMRT_jfile\n\t$ill_jfile\n\t$ann_file\n";
}
close(OUT);