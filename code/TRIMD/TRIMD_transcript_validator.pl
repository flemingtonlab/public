#!usr/bin/perl

# TRIMD_end_validator.pl
# Copyright (C) 2016 Flemington Lab

#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Takes a sam file of Iso-Seq fl isoforms and compares them to a list of validated 5' ends, 3' ends and introns to create a list of validated transcript structures, which are compared to an annotation file.

#USAGE:
# perl <PATH/TRIMD_transcript_validator.pl> </PATH/Iso-Seq_sam_file> </PATH/validated_starts_file> </PATH/validated_ends_file> </PATH/validated_introns_file> </PATH/Annotation_bed_file>

use warnings;
use strict;

die "USAGE: 'perl <PATH/TRIMD_transcript_validator.pl> </PATH/Iso-Seq_sam_file> </PATH/validated_starts_file> </PATH/validated_ends_file> </PATH/validated_introns_file> </PATH/Annotation_bed_file>'" unless @ARGV == 5;

my ($test_file, $valid_starts_file, $valid_ends_file, $valid_introns_file, $ann_file) = (@ARGV);

print "Enter maximum distance from an annotated 5' start to be called annotated (e.g. 10): ";
my $start_dist = <STDIN>;
chomp $start_dist;

print "Enter maximum distance from an annotated 3' end to be called annotated (e.g. 10): ";
my $end_dist = <STDIN>;
chomp $end_dist;

#Convert SMRT sam file to bed:

print "------------------------------------------------\nReformatting Iso-Seq file...\n";

open(INF, "<$test_file") or die "couldn't open input file";
open(OUT, ">$test_file.bed") or die "couldn't open output file";

while (my $line = <INF>) {
    $line =~ s/\r//g;
    chomp($line);
    next if ($line =~ m/\@/); #skips SAM header lines
    my @cols = split("\t", $line);
    my @split_id = split("\/", $cols[0]);
    my $strand;
    my $chr = $cols[2];
    my $chr_start = $cols[3] - 1;
    my $chr_end = 0;
    my $feature_name = $cols[0];
    my $score = $split_id[1];
    my $color = "133,0,33";
    if ($cols[1] == 0) {
        $strand = "+";
    }
    elsif ($cols[1] == 16) {
        $strand = "-";
    }
    else {
        next; #skips isoforms that aren't mapped
    }
    my @split_CIGAR_temp = split(/(\d+\D)/, $cols[5]); #splits CIGAR code into segments and puts segments into an array (but also the empty values between them)
    my @split_CIGAR;
    foreach my $temporary(@split_CIGAR_temp) { #removes empty values from the array
        if ($temporary =~ m/\d+\D/) {
            push(@split_CIGAR, $temporary);
        }
    }
    my $exon_sum = 0;
    my @exon_lengths = ();
    my @block_starts = (0);
    my $count = 0;
    foreach my $split_CIGAR(@split_CIGAR) {
        $count++;
        if (($count == 1) && (my ($five_prime_clipped_bases) = $split_CIGAR =~ m/(\d+)S$/)) { #ignores soft clipping at the beginning
        }
        
        elsif (($count > 1) && (my ($three_prime_clipped_bases) = $split_CIGAR =~ m/(\d+)S$/)) { #ignores soft clipping at the end
        }
        elsif ($split_CIGAR =~ m/N$/) { #if element is an intron...
            push(@exon_lengths, $exon_sum); #...adds the last value to the exon sum...
            my ($intron_length) = $split_CIGAR =~ m/(\d+)/;#...gets intron length...
            my $new_block_start = $exon_lengths[-1] + $block_starts[-1] + $intron_length;#...calculates new blockStart...
            push(@block_starts, $new_block_start);#...adds new blockStart to array...
            $exon_sum = 0;#...and resets the exon sum
        }
        else {
            my ($value) = $split_CIGAR =~ m/(\d+)/;
            if ($split_CIGAR =~ m/I$/) {#ignores insertions
                $exon_sum = $exon_sum - 0;
            }
            else {#adds matches, mismatches and deletions to the exon sum
                $exon_sum = $exon_sum + $value;
            }
        }
    }
    push(@exon_lengths, $exon_sum); #at the end of the CIGAR array, push the last exon sum into the exon_lengths array
    $chr_end = $chr_start + $block_starts[-1] + $exon_lengths[-1];
    my $exon_number = @exon_lengths;
    print OUT $chr, "\t", $chr_start, "\t", $chr_end, "\t", $feature_name, "\t", $score, "\t", $strand, "\t", $chr_start, "\t", $chr_end, "\t", $color, "\t", $exon_number, "\t", join("\,", @exon_lengths), "\t", join("\,", @block_starts), "\n";
}
close(INF);
close(OUT);

#Create an array of validated start sites from the start sites input file:
open(INF, "<$valid_starts_file") or die "couldn't open file";

my @valid_start;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@valid_start, $line); #puts each line of the start sites file into an array to be checked later
}
close(INF);

#Check each start site in the SMRT reads file against the array of validated start sites:
open(INF, "<$test_file.bed") or die "couldn't open file";
#open(OUT, ">$test_file.valid_start.bed.temp"); #uncomment to print out file of isoforms with validated starts

print "Checking start sites...\n";

my @good_start;
my $new_start_line;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts)  = split("\t", $line);
	if ($strand eq "+") {
		foreach my $valid_start (@valid_start) { #checks to see if the 5' end of the (plus strand) SMRT transcript matches a range of possible start site values from the list of validated start sites
			my @start_cols = split("\t", $valid_start);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $start_cols[6]);
			if (($chrom eq $start_cols[0]) and ($strand eq $start_cols[5]) and ($chromStart >= $range_start) and ($chromStart <= $range_end)) {
				$new_start_line = "$line\t$start_cols[1]"; #creates a line for the read, changing the start site to the consensus start site and adding an extra field with the original start site
				push (@good_start, $new_start_line); #if the start site matches, pushes the line into a new array of SMRT transcripts with validated 5' ends
                #print OUT $new_start_line, "\n"; #uncomment to print out file of isoforms with validated starts
                last;
			}
		}
	}
	elsif ($strand eq "-"){
		foreach my $valid_start (@valid_start) { #checks to see if the 5' end of the (minus strand) SMRT transcript matches a range of possible start site values from the list of validated start sites
			my @start_cols = split("\t", $valid_start);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $start_cols[6]);
			if (($chrom eq $start_cols[0]) and ($strand eq $start_cols[5]) and ($chromEnd >= $range_start) and ($chromEnd <= $range_end)) {
				$new_start_line = "$line\t$start_cols[2]"; #creates a line for the read, changing the start site to the consensus start site and adding an extra field with the original start site
				push (@good_start, $new_start_line); #if the start site matches, pushes the line into a new array of SMRT transcripts with validated 5' ends
                #print OUT $new_start_line, "\n"; #uncomment to print out file of isoforms with validated starts
                last;
			}
		}
	}
}

my $good_start_number = scalar @good_start;

#close(OUT); #uncomment to print out file of isoforms with validated starts and ends
close(INF); #have an array in memory of reads that have validated 5' ends, and their newly estimated 5' ends. Can uncomment lines to have an output a file of reads (in their original form) that have validated 5' ends.

#Create an array of validated end sites from the end sites input file:
open(INF, "<$valid_ends_file") or die "couldn't open file";

my @valid_end;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@valid_end, $line); #puts each line of the end sites file into an array to be checked later
}

close(INF);

#open(OUT, ">$test_file.valid_start_and_end.bed.temp"); #uncomment to print out file of isoforms with validated starts and ends

print "Checking end sites...\n";

my @good_start_and_end;
my $new_end_line;

foreach my $good_start (@good_start) { #starts with the array of SMRT transcripts with validated 5' ends
	my ($chrom, $chromStart, $chromEnd, $name, $score, $strand, $thickStart, $thickEnd, $itemRgb, $blockCount, $blockSizes, $blockStarts, $new_coord) = split("\t", $good_start);
	if ($strand eq "+") { #determines 3' end of the SMRT transcript
		foreach my $valid_end (@valid_end) { #checks to see if the 3' end of the SMRT transcript matches a range of possible 3' end values from the list of validated start sites
            my @end_cols = split("\t", $valid_end);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $end_cols[6]);
            if (($chrom eq $end_cols[0]) and ($strand eq $end_cols[5]) and ($chromEnd >= $range_start) and ($chromEnd <= $range_end)) {
                $new_end_line = "$good_start\t$end_cols[2]";
                push (@good_start_and_end, $new_end_line);
                #print OUT $new_end_line, "\n"; #uncomment to print out file of isoforms with validated starts and ends
                last;
            }
        }
    }
	elsif ($strand eq "-") {
        foreach my $valid_end (@valid_end) { #checks to see if the 3' end of the SMRT transcript matches a range of possible 3' end values from the list of validated start sites
            my @end_cols = split("\t", $valid_end);
            my ($range_start, $range_end, $SMRT_depth) = split(":", $end_cols[6]);
            if (($chrom eq $end_cols[0]) and ($strand eq $end_cols[5]) and ($chromStart >= $range_start) and ($chromStart <= $range_end)) {
                $new_end_line = "$chrom\t$chromStart\t$chromEnd\t$name\t$score\t$strand\t$thickStart\t$thickEnd\t$itemRgb\t$blockCount\t$blockSizes\t$blockStarts\t$end_cols[1]\t$new_coord";
                push (@good_start_and_end, $new_end_line);
                #print OUT $new_end_line, "\n"; #uncomment to print out file of isoforms with validated starts and ends
                last;
            }
        }
    }
}

my $good_start_end_number = scalar @good_start_and_end;

#close(OUT); #uncomment to print out file of isoforms with validated starts and ends

open(INF, "<$valid_introns_file") or die "couldn't open file";

my @valid_intron;

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@valid_intron, $line); #creates an array of valid splice junctions
}

close(INF);

open(OUT, ">$test_file.validated_unrefined.bed.temp");

print "Checking splice junctions...\n";

my $start;
my $end;
my @intron_start;
my @intron_end;
my @intron_coord_pair;
my @good_intron_counter;

foreach my $good_start_and_end (@good_start_and_end) { #starts with the array of SMRT transcripts with validated 5' and 3' ends
	my @cols = split("\t", $good_start_and_end);
	my $intron_strand = $cols[5];
	my $intron_number = $cols[9] - 1;
	if ($intron_number == 0) { #if a SMRT transcript has validated 5' and 3' ends and no introns, it is fully validated
		print OUT $good_start_and_end, "\n";
	}
	else {
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
		for (my $i3 = 0; $i3 < $intron_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create an array of complete intron coordinates relative to the genome
			my $intron_coords = "$intron_start[$i3]:$intron_end[$i3]";
			push (@intron_coord_pair, $intron_coords);
		} 
		@intron_start = ();
		@intron_end = (); #intron starts and ends have been assigned to the @intron_coords array; empty them for the next transcript
		foreach my $intron_coord_pair (@intron_coord_pair) { #goes through each intron in the SMRT transcript
			my @coords = split(":", $intron_coord_pair); #allows extraction of the start and end coordinates from each intron in the SMRT transcript
			foreach my $valid_intron (@valid_intron) { #goes through each intron in the array of validated introns
				my @valid_coords = split("\t", $valid_intron); #allows extraction of the start and end coordinates from each validated intron
				next if $cols[0] ne $valid_coords[0]; #enforces chromosome matching
                next if $intron_strand ne $valid_coords[5]; #enforces strand matching
				if (($coords[0] == $valid_coords[1]) and ($coords[1] == $valid_coords[2])) {
					push(@good_intron_counter, $intron_coord_pair);	#puts introns that are validated for this transcript into an array (this really just functions as a counter)
				}			
			}
		}
		@intron_coord_pair = (); #once each intron in the SMRT transcript has been examined, empty the array for the next transcript
		if (@good_intron_counter == $intron_number) { #check to see if all of the introns in the transcript are validated
			print OUT $good_start_and_end, "\n";
		}
		@good_intron_counter = (); #after checking to see if all the introns in the transcript are validated, empties this array for the next transcript

	}
}

close(OUT);

#system("sort -k 2,2n -k 3,3n \Q$test_file\E.valid_start.bed.temp > \Q$test_file\E.valid_start.bed"); #uncomment to print out file of isoforms with validated starts
#system("rm \Q$test_file\E.valid_start.bed.temp"); #uncomment to print out file of isoforms with validated starts

#system("sort -k 2,2n -k 3,3n \Q$test_file\E.valid_start_and_end.bed.temp > \Q$test_file\E.valid_start_and_end.bed"); #uncomment to print out file of isoforms with validated starts and ends
#system("rm \Q$test_file\E.valid_start_and_end.bed.temp"); #uncomment to print out file of isoforms with validated starts and ends

system("sort -k2,2n -k3,3n \Q$test_file\E.validated_unrefined.bed.temp > \Q$test_file\E.validated_unrefined.bed");
system("rm \Q$test_file\E.validated_unrefined.bed.temp");

open(INF, "<$test_file.validated_unrefined.bed");
open(OUT, ">$test_file.validated_refined.temp");

my $new_block_size;
my @exon_start;
my @exon_end;
my @exon_coord_pair;
my $validated_count = 0;

while (my $line = <INF>) {
    chomp($line);
    my @cols = split("\t", $line);
    my $exon_number = $cols[9];
    if ($exon_number == 1) { #if a SMRT transcript has validated 5' and 3' ends and no introns, it is fully validated. Just adjust the start and end to the consensus sites and fix the BlockSize accordingly
        $new_block_size = $cols[13]-$cols[12];
		print OUT "$cols[0]\t$cols[12]\t$cols[13]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[12]\t$cols[13]\t$cols[8]\t$cols[9]\t$new_block_size\t$cols[11]\n";
        $validated_count++;
	}
    else { #need to adjust the start and end sites and also the blockStarts and blockSizes
        my @block_sizes = split(",", $cols[10]);
		my @block_starts = split(",", $cols[11]);
		for (my $i = 0; $i < $exon_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of exon start sites relative to the genome
			$start = $cols[1] + $block_starts[$i];
			push(@exon_start, $start);
		}
		for (my $i2 = 0; $i2 < $exon_number; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
			$end = $cols[1] + $block_starts[$i2] + $block_sizes[$i2];
			push(@exon_end, $end);
		}
        shift(@exon_start); #removes the first exon start
        unshift(@exon_start, $cols[12]); #replaces the first exon start with the adjust chrStart value
        pop(@exon_end); #removes the last exon end
        push(@exon_end, $cols[13]); #replaces the last exon end with the adjusted chrEnd value
		for (my $i3 = 0; $i3 < $exon_number; $i3 = $i3 + 1) { #for the transcript currently in the "while" loop, matches up intron start and end sites to create an array of complete intron coordinates relative to the genome
			my $exon_coords = "$exon_start[$i3]:$exon_end[$i3]";
			push (@exon_coord_pair, $exon_coords);
		} 
		@exon_start = ();
		@exon_end = (); #intron starts and ends have been assigned to the @intron_coords array; empty them for the next transcript
		my @new_starts;
		my @new_sizes;
        #print $cols[3], "\t", @exon_coord_pair, "\n";
		foreach my $exon_coord_pair (@exon_coord_pair) { #goes through each exon in the SMRT transcript
			my @coords = split(":", $exon_coord_pair);
			my $blockStart = $coords[0] - $cols[12];
			my $blockSize = $coords[1] - $coords[0];
			push (@new_starts, $blockStart);
			push (@new_sizes, $blockSize);
		}
        @exon_coord_pair = ();
        shift(@new_starts); #removes the first value of the new_starts array, so we can replace it with 0 (it won't be 0 already if chrStart has been updated)
		my $assembled_starts = join(",", 0, @new_starts);
		my $assembled_sizes = join(",", @new_sizes);
		print OUT "$cols[0]\t$cols[12]\t$cols[13]\t$cols[3]\t$cols[4]\t$cols[5]\t$cols[12]\t$cols[13]\t$cols[8]\t$cols[9]\t$assembled_sizes\t$assembled_starts\n";
        $validated_count++;
    }
}

#then need to add code to collapse the transcripts with identical structure into a single feature

close(INF);
close(OUT);

system("sort -k 2,2n -k 3,3n -k11,11 -k12,12 -k5,5n \Q$test_file\E.validated_refined.temp > \Q$test_file\E.validated_refined.bed");
system("rm \Q$test_file\E.validated_refined.temp");

#Collapsing matching transcripts into single isoforms

open(INF, "<$test_file.validated_refined.bed");
open(OUT, ">$test_file.isoforms.bed");

print "Collapsing matching transcripts into isoforms...\n";

my $prev_chr_plus = "start";
my $prev_chrStart_plus = 0;
my $prev_chrEnd_plus = 0;
my $prev_name_plus;
my $count_plus = 0;
my $prev_rgb_plus;
my $prev_blocks_plus;
my $prev_blockSizes_plus = "1,1";
my $prev_blockStarts_plus = "0,0";
my $prev_chr_minus = "start";
my $prev_chrStart_minus = 0;
my $prev_chrEnd_minus = 0;
my $prev_name_minus;
my $count_minus = 0;
my $prev_rgb_minus;
my $prev_blocks_minus;
my $prev_blockSizes_minus = "1,1";
my $prev_blockStarts_minus = "0,0";
my $iso_count = 0;

while(my $line = <INF>) {
	chomp($line);
	my @cols = split("\t", $line);
    if ($cols[5] eq "+") {
        if (($cols[0] eq $prev_chr_plus) and ($cols[1] == $prev_chrStart_plus) and ($cols[2] == $prev_chrEnd_plus) and ($cols[10] eq $prev_blockSizes_plus) and ($cols[11] eq $prev_blockStarts_plus)) {
            $count_plus = $count_plus + $cols[4];
            $prev_chr_plus = $cols[0];
            $prev_name_plus = $cols[3];
            $prev_rgb_plus = $cols[8];
            $prev_blocks_plus = $cols[9];
            
        }
        else {
            if ($count_plus == 0) {
                $prev_chrStart_plus = $cols[1];
                $prev_chrEnd_plus = $cols[2];
                $prev_blockSizes_plus = $cols[10];
                $prev_blockStarts_plus = $cols[11];
                $count_plus = $cols[4];
                $prev_chr_plus = $cols[0];
                $prev_name_plus = $cols[3];
                $prev_rgb_plus = $cols[8];
                $prev_blocks_plus = $cols[9];
            }
            else {
                print OUT "$prev_chr_plus\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_name_plus\t$count_plus\t\+\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_rgb_plus\t$prev_blocks_plus\t$prev_blockSizes_plus\t$prev_blockStarts_plus\n";
                $iso_count++;
                $prev_chrStart_plus = $cols[1];
                $prev_chrEnd_plus = $cols[2];
                $prev_blockSizes_plus = $cols[10];
                $prev_blockStarts_plus = $cols[11];
                $count_plus = $cols[4];
                $prev_chr_plus = $cols[0];
                $prev_name_plus = $cols[3];
                $prev_rgb_plus = $cols[8];
                $prev_blocks_plus = $cols[9];
            }
        }
    }
    elsif ($cols[5] eq "-") {
        if (($cols[0] eq $prev_chr_minus) and ($cols[1] == $prev_chrStart_minus) and ($cols[2] == $prev_chrEnd_minus) and ($cols[10] eq $prev_blockSizes_minus) and ($cols[11] eq $prev_blockStarts_minus)) {
            $count_minus = $count_minus + $cols[4];
            $prev_chr_minus = $cols[0];
            $prev_name_minus = $cols[3];
            $prev_rgb_minus = $cols[8];
            $prev_blocks_minus = $cols[9];
            
        }
        else {
            if ($count_minus == 0) {
                $prev_chrStart_minus = $cols[1];
                $prev_chrEnd_minus = $cols[2];
                $prev_blockSizes_minus = $cols[10];
                $prev_blockStarts_minus = $cols[11];
                $count_minus = $cols[4];
                $prev_chr_minus = $cols[0];
                $prev_name_minus = $cols[3];
                $prev_rgb_minus = $cols[8];
                $prev_blocks_minus = $cols[9];
            }
            else {
                print OUT "$prev_chr_minus\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_name_minus\t$count_minus\t\-\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_rgb_minus\t$prev_blocks_minus\t$prev_blockSizes_minus\t$prev_blockStarts_minus\n";
                $iso_count++;
                $prev_chrStart_minus = $cols[1];
                $prev_chrEnd_minus = $cols[2];
                $prev_blockSizes_minus = $cols[10];
                $prev_blockStarts_minus = $cols[11];
                $count_minus = $cols[4];
                $prev_chr_minus = $cols[0];
                $prev_name_minus = $cols[3];
                $prev_rgb_minus = $cols[8];
                $prev_blocks_minus = $cols[9];
            }
        }
    }
}
if ($count_plus > 0) {#prints out the last feature (plus strand)
    print OUT "$prev_chr_plus\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_name_plus\t$count_plus\t\+\t$prev_chrStart_plus\t$prev_chrEnd_plus\t$prev_rgb_plus\t$prev_blocks_plus\t$prev_blockSizes_plus\t$prev_blockStarts_plus\n";
    $iso_count++;
}

if ($count_minus > 0) {#prints out the last feature (minus strand)
    print OUT "$prev_chr_minus\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_name_minus\t$count_minus\t\-\t$prev_chrStart_minus\t$prev_chrEnd_minus\t$prev_rgb_minus\t$prev_blocks_minus\t$prev_blockSizes_minus\t$prev_blockStarts_minus\n";
    $iso_count++;
}

close(INF);
close(OUT);

my @ann;

open(INF, "<$ann_file") or die "couldn't open file";
while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
	push (@ann, $line); #puts each line of the annotation file into an array to be checked later
}
close(INF);

open(INF, "<$test_file.isoforms.bed") or die "couldn't open file";
open(OUT, ">$test_file.validated_transcripts.bed");

print "Checking for annotated isoforms...\n";

my $upper_limit_s;
my $lower_limit_s;
my $upper_limit_e;
my $lower_limit_e;
my $ann_count = 0;

print OUT "track type=bed name=\"$test_file.validated_transcripts.bed\" description=\"validated transcript structures from transcript_structure_validator.pl\"\n";

while (my $line = <INF>) {
	chomp($line);
    next if ($line =~ /^track/); #skips the track definition line
    my @val_cols = split("\t", $line);
    my $found_flag=0;
    foreach my $ann (@ann) {
        my $val_introns = 0;
        my $ann_introns = 0;
        my @ann_cols = split("\t", $ann);
        next if ($val_cols[5] ne $ann_cols[5]);
        next if ($val_cols[9] ne $ann_cols[9]);
        if ($val_cols[5] eq "+") {
            $upper_limit_s = $ann_cols[1] + $start_dist;
            $lower_limit_s = $ann_cols[1] - $start_dist;
            $upper_limit_e = $ann_cols[2] + $end_dist;
            $lower_limit_e = $ann_cols[2] - $end_dist;
        }
        if ($val_cols[5] eq "-") {
            $upper_limit_s = $ann_cols[1] + $end_dist;
            $lower_limit_s = $ann_cols[1] - $end_dist;
            $upper_limit_e = $ann_cols[2] + $start_dist;
            $lower_limit_e = $ann_cols[2] - $start_dist;
        }
        if (($val_cols[1] >= $lower_limit_s) and ($val_cols[1] <= $upper_limit_s)) {
            if (($val_cols[2] >= $lower_limit_e) and ($val_cols[2] <= $upper_limit_e)) {
                if ($val_cols[9] == 1) {
                    print OUT $val_cols[0], "\t", $val_cols[1], "\t", $val_cols[2], "\t", $ann_cols[3], "_", $val_cols[3], "\t", $val_cols[4], "\t", $val_cols[5], "\t", $val_cols[6], "\t", $val_cols[7], "\t", $ann_cols[8], "\t", $val_cols[9], "\t", $val_cols[10], "\t", $val_cols[11], "\n";
                    $ann_count++;
                    $found_flag=1;
                }
                else {
                    my $val_intron_number = $val_cols[9] - 1;
                    my @val_block_sizes = split(",", $val_cols[10]);
                    my @val_block_starts = split(",", $val_cols[11]);
                    for (my $i = 0; $i < $val_intron_number; $i = $i + 1) { #for the transcript currently in the "while" loop, creates an array of intron start sites relative to the genome
                        $start = $val_cols[1] + $val_block_sizes[$i] + $val_block_starts[$i];
                        $val_introns = "$val_introns:$start";
                    }
                    for (my $i2 = 1; $i2 < $val_cols[9]; $i2 = $i2 + 1) { #for the transcript currently in the "while" loop, creates an array of intron end sites relative to the genome
                        $end = $val_cols[1] + $val_block_starts[$i2];
                        $val_introns = "$val_introns:$end";
                    }
                    my $ann_intron_number = $ann_cols[9] - 1;
                    my @ann_block_sizes = split(",", $ann_cols[10]);
                    my @ann_block_starts = split(",", $ann_cols[11]);
                    for (my $i = 0; $i < $ann_intron_number; $i = $i + 1) { #for the annotation currently in the "foreach" loop, creates an array of intron start sites relative to the genome
                        $start = $ann_cols[1] + $ann_block_sizes[$i] + $ann_block_starts[$i];
                        $ann_introns = "$ann_introns:$start";
                    }
                    for (my $i2 = 1; $i2 < $ann_cols[9]; $i2 = $i2 + 1) { #for the annotation currently in the "foreach" loop, creates an array of intron end sites relative to the genome
                        $end = $ann_cols[1] + $ann_block_starts[$i2];
                        $ann_introns = "$ann_introns:$end";
                    }
                    if ($val_introns eq $ann_introns) {
                        print OUT $val_cols[0], "\t", $val_cols[1], "\t", $val_cols[2], "\t", $ann_cols[3], "_", $val_cols[3], "\t", $val_cols[4], "\t", $val_cols[5], "\t", $val_cols[6], "\t", $val_cols[7], "\t", $ann_cols[8], "\t", $val_cols[9], "\t", $val_cols[10], "\t", $val_cols[11], "\n";
                        $ann_count++;
                        $found_flag=1;
                    }
                }
            }
        }
    }
    if ($found_flag == 0){
        print OUT $line, "\n";
    }
}

close(INF);
close(OUT);

open(OUT, ">validated_isoforms_stats.txt");

my $novel_count = $iso_count - $ann_count;
print OUT "$iso_count validated transcripts\n\t$novel_count novel\n\t$ann_count annotated\n";

close(OUT);

print "------------------------------------------------\n$good_start_number sequences have validated start sites.\n";
print "$good_start_end_number sequences have validated start and end sites.\n";
print "$validated_count fully validated sequences collapse into $iso_count distinct isoforms.\n";
print "$ann_count isoforms match annotated transcripts.\n";

system("rm \Q$test_file\E.validated_refined.bed");
system("rm \Q$test_file\E.validated_unrefined.bed");
system("rm \Q$test_file\E.isoforms.bed");
system("rm \Q$test_file\E.bed");