#!/usr/bin/perl

# TRIMD_start_validator.pl
# Copyright (C) 2016 Flemington Lab (except Paraclu portion)

#This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

#This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

#Accepts a SAM file of Iso-Seq fl data, a SAM file of CAGE data, and a bed file of annotated polyadenylated transcripts. Counts the number of non-clipped Iso-Seq reads with 5' starts at each genomic position and estimates consensus locations of clusters of 5' starts. Uses Paraclu to identify clusters of 5' starts in the CAGE data. Output includes a bedgraph file of Iso-Seq 5' starts, a bed file of the weighted centers of Iso-Seq start clusters, a bedgraph file of CAGE tag 5' starts, a bed file of the weighted centers of Paraclu-identified CAGE 5' start clusters, and a bed file of Iso_seq 5' starts supported by the CAGE data, with their annotation status noted.

#USAGE:
# perl <PATH/TRIMD_start_validator.pl> </PATH/Iso-Seq_sam_file> </PATH/CAGE_file> </PATH/Annotation_bed_file>

use warnings;
use strict;

die "USAGE: 'perl <PATH/TRIMD_start_validator.pl> </PATH/Iso-Seq_sam_file> </PATH/CAGE_file> </PATH/Annotation_bed_file>'" unless @ARGV == 3;

my ($SMRT_file, $CAGE_file, $ann_file) = @ARGV;

#print "Enter name of viral chromosome (e.g. chrEBV_Akata_inverted): ";
#my $viral_chr = <STDIN>;
#chomp $viral_chr;

my $distance_between_SMRT_peaks;
my $min_tags;
my $min_dens;
my $min_length;
my $max_length;
my $dist_SMRT_CAGE;
my $min_SMRT;
my $ann_dist;

print "Use default parameters [y/n]? ";
my $answer = <STDIN>;
chomp $answer;

if ($answer eq "y") {
    $distance_between_SMRT_peaks = 8;
    $min_tags = 15;
    $min_dens = 2;
    $min_length = 1;
    $max_length = 20;
    $dist_SMRT_CAGE = 3;
    $min_SMRT = 1;
    $ann_dist = 10;
}
else {
    print "Enter desired window for collapsing Iso-Seq 5' starts (e.g. 8): ";
    $distance_between_SMRT_peaks = <STDIN>;
    chomp $distance_between_SMRT_peaks;
    
    print "Enter minimum tags per CAGE cluster (e.g. 15): ";
    $min_tags = <STDIN>;
    chomp $min_tags;
    
    print "Enter minimum relative density for CAGE clusters (e.g. 2): ";
    $min_dens = <STDIN>;
    chomp $min_dens;
    
    print "Enter minimum CAGE cluster length (e.g. 1): ";
    $min_length = <STDIN>;
    chomp $min_length;
    
    print "Enter maximum CAGE cluster length (e.g. 20): ";
    $max_length = <STDIN>;
    chomp $max_length;
    
    print "Enter desired maximum allowable distance between Iso-Seq and CAGE 5' starts (e.g. 3): ";
    $dist_SMRT_CAGE = <STDIN>;
    chomp $dist_SMRT_CAGE;
    
    print "Enter minimum number of SMRT reads to report a 5' start (e.g. 1): ";
    $min_SMRT = <STDIN>;
    chomp $min_SMRT;
    
    print "Enter maximum distance in bp from an annotated start to be called as 'annotated' (e.g. 10): ";
    $ann_dist = <STDIN>;
    chomp $ann_dist;
}

print "================================================\n";

#####----------SMRT FILE CHROMOSOME EXTRACTION-------------######

open(INF, "<$SMRT_file") or die "couldn't open file";

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

my @total_ends_found;
my @novel_ends_found;
my @SMRT_ends_ann;
my @ann_ends_found;
my @ann_ends;

foreach my $chrom (sort keys %chroms) {
    print "Processing $chrom:\n";
    
    
    #####----------SMRT FILE PROCESSING-------------######
    system("awk '\$3==\"$chrom\"' \Q$SMRT_file\E \| sort -k 4,4n > \Q$SMRT_file\E.sorted.temp");
    system("awk '\$2==0' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.plus.sam.temp");
    system("awk '\$2==16' \Q$SMRT_file\E.sorted.temp > \Q$SMRT_file\E.sorted.minus.sam.temp");
    system("rm \Q$SMRT_file\E.sorted.temp");
    
    #processing of PLUS SMRT sam file
    open(INF, "<$SMRT_file.sorted.plus.sam.temp") or die "couldn't open file";
    open(OUT, ">$SMRT_file.sorted.plus.sam.read_starts.bedgraph") or die "couldn't open file";
    
    my $previous_coordinate=1;
    my $count=0;
    my $previous_chr = "start";
    print "Processing Iso-Seq plus strand reads...\n";
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        next if ($cols[5] =~ m/^\d+S/); #skips reads clipped at the 5' end
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
                print OUT $previous_chr, "\t", $previous_coordinate-1, "\t", $previous_coordinate, "\t", $count, "\n"; #prints to output file, converting to chrStart 0-based bedgraph coordinate
                $previous_chr = $cols[2];
                $previous_coordinate = $cols[3];
                $count = $split_id[1];
            }
        }
    }
    
    print OUT $previous_chr, "\t", $previous_coordinate-1, "\t", $previous_coordinate, "\t", $count, "\n"; #prints the last start coordinates to output file
    close(INF);
    close(OUT);
    
    system("rm \Q$SMRT_file\E.sorted.plus.sam.temp");
    
    #processing of MINUS SMRT sam file
    open(INF, "<$SMRT_file.sorted.minus.sam.temp") or die "couldn't open file";
    open(OUT, ">$SMRT_file.sorted.minus.sam.read_starts.bedgraph.temp") or die "couldn't open file";
    
    my @CIGAR_dist;
    my $sum;
    my %minus_starts;
    print "Processing Iso-Seq minus strand reads...\n";
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        next if ($cols[5] =~ m/\d+S$/); #skips reads soft-clipped at the 5' end
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
            push (@CIGAR_dist, $1);
        }
        $sum += $_ for @CIGAR_dist;
        my $start_coord = $cols[3] + $sum - 1; #subtract 1 to account for start/end inclusion
        my $chr_start_coord = "$cols[2]\:$start_coord"; #combines the chromosome and 5' end coordinate into a key to use for the hash
        $sum = 0;
        @CIGAR_dist = ();
        my @split_id = split("\/", $cols[0]); #extracts the read depth for this putative isoform from its id
        if (exists $minus_starts{$chr_start_coord}) { #if the key is already in the hash, increases the value (count) by the read depth for that putative isoform
            $minus_starts{$chr_start_coord} = $minus_starts{$chr_start_coord} + $split_id[1];
        }
        else {
            $minus_starts{$chr_start_coord} = $split_id[1]; #if the key is not already in the hash, adds it with a value (count) of the read depth for that putative isoform
        }
    }
    
    foreach my $chr_start_coord (sort keys %minus_starts) { #prints out a(n inadequately) sorted temporary bedgraph file
        my @split_keys = split("\:", $chr_start_coord);
        print OUT $split_keys[0], "\t", $split_keys[1]-1, "\t", $split_keys[1], "\t-", $minus_starts{$chr_start_coord}, "\n"; #prints to output file, converting chrStart to 0-based bedgraph coordinates
    }
    close(INF);
    close(OUT);
    
    system("sort -k 1,1 -k 2,2n \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph.temp > \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph");
    
    system("cat \Q$SMRT_file\E.sorted.plus.sam.read_starts.bedgraph \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph.temp | sort -k2,3n > \Q$SMRT_file\E.\Q$chrom\E.read_starts.bedgraph.noheader");
    
    system("rm \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph.temp");
    system("rm \Q$SMRT_file\E.sorted.minus.sam.read_starts.bedgraph");
    system("rm \Q$SMRT_file\E.sorted.minus.sam.temp");
    system("rm \Q$SMRT_file\E.sorted.plus.sam.read_starts.bedgraph");
    
    #add header to bedgraph file
    open(INF, "<$SMRT_file.$chrom.read_starts.bedgraph.noheader") or die "couldn't open file";
    open(OUT, ">$SMRT_file.$chrom.read_starts.bedgraph") or die "couldn't open file";
    
    print OUT "track type=bedgraph name=\"$SMRT_file.$chrom.read_starts.bedgraph\" description=\"5' starts of SMRT reads from start_finder_sam_to_bed.pl\"\n";
    while (my $line = <INF>) {
        print OUT $line;
    }
    close(OUT);
    close(INF);
    
    system("rm \Q$SMRT_file\E.\Q$chrom\E.read_starts.bedgraph.noheader");
    
    #make a bed file from the SMRT bedgraph file:
    open(INF, "<$SMRT_file.$chrom.read_starts.bedgraph") or die "couldn't open file";
    open(OUT, ">$SMRT_file.starts.temp.bed") or die "couldn't open file";
    
    print "Combining Iso-Seq 5' starts within $distance_between_SMRT_peaks of each other and calculating consensus 5' starts...\n";
    collapse_bedgraph($chrom, $distance_between_SMRT_peaks);
    
    close(INF);
    close(OUT);
    
    system("sort -k 1,1 -k 2,2n \Q$SMRT_file\E.starts.temp.bed > \Q$SMRT_file\E.starts.bed.noheader");
    system("rm \Q$SMRT_file.starts.temp.bed\E");
    
    #add header to bed file
    open(INF, "<$SMRT_file.starts.bed.noheader") or die "couldn't open file";
    open(OUT, ">$SMRT_file.$chrom.SMRT_starts.bed") or die "couldn't open file";
    
    print OUT "track type=bed name=\"$SMRT_file.$chrom.SMRT_starts.bed\" description=\"consensus 5' starts of Iso-Seq reads within $distance_between_SMRT_peaks bp collapsed to weighted center from start_finder_sam_to_bed.pl\"\n";
    while (my $line = <INF>) {
        print OUT $line;
    }
    close(OUT);
    close(INF);
    
    system("rm \Q$SMRT_file\E.starts.bed.noheader");
    
    #####----------PROCESSING CAGE DATA-------------######
    
    print "Preparing CAGE file...\n";
    
    system("awk '\$3==\"$chrom\"' \Q$CAGE_file\E \| sort -k 4,4n > \Q$CAGE_file\E.sorted.temp");
    system("awk '\$2==0 \|\| \$2==81 \|\| \$2==83 \|\| \$2==89 \|\| \$2==137 \|\| \$2==161 \|\| \$2==163' \Q$CAGE_file\E.sorted.temp > \Q$CAGE_file\E.sorted.plus.sam.temp");
    system("awk '\$2==16 \|\| \$2==73 \|\| \$2==97 \|\| \$2==99 \|\| \$2==145 \|\| \$2==147 \|\| \$2==153' \Q$CAGE_file\E.sorted.temp > \Q$CAGE_file\E.sorted.minus.sam.temp");
    
    #processing of plus CAGE sam file
    
    open(INF, "<$CAGE_file.sorted.plus.sam.temp") or die "couldn't open file";
    open(OUT, ">$CAGE_file.read_starts.txt") or die "couldn't open file";
    open(OUT2, ">$CAGE_file.starts.bedgraph.temp");
    
    my $prev_coord=0.5;
    my $start_count=0;
    
    while (my $line = <INF>) {
        chomp($line);
        next if ($line =~ m/^@/); #skips header lines
        my @cols = split("\t", $line);
        if ($cols[3] == $prev_coord) {
            $start_count++; #increases the count by 1
        }
        
        else {
            if ($prev_coord == 0.5) { #doesn't print out the placeholder first line.
                $prev_coord = $cols[3];
                $start_count = 1;
            }
            
            else {
                print OUT $chrom, "\t+\t", $prev_coord, "\t", $start_count, "\n"; #prints to output txt file for Paraclu
                print OUT2 $chrom, "\t", $prev_coord-1, "\t", $prev_coord, "\t", $start_count, "\n"; #prints to output bedgraph file
                $prev_coord = $cols[3];
                $start_count = 1;
            }
        }
    }
    print OUT "$chrom\t+\t$prev_coord\t$start_count\n"; #prints the last start coordinates to output file
    close(INF);
    
    system("rm \Q$CAGE_file\E.sorted.plus.sam.temp");
    
    #processing of MINUS CAGE sam file
    open(INF, "<$CAGE_file.sorted.minus.sam.temp") or die "couldn't open file";
    
    my @read_dist;
    my $dist_sum;
    my %minus_start;
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        while ($cols[5] =~ /(\d+)[DMNX=]/g) { #these lines use the CIGAR string to determine the downstream coordinate
            push (@read_dist, $1);
        }
        $dist_sum += $_ for @read_dist;
        my $start_coord = $cols[3] + $dist_sum - 1; #subtract one to account for start/end inclusion
        $dist_sum = 0;
        @read_dist = ();
        if (exists $minus_start{$start_coord}) { #if the key is already in the hash, increases the value (count) by the read depth for that putative isoform
            $minus_start{$start_coord} = $minus_start{$start_coord} + 1;
        }
        else {
            $minus_start{$start_coord} = 1; #if the key is not already in the hash, adds it with a value (count) of the read depth for that putative isoform
        }
    }
    foreach my $start_coord (sort keys %minus_start) { #prints out a(n inadequately) sorted file. Doesn't need to be sorted because Paraclu will do that anyways.
        print OUT "$chrom\t-\t$start_coord\t$minus_start{$start_coord}\n";
        print OUT2 $chrom, "\t", $start_coord-1, "\t", $start_coord, "\t-", $minus_start{$start_coord}, "\n"; #prints to output bedgraph file
    }
    close(INF);
    close(OUT);
    close(OUT2);
    
    system("rm \Q$CAGE_file\E.sorted.minus.sam.temp");
    system("rm \Q$CAGE_file\E.sorted.temp");
    system("sort -k2,3n \Q$CAGE_file\E.starts.bedgraph.temp > \Q$CAGE_file\E.starts.bedgraph.noheader");
    system("rm \Q$CAGE_file\E.starts.bedgraph.temp");
    
    #add header to bedgraph file
    open(INF, "<$CAGE_file.starts.bedgraph.noheader") or die "couldn't open file";
    open(OUT, ">$CAGE_file.$chrom.starts.bedgraph") or die "couldn't open file";
    
    print OUT "track type=bedgraph name=\"$CAGE_file.$chrom.starts.bedgraph\" description=\"5' starts of CAGE tags from start_finder_sam_to_bed.pl\"\n";
    while (my $line = <INF>) {
        print OUT $line;
    }
    close(INF);
    close(OUT);
    
    system("rm \Q$CAGE_file\E.starts.bedgraph.noheader");
    
    #Running Paraclu to define clusters
    
    open(INF, "<$CAGE_file.read_starts.txt") or die "couldn't open file";
    open(OUT, ">$CAGE_file.paraclu.txt.temp");
    
    # paraclu.pl: perform parametric clustering of data attached to sequences
    
    # Written by Martin C Frith 2006
    # Genome Exploration Research Group, RIKEN GSC and
    # Institute for Molecular Bioscience, University of Queensland
    
    # This program reads in a list of numeric values attached to positions
    # in sequences. The list should have four tab- (or space-) separated
    # columns containing: the sequence name, the strand, the position, and
    # the value. (Multiple values for the same sequence/strand/position
    # will be summed.) It outputs the clusters as eight tab-separated
    # columns: sequence name, strand, start, end, number of values, sum of
    # values, min d, max d. See below for the meaning of "d".
    
    # An example line of input:
    # chr1    +       17689   3
    # Clustering is performed separately for different strands (as if each
    # strand were a completely different sequence).  It does not matter
    # whether the position uses 0-based or 1-based coordinates: the
    # program does not care, and the output will be consistent with the
    # input.
    
    # The clusters are defined as follows. A cluster is a maximal scoring
    # segment, where the score of any segment is: the sum of the values in
    # the segment minus d times the size of the segment. Large values of d
    # give smaller, tighter clusters and small values of d give larger,
    # looser clusters. The program finds all possible clusters for any
    # value of d, and annotates each cluster with the maximum and minimum
    # values of d that produce it. The ratio max d / min d provides a
    # measure of the cluster's "stability".
    
    # The output will include two types of obvious/trivial/degenerate
    # clusters: those that cover single positions, and those that cover
    # all of the positions in a sequence.  For many purposes, it would be
    # best to ignore these cases.
    
    use strict;
    use List::Util qw(min max);
    
    my %data;
    
    #warn "reading...\n";
    
    while (<INF>) {
        chomp;
        s/#.*//;  # ignore comments
        next unless /\S/;  # skip blank lines
        
        my ($seq, $strand, $pos, $value) = split;
        my $key = "$seq $strand";
        push @{$data{$key}}, [ $pos, $value ];
    }
    
    warn "Clustering CAGE data...\n";
    
    print OUT "# sequence, strand, start, end, sites, sum of values, min d, max d\n";
    
    for my $key (sort keys %data) {  # iterate over sequences / strands
        my ($seq, $strand) = split " ", $key;
        my $sites = $data{$key};
        
        @$sites = sort { $$a[0] <=> $$b[0] } @$sites;  # sort by position
        
        my $clusters = all_clusters($sites);
        
        for my $c (@$clusters) {
            my ($beg, $end, $tot, $sit, $min, $max) = @$c;
            my $beg_pos = $$sites[$beg][0];
            my $end_pos = $$sites[$end][0];
            printf OUT "$seq\t$strand\t$beg_pos\t$end_pos\t$sit\t$tot\t%.3g\t%.3g\n",
            $min, $max;
        }
    }
    
    ### Generic code to find clusters in a sparse sequence of values: ###
    
    sub all_clusters {
        our $inf = 1e100;  # hopefully much bigger than any value in the input
        our $sites = shift;  # input: reference to array of site locations & values
        our $clusters = [];  # output: reference to array of clusters
        get_clusters(0, $#$sites, -$inf);
        return $clusters;
    }
    
    # get clusters of sites between beg and end with density > min_density
    sub get_clusters {
        our ($clusters, $inf);
        my ($beg, $end, $min_density) = @_;
        
        my ($prefix, $pmin, $ptot, $psit) = weakest_prefix($beg, $end);
        my ($suffix, $smin, $stot, $ssit) = weakest_suffix($beg, $end);
        $ptot == $stot and $psit == $ssit or die "internal error!";
        my $max_density = min $pmin, $smin;
        
        unless ($max_density == $inf) {
            my $break = $pmin < $smin ? $prefix + 1 : $suffix;
            my $new_min = max $min_density, $max_density;
            get_clusters($beg, $break-1, $new_min);
            get_clusters($break, $end, $new_min);
        }
        
        push @$clusters, [ $beg, $end, $ptot, $psit, $min_density, $max_density ]
        if $max_density > $min_density;
    }
    
    # get least dense prefix (and total of values & sites)
    sub weakest_prefix {
        our ($sites, $inf);
        my ($beg, $end) = @_;
        
        my $beg_pos = $$sites[$beg][0];
        my $min_density = $inf;
        my $min_prefix = $end;
        my $tot = 0;
        my $sit = 0;
        
        for (my $i = $beg; $i < $end; ++$i) {
            $tot += $$sites[$i][1];
            next if $$sites[$i][0] == $$sites[$i+1][0];  # idiot-proofing
            ++$sit;
            my $dist = $$sites[$i+1][0] - $beg_pos;
            my $density = $tot / $dist;
            if ($density < $min_density) {
                $min_prefix = $i;
                $min_density = $density;
            }
        }
        
        $tot += $$sites[$end][1];
        ++$sit;
        return ($min_prefix, $min_density, $tot, $sit);
    }
    
    # get least dense suffix (and total of values & sites)
    sub weakest_suffix {
        our ($sites, $inf);
        my ($beg, $end) = @_;
        
        my $end_pos = $$sites[$end][0];
        my $min_density = $inf;
        my $min_suffix = $beg;
        my $tot = 0;
        my $sit = 0;
        
        for (my $i = $end; $i > $beg; --$i) {
            $tot += $$sites[$i][1];
            next if $$sites[$i][0] == $$sites[$i-1][0];  # idiot-proofing
            ++$sit;
            my $dist = $end_pos - $$sites[$i-1][0];
            my $density = $tot / $dist;
            if ($density < $min_density) {
                $min_suffix = $i;
                $min_density = $density;
            }
        }
        
        $tot += $$sites[$beg][1];
        ++$sit;
        return ($min_suffix, $min_density, $tot, $sit);
    }
    
    close(INF);
    close(OUT);
    
    system("sort -k2,2 -k3,3n -k4,4rn \Q$CAGE_file\E.paraclu.txt.temp > \Q$CAGE_file\E.paraclu.txt");
    system("rm \Q$CAGE_file\E.paraclu.txt.temp");
    
    #filtering clusters:
    print "Extracting CAGE clusters containing $min_tags tags, density fold change at least $min_dens, from $min_length to $max_length bp long\n";
    
    my $length;
    my $dens;
    my $prev_start = 0;
    my $prev_end = 0;
    
    open(INF, "<$CAGE_file.paraclu.txt") or die "couldn't open file";
    open(OUT, ">$CAGE_file.clusters.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";
    
    while (my $line = <INF>) { #extracts clusters meeting the criteria. Excludes subclusters.
        chomp($line);
        next if ($line =~ /^#/); #skips the header line
            my @cols = split("\t", $line);
        next if ($cols[5] < $min_tags); #checks number of tags
        $length = $cols[3] - $cols[2] + 1; #gets length of cluster
        if (($length >= $min_length) and ($length <= $max_length)) { #checks if length is in specified range
            $dens = $cols[7] / $cols[6]; #gets relative density
            if ($dens >= $min_dens) { #if relative density is high enough:
                next if (($cols[2] >= $prev_start) and ($cols[2] <= $prev_end)); #excludes subclusters
                if ($dens < 100) { #keep everything one-based (like sam) for now; will convert to zero-based in next step
                    printf OUT "%s\t%d\t%d\t%d%s%.1f\t%d\t%s\n", $cols[0], $cols[2], $cols[3], $cols[5], ":", $dens, $cols[5], $cols[1];   #limits the density output to 1 decimal place, but doesn't change huge numbers to exponents
                }
                else {
                    printf OUT "%s\t%d\t%d\t%d%s%.1e\t%d\t%s\n", $cols[0], $cols[2], $cols[3], $cols[5], ":", $dens, $cols[5], $cols[1]; #changes large numbers to exponents
                }
                $prev_start = $cols[2];
                $prev_end = $cols[3];
            }
        }
    }
    close(INF);
    close(OUT);
    
    system("rm \Q$CAGE_file\E.paraclu.txt");
    
    #getting weighted averages of Paraclu clusters:
    
    my $rangeStart_CAGE;
    my $rangeEnd_CAGE;
    my $strand_CAGE;
    my $CAGE_weighted_sum = 0;
    my $CAGE_weighted_average;
    
    open(INF, "<$CAGE_file.clusters.$min_tags.$min_dens.$min_length.$max_length.bed") or die "couldn't open file";
    open(OUT, ">$CAGE_file.$chrom.CAGE_starts.temp") or die "couldn't open file";
    
    while (my $line = <INF>) {
        chomp($line);
        my @cols = split("\t", $line);
        $rangeStart_CAGE = $cols[1];
        $rangeEnd_CAGE = $cols[2];
        $strand_CAGE = $cols[5];
        open(INF2, "<$CAGE_file.read_starts.txt") or die "couldn't open file";
        while (my $line2 = <INF2>) {
            chomp($line2);
            my @cols2 = split("\t", $line2);
            if ((($cols2[2]) >= $rangeStart_CAGE) and (($cols2[2]) <= $rangeEnd_CAGE) and ($cols2[1] eq $strand_CAGE)) {
                $CAGE_weighted_sum = $CAGE_weighted_sum + ($cols2[2]*$cols2[3]);
            }
        }
        $CAGE_weighted_average = sprintf("%1.0f", ($CAGE_weighted_sum/$cols[4]));
        if ($strand_CAGE eq "+") {
            print OUT $cols[0], "\t", $CAGE_weighted_average-1, "\t", $CAGE_weighted_average, "\t", $rangeStart_CAGE-1,  ":", $rangeEnd_CAGE-1, ":", $cols[3], "\t", $cols[4], "\t", $strand_CAGE, "\n"; #prints output, converting chrStart and range to 0-based
        }
        elsif ($strand_CAGE eq "-") {
            print OUT $cols[0], "\t", $CAGE_weighted_average-1, "\t", $CAGE_weighted_average, "\t", $rangeStart_CAGE,  ":", $rangeEnd_CAGE, ":-", $cols[3], "\t", $cols[4], "\t", $strand_CAGE, "\n"; #prints output, converting chrStart to 0-based but keeping range 1-based (because these are all chrEnds)
        }
        $CAGE_weighted_sum = 0;
        close(INF2);
    }
    
    close(INF);
    close(OUT);
    
    system("sort -k2,2n -k 3,3n \Q$CAGE_file\E.\Q$chrom\E.CAGE_starts.temp > \Q$CAGE_file\E.\Q$chrom\E.CAGE_starts.noheader");
    system("rm \Q$CAGE_file\E.\Q$chrom\E.CAGE_starts.temp");
    
    #add header to bed file
    
    open(INF, "<$CAGE_file.$chrom.CAGE_starts.noheader") or die "couldn't open file";
    open(OUT, ">$CAGE_file.$chrom.CAGE_starts.bed");
    
    print OUT "track type=bed name=\"$CAGE_file.$chrom.CAGE_starts.bed\" description=\"weighted averages of CAGE clusters between $min_length and $max_length bases long with at least $min_tags tags and relative density of at least $min_dens from start_finder_sam_to_bed.pl and paraclu\"\n";
    while (my $line = <INF>) {
        print OUT $line;
    }
    close(OUT);
    close(INF);
    
    system("rm \Q$CAGE_file\E.\Q$chrom\E.CAGE_starts.noheader");
    system("rm \Q$CAGE_file\E.read_starts.txt");
    system("rm \Q$CAGE_file\E.clusters.$min_tags.$min_dens.$min_length.$max_length.bed");
    
    #####----------SEEKING CAGE SUPPORT FOR SMRT STARTS-------------######
    
    open(INF, "<$CAGE_file.$chrom.CAGE_starts.bed" ) or die "couldn't open file";
    
    print "Extracting Iso-Seq 5' starts within $dist_SMRT_CAGE bases of CAGE clusters...\n";
    
    my %features_CAGE;
    my $key_combo_CAGE;
    
    while(my $line = <INF> ) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @cols = split("\t", $line);
        if ($cols[5] eq "+") { #for each line in the CAGE bed file, creates a key for the hash combining coordinate and strand. Selects chrStart for starts on the plus strand and chrEnd for starts on the minus strand.
            $key_combo_CAGE = "$cols[1]:$cols[5]";
        }
        if ($cols[5] eq "-") {
            $key_combo_CAGE = "$cols[2]:$cols[5]";
        }
        $features_CAGE{$key_combo_CAGE} = $cols[4]; #enters a count value for the key into the hash
    }
    
    close(INF);
    
    open(INF, "<$SMRT_file.$chrom.SMRT_starts.bed" ) or die "couldn't open file";
    open(OUT, ">$SMRT_file.$chrom.SMRT_starts.bed.CAGE_support.bed.temp");
    
    my $match_count;
    my $lower_limit;
    my $upper_limit;
    
    while(my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @SMRT_cols = split("\t", $line);
        next if (abs $SMRT_cols[4] < $min_SMRT); #skips starts without enough SMRT support
        foreach my $key_combo_CAGE (keys %features_CAGE) {
            my @CAGE_cols = split(":", $key_combo_CAGE);
            if ($SMRT_cols[5] eq "+") {
                $lower_limit = $SMRT_cols[1]-$dist_SMRT_CAGE;
                $upper_limit = $SMRT_cols[1]+$dist_SMRT_CAGE;
            }
            if ($SMRT_cols[5] eq "-") {
                $lower_limit = $SMRT_cols[2]-$dist_SMRT_CAGE;
                $upper_limit = $SMRT_cols[2]+$dist_SMRT_CAGE;
            }
            if (($SMRT_cols[5] eq $CAGE_cols[1]) and ($CAGE_cols[0] >= $lower_limit) and ($CAGE_cols[0] <= $upper_limit)) {
                if ($match_count) { #if more than one CAGE start matches the SMRT start, selects the CAGE end with the most tags
                    if ($features_CAGE{$key_combo_CAGE} > $match_count) {
                        $match_count = $features_CAGE{$key_combo_CAGE};
                    }
                }
                else {
                    $match_count = $features_CAGE{$key_combo_CAGE};
                }
            }
        }
        if ($match_count) {
            my $name = "$SMRT_cols[4].IsoSeq_$match_count.CAGE";
            my $count = $match_count + $SMRT_cols[4];
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\t$name\t$count\t$SMRT_cols[5]\t$SMRT_cols[3]\n";
            undef($match_count);
        }
        else {
            my @range_cols = split (":", $SMRT_cols[3]);
            print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\t$range_cols[2].IsoSeq\t$range_cols[2]\t$SMRT_cols[5]\t$SMRT_cols[3]\n";
        }
    }
    
    close(OUT);
    close(INF);
    
    
    #####----------COMPARING TO ANNOTATED STARTS-------------######
    open(INF, "<$ann_file" ) or die "couldn't open file";
    
    print "Processing annotation file...\n";
    
    #extract 5' starts from the annotation file:
    #annotation file must be sorted by chrStart then chrEnd!
    my @annotated_starts;
    my $plus_prev_coord = 0;
    my $minus_prev_coord = 0;
    
    while(my $line = <INF>) {
        chomp($line);
        next if ($line =~ /^track/); #skips the track definition line
        my @ann_cols = split("\t", $line);
        next if $ann_cols[0] ne $chrom; #skip lines that aren't viral
        if ($ann_cols[5] eq "+") {
            if ($ann_cols[1] != $plus_prev_coord) {
                push (@annotated_starts, "$ann_cols[1]:$ann_cols[5]");
                $plus_prev_coord = $ann_cols[1];
            }
        }
        elsif ($ann_cols[5] eq "-"){
            if ($ann_cols[2] != $minus_prev_coord) {
                push (@annotated_starts, "$ann_cols[2]:$ann_cols[5]");
                $minus_prev_coord = $ann_cols[2];
            }
        }
    }
    
    my $annotated = scalar @annotated_starts;
    
    close(INF);
    
    #compare starts in the altered SMRT starts file (that already has info about CAGE starts) with annotated starts
    
    open(INF, "<$SMRT_file.$chrom.SMRT_starts.bed.CAGE_support.bed.temp" ) or die "couldn't open file";
    open(OUT, ">$SMRT_file.$chrom.validated_starts.bed");
    
    print "Comparing Iso-seq starts to annotated starts...\n";
    
    print OUT "track type=bedDetail name=\"$SMRT_file.$chrom.validated_starts.bed\" description=\"consensus Iso-Seq 5' starts of collapse value 8 supported by at least $min_SMRT read(s) within $dist_SMRT_CAGE bp of CAGE clusters or within $ann_dist bp of annotated starts. From start_finder_sam_to_bed.pl\"\n";
    
    my $annotated_found_by_SMRT = 0;
    my $novel_found_by_SMRT_CAGE = 0;
    my $SMRT_annotated = 0; #this is different than $annotated_found_by_SMRT because depending on input parameters two SMRT starts may correspond to a single annotated start or vice versa.
    
    while(my $line = <INF>) {
        chomp($line);
        my @SMRT_cols = split("\t", $line);
        my $found_flag=0;
        foreach my $ann_start (@annotated_starts) {
            my @ann_cols = split(":", $ann_start);
            my $lower_limit = $ann_cols[0]-$ann_dist;
            my $upper_limit = $ann_cols[0]+$ann_dist;
            if ($SMRT_cols[5] eq "+") {
                if (($SMRT_cols[5] eq $ann_cols[1]) and ($SMRT_cols[1]>=$lower_limit) and ($SMRT_cols[1]<=$upper_limit)) {
                    if ($found_flag == 0) {
                        print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tann_$SMRT_cols[5]_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\t$SMRT_cols[6]\n";
                        $found_flag = 1;
                        $annotated_found_by_SMRT++; #counts multiple annotated starts near SMRT starts
                        $SMRT_annotated++; #only counts one annotated start per SMRT start
                    }
                    elsif ($found_flag == 1) {
                        $annotated_found_by_SMRT++;
                    }
                }
            }
            if ($SMRT_cols[5] eq "-") {
                if (($SMRT_cols[5] eq $ann_cols[1]) and ($SMRT_cols[2]>=$lower_limit) and ($SMRT_cols[2]<=$upper_limit)) {
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
            if ($SMRT_cols[3] =~ /.+IsoSeq_.+CAGE/) {
                print OUT "$SMRT_cols[0]\t$SMRT_cols[1]\t$SMRT_cols[2]\tnov_$SMRT_cols[5]_$SMRT_cols[3]\t$SMRT_cols[4]\t$SMRT_cols[5]\t$SMRT_cols[6]\n";
                $novel_found_by_SMRT_CAGE++;
            }
        }
    }
    
    my $total_found = $SMRT_annotated + $novel_found_by_SMRT_CAGE;
    
    close(INF);
    close(OUT);
    
    print "------------------------------------------------\n";
    
    push (@total_ends_found, $total_found);
    push (@novel_ends_found, $novel_found_by_SMRT_CAGE);
    push (@SMRT_ends_ann, $SMRT_annotated);
    push (@ann_ends_found, $annotated_found_by_SMRT);
    push (@ann_ends, $annotated);
    
    #open(OUT, ">${chrom}_validated_starts_stats.txt");
    
    if ($total_found > 0) {
        if ($SMRT_annotated != $annotated_found_by_SMRT) {
            print "$total_found 5' starts found. $novel_found_by_SMRT_CAGE are novel, $SMRT_annotated are annotated.  $annotated_found_by_SMRT out of $annotated total annotated 5' starts are found.\nNote that two annotated starts may be within $ann_dist bp of a single Iso-Seq start or vice versa.\n";
            #print OUT "$chrom\n$total_found 5' starts\n\t$novel_found_by_SMRT_CAGE novel\n\t$SMRT_annotated annotated\n$annotated starts in annotation file\n\t$annotated_found_by_SMRT detected by Iso-Seq\n\ninput files:\n\t$SMRT_file\n\t$CAGE_file\n\t$ann_file\n";
        }
        else {
            print "$total_found 5' starts found. $novel_found_by_SMRT_CAGE are novel, $SMRT_annotated are annotated (out of a total of $annotated annotated 5' starts).\n";
            #print OUT "$chrom\n\n$total_found 5' starts\n\t$novel_found_by_SMRT_CAGE novel\n\t$SMRT_annotated annotated\n$annotated starts in annotation file\n\ninput files:\n\t$SMRT_file\n\t$CAGE_file\n\t$ann_file\n";
        }
    }
    else {
        print "No validated starts found.\n";
        #print OUT "No validated starts found.\n\ninput files:\n\t$SMRT_file\n\t$CAGE_file\n\t$ann_file\n";
    }
    
    #close(OUT);
    
    system("rm \Q$SMRT_file\E.\Q$chrom\E.SMRT_starts.bed.CAGE_support.bed.temp");
    
    print "================================================\n";
    
}



system("cat $CAGE_file.*.CAGE_starts.bed > ${CAGE_file}_CAGE_starts.bed");
system("cat $CAGE_file.*.starts.bedgraph > ${CAGE_file}_starts.bedgraph");
system("cat $SMRT_file.*.SMRT_starts.bed > ${SMRT_file}_SMRT_starts.bed");
system("cat $SMRT_file.*.read_starts.bedgraph > ${SMRT_file}_read_starts.bedgraph");
system("cat $SMRT_file.*.validated_starts.bed > ${SMRT_file}_validated_starts.bed");
system("rm $CAGE_file.*.CAGE_starts.bed");
system("rm $CAGE_file.*.starts.bedgraph");
system("rm $SMRT_file.*.SMRT_starts.bed");
system("rm $SMRT_file.*.read_starts.bedgraph");
system("rm $SMRT_file.*.validated_starts.bed");

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

open(OUT, ">validated_starts.txt");

print OUT "$sum_total_found 5' starts\n\t$sum_novel_found novel\n\t$sum_SMRT_ann annotated\n$sum_total_ann starts in annotation file\n\t$sum_ann_found detected by Iso-Seq\n\ninput files:\n\t$SMRT_file\n\t$CAGE_file\n\t$ann_file\n";

close(OUT);

#########################
sub collapse_bedgraph {
    my ($chrom, $distance_between_peaks) = @_;
    my $prev_coord_plus = 1;
    my $prev_coord_minus = 1;
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
            if ($cols[1] <= $prev_coord_plus + ($distance_between_peaks)) { #if the coordinate is within the specified number of bp of the previous coordinate
                $count_sum_plus = $count_sum_plus + $cols[3]; #adds to the sums to eventually calculate the weighted average
                $weighted_coordinate_sum_plus = $weighted_coordinate_sum_plus + ($cols[1]*$cols[3]);
                push (@coords_plus, $cols[1]);
                $prev_coord_plus = $cols[1]; #sets the current coordinate as the "previous coordinate" before moving on
            }
            else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
                if ($first_plus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                    $count_sum_plus = $cols[3];
                    $weighted_coordinate_sum_plus = $cols[1]*$cols[3];
                    $prev_coord_plus = $cols[1];
                    push (@coords_plus, $cols[1]);
                    $first_plus = 0;
                }
                else {
                    $weighted_average_plus = sprintf("%1.0f", ($weighted_coordinate_sum_plus/$count_sum_plus)); #calculates weighted average
                    $chrStart_plus = $coords_plus[0];
                    $chrEnd_plus = pop(@coords_plus);
                    print OUT $chrom, "\t", $weighted_average_plus, "\t", $weighted_average_plus+1,  "\t", $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, "\t", $count_sum_plus, "\t+\n"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
                    @coords_plus = ($cols[1]);
                    $count_sum_plus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
                    $weighted_coordinate_sum_plus = $cols[1]*$cols[3];
                    $prev_coord_plus = $cols[1];
                }
            }
        }
        elsif ($cols[3] < 0) { #if this coordinate has a negative count...
            if ($cols[2] <= $prev_coord_minus + ($distance_between_peaks)) { #if the coordinate is within the specified number of bp of the previous coordinate
                $count_sum_minus = $count_sum_minus + $cols[3]; #adds to the sums to eventually calculate the weighted average
                $weighted_coordinate_sum_minus = $weighted_coordinate_sum_minus + ($cols[2]*$cols[3]);
                push (@coords_minus, $cols[2]);
                $prev_coord_minus = $cols[2]; #sets the current coordinate as the "previous coordinate" before moving on
            }
            else { #if the present coordinate is not within the specified number of bp of the previous coordinate, need to print out a feature
                if ($first_minus == 1) { #"first" flag avoids wonkiness if the first coordinate is far from coordinate 1 (don't need to print out a feature yet)
                    $count_sum_minus = $cols[3];
                    $weighted_coordinate_sum_minus = $cols[2]*$cols[3];
                    $prev_coord_minus = $cols[2];
                    push (@coords_minus, $cols[2]);
                    $first_minus = 0;
                }
                else {
                    $weighted_average_minus = sprintf("%1.0f", ($weighted_coordinate_sum_minus/$count_sum_minus)); #calculates weighted average.
                    $chrStart_minus = $coords_minus[0];
                    $chrEnd_minus = pop(@coords_minus);
                    print OUT $chrom, "\t", $weighted_average_minus-1, "\t", $weighted_average_minus, "\t", $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, "\t", abs($count_sum_minus), "\t-\n";
                    @coords_minus = ($cols[2]);
                    @coords_minus = ($cols[2]);
                    $count_sum_minus = $cols[3]; #sets "previous coordinate", count and sum of counts for the current coordinate
                    $weighted_coordinate_sum_minus = $cols[2]*$cols[3];
                    $prev_coord_minus = $cols[2];
                }
            }
        }
    }
    
    if ($count_sum_plus > 0) {#calculates and prints out weighted average for the last feature (plus strand)
        $weighted_average_plus = sprintf("%1.0f", ($weighted_coordinate_sum_plus/$count_sum_plus));
        $chrStart_plus = $coords_plus[0];
        $chrEnd_plus = pop(@coords_plus);
        print OUT $chrom, "\t", $weighted_average_plus, "\t", $weighted_average_plus+1,  "\t", $chrStart_plus, ":", $chrEnd_plus, ":", $count_sum_plus, "\t", $count_sum_plus, "\t+\n"; #prints out weighted average for plus strand features. Use printf to round the weighted average.
    }
    
    if ($count_sum_minus < 0) {#calculates and prints out weighted average for the last feature (minus strand)
        $weighted_average_minus = sprintf("%1.0f", ($weighted_coordinate_sum_minus/$count_sum_minus));
        $chrStart_minus = $coords_minus[0];
        $chrEnd_minus = pop(@coords_minus);
        print OUT $chrom, "\t", $weighted_average_minus-1, "\t", $weighted_average_minus, "\t", $chrStart_minus, ":", $chrEnd_minus, ":", $count_sum_minus, "\t", abs($count_sum_minus), "\t-\n";
    }
}