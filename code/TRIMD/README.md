# TRIMD
Transcriptome Resolution by Integration of Multi-platform Data
## As of 2019, all TRIMD code will be maintained in the repository  https://github.com/flemingtonlab/TRIMD. 

Scripts included:
* TRIMD_start_validator.pl
* TRIMD_junction_validator.pl
* TRIMD_end_validator.pl
* TRIMD_transcript_validator.pl

### Notes for all scripts:

Defaults are set using the Epstein-Barr virus Akata strain as a model.

Annotation file should contain only features for polyadenylated transcripts and must be sorted by chrStart, then chrEnd.

Fasta files of Iso-Seq data must have names formatted as putative_isoform_id/number_of_SMRT_reads/length, as from the Iso-Seq pipeline.

Circular genomes present problems for most aligners: they are generally depicted in a linear form in aligner reference files. Reads that span the locus where the genome is linearized do not map properly. To overcome this limitation you may run the alignments and TRIMD twice: once to the standard genome and once to a version of the genome that is linearized at its opposite point (position N/2, where N is the length of the genome). This will allow detection of transcripts that span the linearization site.


## TRIMD_start_validator.pl
```
USAGE: perl /PATH/TRIMD_start_validator.pl </PATH/SMRT_sam_file> </PATH/CAGE_file> </PATH/Annotation_bed_file>
```
Accepts a SAM file of Iso-Seq fl data, a SAM file of CAGE data, and a bed file of annotated polyadenylated transcripts. Counts the number of non-clipped SMRT reads with 5' starts at each genomic position and estimates consensus locations of clusters of 5' starts. Uses Paraclu to identify clusters of 5' starts in the CAGE data. Output includes BEDGRAPH files of all 5' starts, BED files of the weighted centers of start clusters and a BED file of Iso-Seq 5' starts supported by either the annotation or the CAGE data.
Paraclu was written by Martin C Frith 2006, Genome Exploration Research Group, RIKEN GSC and Institute for Molecular Bioscience, University of Queensland and distributed under the GNU General Public License.

INPUT

1. SAM file of Iso-Seq fl isoforms: this script was developed using data aligned with GMAP (-f samse option). Other aligners may also be appropriate.

2. SAM file of CAGE reads: this script was developed using data aligned with STAR. Other aligners may also be appropriate.

3. BED file of annotated polyadenylated transcripts: the annotation file MUST be sorted by chrStart, then chrEnd. If your annotation file contains non-polyadenylated transcripts or other features (e.g. repeat regions, promoters) these should be removed to avoid false positives.


PARAMETERS

1. (Viral) chromosome name: the name of the chromosome under investigation. This must match between both SAM files (field 3) and the BED annotation file (field 1). Does not have to be viral, but for organisms with multiple chromosomes only one chromosome can be examined at a time.
No default: must be entered at prompt

2. Window for collapsing Iso-Seq 5' starts: Iso-Seq 5' starts within this number of bases of each other will be considered to represent the same transcription start site. The consensus transcription start site is determined by calculating an average of the coordinates in the cluster, weighted by read depth at each coordinate.
Default: 8

3. Minimum tags per CAGE cluster: the minimum number of CAGE tags in a cluster to be considered a potential transcription start site.
Default: 15

4. Minimum relative density for CAGE clusters: a measure of change in tag density between the cluster and its surroundings. Higher number = bigger change. For more information, see the Paraclu paper referenced
Default: 2

5. Minimum CAGE cluster length: the minimum cluster length, in base pairs, to be considered a potential transcription start site.
Default: 1

6. Maximum CAGE cluster length: the maximum cluster length, in base pairs, to be considered a potential transcription start site.
Default: 20

7. Maximum allowable distance between Iso-Seq and CAGE 5' starts: Maximum distance of a CAGE consensus start site from an Iso-Seq consensus start site to consider the start site validated.
Default: 3

8. Minimum number of SMRT reads to report a 5' start: the minimum number of Iso-Seq 5' starts in a cluster to be considered a potential transcription start site
Default: 1

9. Maximum distance in bp from an annotated start site to be called as "annotated"
Default: 10


OUTPUT

1. A BED file of validated 5' starts: this file is in bedDetail format, using the first 6 standard bed fields and an additional field that is necessary for the TRIMD_transcript_validator.pl script. The coordinates are those of the Iso-Seq consensus 5' start. The name field is in the format [nov|ann]_[+|-]_123.IsoSeq_456.CAGE and indicates whether the start site is novel or annotated, strand, the number of SMRT reads supporting it and the number of CAGE tags supporting it. The score is the number of SMRT reads and CAGE tags added together. The additional field indicates the range of the SMRT start cluster and number of supporting SMRT reads. Note that the validated starts file includes Iso-Seq starts that are supported by CAGE and/or annotation data. You may wish to filter the results to contain only starts that are supported by CAGE.

2. A text file with the number of total, novel and annotated 5' start sites validated, and a record of the input files.

3. BED file of Iso-Seq consensus 5' starts

4. BEDGRAPH file of all nonclipped Iso-Seq 5' starts

5. BED file of CAGE consensus 5' starts

6. BEDGRAPH file of all nonclipped CAGE 5' starts


## TRIMD_junction_matcher.pl
```
USAGE: perl /PATH/TRIMD_junction_matcher.pl </PATH/SMRT_introns_file> </PATH/Illumina_SJ.out.tab_file> </PATH/transcript_annotation_bed_file> <coordinates_to_ignore_bed_file(optional)>
```

Accepts a junctions files from GMAP/SMRT (generated with the -f introns argument) and an SJ.out.tab files from STAR/Illumina.
Returns 3 bed files: one of SMRT splice junctions, one of Illumina splice junctions and one of junctions detected by both methods. The coordinates in the output bed files correspond to the first and last bases of the introns. Optionally, splice junctions in repeat regions can be ignored by providing a bed file with the coordinates of those junctions.


INPUT

1. A file of Iso-Seq splice junctions data generated by GMAP (-f introns argument)

2. A file of Illumina splice junctions data generated by STAR (SJ.out.tab file in default STAR output).

3. BED file of annotated polyadenylated transcripts

4. (optional) BED file of genomic regions to ignore (e.g. repeat regions)


PARAMETERS

1. (Viral) chromosome name: the name of the chromosome under investigation. This must match between both junction files and the BED annotation file (field 1). Does not have to be viral, but for organisms with multiple chromosomes only one chromosome can be examined at a time.
No default: must be entered at prompt

2. Minimum SMRT read depth to report a splice junction

3. Minimum Illumina RNA-seq read depth to report a splice junction


OUTPUT

1. BED file of validated splice junctions: this file is in BED format, using the first 6 standard bed fields. The coordinates correspond to the first and last base of the excised intron. The name field is in the format [nov|ann]_123.IsoSeq_456.CAGE and indicates whether the junction is novel or annotated, the number of SMRT reads supporting it and the number of Illumina reads supporting it. The score is the number of SMRT reads and Illumina reads added together.

2. A text file with the number of total, novel and annotated splice junctions validated, and a record of the input files.

3. BED file of splice junctions detected on the specified chromosome in the Iso-Seq data

4. BED file of splice junctions detected on the specified chromosome in the Illumina data

## TRIMD_end_validator.pl
```
USAGE: perl /PATH/TRIMD_end_validator.pl </PATH/SMRT_sam_file> </PATH/Illumina_sam_file> </PATH/Annotation_bed_file>
```
Accepts a SAM file using Iso-Seq fl data, a SAM file using Illumina data, and a BED file of annotated polyadenylated transcripts. Counts the number of non-clipped SMRT reads with 3' ends at each genomic position and estimates consensus locations of clusters of 3' ends. Extracts Illumina reads containing apparent poly(A) tails and estimates consensus locations of clusters of polyadenylation sites. Output includes BEDGRAPH files of all 3' ends, BED files of the weighted centers of end clusters, a sam file of reads with polyA tails and a BED file of Iso-Seq 3' ends supported by either the annotation or the Illumina data.

SMRT fl read names must be formatted as putative_isoform_id/number_of_SMRT_reads/length.

INPUT

1. SAM file of Iso-Seq fl isoforms: this script was developed using data aligned with GMAP (-f samse option). Other aligners may also be appropriate.

2. SAM file of Illumina RNA-seq data: Illumina libraries should have been prepared with stranded TruSeq or a similar protocol. Sequence data can be paired-end or single-end. This script was developed using data aligned with STAR. Other aligners may also be appropriate.

3. BED file of annotated polyadenylated transcripts: the annotation file MUST be sorted by chrStart, then chrEnd. If your annotation file contains non-polyadenylated transcripts or other features (e.g. repeat regions, promoters) these should be removed to avoid false positives.

PARAMETERS

1. (Viral) chromosome name: the name of the chromosome under investigation. This must match between both SAM files (field 3) and the BED annotation file (field 1). Does not have to be viral, but for organisms with multiple chromosomes only one chromosome can be examined at a time.
No default: must be entered at prompt

2. Window for collapsing Iso-Seq 3' ends: Iso-Seq 3' ends within this number of bases of each other will be considered to represent the same polyadenlyation site. The consensus Iso-Seq polyadenylation site is determined by calculating an average of the coordinates in the cluster, weighted by read depth at each coordinate.
Default: 8

3. Minimum number of As for Illumina poly(A) tails: the number of As (or Ts, as appropriate) required in a read to indicate the presence of a poly(A) tail.
Default: 5

4. Minimum number of mismatches for Illumina poly(A) tails: the number of terminal mismatches in a read relative to the genome sequence to indicate the presence of a poly(A) tail.
Default: 2

5. Window for collapsing Illumina 3' ends: Illumina reads containing poly(A) tails within this number of bases of each other will be considered to represent the same polyadenylation site. The consensus Illumina polyadenylation site is determined by calculating an average of the coordinates in the cluster, weighted by read depth at each coordinate.
Default: 8

6. Number of bases downstream of Iso-Seq consensus 3' ends to look for Illumina support
Default: 10

7. Number of bases upstream of Iso-Seq consensus 3' ends to look for Illumina support
Default: 4

8. Minimum number of SMRT reads to report a 3' end
Default: 5

9. Minimum number of Illumina poly(A) tail reads to report a 3' end
Default: 1

10. Maximum distance in bp from an annotated end to be called as "annotated"
Default: 25

OUTPUT

1. BED file of validated 3' ends
This file is in bedDetail format, using the first 6 standard bed fields and an additional field that is necessary for the TRIMD_transcript_validator.pl script. The coordinates are those of the Illumina consensus 3' end. The name field is in the format [nov|ann]_[+|-]_123.IsoSeq_456.CAGE and indicates whether the end is novel or annotated, strand, the number of SMRT reads supporting it and the number of Illumina poly(A) reads supporting it. The score is the number of SMRT reads and Illumina poly(A) reads added together. The additional field indicates the range of the SMRT end cluster and number of supporting SMRT reads. Note that the validated starts file includes SMRT ends that are supported by Illumina and/or annotation data. You may wish to filter the results to contain only starts that are supported by Illumina.

2. A text file with the number of total, novel and annotated 3' end sites validated, and a record of the input files.

3. BED file of Iso-Seq consensus 3' ends

4. BEDGRAPH file of nonclipped Iso-Seq 3' ends

5. BED file of Illumina consensus 3' ends

6. BEDGRAPH file of nonclipped Illumina 3' ends

7. SAM file of Illumina reads containing putative poly(A) tails

## TRIMD_transcript_validator.pl
```
USAGE: perl /PATH/TRIMD_transcript_validator.pl </PATH/SMRT_sam_file> </PATH/validated_starts_file> </PATH/validated_ends_file> </PATH/validated_introns_file> </PATH/Annotation_bed_file>
```
Takes a SAM file of Iso-Seq fl isoforms and compares them to a list of validated 5' ends, 3' ends and introns to create a list of validated isoform structures, which are compared to an annotation file.

INPUT

1. SAM file of Iso-Seq fl isoforms

2. BedDetail file of validated starts: the output file from TRIMD_start_validator.pl that ends ".validated_starts.bed"

3. BedDetail file of validated ends: the output file from TRIMD_end_validator.pl that ends ".validated_ends.bed"

4. BedDetail file of validated splice junctions: the output file from TRIMD_junction_matcher.pl that ends ".validated_introns.bed"

5. BED file of annotated polyadenylated transcripts: the annotation file MUST be sorted by chrStart, then chrEnd. If your annotation file contains non-polyadenylated transcripts or other features (e.g. repeat regions, promoters) these should be removed to avoid false positives.

Note that the chromosome names must match exactly between all of the input files.

PARAMETERS

1. Maximum distance from an annotated 5' start to be called annotated
Default: 10

2. Maximum distance from an annotated 3' end to be called annotated
Default: 10

Note that if there are overlapping transcripts, more than one may be called as the same "annotated" transcript, depending on how the distance parameters are set.

OUTPUT

1. BED file of validated isoform structures: this file is in bed format, using all 12 standard fields. The coordinate for the 5' start is taken from the Iso-Seq consensus 5' start site and the coordinate for the 3' end is taken from the Illumina poly(A) read consensus 3' end. The name is that of one of the Iso-Seq isoforms representing that transcript, prefixed by any matching annotated transcript. The score is the number of SMRT reads that support that transcript. ThickStart and thickEnd (fields 7 and 8) match chrStart and chrEnd (fields 2 and 3): no information about ORFs is inferred. If the transcript is called as annotated, the color (field 9) is imported from the annotation file.

2. A text file with the number of total, novel and annotated isoforms validated, and a record of the input files.


## LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program as License.txt.  If not, see (http://www.gnu.org/licenses/).


## PLEASE CITE

O’Grady T, Wang X, Höner zu Bentrup K, Baddoo M, Concha M and Flemington EK (2016) Global transcript structure resolution of high gene density genomes through multi-platform data integration. Nucleic Acids Research  44(18), e145-e145.

and

Frith MC, Valen E, Krogh A, Hayashizaki Y, Carninci P, Sandelin A (2008) A code for transcription initiation in mammalian genomes. Genome Research 18(1):1-12.

## CONTACT

Erik Flemington (erik@tulane.edu) or Tina O'Grady (tmogrady@gmail.com)
