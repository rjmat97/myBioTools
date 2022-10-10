# myBioTools
some basic scripts that can be used to simplify repetitive biological data handling

[compFastas](#1compfastas)<br/>
[vcfAnotator](#2vcfAnotator)<br/>
[motifFinder](#3motifFinder)<br/>

## 1.compFastas
    This script can be used to compare two multi FASTA files and determine the common entries based on the sequence rather than the ID itself.<br/>
    The script can also be used to compare multi FASTA files to JSON format.
    use -h when running the scriot to get the necessary documentation.

## 2.vcfAnotator
    This script can be used to annotate VCF files to identify which genes contain the mutations.
    Input files: vcf file and gtf file

## 3.motifFinder
    MotifFinder can be used to generate a bedfile conatining positions where a given motif is present. It was created with the intention to identify endocuvlease cut sites, hence by default it generates the a bed file with atack sites.
