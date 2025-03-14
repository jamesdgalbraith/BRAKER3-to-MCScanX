# BRAKER3 to MCScanX
A collection of scripts to: 
1) Transform BRAKER3 outputs into MCScanX format and compile data
2) Run the blastp search needed for MCScanX
3) Run MCScanX

Has been successfully run on both on Linux and macOS

### Requirements
NCBI BLAST+, MCScanX (https://github.com/wyp1125/MCScanX)

__R packages__: tidyverse, Biostrings, plyranges, optparse

### Scripts
__BRAKER3-to-MCScanX.R__
Used to convert BRAKER3 gff and fastas into MCScanX approriate format. Requires BRAKER3 `.aa`, `.codingseq` and `.gff3` along with a two letter prefix for the genome. The longest transcript of each gene will be used as its representative sequence.

__BRAKER3-to-MCScanX.sh__
Takes BRAKER3 outputs, uses above Rscript to reformat for MCScanX, runs BLASTP searches, compiles data and runs MCScanX.
As written it requires a tsv file `species_list.tsv` containing the species name in column 1 and desired 2 letter abbreviation in column 2; and BRAKER3 output to be in the directory `data/braker/`, with the file names corresponding to the species name in the `species_list.tsv`.
