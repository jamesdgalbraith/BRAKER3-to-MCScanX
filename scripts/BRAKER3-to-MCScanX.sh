#!/bin/bash

set -e

# Make directory for MCScanX formatted data
mkdir -p data/MCScanX/raw/

### Part 1 - Reformat BRAKER3 ###
# Convert BRAKER3 output for each species into format necessary for MCScanX
# At this stage the output of BRAKER3 for each species is in a subdirectory of data/braker
# and has had the species name added to the names of the BRAKER3 files
# The output will be the amino acid, nucleotide and MCScanX gff of the longest transcript
# of each gene correctly prefixed

while read SPECIES ABBREV; do
    echo $SPECIES
    Rscript scripts/BRAKER3-to-MCScanX.R \
        -a data/braker/${SPECIES}.braker.aa \
        -c data/braker/${SPECIES}.braker.codingseq \
        -g data/braker/${SPECIES}.braker.gff3 \
        -p ${ABBREV} \
        -o data/MCScanX/raw/
done < "species_list.tsv"

# Compile data
cat data/MCScanX/raw/*.aa > data/MCScanX/all.aa
cat data/MCScanX/raw/*.cds > data/MCScanX/all.cds
cat data/MCScanX/raw/*.gff > data/MCScanX/all.gff

### Part 2 - Run blastp ###
# For each species perform blastp of braker.aa against all other species braker.aa
# and compile into single file
mkdir -p data/MCScanX/blast
while read SPECIES_A ABBREV_A; do
    makeblastdb -in data/MCScanX/raw/MCScanX_${SPECIES_A}.braker.aa -dbtype prot
    echo $SPECIES_A
    while read SPECIES_B ABBREV_B; do
        if [[ ! ${SPECIES_A} = ${SPECIES_B} ]]; then
            echo "\t vs " $SPECIES_B
            blastp -db data/MCScanX/raw/MCScanX_${SPECIES_A}.braker.aa \
                -query data/MCScanX/raw/MCScanX_${SPECIES_B}.braker.aa \
                -num_threads 8 -evalue 1e-10 -outfmt 6 -max_target_seqs 5 -seg yes \
                -out data/MCScanX/blast/${SPECIES_B}_in_${SPECIES_A}.blast 
        fi
    done < "species_list.tsv"
    rm data/MCScanX/raw/MCScanX_${SPECIES_A}.braker.aa.p*
done < "species_list.tsv"

# Compile blast data for all vs all
cat data/MCScanX/blast/*.blast > data/MCScanX/all.blast

### PART 3 - Run MCScanX ###
# Downloaded from https://github.com/wyp1125/MCScanX
DATASET=data/MCScanX/all
MCScanX -b 2 ${DATASET}
