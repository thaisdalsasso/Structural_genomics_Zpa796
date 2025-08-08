#!/bin/bash

# Input file
input_file="G.12_mature.fasta"

baseline=$(basename "$input_file" .fasta)

# Run MAFFT
mafft-linsi "$input_file" > "${baseline}_mafft-linsi.fasta"

# Run TrimAl
trimal -in "${baseline}_mafft-linsi.fasta" -out "${baseline}_mafft-linsi_trimmed-gappyout.fasta" -gappyout
