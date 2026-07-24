#!/usr/bin/env bash

root_dir="./phylogeny_homologs"
input_multifasta="${root_dir}/G.12/G.12_homologs_Zymospp.fasta"
outdir="${root_dir}/phylogeny_Zymospp/G.12"

mkdir -p "${outdir}"

alignment="${outdir}/G.12_Zymospp_proteins_mafft_linsi.fasta"
trimmed_alignment="${outdir}/G.12_Zymospp_homologs_mafft_linsi_trimal_gt0.9_cons60_w3.fasta"
iqtree_prefix="${outdir}/G.12_Zymospp_homologs_iqtree"


# Align sequences
mafft-linsi "${input_multifasta}" > "${alignment}"

# Trim alignment
trimal -in "${alignment}" -out "${trimmed_alignment}" -gt 0.9 -cons 60 -w 3

# Run ML phylogeny
iqtree2 -s "${trimmed_alignment}" -st AA --runs 10 --ufboot 1000 -m TEST -pre "${iqtree_prefix}"

echo "Done!"
