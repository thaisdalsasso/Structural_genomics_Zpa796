#!/usr/bin/env bash

root_dir="/path2/phylogeny_homologs"
work_dir="${root_dir}/Zpa796_copies"


# Prepare data input for HyPhy
python3.9 "${root_dir}/prepare_hyphy_inputs_Zpa796.py"

# Run phylogeny and selection models
for cluster in "G.04" "G.12"; do
  cluster_dir="${work_dir}/${cluster}"
  alignment="${cluster_dir}/${cluster}_Zpa796_mature_codon_alignment_trimmed_gappyout.fasta"
  tree_prefix="${cluster_dir}/${cluster}_Zpa796_mature_codon_iqtree"

  echo "Inferring nucleotide ML tree for ${cluster}"
  iqtree2 \
    -s "${alignment}" \
    -st DNA \
    -m TEST \
    --runs 10 \
    -pre "${tree_prefix}"

  echo "Running HyPhy MEME for ${cluster}"
  mkdir -p "${cluster_dir}/hyphy_MEME"
  hyphy meme \
    --alignment "${alignment}" \
    --tree "${tree_prefix}.treefile" \
    --branches All \
    --pvalue 0.05 \
    --output "${cluster_dir}/hyphy_MEME/${cluster}_Zpa796_MEME.json" \
    2>&1 | tee "${cluster_dir}/hyphy_MEME/${cluster}_Zpa796_MEME.log"

  echo "Running HyPhy FEL for ${cluster}"
  mkdir -p "${cluster_dir}/hyphy_FEL"
  hyphy fel \
    --alignment "${alignment}" \
    --tree "${tree_prefix}.treefile" \
    --branches All \
    --pvalue 0.05 \
    --output "${cluster_dir}/hyphy_FEL/${cluster}_Zpa796_FEL.json" \
    2>&1 | tee "${cluster_dir}/hyphy_FEL/${cluster}_Zpa796_FEL.log"
done

echo "Done!"
