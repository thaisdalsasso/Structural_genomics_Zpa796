#!/usr/bin/env bash

ROOT="/path2/phylogeny_homologs"
WORKDIR="${ROOT}/Zpa796_copies"


# Prepare data input for HyPhy
python3.9 "${ROOT}/prepare_hyphy_inputs_Zpa796.py"

# Run phylogeny and selection models
for CLUSTER in "G.04" "G.12"; do
  CLUSTER_DIR="${WORKDIR}/${CLUSTER}"
  ALIGNMENT="${CLUSTER_DIR}/${CLUSTER}_Zpa796_mature_codon_alignment_trimmed_gappyout.fasta"
  TREE_PREFIX="${CLUSTER_DIR}/${CLUSTER}_Zpa796_mature_codon_iqtree"

  echo "Inferring nucleotide ML tree for ${CLUSTER}"
  iqtree2 \
    -s "${ALIGNMENT}" \
    -st DNA \
    -m TEST \
    --runs 10 \
    -pre "${TREE_PREFIX}"

  echo "Running HyPhy MEME for ${CLUSTER}"
  mkdir -p "${CLUSTER_DIR}/hyphy_MEME"
  hyphy meme \
    --alignment "${ALIGNMENT}" \
    --tree "${TREE_PREFIX}.treefile" \
    --branches All \
    --pvalue 0.05 \
    --output "${CLUSTER_DIR}/hyphy_MEME/${CLUSTER}_Zpa796_MEME.json" \
    2>&1 | tee "${CLUSTER_DIR}/hyphy_MEME/${CLUSTER}_Zpa796_MEME.log"

  echo "Running HyPhy FEL for ${CLUSTER}"
  mkdir -p "${CLUSTER_DIR}/hyphy_FEL"
  hyphy fel \
    --alignment "${ALIGNMENT}" \
    --tree "${TREE_PREFIX}.treefile" \
    --branches All \
    --pvalue 0.05 \
    --output "${CLUSTER_DIR}/hyphy_FEL/${CLUSTER}_Zpa796_FEL.json" \
    2>&1 | tee "${CLUSTER_DIR}/hyphy_FEL/${CLUSTER}_Zpa796_FEL.log"
done

echo "Done!"
