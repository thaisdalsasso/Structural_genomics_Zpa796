#!/bin/bash
#SBATCH --job-name=TMalign
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

# Load modules
module load gcc12-env/12.3.0

# Directory containing PDB files for all-vs-all alignment
pdb_dir="/gxfs_work/cau/sunbo511/alphafold/Zpa796_AF2bestmodels"

# Directory for outputs
work_dir="/gxfs_work/cau/sunbo511/tmalign/Zpa796_all-vs-all"

mkdir -p "$work_dir"

# Iterate over all PDB files
for protein1 in "$pdb_dir"/*.pdb; do
  for protein2 in "$pdb_dir"/*.pdb; do
    if [ "$protein1" != "$protein2" ]; then
      protein1_id=$(basename "$protein1" .pdb)
      protein2_id=$(basename "$protein2" .pdb)

      # Set output files and run TMalign
      output_file="${protein1_id}_vs_${protein2_id}.tmalign.out"

      ~/TMalign/TMalign "$protein1" "$protein2" > "${work_dir}/$output_file"

    fi
  done
done

