#!/bin/bash
#SBATCH --job-name=TMalign
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err

# Load modules
module load gcc12-env/12.3.0

# Define species name and directory paths
species="Zpa796"
ALPHAFOLD_BASE="/gxfs_work/cau/sunbo511/alphafold/${species}"
ESMFOLD_BASE="/gxfs_work/cau/sunbo511/esmfold/${species}"
PROTEIN_IDS_FILE="/gxfs_home/cau/sunbo511/data/References/mature_secretomes/Zpa796.mature_secretome.ids"
OUTPUT_BASE="./${species}"  # Base directory for output


# Iterate over each range 
# In this step, protein structures were copied and stored in folders according to their protein lenth categories
for range in "1-100" "101-200" "201-300" "301-400" "401-500" "501-600" "601-700" "701-800" "801-900" "901-1000" "1001_above"; do
    alpha_range_dir="${ALPHAFOLD_BASE}/${range}"
    esm_range_dir="${ESMFOLD_BASE}/${range}"
    output_range_dir="${OUTPUT_BASE}/${range}"  

    echo "Checking range: ${range}"
    echo "Alpha range dir: ${alpha_range_dir}"
    echo "ESM range dir: ${esm_range_dir}"

    mkdir -p "$output_range_dir"

    # Check if the range directories exist
    if [ -d "$alpha_range_dir" ] && [ -d "$esm_range_dir" ]; then
        while IFS= read -r protein_id; do
            alpha_file="${alpha_range_dir}/${protein_id}/ranked_0.pdb"
            esm_file="${esm_range_dir}/${protein_id}.pdb"

            if [ -f "$alpha_file" ] && [ -f "$esm_file" ]; then

        	# Perform TMalign and save the result in the output directory
                output_file="${output_range_dir}/${protein_id}.tmalign.out"
                matrix_file="${output_range_dir}/${protein_id}.TMmatrix.txt"
                supplementary_files="${output_range_dir}/${protein_id}.TMsup"

                ~/TMalign/TMalign "$alpha_file" "$esm_file" -m "$matrix_file" -o "$supplementary_files"> "$output_file"

		echo "TMalign completed for ${species}/${range}/${protein_id}"
            else
                echo "Error: ${protein_id} structure not predicted by Alphafold and/or ESMFold"
            fi
        done < "$PROTEIN_IDS_FILE"
    else
        echo "Error: One or both structure directories not found for ${species}/${range}"
    fi
done