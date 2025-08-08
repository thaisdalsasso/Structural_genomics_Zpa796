#!/bin/bash
#SBATCH --job-name=Zpa796_AF2_batch_submission
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --mem=10G
#SBATCH --time=48:00:00
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --partition=gpu


# Directory containing protein fasta files (one sequence per file)
fasta_dir="/gxfs_home/cau/sunbo511/data/References/mature_secretomes/Zpa796.mature_secretome-individual_fasta_files"

# Directory for output files
out_dir="./Zpa796_AF2"

# Directory to store submission scripts
submission_dir="./Zpa796_AF2_submission_scripts"

mkdir -p "$submission_dir"
mkdir -p "$out_dir"

# Loop through fasta and generate a submission script for each protein
for fasta_file in "$fasta_dir"/*.fasta; do
    protein_name=$(basename "$fasta_file" .fasta)
    submission_script="$submission_dir/${protein_name}_submission.sh"
    
    cat <<EOT > "$submission_script"
#!/bin/bash
#SBATCH --job-name=${protein_name}
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --mem=80G
#SBATCH --time=48:00:00
#SBATCH --output=job.%J.${protein_name}.out
#SBATCH --error=job.%J.${protein_name}.err
#SBATCH --partition=gpu
#SBATCH --constraint=V100

module load gpu-env/default
module load alphafold/2.3.1
source activate
conda activate alphafold2.3.1

export OMP_NUM_THREADS=8

run_alphafold.sh -d \$ALPHAFOLD_DATA -a \$CUDA_VISIBLE_DEVICES \
        -m monomer -t 2020-05-1  \
        -o "$out_dir" -f "$fasta_file"

conda deactivate
jobinfo
EOT

    # Submit script
    sbatch "$submission_script"
done
