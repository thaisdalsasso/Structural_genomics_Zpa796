#!/bin/bash
#SBATCH --job-name=Zpa796_esmfold
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --qos=long
#SBATCH --time=10-00:00:00
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out


#Load modules:
module load gcc12-env/12.3.0
module load miniconda3/23.5.2

#Activate environment
conda activate  my_esmfold_env

# Directory containing protein fasta files (one sequence per file)
fasta_dir="/gxfs_home/cau/sunbo511/data/References/mature_secretomes/Zpa796.mature_secretome-individual_fasta_files"

# Directory for output files
output_dir="./Zpa796_ESMFold"


mkdir -p ./"$output_dir"

for fasta_file in "$fasta_dir"/*.fasta; do
    esm-fold -i "$fasta_file" -o "$output_dir"/ --cpu-only --chunk-size 32
done


conda deactivate
jobinfo