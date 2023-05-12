#!/bin/bash
#SBATCH --job-name=array
#SBATCH -e logs/create_individual_features_%A_%a_.e
#SBATCH -o logs/create_individual_features_%A_%a_.o
#SBATCH -p normal	        # Queue (partition) name
#SBATCH -N 1                    # Total # of nodes (must be 1 for serial)
#SBATCH -n 8                    # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 02:30:00              # Run time (hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH -A #######	        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=########

cd /scratch/08789/rm57578/protein_interactions
source /work/08789/rm57578/miniconda3/etc/profile.d/conda.sh
conda activate alphapulldown2

create_individual_features.py \
  --fasta_paths=ATP8_protein.fa,GCNT1_protein_network.fa \
  --data_dir=/scratch/tacc/apps/bio/alphafold/data/ \
  --save_msa_files=False \
  --output_dir=/scratch/08789/rm57578/protein_interactions/features/ \
  --use_precomputed_msas=False \
  --max_template_date=2050-01-01 \
  --skip_existing=True \
  --seq_index=$SLURM_ARRAY_TASK_ID

