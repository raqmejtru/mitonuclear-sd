#!/bin/bash
#SBATCH --job-name=array
#SBATCH -e logs/atp8_predict_interaction_%A_%a_.e
#SBATCH -o logs/atp8_predict_interaction_%A_%a_.o
#SBATCH -p gpu-a100	        # Queue (partition) name
#SBATCH -N 1                    # Total # of nodes (must be 1 for serial)
#SBATCH -n 8                    # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 04:00:00             # Run time (hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH -A #######	        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=########

cd /scratch/08789/rm57578/protein_interactions
source /work/08789/rm57578/miniconda3/etc/profile.d/conda.sh
conda activate alphapulldown2

run_multimer_jobs.py --mode=pulldown \
    --num_cycle=3 \
    --num_predictions_per_model=1 \
    --output_path=/scratch/08789/rm57578/protein_interactions/models/ \
    --data_dir=/scratch/tacc/apps/bio/alphafold/data/ \
    --protein_lists=ATP8_protein.txt,GCNT1_protein_network.txt \
    --monomer_objects_dir=/scratch/08789/rm57578/protein_interactions/features/ \
    --job_index=$SLURM_ARRAY_TASK_ID