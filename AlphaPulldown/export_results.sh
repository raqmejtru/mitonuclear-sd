#!/bin/bash
#SBATCH -J export_results
#SBATCH -e logs/export_results.e
#SBATCH -o logs/export_results.o
#SBATCH -p vm-small	        # Queue (partition) name
#SBATCH -N 1                    # Total # of nodes (must be 1 for serial)
#SBATCH -n 1                    # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 00:45:00             # Run time (hh:mm:ss)
#SBATCH --mail-type=all         # Send email at begin and end of job
#SBATCH -A #######	        # Project/Allocation name (req'd if you have more than 1)
#SBATCH --mail-user=########

source /work/08789/rm57578/miniconda3/etc/profile.d/conda.sh
conda activate alphapulldown2

module load tacc-singularity
singularity exec \
    --no-home \
    --bind /scratch/08789/rm57578/protein_interactions/models:/mnt \
    /scratch/08789/rm57578/protein_interactions/alpha-analysis.sif \
    run_get_good_pae.sh \
    --output_dir=/mnt \
    --cutoff=500
