___
## **AlphaPulldown**
*Install and Run on TACC HPC*  
April 2023
___


### **Description**
In `pulldown` mode, AlphaPulldown screens a list of "bait" proteins against a list or lists of other proteins. 
It then returns models of all interactions, summarizes results in a score table, and provides PAE plots and 3D displays of models colored by chain and pLDDT score.

Note: current AlphaPulldown version not compatible with AlphaFold v 2.2.0 on TACC. Instead, use older version: `alphapulldown==0.22.3`. 



### **Install** 
```bash
mamba create -n alphapulldown2
mamba activate alphapulldown2
mamba install -c omnia -c bioconda -c conda-forge -c anaconda \
    -c "nvidia/label/cuda-11.3.0" \
	python==3.7 openmm=7.5.1 pdbfixer kalign2=2.04 cctbx-base \
    hmmer hhsuite threadpoolctl \
	docopt xlrd==1.2.0 xlsxwriter cuda-nvcc
python3 -m pip install alphapulldown==0.22.3
pip install -q "jax[cuda]>=0.3.8,<0.3.10"  \
    -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html

# Download singularity container to summarize prediction scores. 
cd /scratch/08789/rm57578/forf_interactions
wget https://www.embl-hamburg.de/AlphaPulldown/downloads/alpha-analysis.sif
```


### **Run** 
```bash
cd /scratch/08789/rm57578/forf_interactions

# array count: total # of proteins (bait + targets)
sbatch --array=1-19 create_individual_features.sh

# array count: total # of interactions (bait * targets)
sbatch --array=1-18 predict_interaction.sh
```


### **Summarize results** 
AlphaPulldown generates a CSV table with structural properties and scores.
```bash
sbatch export_results.sh
```


### **Interpreting models**
DockQ score cutoffs:

From https://github.com/bjornwallner/DockQ:
```
***********************************************************
*                       DockQ                             *
*   Scoring function for protein-protein docking models   *
*   Statistics on CAPRI data:                             *
*    0    <  DockQ <  0.23 - Incorrect                    *
*    0.23 <= DockQ <  0.49 - Acceptable quality           *
*    0.49 <= DockQ <  0.80 - Medium quality               *
*            DockQ >= 0.80 - High quality                 *
*   Reference: Sankar Basu and Bjorn Wallner, DockQ:...   *
*   For comments, please email: bjornw@ifm.liu.se         *
***********************************************************
```

From https://doi.org/10.1038/s41467-021-23692-x:  
Additionally, PI_scores > 1 indicate good likelihood of protein interface. 

### **Downstream Visualization**
In ChimeraX, use following command to color residues by confidence:
``` 
color bfactor #* palette alphafold
```
