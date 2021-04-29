#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:10:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - refine seed"
#SBATCH --dependency=afterok:43670706

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

mpirun -np ${SLURM_NTASKS} ./refine_Helmholtz_seed.x \
    --coarse_file ./coarse_toroidal_projection.nc \
    --fine_file ./velocity_sample.nc \
    --output_file toroidal_seed.nc \
    --var_in_coarse F \
    --var_in_output seed

mpirun -np ${SLURM_NTASKS} ./refine_Helmholtz_seed.x \
    --coarse_file ./coarse_potential_projection.nc \
    --fine_file ./velocity_sample.nc \
    --output_file potential_seed.nc \
    --var_in_coarse F \
    --var_in_output seed
