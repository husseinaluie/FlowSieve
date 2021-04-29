#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:05:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - coarsen"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

mpirun -np ${SLURM_NTASKS} ./coarsen_grid.x \
    --input_file ./velocity_sample.nc \
    --output_file ./coarsened_sample.nc \
    --Nlat_reduce_factor 4 \
    --Nlon_reduce_factor 4
