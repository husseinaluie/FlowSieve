#!/bin/bash
#SBATCH --output=sim-%j.log
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:30:00                              # time (DD-HH:MM:SS)
#SBATCH --ntasks=1                                      # can only use one MPI process
#SBATCH --cpus-per-task=6                               # but we can use multiple openmp threads
#SBATCH --job-name="Coarse-grain : Postprocess Demo"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}  # SLURM_CPUS_PER_TASK is the number of Openmp threads per MPI process
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

FILTER_SCALES="5.e+04 7.53e+04 1.13e+05 1.71e+05 2.58e+05 3.88e+05 5.85e+05 8.81e+05 1.33e+06 2.00e+06"

# SLURM_NTASKS is the number of MPI processes
mpirun -n ${SLURM_NTASKS} ./coarse_grain.x              \
    --input_file ./velocity_sample.nc                   \
    --region_definitions_file ./region_definitions.nc   \
    --filter_scales "${FILTER_SCALES}"                  
