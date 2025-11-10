#!/bin/bash
#SBATCH --output=sim-%j.log
#SBATCH --error=sim-%j.err
#SBATCH --time=00-03:00:00         # time (DD-HH:MM:SS)
#SBATCH -p spr
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --job-name="tutorial - Streamlines"

set -ex

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

FINAL_TIME="26280" # three years in hours [unit of 'time_unit']
OUTPUT_FREQUENCY="28800" # eight hours in seconds [always in seconds]

NUM_PARTICLES=1000  # this is the number of particles *per MPI rank*

./particles.x --zonal_vel "uo" \
              --merid_vel "vo" \
              --input_file "velocity_sample.nc" \
              --output_file "particles.nc" \
              --time_unit "hours" \
              --particle_per_mpi ${NUM_PARTICLES} \
              --output_frequency ${OUTPUT_FREQUENCY} \
              --final_time ${FINAL_TIME} \
              --particle_lifespan -1
