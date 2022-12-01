#!/bin/bash
#SBATCH --output=sim-%j.log
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:15:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="test"

set -ex

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

##
## Generate sample data
##
python generate_data_sphere.py

##
##  Apply coarse-graining
##
./coarse_grain_scalars.x \
    --input_file ./density_sample.nc \
    --filter_scales "100e3 500e3 2000e3" \
    --variables "rho"

