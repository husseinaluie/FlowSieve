#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:20:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - coarse projection"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

mpirun -np ${SLURM_NTASKS} ./toroidal_projection.x \
    --input_file ./coarsened_sample.nc \
    --seed_file zero \
    --output_file coarse_toroidal_projection.nc \
    --max_iterations 500000 \
    --tolerance 1e-12

mpirun -np ${SLURM_NTASKS} ./potential_projection.x \
    --input_file ./coarsened_sample.nc \
    --seed_file zero \
    --output_file coarse_potential_projection.nc \
    --max_iterations 500000 \
    --tolerance 1e-12
