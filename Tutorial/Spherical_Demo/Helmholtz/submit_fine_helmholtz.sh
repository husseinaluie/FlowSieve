#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-01:00:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - fine projection"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

mpirun -np ${SLURM_NTASKS} ./toroidal_projection.x \
    --input_file ./velocity_sample.nc \
    --seed_file ./toroidal_seed.nc \
    --output_file toroidal_projection.nc \
    --max_iterations 50000 \
    --tolerance 1e-8

mpirun -np ${SLURM_NTASKS} ./potential_projection.x \
    --input_file ./velocity_sample.nc \
    --seed_file ./potential_seed.nc \
    --output_file potential_projection.nc \
    --max_iterations 50000 \
    --tolerance 1e-8
