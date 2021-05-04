#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:40:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - Helmholtz - full"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

###
##      First, create the dataset
###

python generate_data_sphere.py

###
##      Then extract a coarsened grid
###

mpirun -np ${SLURM_NTASKS} ./coarsen_grid.x \
    --input_file ./velocity_sample.nc \
    --output_file ./coarsened_sample.nc \
    --Nlat_reduce_factor 4 \
    --Nlon_reduce_factor 4


###
##      Now do the Helmholtz projection on the coarse grid
###

mpirun -np ${SLURM_NTASKS} ./toroidal_projection.x \
    --input_file ./coarsened_sample.nc \
    --seed_file zero \
    --output_file coarse_toroidal_projection.nc \
    --max_iterations 50000 \
    --tolerance 1e-12

mpirun -np ${SLURM_NTASKS} ./potential_projection.x \
    --input_file ./coarsened_sample.nc \
    --seed_file zero \
    --output_file coarse_potential_projection.nc \
    --max_iterations 50000 \
    --tolerance 1e-12


###
##      Next, refine the coarse results to produce a seed for the fine grid
###

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


###
##      Finally, do the Helmholtz projection on the fine grid
###

mpirun -np ${SLURM_NTASKS} ./toroidal_projection.x \
    --input_file ./velocity_sample.nc \
    --seed_file ./toroidal_seed.nc \
    --output_file toroidal_projection.nc \
    --max_iterations 5000 \
    --tolerance 1e-8

mpirun -np ${SLURM_NTASKS} ./potential_projection.x \
    --input_file ./velocity_sample.nc \
    --seed_file ./potential_seed.nc \
    --output_file potential_projection.nc \
    --max_iterations 5000 \
    --tolerance 1e-8


###
##      Once that's all done, run the final analysis / plotting script 
###

python analyze_projection.py
