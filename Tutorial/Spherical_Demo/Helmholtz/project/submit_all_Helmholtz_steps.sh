#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-02:45:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - Helmholtz : project"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

###
##      First, create the dataset
###

cd ..
python generate_data_sphere.py
cd project/

###
##      Then extract the data a coarsened grid
###

#python generate_coarse_grid.py

echo "\n\n COARSEN GRID \n\n"
mpirun -np ${SLURM_NTASKS} ./coarsen_grid_linear.x \
    --fine_file ../velocity_sample.nc \
    --coarse_file coarse_grid.nc \
    --output_file coarsened_sample.nc \
    --input_variables "uo vo" \
    --output_variables "uo vo"


###
##      Now do the Helmholtz projection on the coarse grid
###

echo "\n\n COARSE HELMHOLTZ Ui\n\n"
mpirun -np ${SLURM_NTASKS} ./Helmholtz_projection.x \
    --input_file ./coarsened_sample.nc \
    --seed_file zero \
    --output_file coarse_projection_ui.nc \
    --max_iterations 50000 \
    --tolerance 1e-12 \
    --Tikhov_Laplace 1




###
##      Next, refine the coarse results to produce a seed for the fine grid
###

echo "\n\n REFINE SEED FOR Ui \n\n"
mpirun -np ${SLURM_NTASKS} ./refine_Helmholtz_seed.x \
    --coarse_file ./coarse_projection_ui.nc \
    --fine_file ../velocity_sample.nc \
    --output_file seed_ui.nc \
    --input_variables "Phi Psi" \
    --output_variables "Phi_seed Psi_seed"




###
##      Finally, do the Helmholtz projection on the fine grid
###

echo "\n\n FINE HELMHOLTZ Ui\n\n"
mpirun -np ${SLURM_NTASKS} ./Helmholtz_projection.x \
    --input_file ../velocity_sample.nc \
    --seed_file ./seed_ui.nc \
    --output_file projection_ui.nc \
    --max_iterations 15000 \
    --tolerance 1e-18 \
    --Tikhov_Laplace 1



###
##      Once that's all done, run the final analysis / plotting script 
###

python analyze_projection.py
