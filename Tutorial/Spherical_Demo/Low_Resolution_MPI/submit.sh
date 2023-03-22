#!/bin/bash
#SBATCH --output=sim-%j.log 
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:15:00         # time (DD-HH:MM:SS)
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - MPI"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

set -ex

###
##      First, create the dataset and regions file
###

python generate_data_sphere.py
python define_geographic_regions.py

###
##      Do the Helmholtz projection on the fine grid
###

mpirun -np ${SLURM_NTASKS} ./Helmholtz_projection.x \
    --input_file ./velocity_sample.nc \
    --seed_file zero \
    --output_file projection_ui.nc \
    --max_iterations 1500 \
    --tolerance 1e-18 \
    --Tikhov_Laplace 1

###
##      Run coarse-graining routine
###

FILTER_SCALES="
1.e5 1.29e5 1.67e5 2.15e5 2.78e5 3.59e5 4.64e5 5.99e5 7.74e5 
1.e6 1.29e6 1.67e6 2.15e6 2.78e6 3.59e6 4.64e6 5.99e6 7.74e6 
1.e7 1.29e7 1.67e7 2.15e7 2.78e7 3.59e7 4.64e7 5.99e7 7.74e7
"

mkdir -p outputs
cd outputs

mpirun -np ${SLURM_NTASKS} ../coarse_grain_helmholtz.x \
    --Helmholtz_input_file ../projection_ui.nc \
    --velocity_input_file ../velocity_sample.nc \
    --tor_field Psi \
    --pot_field Phi \
    --vel_field uo \
    --region_definitions_file ../region_definitions.nc \
    --filter_scales "${FILTER_SCALES}"

cd ..

###
##      Now merge the postprocess outputs into single files
###

python ../../../PythonTools/OutputHandling/merge_postprocess_results.py \
        --file_pattern "outputs/postprocess_full_*nc" \
        --output_filename "RESULTS_full.nc"
python ../../../PythonTools/OutputHandling/merge_postprocess_results.py \
        --file_pattern "outputs/postprocess_toroidal_*nc" \
        --output_filename "RESULTS_toroidal.nc"
python ../../../PythonTools/OutputHandling/merge_postprocess_results.py \
        --file_pattern "outputs/postprocess_potential_*nc" \
        --output_filename "RESULTS_potential.nc"

##

echo "Process successfully completed!"
