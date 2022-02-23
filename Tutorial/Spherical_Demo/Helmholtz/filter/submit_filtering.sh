#!/bin/bash
#SBATCH --output=sim-%j.log
#SBATCH --error=sim-%j.err
#SBATCH --time=00-00:30:00         # time (DD-HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --job-name="tutorial - Helmholtz : filter"

echo "${SLURM_NTASKS} MPI processors with ${SLURM_CPUS_PER_TASK} threads each"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

# Recall that eddy scales are [ 250e3, 750e3, 3500e3 ]

FILTER_SCALES="100e3 500e3 2000e3 5000e3"

mpirun -n ${SLURM_NTASKS} ./coarse_grain_helmholtz.x \
        --Helmholtz_input_file ../project/projection_ui.nc \
        --velocity_input_file ../velocity_sample.nc \
        --tor_field Psi \
        --pot_field Phi \
        --vel_field uo \
        --filter_scales "${FILTER_SCALES}"

#mpirun -n ${SLURM_NTASKS} ./coarse_grain_helmholtz.x \
#        --toroidal_input_file ../project/projection_ui.nc \
#        --potential_input_file ../project/projection_ui.nc \
#        --tor_field Psi \
#        --pot_field Phi \
#        --velocity_input_file ../velocity_sample.nc \
#        --uiuj_Helmholtz_input_file ../project/projection_uiuj.nc \
#        --uiuj_F_r v_r \
#        --uiuj_F_lon v_lon \
#        --uiuj_F_lat v_lat \
#        --vel_field uo \
#        --filter_scales "${FILTER_SCALES}"
