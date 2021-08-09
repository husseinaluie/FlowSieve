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

FILTER_SCALES="10e3 20e3 30e3 4.03e+04 5.33e+04 7.04e+04 9.31e+04 1.23e+05 1.63e+05 2.15e+05 2.84e+05 3.75e+05 4.96e+05 6.56e+05 8.66e+05 1.15e+06 1.51e+06 2.00e+06"

mpirun -n ${SLURM_NTASKS} ./coarse_grain_helmholtz.x \
        --toroidal_input_file ../project/projection_ui.nc \
        --potential_input_file ../project/projection_ui.nc \
        --tor_field Psi \
        --pot_field Phi \
        --velocity_input_file ../velocity_sample.nc \
        --uiuj_Helmholtz_input_file ../project/projection_uiuj.nc \
        --uiuj_F_r v_r \
        --uiuj_F_lon v_lon \
        --uiuj_F_lat v_lat \
        --vel_field uo \
        --filter_scales "${FILTER_SCALES}"
