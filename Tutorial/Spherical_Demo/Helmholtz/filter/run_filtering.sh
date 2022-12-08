##
## This section can take ~1 hour to run on one thread
##      so running one multiple threads (if available)
##      is recommended. Alternatively, the job can be
##      submitted to a SLURM schedule using the submit
##      script.
##

export OMP_NUM_THREADS=1

set -ex

# Recall that eddy scales are [ 250e3, 750e3, 3500e3 ]

FILTER_SCALES="100e3 500e3 2000e3 5000e3"

./coarse_grain_helmholtz.x \
    --Helmholtz_input_file ../project/projection_ui.nc \
    --velocity_input_file ../velocity_sample.nc \
    --tor_field Psi \
    --pot_field Phi \
    --vel_field uo \
    --filter_scales "${FILTER_SCALES}"


##
echo "Process successfully completed!"
