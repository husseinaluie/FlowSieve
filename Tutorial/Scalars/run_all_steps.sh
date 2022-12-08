export OMP_NUM_THREADS=1
export KMP_AFFINITY="compact"
export I_MPI_PIN_DOMAIN="auto"

set -ex

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


##
echo "Process successfully completed!"
