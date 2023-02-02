set -ex

##
## Generate our sample data
##
python generate_data.py


##
## Run coarse-graining
##

mpirun -n 1 ./coarse_grain.x --input_file velocity_sample.nc --filter_scales "1e3 15e3 50e3 100e3"

##
## Make figures from the outputs
##

python process_results.py

##
echo "Process successfully completed!"
