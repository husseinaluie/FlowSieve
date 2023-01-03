##
## This section can take ~4.5 hours to run on one thread
##      so running one multiple threads (if available)
##      is recommended. Alternatively, the job can be
##      submitted to a SLURM schedule using the submit
##      script.
##

export OMP_NUM_THREADS=1

set -ex

###
##      Define the geographic regions
###

python define_geographic_regions.py


###
##      Run the filtering routine
###

FILTER_SCALES="
1.e4 1.29e4 1.67e4 2.15e4 2.78e4 3.59e4 4.64e4 5.99e4 7.74e4 
1.e5 1.29e5 1.67e5 2.15e5 2.78e5 3.59e5 4.64e5 5.99e5 7.74e5 
1.e6 1.29e6 1.67e6 2.15e6 2.78e6 3.59e6 4.64e6 5.99e6 7.74e6 
1.e7" 

REGIONS_FILE="region_definitions.nc"

./coarse_grain_helmholtz.x \
    --Helmholtz_input_file ../project/projection_ui.nc \
    --velocity_input_file ../velocity_sample.nc \
    --tor_field Psi \
    --pot_field Phi \
    --vel_field uo \
    --region_definitions_file ${REGIONS_FILE} \
    --filter_scales "${FILTER_SCALES}"


###
##      Now merge the postprocess outputs into single files
###

python ../../../../PythonTools/OutputHandling/merge_postprocess_results.py \
    --file_pattern "postprocess_full_*nc" \
    --output_filename "RESULTS_full.nc"
python ../../../../PythonTools/OutputHandling/merge_postprocess_results.py \
    --file_pattern "postprocess_toroidal_*nc" \
    --output_filename "RESULTS_toroidal.nc"
python ../../../../PythonTools/OutputHandling/merge_postprocess_results.py \
    --file_pattern "postprocess_potential_*nc" \
    --output_filename "RESULTS_potential.nc"

##
echo "Process successfully completed!"
