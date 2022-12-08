export OMP_NUM_THREADS=1
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

./Helmholtz_projection.x \
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

../coarse_grain_helmholtz.x \
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


KINDS=("full" "toroidal" "potential")

for KIND in "${KINDS[@]}"
do
    python ../../../PythonTools/OutputHandling/merge_postprocess_results.py \
        --file_pattern "outputs/postprocess_${KIND}_*nc" \
        --output_filename "RESULTS_${KIND}.nc"
done

##

echo "Process successfully completed!"
