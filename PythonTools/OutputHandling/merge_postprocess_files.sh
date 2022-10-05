
# This sample assumes that the postprocessing files are stored in
#   a directory called 'outputs' in the current working directory.
# Structure assume Helmholtz-type outputs.

KINDS=("full" "toroidal" "potential")

for KIND in "${KINDS[@]}"
do
    python merge_postprocess_results.py \
        --file_pattern "outputs/postprocess_${KIND}_*nc" \
        --output_filename "postprocess_${KIND}.nc" \
        --exclude_time_means \
        --exclude_OkuboWeiss \
        --exclude_zonal_means
done
