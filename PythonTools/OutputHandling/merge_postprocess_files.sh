
# This sample assumes that the postprocessing files are stored in
#   a directory called 'outputs' in the current working directory.
# Structure assume Helmholtz-type outputs.

python merge_postprocess_results.py --file_pattern "outputs/postprocess_full_*nc" --output_filename "postprocess_full.nc"
python merge_postprocess_results.py --file_pattern "outputs/postprocess_toroidal_*nc" --output_filename "postprocess_toroidal.nc"
python merge_postprocess_results.py --file_pattern "outputs/postprocess_potential_*nc" --output_filename "postprocess_potential.nc"
