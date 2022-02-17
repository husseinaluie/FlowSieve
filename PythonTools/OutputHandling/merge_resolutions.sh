# Merge the full-field results
python merge_downsampling_results.py \
    --input_files postprocess_full_1.nc       postprocess_full_4.nc       postprocess_full_12.nc \
    --output_filename results_full.nc \
    --print_level 0

# Merge the toroidal results
python merge_downsampling_results.py \
    --input_files postprocess_toroidal_1.nc   postprocess_toroidal_4.nc   postprocess_toroidal_12.nc \
    --output_filename results_toroidal.nc  \
    --print_level 0

# Merge the potential results
python merge_downsampling_results.py \
    --input_files postprocess_potential_1.nc  postprocess_potential_4.nc  postprocess_potential_12.nc \
    --output_filename results_potential.nc  \
    --print_level 0
