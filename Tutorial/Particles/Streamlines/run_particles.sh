set -ex

FINAL_TIME="26280" # three years in hours [unit of 'time_unit']
OUTPUT_FREQUENCY="28800" # eight hours in seconds [always in seconds]

NUM_PARTICLES="1000"  # this is the number of particles *per MPI rank*

./particles.x --zonal_vel "uo" \
              --merid_vel "vo" \
              --input_file "velocity_sample.nc" \
              --output_file "particles.nc" \
              --time_unit "hours" \
              --particle_per_mpi ${NUM_PARTICLES} \
              --output_frequency ${OUTPUT_FREQUENCY} \
              --final_time ${FINAL_TIME} \
              --particle_lifespan -1
