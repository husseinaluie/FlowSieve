FINAL_TIME="8760" # one year in hours [unit of 'time_unit']
OUTPUT_FREQUENCY="28800" # eight hours in seconds [always in seconds]

NUM_PARTICLES=1000  # this is the number of particles *per MPI rank*

./particles.x --zonal_vel "uo" \
              --merid_vel "vo" \
              --input_file "velocity_sample.nc" \
              --output_file "particles.nc" \
              --time_unit "hours" \
              --particles_per_mpi ${NUM_PARTICLE} \
              --output_frequency ${OUTPUT_FREQUENCE} \
              --final_time ${FINAL_TIME} \
              --particle_lifespan -1
