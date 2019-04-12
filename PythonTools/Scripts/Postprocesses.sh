#!/bin/sh
COL='\033[1;37m'
NC='\033[0m' # No Color

mpi=mpirun-openmpi-gcc8

echo -e "${COL}Preparing output directories${NC}"
python $PythonPlotScripts/prepare_directories.py

echo -e "${COL}Plotting KE${NC}"
${mpi} -n 2 python $PythonPlotScripts/plot_KE_dev.py
${mpi} -n 2 python $PythonPlotScripts/plot_KE_from_vels.py
echo -e "${COL}    done plotting KE${NC}"

#echo -e "${COL}Plotting various scatters${NC}"
#${mpi} -n 2 python $PythonPlotScripts/plot_scatters_hist.py
#echo -e "${COL}    done plotting various scatters${NC}"

echo -e "${COL}Plotting integrated fluxes${NC}"
${mpi} -n 2 python $PythonPlotScripts/plot_integrated_fluxes.py
${mpi} -n 2 python $PythonPlotScripts/compare_integrated_fluxes.py
echo -e "${COL}    done plotting KE / EN fluxes${NC}"
