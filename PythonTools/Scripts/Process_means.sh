#!/bin/sh
COL='\033[1;37m'
NC='\033[0m' # No Color

mpi=mpirun-openmpi-gcc8

echo -e "${COL}Preparing output directories${NC}"
python $PythonPlotScripts/prepare_directories.py

echo -e "${COL}Plotting mean velocities${NC}"
${mpi} -n 2 python $PythonPlotScripts/Means/plot_vels.py
echo -e "${COL}    done plotting velocities${NC}"

echo -e "${COL}Plotting mean vorticities${NC}"
${mpi} -n 2 python $PythonPlotScripts/Means/plot_vorts.py
echo -e "${COL}    done plotting vorticities${NC}"

echo -e "${COL}Plotting mean KE${NC}"
${mpi} -n 2 python $PythonPlotScripts/Means/plot_KE_dev.py
${mpi} -n 2 python $PythonPlotScripts/Means/plot_KE_from_vels.py
echo -e "${COL}    done plotting KE${NC}"

echo -e "${COL}Plotting mean EN${NC}"
${mpi} -n 2 python $PythonPlotScripts/Means/plot_EN_from_vorts.py
echo -e "${COL}    done plotting EN${NC}"

echo -e "${COL}Plotting mean energy transfers${NC}"
${mpi} -n 2 python $PythonPlotScripts/Means/plot_transfers.py
echo -e "${COL}    done plotting energy transfers${NC}"
