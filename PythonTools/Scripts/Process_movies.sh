#!/bin/sh
COL='\033[1;37m'
NC='\033[0m' # No Color

mpi=mpirun-openmpi-gcc8

echo -e "${COL}Plotting original fields${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_originals.py
echo -e "${COL}    done plotting original fields${NC}"

echo -e "${COL}Plotting velocities${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_vels.py
echo -e "${COL}    done plotting velocities${NC}"

echo -e "${COL}Plotting vorticities${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_vorts.py
echo -e "${COL}    done plotting vorticities${NC}"

echo -e "${COL}Plotting KE${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_KE_dev.py
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_KE_from_vels.py
echo -e "${COL}    done plotting KE${NC}"

echo -e "${COL}Plotting EN${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_EN_from_vorts.py
echo -e "${COL}    done plotting EN${NC}"

echo -e "${COL}Plotting energy transfers${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_transfers.py
echo -e "${COL}    done plotting energy transfers${NC}"

echo -e "${COL}Plotting KE / EN fluxes${NC}"
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_KE_vel_flux.py
${mpi} -n 2 python $PythonPlotScripts/Movies/plot_EN_vort_flux.py
echo -e "${COL}    done plotting KE / EN fluxes${NC}"
