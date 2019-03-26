#!/bin/sh
COL='\033[1;37m'
NC='\033[0m' # No Color

echo -e "${COL}Preparing output directories${NC}"
python $PythonPlotScripts/prepare_directories.py

#echo -e "${COL}Plotting original fields${NC}"
#python $PythonPlotScripts/plot_originals.py
#echo -e "${COL}    done plotting original fields${NC}"

echo -e "${COL}Plotting velocities${NC}"
python $PythonPlotScripts/plot_vels.py
echo -e "${COL}    done plotting velocities${NC}"

echo -e "${COL}Plotting vorticities${NC}"
python $PythonPlotScripts/plot_vorts.py
echo -e "${COL}    done plotting vorticities${NC}"

echo -e "${COL}Plotting KE${NC}"
python $PythonPlotScripts/plot_KE_dev.py
python $PythonPlotScripts/plot_KE_from_vels.py
echo -e "${COL}    done plotting KE${NC}"

#echo -e "${COL}Plotting spectra${NC}"
#python $PythonPlotScripts/plot_spectra.py
#echo -e "${COL}    done plotting spectra${NC}"

echo -e "${COL}Plotting energy transfers${NC}"
python $PythonPlotScripts/plot_transfers.py
echo -e "${COL}    done plotting energy transfers${NC}"

echo -e "${COL}Plotting various scatters${NC}"
#python $PythonPlotScripts/plot_scatters.py
python $PythonPlotScripts/plot_scatters_hist.py
echo -e "${COL}    done plotting various scatters${NC}"

echo -e "${COL}Plotting KE / EN fluxes${NC}"
python $PythonPlotScripts/plot_KE_vel_flux.py
python $PythonPlotScripts/plot_EN_vort_flux.py
echo -e "${COL}    done plotting KE / EN fluxes${NC}"
