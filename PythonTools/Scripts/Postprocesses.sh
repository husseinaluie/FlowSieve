COL='\033[1;37m'
NC='\033[0m' # No Color

# Plot KE
echo -e "${COL}Plotting velocities${NC}"
python $PythonPlotScripts/plot_vels.py
echo -e "${COL}    done plotting velocities${NC}"

echo -e "${COL}Plotting vorticities${NC}"
python $PythonPlotScripts/plot_vorts.py
echo -e "${COL}    done plotting vorticities${NC}"

echo -e "${COL}Plotting KE${NC}"
python $PythonPlotScripts/plot_KE.py
echo -e "${COL}    done plotting KE${NC}"

echo -e "${COL}Plotting spectra${NC}"
python $PythonPlotScripts/plot_spectra.py
echo -e "${COL}    done plotting spectra${NC}"

echo -e "${COL}Plotting energy transfers${NC}"
python $PythonPlotScripts/plot_transfers.py
echo -e "${COL}    done plotting energy transfers${NC}"
