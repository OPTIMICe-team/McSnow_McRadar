# McSnow_McRadar
this runs McSnow and McRadar, aswell as some plotting routines. The main routine that controls everything is McSnow_McRadar.sh. 
Run with bash McSnow_McRadar.sh 1 1 1 1 1 in order to recompile McSnow, run McSnow with the specified setup, convert McSnow outputfile from .dat to .nc, 
run McRadar, plot output

The McSnow runskript is in McSnow_runskripts/1d_habit

The routine needed to plot McSnow output is plot_output.py, which uses the functions in plotRoutines.py

The version of plot_output.py in BIMOD can run without McRadar. Therefore, if McRadar is not installed set that to zero when running (so run bash McSnow_McRadar.sh 1 1 1 0 1 to recompile McSnow, run McSnow, convert mass2fr.dat to mass2fr.nc and then plot the properties)
