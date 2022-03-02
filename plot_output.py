import mcradar as mcr
from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt

freq = np.array([float(os.environ['freq'])*1e9])
print(freq)
experimentID = os.environ['experiment']
#print(experimentID)
inputPath = os.environ['MCexp']+'experiments/'+experimentID+'/'
scatMode = os.environ['scatMode']
print(inputPath)

print('loading the settings')
splitPath = inputPath.split('domtop')[1]
domTop = splitPath[0:4]
lutPath = '/work/lvonterz/SSRGA/snowScatt/ssrga_LUT/'
#-- load the settings of McSnow domain, as well as elevation you want to plot:
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',
                               elv=30, freq=freq,gridBaseArea=5.0,maxHeight=int(domTop),
                               ndgsVal=50,heightRes=50,scatSet={'mode':scatMode, 'lutPath':lutPath,'particle_name':'vonTerzi_dendrite','safeTmatrix':False})

print('loading the McSnow output')
# now generate a table from the McSnow output. You can specify xi0, if it is not stored in the table (like in my cases)
mcTable = mcr.getMcSnowTable(dicSettings['dataPath'])
print('selecting time step = 600 min  ')
#-- now select time step to use (600s is usually used)
selTime = 600.
times = mcTable['time']
mcTable = mcTable[times==selTime]
mcTable = mcTable.sort_values('sHeight')
print(mcTable)

#- now reading in McSnow output. If we did not run McSnow, need to comment that out TODO: make that automatic
McRadar_Outname = os.environ['freq']+'GHz_output_{mode}.nc'.format(mode=dicSettings['scatSet']['mode'])
output = xr.open_dataset(inputPath+McRadar_Outname)
print('plotting spectra')
plot.plotSpectra(dicSettings,output,inputPath)#,convoluted=True)
print('plotting moments')
plot.plotMoments(dicSettings,output,inputPath)
if len(freq)==2:
  print('plotting DWR')
  wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
  wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
  plot.plotDWR(dicSettings,wlStr1,wlStr2,output,inputPath)
  plot.plotDWRspectra(dicSettings,wlStr1,wlStr2,output,inputPath)
elif len(freq)==3:
  print('plotting DWR')
  wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
  wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
  wlStr3 = '{:.2e}'.format(dicSettings['wl'][2])
  plot.plotDWR(dicSettings,wlStr1,wlStr2,output,inputPath)
  plot.plotDWR(dicSettings,wlStr2,wlStr3,output,inputPath)
  plot.plotDWRspectra(dicSettings,wlStr1,wlStr2,output,inputPath)
  plot.plotDWRspectra(dicSettings,wlStr2,wlStr3,output,inputPath)


# now plotting stuff directly from McSnow output but in the shape of a velocity spectrum:
print('plotting aspect ratios')
velBins = np.linspace(-3,0,100)
dBins = 10**(np.linspace(-3,0,100))
#print(dBins)
print('plotting sizes')
plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'dia_cm')
print('plotting masses')
plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'mTot_g')
print('plotting sNmono')
plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'sNmono')
print('plotting number_conc')
plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'number_conc')
