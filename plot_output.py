import mcradar as mcr
from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
freq = np.array([35.5e9])
experimentID = os.environ['experiment']
print(experimentID)
inputPath = os.environ['MCexp']+'experiments/'+experimentID+'/'
print('loading the settings')
#-- load the settings of McSnow domain, as well as elevation you want to plot:
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',
                               elv=30, freq=freq,gridBaseArea=5.0,maxHeight=3850,ndgsVal=50,heightRes=50,scatSet={'mode':'full', 'safeTmatrix':False})

print('loading the McSnow output')
# now generate a table from the McSnow output. You can specify xi0, if it is not stored in the table (like in my cases)
mcTable = mcr.getMcSnowTable(dicSettings['dataPath'])
print('selecting time step = 600 min  ')
#-- now select time step to use (600s is usually used)
selTime = 600.
times = mcTable['time']
mcTable = mcTable[times==selTime]
mcTable = mcTable.sort_values('sHeight')

output = xr.open_dataset(inputPath+'KaBand_output.nc')
print('plotting spectra')
plot.plotSpectra(dicSettings,output,inputPath)
print('plotting moments')
plot.plotMoments(dicSettings,output,inputPath)
print('plotting aspect ratios')
velBins = np.linspace(-3,0,100)
dBins = 10**(np.linspace(-3,0,100))
plot.plotArSpec(dicSettings,mcTable,velBins,inputPath)

