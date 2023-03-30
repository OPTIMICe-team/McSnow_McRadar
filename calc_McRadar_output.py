#-*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Author: Leonie von Terzi


# this calculates the polarimetric variables at Wband for McSnow output. 
# It is intended to test habit prediction, aggregation has not been implemented in this McSnow run.
# The McSnow data was produced by Jan-Niklas Welß

import numpy as np
import mcradar as mcr
from scipy import constants
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib as mpl
import pandas as pd
import xarray as xr
import os
def str2bool(v):
  return v.lower() in ("yes", "True", "t", "1","true")

#- get all variables necessary to calculate scattering from environment
allTimes=False
single_particle=str2bool(os.environ['singleParticle'])
convolute=str2bool(os.environ['convolute'])
freqEnv = os.environ['freq'].split('_')
elvEnv = os.environ['elv'].split('_')
#elv = float(os.environ['elv'])
elv = np.array([float(e) for e in elvEnv])
particle_name = os.environ['particle']
outName = os.environ['McRadarfileName']
freq = np.array([float(f)*1e9 for f in freqEnv])
experimentID = os.environ['experiment']
inputPath = os.environ['MCexp']+'/'+experimentID+'/'
scatMode = os.environ['scatMode']
lutPath = os.environ['LUT_dir'] #'/project/meteo/work/L.Terzi/McRadar/LUT/' #'/work/lvonterz/SSRGA/snowScatt/ssrga_LUT/' #'/data/optimice/McRadarLUTs/'
if 'DDA' in scatMode:
	lutPath = lutPath + 'DDA/'
elif 'SSRGA' in scatMode:
	lutPath = lutPath+ 'SSRGA/'
elif scatMode == 'wisdom' or scatMode == 'table':
	lutPath = lutPath + 'Tmatrix'
domTop = inputPath.split('domtop')[1].split('_')[0].split('.')[0]

try:
  minmax = os.environ['minmax'] 
  vmin=int(minmax.split('_')[0]); vmax=int(minmax.split('_')[1]) 
  print(minmax)
except: 
  minmax=False
  print('no minmax')
print('loading the settings')
#minmax=True
#vmin= 180; vmax=350

#-- load the settings of McSnow domain, as well as elevation you want to plot: 
if ('trajectories' not in experimentID) and ('trajectories' not in inputPath):
	heightRes = 50
else:
	heightRes = 2
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',
                               elv=elv, freq=freq,gridBaseArea=5.0,maxHeight=int(domTop),
                               ndgsVal=50,heightRes=heightRes,convolute=convolute,scatSet={'mode':scatMode,'lutPath':lutPath,'particle_name':particle_name,'safeTmatrix':True})

print('loading the McSnow output')
# now generate a table from the McSnow output.
mcTable = mcr.getMcSnowTable(dicSettings['dataPath'])

# for trajectories:
if minmax:
  mcTable = mcTable[(mcTable['sMult']>vmin) & (mcTable['sMult']<=vmax)] 

#quit()

times = mcTable['time']
if ('trajectories' not in experimentID) and ('trajectories' not in inputPath):
	selTime = mcTable['time'].max()
	
	mcTableTmp = mcTable.where(times==selTime,drop=True)	#mcTable[times==selTime]#
	
	#quit()
else: # if we have trajectories it makes sense to have single_particles = True!!!
	if single_particle == False:
		mcTable['sMult'] = 1.0 
	mcTableTmp = mcTable
print('getting things done :) -> calculating radar variables for '+str(freq)+'GHz')

if single_particle == False:
	output = mcr.fullRadar(dicSettings, mcTableTmp)
	print(output)
	if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
		output['Ze_H'] = output['spec_H'].sum(dim='vel')
		output['MDV_H'] = (output['spec_H']*output['vel']).sum(dim='vel')/output['Ze_H']
	else:	
		output['Ze_H'] = output['spec_H'].sum(dim='vel')
		output['Ze_V'] = output['spec_V'].sum(dim='vel')
		output['Ze_HV'] = output['spec_HV'].sum(dim='vel')
		output['LDR'] = mcr.lin2db(output['Ze_HV']/output['Ze_H'])
		output['MDV_H'] = (output['spec_H']*output['vel']).sum(dim='vel')/output['Ze_H']
		output['MDV_V'] = (output['spec_V']*output['vel']).sum(dim='vel')/output['Ze_V']

else:

	output = mcr.singleParticleTrajectories(dicSettings, mcTableTmp)
	print(output)

print('saving the output file at: '+inputPath+outName)
#-- now save it
output.to_netcdf(inputPath+outName)

