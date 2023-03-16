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
elv = float(os.environ['elv'])
particle_name = os.environ['particle']
outName = os.environ['McRadarfileName']
freq = np.array([float(f)*1e9 for f in freqEnv])
experimentID = os.environ['experiment']
inputPath = os.environ['MCexp']+'/'+experimentID+'/'
scatMode = os.environ['scatMode']
lutPath = '/project/meteo/work/L.Terzi/McRadar/LUT/' #'/work/lvonterz/SSRGA/snowScatt/ssrga_LUT/' #'/data/optimice/McRadarLUTs/'
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
	heightRes = 36
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
if allTimes==True:
	for i,selTime in enumerate(times.unique()):
		mcTableTmp = mcTable[times==selTime]
		mcTableTmp = mcTableTmp.sort_values('sHeight')
		print('getting things done :) -> calculating radar variables for '+str(freq)+'GHz')
		outputTime = mcr.fullRadar(dicSettings, mcTableTmp)
		print(output)
		for wl in dicSettings['wl']:
			wlStr = '{:.2e}'.format(wl)
			if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
				outputTime['Ze_H_{0}'.format(wlStr)] = outputTime['spec_H_{0}'.format(wlStr)].sum(dim='vel')
				outputTime['MDV_H_{0}'.format(wlStr)] = (outputTime['spec_H_{0}'.format(wlStr)]*outputTime['vel']).sum(dim='vel')/outputTime['Ze_H_{0}'.format(wlStr)]
				      
			else:
				outputTime['Ze_H_{0}'.format(wlStr)] = outputTime['spec_H_{0}'.format(wlStr)].sum(dim='vel')
				outputTime['Ze_V_{0}'.format(wlStr)] = outputTime['spec_V_{0}'.format(wlStr)].sum(dim='vel')
				outputTime['MDV_H_{0}'.format(wlStr)] = (outputTime['spec_H_{0}'.format(wlStr)]*outputTime['vel']).sum(dim='vel')/outputTime['Ze_H_{0}'.format(wlStr)]
				outputTime['MDV_V_{0}'.format(wlStr)] = (outputTime['spec_V_{0}'.format(wlStr)]*outputTime['vel']).sum(dim='vel')/outputTime['Ze_V_{0}'.format(wlStr)]
			print(output)
			print(selTime)
		outputTime = outputTime.expand_dims(dim={'time':1})
		outputTime = outputTime.assign_coords(time=np.asarray(selTime).reshape(1))
		if i == 0:
			output = outputTime
		else:
			output = xr.merge([outputTime,output])  
	outName = os.environ['freq']+'GHz_output_{mode}_alltimes.nc'.format(mode=dicSettings['scatSet']['mode'])
	print(outputTime)
else:
	if ('trajectories' not in experimentID) and ('trajectories' not in inputPath):
		selTime = mcTable['time'].max()
		
		mcTableTmp = mcTable[times==selTime]	
	else: # if we have trajectories it makes sense to have single_particles = True!!!
		if single_particle == False:
			mcTable['sMult'] = 1.0 
		mcTableTmp = mcTable
	print(mcTableTmp.time)
	
	print('getting things done :) -> calculating radar variables for '+str(freq)+'GHz')
	print(single_particle)
	
	if single_particle == False:
		#print(mcTableTmp)
		output = mcr.fullRadar(dicSettings, mcTableTmp)
		print(output)
		for wl in dicSettings['wl']:
			wlStr = '{:.2e}'.format(wl)
			if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
				output['Ze_H_{0}'.format(wlStr)] = output['spec_H_{0}'.format(wlStr)].sum(dim='vel')
				output['MDV_H_{0}'.format(wlStr)] = (output['spec_H_{0}'.format(wlStr)]*output['vel']).sum(dim='vel')/output['Ze_H_{0}'.format(wlStr)]
				      
			else:
				output['Ze_H_{0}'.format(wlStr)] = output['spec_H_{0}'.format(wlStr)].sum(dim='vel')
				output['Ze_V_{0}'.format(wlStr)] = output['spec_V_{0}'.format(wlStr)].sum(dim='vel')
				output['MDV_H_{0}'.format(wlStr)] = (output['spec_H_{0}'.format(wlStr)]*output['vel']).sum(dim='vel')/output['Ze_H_{0}'.format(wlStr)]
				output['MDV_V_{0}'.format(wlStr)] = (output['spec_V_{0}'.format(wlStr)]*output['vel']).sum(dim='vel')/output['Ze_V_{0}'.format(wlStr)]
	
	else:
	
		output = mcr.singleParticleTrajectories(dicSettings, mcTableTmp)
		print(output)

print('saving the output file at: '+inputPath+outName)
#-- now save it
output.to_netcdf(inputPath+outName)

