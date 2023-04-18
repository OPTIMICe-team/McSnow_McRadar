#-*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
# Author: Leonie von Terzi
# This uses McRadar to calculate the scattering properties from the McSnow output. 
# The needed variables are provided in McSnow_McRadar.sh


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
freqEnv = os.environ['freq'].split('_')
elvEnv = os.environ['elv'].split('_')
freq = np.array([float(f)*1e9 for f in freqEnv])
elv = np.array([float(e) for e in elvEnv])

convolute=str2bool(os.environ['convolute'])
particle_name = os.environ['particle']
outName = os.environ['McRadarfileName']
experimentID = os.environ['experiment']
inputPath = os.environ['MCexp']+'/'+experimentID+'/'
scatMode = os.environ['scatMode']
lutPath = os.environ['LUT_dir'] 
if 'DDA' in scatMode:
	lutPath = lutPath + 'DDA/'
elif 'SSRGA' in scatMode:
	lutPath = lutPath+ 'SSRGA/'
elif scatMode == 'wisdom' or scatMode == 'table':
	lutPath = lutPath + 'Tmatrix'

domTop = inputPath.split('domtop')[1].split('_')[0].split('.')[0]
print('loading the settings')
#-- load the settings of McSnow domain, as well as elevation you want to plot: 
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',#'mass2fr.nc',#inputPath+'mass2fr.nc',
                               elv=elv, freq=freq,gridBaseArea=5.0,maxHeight=int(domTop),
                               ndgsVal=50,heightRes=36,convolute=convolute,#k_theta=0,k_phi=0,k_r=0,shear_height0=0,shear_height1=0,
                               scatSet={'mode':scatMode,'lutPath':lutPath,'particle_name':particle_name,'safeTmatrix':True})

print('loading the McSnow output')
# now generate a table from the McSnow output.
mcTable = mcr.getMcSnowTable(dicSettings['dataPath'])

#only calculate for last output step:
times = mcTable['time']
selTime = mcTable['time'].max()
mcTableTmp = mcTable.where(times==selTime,drop=True)	#mcTable[times==selTime]#
	
print('getting things done :) -> calculating radar variables for '+str(freq)+'Hz and '+str(elv)+'° elevation')

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


print('saving the output file at: '+inputPath+outName)
#-- now save it
output.to_netcdf(inputPath+outName)#inputPath+outName)

