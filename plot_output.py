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
'''
#- now reading in McSnow output. If we did not run McSnow, need to comment that out TODO: make that automatic
McRadar_Outname = os.environ['freq']+'GHz_output_{mode}.nc'.format(mode=dicSettings['scatSet']['mode'])
output = xr.open_dataset(inputPath+McRadar_Outname)
print('plotting spectra')
plot.plotSpectra(dicSettings,output,inputPath)#,convoluted=True)
print('plotting moments')
plot.plotMoments(dicSettings,output,inputPath)
#quit()

#- read in atmo file to get Temperature information if you want to plot it with that
atmoFile = np.loadtxt(inputPath+'atmo.dat')

plot.plotAtmo(atmoFile,inputPath)
#quit()
'''
'''
# now plotting stuff directly from McSnow output but in the shape of a velocity spectrum:
print('plotting aspect ratios')
velBins = np.linspace(-3,0,100)
dBins = 10**(np.linspace(-3,0,100))
#print(dBins)
plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'sPhi')
#quit()
print('plotting sizes')
plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'dia_cm')
print('plotting masses')
#plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'mTot_g')
print('plotting sNmono')
#plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'sNmono')
#print('plotting sNmono min')
#plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'sNmono_min')
#print('plotting sNmono max')
#plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'sNmono_max')
print('plotting number concentration')
#plot.plotPropSpec(dicSettings,mcTable,velBins,inputPath,'number_conc')
#-- plot PSD
#plot.plotPSD(dicSettings,mcTable,dBins,inputPath)
#quit()

#-- now plotting PSD if you want to do that
dBins = 10**(np.linspace(-4,0,100))*10
phiBins = 10**(np.linspace(-3,0,100))
for i, heightEdge0 in enumerate(dicSettings['heightRange'][::-1]):
  heightEdge1 = heightEdge0 + dicSettings['heightRes']
  height = heightEdge0+dicSettings['heightRes']/2 
  mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &(mcTable['sHeight']<=heightEdge1)].copy()
  
  fig,ax = plt.subplots()
  ax.hist(mcTableTmp.dia_cm,bins=dBins,weights=mcTableTmp.sMult)
  #ax.set_xscale('log')
  #ax.set_yscale('log')
  ax.grid()
  ax.set_xlabel('dia [cm]')
  ax.set_ylabel('#')
  ax.set_xlim(0,0.1)#[10**-3*10,0.5*10**0])
  plt.tight_layout()
  plt.savefig(inputPath+'{height}_dia_hist_weigth.png'.format(height=height))
  plt.close()
  print(height,' done')
  quit()  
'''
'''
#-- plot proposal
atmoFile = pd.read_csv(inputPath+'atmo.dat')#,delimiter=' ')#,header=None,names=['height','1','temp','3','4','5','rh','S_i','8''9'])
atmoPD = pd.DataFrame(atmoFile)
atmoPD = atmoPD[[0,2,6]]
atmoPD.columns = ['height','temp','RH']
atmoPD['temp'] = atmoPD['temp'] - 273.15
atmoPD = atmoPD.set_index('height')
#Temp = xr.DataArray(np.empty(len(output.range)),dims='range',coords={'range':output.range.values})
Temp = np.empty(len(output.range))

concMono = Temp.copy()
concAgg = Temp.copy()
for i, heightEdge0 in enumerate(dicSettings['heightRange']):
  heightEdge1 = heightEdge0 + dicSettings['heightRes']
  height = heightEdge0+dicSettings['heightRes']/2
  atmoTmp = atmoPD[(atmoPD['height']>heightEdge0) & (atmoPD['height']<=heightEdge1)].copy()
  Temp[i] = atmoTmp.temp.mean()

  mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &(mcTable['sHeight']<=heightEdge1)].copy()
  mcTableMono = mcTableTmp[mcTableTmp['sNmono'] == 1]
  mcTableAgg = mcTableTmp[mcTableTmp['sNmono'] > 1]
  concMono[i] = mcTableMono['sMult'].sum()
  concAgg[i] = mcTableAgg['sMult'].sum()

ConcMono = xr.DataArray(concMono,name='conc_mono',dims='temperature [°C]',coords={'temperature [°C]':Temp})
ConcAgg = xr.DataArray(concAgg,name='conc_agg',dims='temperature [°C]',coords={'temperature [°C]':Temp})

#TempXR = xr.DataArray(Temp,name='Temp',dims='range',coords={'range':output.range.values})
#output = xr.merge([output,TempXR])
outputNew = output.assign({'range':Temp})
outputNew = outputNew.rename({'range':'temperature [°C]'})
#read in concoluted file
outputConv = xr.open_dataset(inputPath+'specXR_convolution_ar_0.008.nc')
outputConv = outputConv.assign({'range':Temp})
outputConv = outputConv.rename({'range':'temperature [°C]'})

outputAll = xr.merge([ConcMono,ConcAgg,outputNew['Ze_H_3.12e+01'],outputNew['Ze_V_3.12e+01'],outputNew['kdpInt_3.12e+01'],outputConv['spec_H_3.12e+01'],outputConv['spec_V_3.12e+01']])
#atmoPD = atmoPD.rename({'height','temp'})

plot.plotProposal(dicSettings,outputAll,inputPath)
quit()
'''

#- plot ar test setup (with bnd_type==3)
atmoFile = np.loadtxt(inputPath+'atmo.dat')

#Temp = xr.DataArray(np.empty(len(output.range)),dims='range',coords={'range':output.range.values})  
# plot only certain sizes
#mcTable = mcTable[times>=0]
#mcTable = mcTable[times<=20]
mcTable = mcTable.sort_values('time')
print(min(mcTable['sPhi']))
quit()
# setup color range for plotting (color is dependent on ID of particle)
plot.plot_var_particle_ID(mcTable,inputPath,'sPhi','time',atmoFile)#,zoom=[0,0.2])
plot.plot_var_particle_ID(mcTable,inputPath,'dia','time',atmoFile)
plot.plot_var_particle_ID(mcTable,inputPath,'dia_mum','time',atmoFile)
plot.plot_var_particle_ID(mcTable,inputPath,'mTot','time',atmoFile)

