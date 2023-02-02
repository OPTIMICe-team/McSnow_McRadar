from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import string

def loadSettings(dataPath=None, freq=np.array([9.5e9, 35e9, 95e9]),
				 elv=90, nfft=512, maxVel=3, minVel=-3, ndgsVal=30,  # setup Doppler spectra
				 convolute=True,nave=19,noise_pow=10**(-40/10), eps_diss=1e-6, theta=0.6 , uwind=10.0 , time_int=2.0 , # for noise convolution and turbulence broadening
                 maxHeight=5500, minHeight=0, heightRes=50, gridBaseArea=1, # setup of McSnow simulation
                 scatSet={'mode':'full',
                          'safeTmatrix':False}):
    #- copied here from McRadar
    
    dicSettings = {'dataPath':dataPath,
                   'elv':elv,
                   'nfft':nfft,
                   'maxVel':maxVel,
                   'minVel':minVel,
                   'velRes':(maxVel - minVel)/nfft,
                   'freq':freq,
                   'wl':(constants.c / freq) * 1e3, #[mm]
                   'ndgsVal':ndgsVal,
                   'maxHeight':maxHeight,
                   'minHeight':minHeight,
                   'heightRes':heightRes,
                   'heightRange':np.arange(minHeight, maxHeight, heightRes),
                   'gridBaseArea':gridBaseArea,
                   'scatSet':scatSet,
                   'convolute':convolute,
                   'nave':nave,
                   'noise_pow':noise_pow,
                   'eps_diss':eps_diss,
                   'theta':theta,
                   'time_int':time_int,
                   'uwind':uwind,
                   }

    velBins = np.arange(minVel, maxVel, dicSettings['velRes'])
    velCenterBin = velBins[0:-1]+np.diff(velBins)/2.

    dicSettings['velBins']=velBins
    dicSettings['velCenterBin']=velCenterBin

	return dicSettings

def str2bool(v):
  return v.lower() in ("yes", "True", "t", "1","true")

experimentID = os.environ['experiment']
print(experimentID)
#quit()
inputPath = os.environ['MCexp']+'/'+experimentID+'/'

splitPath = inputPath.split('domtop')[1]
domTop = splitPath[0:4]

# decide what you want to plot
plot_initemp = True
plot_inidia = False
plot_thesis_particle_evolution = False
plotTemp=True

#-- load the settings of McSnow domain, as well as elevation you want to plot:
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation

dicSettings = loadSettings(dataPath=inputPath+'mass2fr.nc',
                               gridBaseArea=5.0,maxHeight=int(domTop),heightRes=36,maxVel=0)

print('loading the McSnow output')
# now generate a table from the McSnow output. You can specify xi0, if it is not stored in the table (like in my cases)
mcTableXR = xr.open_dataset(dicSettings['dataPath'])
#print(mcTableXR)
mcTableXR = mcTableXR.astype('float64')
#change to pandas dataframe, since McRadar has been working with that
mcTable = mcTableXR.to_dataframe() 
mcTable['vel'] = -1. * mcTable['vel']
mcTable['radii_mm'] = mcTable['dia'] * 1e3 / 2.
mcTable['dia_mum'] = mcTable['dia'] * 1e6 
mcTable['mTot_g'] = mcTable['mTot'] * 1e3
mcTable['dia_cm'] = mcTable['dia'] * 1e2

mcTable = mcTable.sort_values('sHeight')

#- read in atmo file to get Temperature information if you want to plot it with that
selTime = mcTable['time'].max()
times = mcTable['time']
mcTable = mcTable.sort_values('sHeight')
mcTable = mcTable[times==selTime]

#- add temperature to the McSnow output file
atmoFile = np.loadtxt(inputPath+'atmo.dat')
plot.plotAtmo(atmoFile,inputPath)
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
mcTableNew=mcTable.set_index('sHeight',drop=False)
mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
mcTableXR = mcTableNew.to_xarray(); atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
mcTableTmp = xr.merge([atmoReindex,mcTableXR])
mcTableTmp = mcTableTmp.to_dataframe()

# now plotting stuff from McSnow output but in the shape of a velocity spectrum:
#print('plotting aspect ratios')
velBins = dicSettings['velBins']#np.linspace(-3,0,100)

print('plotting particle properties')
fig,ax=plt.subplots(ncols=3,nrows=2,figsize=(15,10))
varVec = ['dia','mTot','sRho_tot']
for i,var in enumerate(varVec):	
	ax[0,i]=plot.plotPropSpecThesis(ax[0,i],dicSettings['heightRange'],dicSettings['heightRes'],mcTableTmp,velBins,var)
	ax[0,i].set_ylim([0,-30])
	ax[0,i].set_xlim([-2,0])
	ax[0,i].tick_params(axis='both',labelsize=16)
	ax[0,i].text(ax[0,i].get_xlim()[0]+0.04*(ax[0,i].get_xlim()[1]-ax[0,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
	ax[0,i].grid(ls='-.')
	ax[0,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[0,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[0,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[0,i].set_ylabel('')
		ax[0,i].set_yticklabels('')
#plot.plotPropSpecThesis(dicSettings,mcTable,velBins,inputPath,'dia_cm')
varVec = ['sNmono','sPhi','number_conc']
for i,var in enumerate(varVec):	
	ax[1,i]=plot.plotPropSpecThesis(ax[1,i],dicSettings['heightRange'],dicSettings['heightRes'],mcTableTmp,velBins,var)
	ax[1,i].set_xlim([-2,0])
	ax[1,i].set_ylim([0,-30])
	ax[1,i].tick_params(axis='both',labelsize=16)
	ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+3]+')',fontsize=18)
	ax[1,i].grid(ls='-.')
	ax[1,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[1,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[1,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[1,i].set_ylabel('')
		ax[1,i].set_yticklabels('')
plt.tight_layout()
plt.savefig(inputPath+'properties.png')
plt.close()
#quit()



	
	
		
	








