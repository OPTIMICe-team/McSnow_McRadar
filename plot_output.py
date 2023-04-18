import mcradar as mcr
from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
def str2bool(v):
  return v.lower() in ("yes", "True", "t", "1","true")
  
freqEnv = os.environ['freq'].split('_')
print(freqEnv)
elvEnv = os.environ['elv'].split('_')
#elv = float(os.environ['elv'])
elv = np.array([float(e) for e in elvEnv])
particle_name = os.environ['particle']
freq = np.array([float(f)*1e9 for f in freqEnv])
experimentID = os.environ['experiment']
print(experimentID)
#quit()
inputPath = os.environ['MCexp']+'/'+experimentID+'/'
scatMode = os.environ['scatMode']
print(inputPath)
#convolute=os.environ['convolute']
print('loading the settings')
splitPath = inputPath.split('domtop')[1]
domTop = splitPath[0:4]
lutPath = os.environ['LUT_dir']
convoluted=str2bool(os.environ['convolute'])

# decide what you want to plot
plotForwardSim = True # plot forward simulations
plotProperties = True # plot properties
plotTemp=True # use T as y-axis

#-- load the settings of McSnow domain, as well as elevation you want to plot:
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',
                               elv=elv, freq=freq,gridBaseArea=5.0,maxHeight=int(domTop),
                               ndgsVal=50,heightRes=36,scatSet={'mode':scatMode, 'lutPath':lutPath,'particle_name':particle_name,'safeTmatrix':False})

print('loading the McSnow output')
# now generate a table from the McSnow output. You can specify xi0, if it is not stored in the table (like in my cases)
mcTable = mcr.getMcSnowTable(dicSettings['dataPath'])

#- now reading in McSnow output. If we did not run McSnow, need to comment that out TODO: make that automatic
McRadar_Outname = os.environ['McRadarfileName']

# plot property spectra from McSnow output. THis is done for all simulations that are NOT a trajectory
#- read in atmo file to get Temperature information if you want to plot it with that
selTime = mcTable['time'].max()
times = mcTable['time']
mcTable = mcTable.where(times==selTime,drop=True)
mcTable = mcTable.sortby('sHeight').set_index(index='sHeight').rename({'index':'sHeight'})

atmoFile = np.loadtxt(inputPath+'atmo.dat')
plot.plotAtmo(atmoFile,inputPath)
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcTable,method='nearest')
mcTableTmp = xr.merge([atmoReindex,mcTable])

# now plotting stuff directly from McSnow output but in the shape of a velocity spectrum:
velBins = np.linspace(-3,0,100)
dBins = 10**(np.linspace(-3,0,100))
#print(dBins)
if plotProperties:
	print('plotting particle properties')
	fig,ax=plt.subplots(ncols=3,nrows=2,figsize=(15,10))
	varVec = ['dia','mTot','sRho_tot']
	for i,var in enumerate(varVec):	
		print(var)
		ax[0,i]=plot.plotPropSpecThesis(ax[0,i],dicSettings['heightRange'],dicSettings['heightRes'],mcTableTmp,velBins,var)
		ax[0,i].set_ylim([mcTableTmp['Temp'].max()+0.1,mcTableTmp['Temp'].min()-0.1])
		ax[0,i].set_xlim([-2,0])
		ax[0,i].tick_params(axis='both',labelsize=16)
		ax[0,i].grid(ls='-.')
		if mcTableTmp['Temp'].min() < -20:
			ax[0,i].axhline(y=-20,ls='--',color='r',linewidth=2)
			ax[0,i].axhline(y=-10,ls='--',color='r',linewidth=2)
			ax[0,i].text(ax[0,i].get_xlim()[0]+0.04*(ax[0,i].get_xlim()[1]-ax[0,i].get_xlim()[0]),mcTableTmp['Temp'].min()+2,'('+string.ascii_lowercase[i]+')',fontsize=18)
		else:
			ax[0,i].text(ax[0,i].get_xlim()[0]+0.04*(ax[0,i].get_xlim()[1]-ax[0,i].get_xlim()[0]),mcTableTmp['Temp'].min()+3,'('+string.ascii_lowercase[i]+')',fontsize=18)
		if i == 0:
			ax[0,i].set_ylabel('T [°C]',fontsize=18)
		else:
			ax[0,i].set_ylabel('')
			ax[0,i].set_yticklabels('')
	#plot.plotPropSpecThesis(dicSettings,mcTable,velBins,inputPath,'dia_cm')
	varVec = ['sNmono','sPhi','number_conc']
	for i,var in enumerate(varVec):	
		print(var)
		ax[1,i]=plot.plotPropSpecThesis(ax[1,i],dicSettings['heightRange'],dicSettings['heightRes'],mcTableTmp,velBins,var)
		ax[1,i].set_xlim([-2,0])
		ax[1,i].set_ylim([mcTableTmp['Temp'].max()+.1,mcTableTmp['Temp'].min()-.1])
		ax[1,i].tick_params(axis='both',labelsize=16)
		ax[1,i].grid(ls='-.')
		if mcTableTmp['Temp'].min() < -20:
			ax[1,i].axhline(y=-20,ls='--',color='r',linewidth=2)
			ax[1,i].axhline(y=-10,ls='--',color='r',linewidth=2)
			ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),mcTableTmp['Temp'].min()+2,'('+string.ascii_lowercase[i+3]+')',fontsize=18)
		
		else:
			ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),mcTableTmp['Temp'].min()+3,'('+string.ascii_lowercase[i+3]+')',fontsize=18)
		
		if i == 0:
			ax[1,i].set_ylabel('T [°C]',fontsize=18)
		else:
			ax[1,i].set_ylabel('')
			ax[1,i].set_yticklabels('')
	plt.tight_layout()
	plt.savefig(inputPath+'properties.png')
	plt.close()
#quit()

#-- plot number of superparticles per grid cell

nz = float(inputPath.split('nz')[1].split('_')[0])
#heightProf = pd.read_csv(inputPath+'hei2massdens.dat',header=0)#,skiprows=1, # header=0
#print(heightProf)
print('plotting number of sp')
plot.plotHeightProf(nz,mcTable,inputPath,dicSettings)
print('plotting PSD')
mBins = 10**(np.linspace(-12,-9,100))
plot.plotPSD(mcTable,dicSettings,inputPath,mBins,'mTot',heightEdge0=1900,unit='[kg]',sepMono=True,yscale='log',xscale='log')
dBins = 10**(np.linspace(-6,-2.5,100))
plot.plotPSD(mcTable,dicSettings,inputPath,dBins,'dia',heightEdge0=1900,unit='[m]',sepMono=True,yscale='log',xscale='log')

# plot McRadar output
if os.path.exists(inputPath+McRadar_Outname) and (plotForwardSim==True):
	output = xr.open_dataset(inputPath+McRadar_Outname)
	print('now plotting McRadar')
	if plotTemp == True:
		height = atmoFile[:,0]
		Temp = atmoFile[:,2] -273.15
		atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
		atmoPD.index.name='range'
		atmoXR = atmoPD.to_xarray()
		atmoReindex = atmoXR.reindex_like(output,method='nearest')
		output = xr.merge([atmoReindex,output])
	if (len(freq)==2) and (len(elv)>1):
		plot.plotOverview(output,dicSettings,inputPath,dicSettings['wl'][0],dicSettings['wl'][1])
	elif (len(freq)==3) and (len(elv)>1):
		plot.plotOverview(output,dicSettings,inputPath,dicSettings['wl'][1],dicSettings['wl'][2])
	
	print('plotting spectra')
	plot.plotSpectra(dicSettings,output,inputPath,plotTemp=plotTemp)

	print('plotting moments')
	plot.plotMoments(dicSettings,output,inputPath,plotTemp=plotTemp)

	if len(freq)==2:
		print('plotting DWR')
		plot.plotDWR(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
		plot.plotDWRspectra(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
	elif len(freq)==3:
		print('plotting DWR')
		plot.plotDWR(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
		plot.plotDWR(dicSettings,dicSettings['wl'][1],dicSettings['wl'][2],output,inputPath,plotTemp=plotTemp)
		plot.plotDWRspectra(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
		plot.plotDWRspectra(dicSettings,dicSettings['wl'][1],dicSettings['wl'][2],output,inputPath,plotTemp=plotTemp)










