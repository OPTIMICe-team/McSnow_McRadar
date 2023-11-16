from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines_noMcRadar as plot
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import string

experimentID = os.environ['experiment']
print(experimentID)
#quit()
inputPath =  os.environ['MCexp']+'/'+experimentID+'/'
inputPath = '1d_habit_soumi_habit1_IGF2_xi100_nz200_iwc3_nugam3.545_mugam0.5_dtc5_nh13000_nh24000_ncl0_nclmass4.8_nuclType1_nrp10_at2_stick2_spkernsig10_timeend36000_dt1d600_domtop8000._atmo2_radiosondes_juelich_20220206_042141_ssat50/'
print(inputPath)
#convolute=os.environ['convolute']
print('loading the settings')
domTop = inputPath.split('domtop')[1].split('_')[0].split('.')[0]
box_area=5

maxHeight = float(inputPath.split('domtop')[1].split('_')[0].split('.')[0])
minHeight=0;  heightRes=36
heightRange = np.arange(minHeight, maxHeight, heightRes)
velBins = np.linspace(-3,0,100)
mcTable = xr.open_dataset(inputPath+'mass2fr.nc')
mcTable['vel'] = -1. * mcTable['vel']
mcTable['radii_mm'] = mcTable['dia'] * 1e3 / 2.
mcTable['dia_mum'] = mcTable['dia'] * 1e6 
mcTable['mTot_g'] = mcTable['mTot'] * 1e3
mcTable['dia_cm'] = mcTable['dia'] * 1e2
selTime = mcTable['time'].max()
times = mcTable['time']
#mcTable = mcTable.sort_values('sHeight')
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
#mcTableTmp = mcTableTmp.to_dataframe()

#mcTableMono = mcTableTmp[mcTableTmp.sNmono==1]
#quit()
# now plotting stuff directly from McSnow output but in the shape of a velocity spectrum:
#print('plotting aspect ratios')
velBins = np.linspace(-3,0,100)
dBins = 10**(np.linspace(-3,0,100))

#print(dBins)
#if plotProp:
print(mcTableTmp)
iwp = (mcTableTmp.mTot*mcTable.sMult).sum()/box_area # iwp in kg/m2

#quit()
if not os.path.isfile(inputPath+'properties.png'):
	print('plotting particle properties')
	fig,ax=plt.subplots(ncols=3,nrows=2,figsize=(15,10))
	varVec = ['dia','mTot','sRho_tot']
	for i,var in enumerate(varVec):	
		print(var)
		ax[0,i]=plot.plotPropSpecThesis(ax[0,i],heightRange,heightRes,mcTableTmp,velBins,var)
		ax[0,i].set_ylim([0,mcTableTmp['Temp'].min()-0.1])#[mcTableTmp['Temp'].max()+0.1,mcTableTmp['Temp'].min()-0.1])
		ax[0,i].set_xlim([-2,0])
		ax[0,i].tick_params(axis='both',labelsize=16)
		#ax[0,i].text(ax[0,i].get_xlim()[0]+0.04*(ax[0,i].get_xlim()[1]-ax[0,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
		ax[0,i].grid(ls='-.')
		if (mcTableTmp['Temp'].min() < -20) and (plotDGLlines==True):
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
		print(var)
		ax[1,i]=plot.plotPropSpecThesis(ax[1,i],heightRange,heightRes,mcTableTmp,velBins,var)
		ax[1,i].set_xlim([-2,0])
		#ax[1,i].set_ylim([mcTableTmp['Temp'].max()+.1,mcTableTmp['Temp'].min()-.1])
		ax[1,i].set_ylim([0,mcTableTmp['Temp'].min()-0.1])
		ax[1,i].tick_params(axis='both',labelsize=16)
		#ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+3]+')',fontsize=18)
		ax[1,i].grid(ls='-.')
		if (mcTableTmp['Temp'].min() < -20) and (plotDGLlines==True):
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
	

		
#-- plot number of superparticles per grid cell
	
if 'nz' in inputPath and not os.path.isfile(inputPath+'number_sParticle.png'):
	nz = float(inputPath.split('nz')[1].split('_')[0])
	print('plotting number of sp')
	plot.plotconcHeightProf('Super',mcTable,inputPath,box_area,maxHeight,nz=nz)
#heightProf = pd.read_csv(inputPath+'hei2massdens.dat',header=0)#,skiprows=1, # header=0
#print(heightProf)
if not os.path.isfile(inputPath+'concentration_realParticle1.png'):
	print('plotting number concentration profile')
	plot.plotconcHeightProf('Real',mcTableTmp,inputPath,box_area,maxHeight,vmax=15000)
print('plotting PSD')
mBins = 10**(np.linspace(-12,-4,100))
plot.plotPSD(mcTable,box_area,heightRes,inputPath,mBins,'mTot',heightEdge0=1900,unit='[kg]',sepMono=True,yscale='log',xscale='log')
plot.plotPSD(mcTable,box_area,heightRes,inputPath,mBins,'mTot',heightEdge0=2300,unit='[kg]',sepMono=True,yscale='log',xscale='log')
plot.plotPSD(mcTable,box_area,heightRes,inputPath,mBins,'mTot',heightEdge0=2600,unit='[kg]',sepMono=True,yscale='log',xscale='log')
plot.plotPSD(mcTable,box_area,heightRes,inputPath,mBins,'mTot',heightEdge0=2900,unit='[kg]',sepMono=True,yscale='log',xscale='log')
#dBins = 10**(np.linspace(-6,-1,100))
dBins = 10**(np.linspace(-3,2,100))
mcTable['dia_mm'] = mcTable['dia'] * 1e3
plot.plotPSD(mcTable,box_area,heightRes,inputPath,dBins,'dia_mm',heightEdge0=1900,unit='[mm]',sepMono=True,yscale='log',xscale='log')
plot.plotPSD(mcTable,box_area,heightRes,inputPath,dBins,'dia_mm',heightEdge0=2300,unit='[mm]',sepMono=True,yscale='log',xscale='log')
plot.plotPSD(mcTable,box_area,heightRes,inputPath,dBins,'dia_mm',heightEdge0=2600,unit='[mm]',sepMono=True,yscale='log',xscale='log')
plot.plotPSD(mcTable,box_area,heightRes,inputPath,dBins,'dia_mm',heightEdge0=2900,unit='[mm]',sepMono=True,yscale='log',xscale='log')











