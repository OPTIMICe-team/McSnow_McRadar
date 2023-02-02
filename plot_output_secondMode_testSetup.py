import mcradar as mcr
from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap
from scipy.optimize import curve_fit
import scipy.interpolate as intp
from scipy import constants
import postProcessSpectra as post
import string


'''
First: only plot one setup to show how variables generally look like.
'''

dataPath = '/data/optimice/McSnowoutput/habit/second_nucleation/'
nugam = 3.5; mugam = 0.5; IWC = 3; nrp=2
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
'''
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass1_ncl0_2nd_prim_testhabit'.format(IWC,nugam,mugam,nrp)
McRadar_Outname = '{freq}GHz_output_{mode}_convolute.nc'.format(freq=freq,mode=mode)
mcTable = xr.open_dataset(dataPath+expID+'/mass2fr.nc')
mcTable = mcTable.to_dataframe()
mcTable = mcTable.sort_values('sHeight')
selTime = mcTable['time'].max()
times = mcTable['time']
print(selTime)
#quit()
mcTable = mcTable[times==selTime]
mcTable['vel'] = -1*mcTable['vel']

atmoFile = np.loadtxt(dataPath+expID+'/atmo.dat')
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
mcTableNew=mcTable.set_index('sHeight',drop=False)
mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
mcTableXR = mcTableNew.to_xarray()
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
mcTableTmp = xr.merge([atmoReindex,mcTableXR])
mcTableTmp = mcTableTmp.to_dataframe()
print(mcTable)

mcRadar = xr.open_dataset(dataPath+expID+'/'+McRadar_Outname)
atmoPD.index.name='range'
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcRadar,method='nearest')
mcRadar = xr.merge([atmoReindex,mcRadar])

fig = plt.figure(figsize=(15,20))
gs0 = mpl.gridspec.GridSpec(2, 1, figure=fig,hspace=0.25)

gs00 = mpl.gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs0[0],hspace=0.25,wspace=0.25)
ax00 = fig.add_subplot(gs00[0,0])
ax01 = fig.add_subplot(gs00[0,1])
ax02 = fig.add_subplot(gs00[0,2])
ax10 = fig.add_subplot(gs00[1,0])
ax11 = fig.add_subplot(gs00[1,1])
ax12 = fig.add_subplot(gs00[1,2])
gs01 = mpl.gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs0[1],hspace=0.25,wspace=0.25)
ax20 = fig.add_subplot(gs01[0,0])
ax21 = fig.add_subplot(gs01[0,1])
ax22 = fig.add_subplot(gs01[0,2])
ax30 = fig.add_subplot(gs01[1,0])
ax31 = fig.add_subplot(gs01[1,1])
ax32 = fig.add_subplot(gs01[1,2])



ax = np.array([[ax00,ax01,ax02],[ax10,ax11,ax12],[ax20,ax21,ax22],[ax30,ax31,ax32]])

#fig,ax = plt.subplots()
varVec = ['number_conc','dia','sNmono']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	print(var)
	ax[0,i] = plot.plotPropSpecThesis(ax[0,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
	#ax[0,i].set_xticks([-15,-10,-5,0])
	ax[0,i].set_xlim([-2,0])
	ax[0,i].set_ylim([0,-30])
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
#ax[0,1].text(-1,-32,r'$S_{wi}$',fontsize=24)
ax[0,1].set_title(r'prim nucl., habit',fontsize=24,y=1.05)
varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	print(var)
	ax[1,i] = plot.plotRadar(mcRadar,var,ax[1,i])
	ax[1,i].set_ylim([0,-30])
	#ax[1,i].set_xlim([-3,0])
	ax[1,i].tick_params(axis='both',labelsize=16)
	ax[1,i].grid(ls='-.')
	ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+3]+')',fontsize=18)
	ax[1,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[1,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[1,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[1,i].set_ylabel('')
		ax[1,i].set_yticklabels('')

# now narrow PSD
#nugam = 10; mugam = 10; IWC = 1; nrp=1
expID ='1d_habit_habit0_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass1_ncl0_2nd_prim_testhabit'.format(IWC,nugam,mugam,nrp)
McRadar_Outname = '{freq}GHz_output_{mode}_convolute.nc'.format(freq=freq,mode=mode)
mcTable = xr.open_dataset(dataPath+expID+'/mass2fr.nc')
mcTable = mcTable.to_dataframe()
mcTable = mcTable.sort_values('sHeight')
selTime = mcTable['time'].max()
times = mcTable['time']
print(selTime)
#quit()
mcTable = mcTable[times==selTime]
mcTable['vel'] = -1*mcTable['vel']

atmoFile = np.loadtxt(dataPath+expID+'/atmo.dat')
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
mcTableNew=mcTable.set_index('sHeight',drop=False)
mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
mcTableXR = mcTableNew.to_xarray()
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
mcTableTmp = xr.merge([atmoReindex,mcTableXR])
mcTableTmp = mcTableTmp.to_dataframe()
print(mcTable)

mcRadar = xr.open_dataset(dataPath+expID+'/'+McRadar_Outname)
atmoPD.index.name='range'
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcRadar,method='nearest')
mcRadar = xr.merge([atmoReindex,mcRadar])

varVec = ['number_conc','dia','sNmono']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	print(var)
	ax[2,i] = plot.plotPropSpecThesis(ax[2,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
	#ax[0,i].set_xticks([-15,-10,-5,0])
	ax[2,i].set_xlim([-2,0])
	ax[2,i].set_ylim([0,-30])
	ax[2,i].tick_params(axis='both',labelsize=16)
	ax[2,i].text(ax[2,i].get_xlim()[0]+0.04*(ax[2,i].get_xlim()[1]-ax[2,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+6]+')',fontsize=18)
	ax[2,i].grid(ls='-.')
	ax[2,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[2,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[2,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[2,i].set_ylabel('')
		ax[2,i].set_yticklabels('')

#ax[2,1].text(-1,-32,r'$S_{na}$',fontsize=24)	
ax[2,1].set_title(r'prim nucl., no habit',fontsize=24,y=1.05)
varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	print(var)
	ax[3,i] = plot.plotRadar(mcRadar,var,ax[3,i])
	ax[3,i].set_ylim([0,-30])
	#ax[1,i].set_xlim([-3,0])
	ax[3,i].tick_params(axis='both',labelsize=16)
	ax[3,i].grid(ls='-.')
	ax[3,i].text(ax[3,i].get_xlim()[0]+0.04*(ax[3,i].get_xlim()[1]-ax[3,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+9]+')',fontsize=18)
	ax[3,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[3,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[3,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[3,i].set_ylabel('')
		ax[3,i].set_yticklabels('')



plt.savefig('plots/plot_second_nucleation/Agg_mode_habit_nohabit.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Wide_narrow_convolute1.pdf',bbox_inches='tight',dpi=1000)
plt.close()
quit()
'''
#- no habit

nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit0_xi500_10xi0_nz200_iwc3_nugam3.5_mugam0.5_dtc5_nrp2_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass1_ncl50_2nd_frag_PMD'
McRadar_Outname = '{freq}GHz_output_{mode}_convolute.nc'.format(freq=freq,mode=mode)
mcTable = xr.open_dataset(dataPath+expID+'/mass2fr.nc')

mcTable = mcTable.to_dataframe()
mcTable = mcTable.sort_values('sHeight')
selTime = mcTable['time'].max()
times = mcTable['time']
mcTable = mcTable[times==selTime]
mcTable['vel'] = -1*mcTable['vel']

atmoFile = np.loadtxt(dataPath+expID+'/atmo.dat')
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


mcRadar = xr.open_dataset(dataPath+expID+'/'+McRadar_Outname)
atmoPD.index.name='range'
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcRadar,method='nearest')
mcRadar = xr.merge([atmoReindex,mcRadar])

#- now second layer
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc3_nugam3.5_mugam0.5_dtc5_nrp2_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass1_ncl50_2nd_PMD'
McRadar_Outname = '{freq}GHz_output_{mode}_convolute.nc'.format(freq=freq,mode=mode)
mcTable = xr.open_dataset(dataPath+expID+'/mass2fr.nc')
mcTable = mcTable.to_dataframe()
mcTable = mcTable.sort_values('sHeight')
selTime = mcTable['time'].max()
times = mcTable['time']
mcTable = mcTable[times==selTime]
mcTable['vel'] = -1*mcTable['vel']

atmoFile = np.loadtxt(dataPath+expID+'/atmo.dat')
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
mcTableNew=mcTable.set_index('sHeight',drop=False)
mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
mcTableXR = mcTableNew.to_xarray(); atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
mcTableTmp1 = xr.merge([atmoReindex,mcTableXR])
mcTableTmp1 = mcTableTmp1.to_dataframe()


mcRadar1 = xr.open_dataset(dataPath+expID+'/'+McRadar_Outname)
atmoPD.index.name='range'
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcRadar1,method='nearest')
mcRadar1 = xr.merge([atmoReindex,mcRadar1])



varVec = ['mTot','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)
fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(10,10))
axes = [ax[0,0],ax[0,1],ax[1,0],ax[1,1]]
for i,var,ax in zip(range(len(varVec)),varVec,axes):	
	if var == 'mTot' or var == 'sNmono':
		ax = plot.plotPropSpecThesis(ax,mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax.set_xlim([-2,0])
		
	else:
		ax = plot.plotRadar(mcRadar,var,ax,diff=True,data1=mcRadar1)
		ax.axvline(x=0,c='k')
	ax.set_ylim([0,-30])
	ax.tick_params(axis='both',labelsize=16)
	ax.text(ax.get_xlim()[0]+0.04*(ax.get_xlim()[1]-ax.get_xlim()[0]),-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
	ax.grid(ls='-.')
	ax.axhline(y=-20,ls='--',color='r',linewidth=2)
	ax.axhline(y=-10,ls='--',color='r',linewidth=2)
	if (i == 0) or (i == 2):
		ax.set_ylabel('T [°C]',fontsize=18)
		ax.set_ylabel('T [°C]',fontsize=18)
	else:
		ax.set_ylabel('')
		ax.set_yticklabels('')
#ax[4,1].text(-1.0,-32,r'$S_{na}-S_{n,frag}$',fontsize=24)
axes[0].set_title(r'Habit - no habit',fontsize=24,x=1.25,y=1.1)

plt.tight_layout()
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_all_convoluted_habit_no_habit.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_all_convoluted_4diss1.pdf',bbox_inches='tight')
plt.close()

