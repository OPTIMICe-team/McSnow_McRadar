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
plt.semilogy(D*1e6,D/L)
plt.semilogy(D2*1e6,D2/L2,c='C1')
plt.semilogy(D3*1e6,D3/D3,c='C1')
plt.xticks(np.arange(0,110,10))
plt.xlabel(r'Dmax [$\mu$m]')
plt.ylabel('aspect ratio')
plt.grid(True,which='both')
plt.savefig('phi_Dmax_10.png')
plt.show()
'''

'''
test 2nd mode mix
'''
dataPath = '/project/meteo/work/L.Terzi/McSnowoutput/habit/second_nucleation/'
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
#- narrow minus frag

nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2nd_mix20mum'.format(iwc,nugam,mugam,nrp0)
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


fig,ax = plt.subplots(nrows=2,ncols=4,figsize=(17,10))
varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[0,i] = plot.plotPropSpecThesis(ax[0,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[0,i].set_xlim([-2,0])
	else:
		ax[0,i] = plot.plotRadar(mcRadar,var,ax[0,i],diff=True,data1=mcRadar1)
		ax[0,i].axvline(x=0,c='k')
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
#ax[4,1].text(-1.0,-32,r'$S_{na}-S_{n,frag}$',fontsize=24)
ax[0,1].set_title(r'$S_{n,mix}-S_{na}$',fontsize=24,x=1.25,y=1.1)

#- now wide
nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2nd_mix20mum'.format(iwc,nugam,mugam,nrp0)
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


#fig,ax = plt.subplots(nrows=5,ncols=4,figsize=(17,22))
varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[1,i] = plot.plotPropSpecThesis(ax[1,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[1,i].set_xlim([-2,0])
	else:
		ax[1,i] = plot.plotRadar(mcRadar,var,ax[1,i],diff=True,data1=mcRadar1)
		ax[1,i].axvline(x=0,c='k')
	ax[1,i].set_ylim([0,-30])
	ax[1,i].tick_params(axis='both',labelsize=16)
	ax[1,i].text(ax[0,i].get_xlim()[0]+0.04*(ax[0,i].get_xlim()[1]-ax[0,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+4]+')',fontsize=18)
	ax[1,i].grid(ls='-.')
	ax[1,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[1,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[1,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[1,i].set_ylabel('')
		ax[1,i].set_yticklabels('')
#ax[4,1].text(-1.0,-32,r'$S_{na}-S_{n,frag}$',fontsize=24)
ax[1,1].set_title(r'$S_{w,mix}-S_{wi}$',fontsize=24,x=1.25,y=1.1)

plt.tight_layout()
plt.savefig('/project/meteo/work/L.Terzi/McSnowoutput/habit/second_nucleation/plots/Diff_all_convoluted_2nd_mix.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_all_convoluted_4diss1.pdf',bbox_inches='tight')
plt.close()
quit()


'''
First: only plot one setup to show how variables generally look like.
'''

dataPath = '/data/optimice/McSnowoutput/habit/second_nucleation/'
nugam = 3.5; mugam = 0.5; IWC = 3; nrp=2
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2nd_prim'.format(IWC,nugam,mugam,nrp)
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


#fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,9))

#fig = plt.figure(figsize=(15,20))
#gs_top1 = mpl.gridspec.GridSpec(1, 3)
#gs.update(wspace=0.25,hspace=0.4)        
#ax00 = fig.add_subplot(gs[0, 0]) 
#ax01 = fig.add_subplot(gs[0, 1]) 
#ax02 = fig.add_subplot(gs[0, 2]) 
#ax10 = fig.add_subplot(gs[1, 0]) 
#ax11 = fig.add_subplot(gs[1, 1]) 
#ax12 = fig.add_subplot(gs[1, 2]) 
#ax20 = fig.add_subplot(gs[2, 0]) 
#ax21 = fig.add_subplot(gs[2, 1]) 
#ax22 = fig.add_subplot(gs[2, 2]) 
#ax30 = fig.add_subplot(gs[3, 0]) 
#ax31 = fig.add_subplot(gs[3, 1]) 
#ax32 = fig.add_subplot(gs[3, 2]) 


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
ax[0,1].set_title(r'$S_{wi}$',fontsize=24,y=1.05)
varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
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
nugam = 10; mugam = 10; IWC = 1; nrp=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2nd_mix20mum'.format(IWC,nugam,mugam,nrp)
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
ax[2,1].set_title(r'$S_{na}$',fontsize=24,y=1.05)
varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
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



plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Wide_narrow_second_mode_mix20mum.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Wide_narrow_convolute1.pdf',bbox_inches='tight',dpi=1000)
plt.close()
quit()


'''
Second: plot difference plot.
'''
'''
dataPath = '/data/optimice/McSnowoutput/habit/second_nucleation/'
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2ndar_Dgamma'.format(iwc,nugam,mugam,nrp0)
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


#fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,9))

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
varVec = ['dia','dia','sNmono']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	print(var)
	ax[0,i] = plot.plotPropSpecThesis(ax[0,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
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
quit()

varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	ax[1,i] = plot.plotRadar(mcRadar,var,ax[1,i],diff=True,data1=mcRadar1)
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


nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2ndar_Dgamma'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['number_conc','dia','sNmono']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	ax[2,i] = plot.plotPropSpecThesis(ax[2,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
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

varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	ax[3,i] = plot.plotRadar(mcRadar,var,ax[3,i],diff=True,data1=mcRadar1)
	ax[3,i].set_ylim([0,-30])
	#ax1,i].set_xlim([-3,0])
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

plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_2nd_mode_frag_convoluted.png',bbox_inches='tight')
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_2nd_mode_frag_convoluted.pdf',bbox_inches='tight')
plt.close()
'''	
'''
Third: plot single difference plot.
'''
'''
dataPath = '/data/optimice/McSnowoutput/habit/second_nucleation/'
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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

#- now wide psd
nugam = 3.5; mugam = 0.5 ; iwc=1; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
'''
'''
Second: plot difference plot.
'''
'''
dataPath = '/data/optimice/McSnowoutput/habit/second_nucleation/'
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2ndar_Dgamma'.format(iwc,nugam,mugam,nrp0)
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


#fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,9))

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
	ax[0,i] = plot.plotPropSpecThesis(ax[0,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
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

varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	ax[1,i] = plot.plotRadar(mcRadar,var,ax[1,i],diff=True,data1=mcRadar1)
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


nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2ndar_Dgamma'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['number_conc','dia','sNmono']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	ax[2,i] = plot.plotPropSpecThesis(ax[2,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
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

varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	ax[3,i] = plot.plotRadar(mcRadar,var,ax[3,i],diff=True,data1=mcRadar1)
	ax[3,i].set_ylim([0,-30])
	#ax1,i].set_xlim([-3,0])
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

plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_2nd_mode_frag_convoluted.png',bbox_inches='tight')
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_2nd_mode_frag_convoluted.pdf',bbox_inches='tight')
plt.close()

#fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,9))

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
	ax[0,i] = plot.plotPropSpecThesis(ax[0,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
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

varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	ax[1,i] = plot.plotRadar(mcRadar,var,ax[1,i],diff=True,data1=mcRadar1)
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


nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['number_conc','dia','sNmono']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	ax[2,i] = plot.plotPropSpecThesis(ax[2,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
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

varVec = ['sZe','DWR','ZeMDV']

for i, var in enumerate(varVec):
	ax[3,i] = plot.plotRadar(mcRadar,var,ax[3,i],diff=True,data1=mcRadar1)
	ax[3,i].set_ylim([0,-30])
	#ax1,i].set_xlim([-3,0])
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

plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_wide_IWC3vs1_nrp2vs1_convoluted.png',bbox_inches='tight')
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_wide_IWC3vs1_nrp2vs1_convoluted.pdf',bbox_inches='tight')
plt.close()
'''	
		
	
	
'''
Fourth: plot difference plot wih only dD, dNmono, dZe, dDWR.
'''

dataPath = '/data/optimice/McSnowoutput/habit/second_nucleation/'
freq='9.6_35.5_94.0'
mode='SSRGA-Rayleigh' 
nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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

#- now narrow
nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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


#fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,9))

fig,ax = plt.subplots(nrows=5,ncols=4,figsize=(17,22))

#fig,ax = plt.subplots()
varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[0,i] = plot.plotPropSpecThesis(ax[0,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[0,i].set_xlim([-2,0])
	else:
		ax[0,i] = plot.plotRadar(mcRadar,var,ax[0,i],diff=True,data1=mcRadar1)
		ax[0,i].axvline(x=0,c='k')
	
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
#ax[0,1].text(-1.0,-32,r'$S_{wi}-S_{na}$',fontsize=24)
ax[0,1].set_title(r'$S_{wi}-S_{na}$',fontsize=24,x=1.25,y=1.1)
	#ax[0,i].set_xlabel('')
	#ax[0,i].set_xticklabels('')

#- wide minus 2nd mode

nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass1_ncl50'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[1,i] = plot.plotPropSpecThesis(ax[1,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[1,i].set_xlim([-2,0])
	else:
		ax[1,i] = plot.plotRadar(mcRadar,var,ax[1,i],diff=True,data1=mcRadar1)
		ax[1,i].axvline(x=0,c='k')
	ax[1,i].set_ylim([0,-30])
	ax[1,i].tick_params(axis='both',labelsize=16)
	ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+4]+')',fontsize=18)
	ax[1,i].grid(ls='-.')
	ax[1,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[1,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[1,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[1,i].set_ylabel('')
		ax[1,i].set_yticklabels('')

#ax[1,1].text(-1.0,-32,r'$S_{wi}-S_{w,2nd}$',fontsize=24)
ax[1,1].set_title(r'$S_{wi}-S_{w,2nd}$',fontsize=24,x=1.25,y=1.1)
#- narrow minus frag

nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass1_ncl50'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[2,i] = plot.plotPropSpecThesis(ax[2,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[2,i].set_xlim([-2,0])
	else:
		ax[2,i] = plot.plotRadar(mcRadar,var,ax[2,i],diff=True,data1=mcRadar1)
		ax[2,i].axvline(x=0,c='k')
	ax[2,i].set_ylim([0,-30])
	ax[2,i].tick_params(axis='both',labelsize=16)
	ax[2,i].text(ax[2,i].get_xlim()[0]+0.04*(ax[2,i].get_xlim()[1]-ax[2,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+8]+')',fontsize=18)
	ax[2,i].grid(ls='-.')
	ax[2,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[2,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[2,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[2,i].set_ylabel('')
		ax[2,i].set_yticklabels('')
#ax[2,1].text(-1.0,-32,r'$S_{na}-S_{n,2nd}$',fontsize=24)
ax[2,1].set_title(r'$S_{na}-S_{n,2nd}$',fontsize=24,x=1.25,y=1.1)
#- wide minus frag

nugam = 3.5; mugam = 0.5 ; iwc=3; nrp0=2
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2ndar_Dgamma'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[3,i] = plot.plotPropSpecThesis(ax[3,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[3,i].set_xlim([-2,0])
	else:
		ax[3,i] = plot.plotRadar(mcRadar,var,ax[3,i],diff=True,data1=mcRadar1)
		ax[3,i].axvline(x=0,c='k')
	ax[3,i].set_ylim([0,-30])
	ax[3,i].tick_params(axis='both',labelsize=16)
	ax[3,i].text(ax[3,i].get_xlim()[0]+0.04*(ax[3,i].get_xlim()[1]-ax[3,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+12]+')',fontsize=18)
	ax[3,i].grid(ls='-.')
	ax[3,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[3,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[3,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[3,i].set_ylabel('')
		ax[3,i].set_yticklabels('')

#ax[3,1].text(-1.0,-32,r'$S_{wi}-S_{w,frag}$',fontsize=24)
ax[3,1].set_title(r'$S_{wi}-S_{w,frag}$',fontsize=24,x=1.25,y=1.1)
#- narrow minus frag

nugam = 10; mugam = 10 ; iwc=1; nrp0=1
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass0_ncl0'.format(iwc,nugam,mugam,nrp0)
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
expID ='1d_habit_habit1_xi500_10xi0_nz200_iwc{0}_nugam{1}_mugam{2}_dtc5_nrp{3}_vt3_coll_kern1_at2_stick2_colleffi1_dt1_bndtype2_ba500_domtop5000._atmo1_ssat50_nuclLayer_nh2300_3000_nclmass100_ncl50_2ndar_Dgamma'.format(iwc,nugam,mugam,nrp0)
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



varVec = ['dia','sNmono','Ze','DWR']
#velBins = np.linspace(-3,0,100)
velBins = np.arange(-3, 3, 6/512)

for i,var in enumerate(varVec):	
	if var == 'dia' or var == 'sNmono':
		ax[4,i] = plot.plotPropSpecThesis(ax[4,i],mcRadar.range,np.diff(mcRadar.range)[0],mcTableTmp,velBins,var,diff=True,mcTable1=mcTableTmp1)#,cblabel=r'conc [m$^{-3}$]')#,zoom=[-30,-10])#,zoom=[0,0.2])
		#ax[0,i].set_xticks([-15,-10,-5,0])
		ax[4,i].set_xlim([-2,0])
	else:
		ax[4,i] = plot.plotRadar(mcRadar,var,ax[4,i],diff=True,data1=mcRadar1)
		ax[4,i].axvline(x=0,c='k')
	ax[4,i].set_ylim([0,-30])
	ax[4,i].tick_params(axis='both',labelsize=16)
	ax[4,i].text(ax[4,i].get_xlim()[0]+0.04*(ax[4,i].get_xlim()[1]-ax[4,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+16]+')',fontsize=18)
	ax[4,i].grid(ls='-.')
	ax[4,i].axhline(y=-20,ls='--',color='r',linewidth=2)
	ax[4,i].axhline(y=-10,ls='--',color='r',linewidth=2)
	if i == 0:
		ax[4,i].set_ylabel('T [°C]',fontsize=18)
	else:
		ax[4,i].set_ylabel('')
		ax[4,i].set_yticklabels('')
#ax[4,1].text(-1.0,-32,r'$S_{na}-S_{n,frag}$',fontsize=24)
ax[4,1].set_title(r'$S_{na}-S_{n,frag}$',fontsize=24,x=1.25,y=1.1)

plt.tight_layout()
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_all_convoluted_4diss1.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/Diff_all_convoluted_4diss1.pdf',bbox_inches='tight')
plt.close()

	
	
	
	
