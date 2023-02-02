from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import string

dataPath = '/data/optimice/McSnowoutput/habit/trajectories/'
exp = '1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo2_RHi105_nucleating_10i_interpPPTrue'
output = xr.open_dataset(dataPath+exp+'/9.6GHz_output_full_0_350_singleParticle.nc')
print(output)
#quit()
atmoFile = np.loadtxt(dataPath+exp+'/atmo.dat')
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
print(atmoPD)

atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(output,method='nearest')
output = xr.merge([atmoReindex,output])
output['sZDR'] = 10*np.log10(output['sZeH_3.12e+01']) - 10*np.log10(output['sZeV_3.12e+01'])
datasZDRmax = pd.read_csv('sZDRmax_median_DWRclass_2.txt',delimiter=' ',header=0)
fig,ax = plt.subplots(ncols=2,figsize=(10,4),sharey=True)
'''
vmin=0;vmax=90
outputSel = output.sel(sMult=slice(vmin,vmax))
outputSel['ZeH'] = outputSel['sZeH_3.12e+01'].sum(dim='sMult')
outputSel['ZeV'] = outputSel['sZeV_3.12e+01'].sum(dim='sMult')
outputSel['ZDR'] = 10*np.log10(outputSel.ZeH)-10*np.log10(outputSel.ZeV)
sZDRmax = outputSel.sZDR.max(dim='sMult')
sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
ZDR = outputSel.ZDR.rolling(sHeight=20,center=True).mean()
ls1=ax[0].plot(ZDR,outputSel.Temp,c='#66CCEE',label='cold',lw=2)#c='#4477AA'
ls2=ax[1].plot(sZDRmax,outputSel.Temp,ls='-',c='#66CCEE',label=r'sZDR$_{\rm max}$ cold',lw=2)#c='#66CCEE',
'''
output['ZeH'] = output['sZeH_3.12e+01'].sum(dim='sMult')
output['ZeV'] = output['sZeV_3.12e+01'].sum(dim='sMult')
output['ZDR'] = 10*np.log10(output.ZeH)-10*np.log10(output.ZeV)
sZDRmax = output.sZDR.max(dim='sMult')
sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
ZDR = output.ZDR.rolling(sHeight=20,center=True).mean()
ls3=ax[0].plot(ZDR,output.Temp,lw=2,c='#228833',label='all')
ls4=ax[1].plot(sZDRmax,output.Temp,lw=2,c='#66CCEE',label=r'all')#c='#CCBB44'


vmin=180;vmax=350
outputSel = output.sel(sMult=slice(vmin,vmax))
outputSel['ZeH'] = outputSel['sZeH_3.12e+01'].sum(dim='sMult')
outputSel['ZeV'] = outputSel['sZeV_3.12e+01'].sum(dim='sMult')
outputSel['ZDR'] = 10*np.log10(outputSel.ZeH)-10*np.log10(outputSel.ZeV)
sZDRmax = outputSel.sZDR.max(dim='sMult')
sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
ZDR = outputSel.ZDR.rolling(sHeight=20,center=True).mean()
ls5=ax[0].plot(ZDR,outputSel.Temp,c='#EE6677',label='DGL',lw=2)
ls6=ax[1].plot(sZDRmax,outputSel.Temp,ls='-',c='#EE6677',label=r'sZDR$_{\rm max}$ DGL',lw=2)#c='#AA3377'

ax1 = ax[1].twiny()
ls7=ax1.plot(datasZDRmax.med,datasZDRmax.Temp,c='k',lw=2,label=r'obs')
ax1.set_xlim([0,2.2])
ax1.xaxis.set_ticks_position("bottom")
ax1.xaxis.set_label_position("bottom")
# Offset the twin axis below the host
ax1.spines["bottom"].set_position(("axes", -0.25))
# Turn on the frame for the twin axis, but then hide all 
# but the bottom spine
ax1.set_frame_on(True)
ax1.patch.set_visible(False)
for sp in ax1.spines.values():
	sp.set_visible(False)
ax1.spines["bottom"].set_visible(True)
ax1.tick_params(axis='both',labelsize=16)
ax1.set_xlabel(r'observed sZDR$_{\rm max}$ [dB]',fontsize=18)

#ls = ls1+ls2+ls3+ls4+ls5+ls6+ls7
#labs = [l.get_label() for l in ls]
#ax[0].legend()
ls = ls4+ls5+ls7
labs = [l.get_label() for l in ls]
ax[1].legend(ls,labs,ncol=1,fontsize=12)
ax[0].set_ylabel('T [°C]',fontsize=18)

ax[0].set_xlabel('simulated ZDR [dB]',fontsize=18)
ax[1].set_xlabel(r'simulated sZDR$_{\rm max}$ [dB]',fontsize=18)
for i,a in enumerate(ax):
	a.set_xlim([0,7])
	a.set_ylim([0,-30])
	a.axhline(y=-10,ls='--',lw=2,c='r')
	a.axhline(y=-20,ls='--',lw=2,c='r')
	a.tick_params(axis='both',labelsize=16)
	a.set_yticks([0,-5,-10,-15,-20,-25,-30])
	a.set_xticks([0,1,2,3,4,5,6,7])
	#a.text(0.1,-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
	a.grid()
#ax[0].set_zorder(1)
#ax[0].legend(ls,labs,bbox_to_anchor=(-0,0.8),loc='lower left',fontsize=12,ncol=3)
#plt.tight_layout()
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/ZDR_int4.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/ZDR_int3.pdf',bbox_inches='tight')
plt.close()
quit()
