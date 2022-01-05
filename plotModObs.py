import mcradar as mcr
from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt

def getNewNipySpectral():

    from matplotlib import cm
    from matplotlib.colors import ListedColormap

    numEnt = 15

    viridis = cm.get_cmap('nipy_spectral', 256)
    newcolors = viridis(np.linspace(0, 1, 256))

    colorSpace = np.linspace(198, 144, numEnt)/256
    colorTest=np.zeros((numEnt,4))
    colorTest[:,3] = 1
    colorTest[:,0]=colorSpace

    newcolors[- numEnt:, :] = colorTest
    newcmp = ListedColormap(newcolors)

    return newcmp

#freq = np.array([94.0e9])
freq = np.array([float(os.environ['freq'])*1e9])
print(freq)
#experimentID = os.environ['experiment']
#print(experimentID)
#inputPath = os.environ['MCexp']+'experiments/'+experimentID+'/'
inputPath = os.environ['MCexp']+'experiments/for_proposal_tests/1d_2_Dendrite_habit1_frag0_xi100_nz200_lwc0_iwc7_nugam2.0_dtc5_nrp2_rm10_rt0_vt3_coll_kern0_at2_stick3_dt1_geo3_h10-_ba500_domtop3500._atmo2_20190122_ar_top_0.05_rhi_105/'
print('loading the settings')
splitPath = inputPath.split('domtop')[1]
domTop = splitPath[0:4]
#inputPath = '/work/lvonterz/McSnow_habit/experiments/Jan_Niklas_frag_Leonie_setup/'
#-- load the settings of McSnow domain, as well as elevation you want to plot:
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',
                               elv=30, freq=freq,gridBaseArea=5.0,maxHeight=int(domTop),
                               ndgsVal=50,heightRes=50,scatSet={'mode':'full', 'safeTmatrix':False})

McRadar_Outname = os.environ['freq']+'GHz_output.nc'
output = xr.open_dataset(inputPath+McRadar_Outname)
#-- plot Prom appraisal
atmoFile = np.loadtxt(inputPath+'atmo.dat')
atmoPD = pd.DataFrame(atmoFile)
atmoPD = atmoPD[[0,2,6]]
atmoPD.columns = ['height','temp','RH']
atmoPD['temp'] = atmoPD['temp'] - 273.15

#atmoPD = atmoPD.set_index('height')
#Temp = xr.DataArray(np.empty(len(output.range)),dims='range',coords={'range':output.range.values})
Temp = np.empty(len(output.range))

for i, heightEdge0 in enumerate(dicSettings['heightRange']):
  heightEdge1 = heightEdge0 + dicSettings['heightRes']
  height = heightEdge0+dicSettings['heightRes']/2
  atmoTmp = atmoPD[(atmoPD['height']>heightEdge0) & (atmoPD['height']<=heightEdge1)].copy()
  Temp[i] = atmoTmp.temp.mean()

#TempXR = xr.DataArray(Temp,name='Temp',dims='range',coords={'range':output.range.values})
#output = xr.merge([output,TempXR])
outputNew = output.assign({'range':Temp})
outputNew = outputNew.rename({'range':'temperature [°C]'})
#read in concoluted file
outputConv = xr.open_dataset(inputPath+'specXR_convolution_ar_0.008.nc')
outputConv = outputConv.assign({'range':Temp})
outputConv = outputConv.rename({'range':'temperature [°C]'})

outputAll = xr.merge([outputNew['Ze_H_3.12e+01'],outputNew['Ze_V_3.12e+01'],outputNew['kdpInt_3.12e+01'],outputConv['spec_H_3.12e+01'],outputConv['spec_V_3.12e+01']])
#atmoPD = atmoPD.rename({'height','temp'})
print(outputAll)

# select observational sZDR to plot:

date = '20190122 15:02:00' 
#date = '20190130 13:30:20'
date2proc = pd.to_datetime(date)
year = date2proc.strftime('%Y'); month = date2proc.strftime('%m'); day = date2proc.strftime('%d'); hour = date2proc.strftime('%H') 
pathProcessed = '/data/obs/campaigns/tripex-pol/processed/'
level2_ID = '_tripex_pol_3fr_L2_mom.nc'
pol_level0_ID = '_tripex_pol_poldata_L0_spec_regridded_dealized.nc'

dataLV2 = xr.open_dataset(pathProcessed+'tripex_pol_level_2/'+year+month+day+level2_ID)
dataPolLV0 = xr.open_dataset(pathProcessed+'tripex_pol_level_0/'+year+'/'+month+'/'+day+'/'+year+month+day+'_'+hour+pol_level0_ID)
print(dataPolLV0)
dataPolsel = dataPolLV0.sel(time=date2proc)
dataPolsel = dataPolsel.where(dataPolsel['sSNR_H'] > 10.0)
dataLV2sel = dataLV2.sel(time=date2proc)
fig,ax = plt.subplots(ncols=2,figsize=(7.5,5.5),sharey=True)
p1 = ax[0].pcolormesh(dataPolsel.Vel2ZeroH.fillna(0),
                 dataLV2sel.ta.values,
                 dataPolsel.sZDR.values,
                 vmin=-0.5,vmax=3,cmap=getNewNipySpectral())
#cb = plt.colorbar(p1,ax=ax[0])
#cb.set_label('sZDR [dB]')
ax[0].grid(True,ls='-.')
ax[0].set_xlabel(r'Doppler vel [ms$^{-1}$]',fontsize=18)
ax[0].set_ylabel('Temp [°C]',fontsize=18)
ax[0].set_xlim(-1,0.05)
ax[0].set_ylim(-5,-20)
ax[0].set_title('Observation',fontsize=18)
ax[0].tick_params(axis='both',labelsize=16)

#plt.tight_layout()
specH = mcr.lin2db(outputAll['spec_H_3.12e+01'])
specV = mcr.lin2db(outputAll['spec_V_3.12e+01'])
specH = specH.where(specH > -40)
specV =  specV.where(specV > -40)
#dataSmooth = specTable.rolling(vel=10,min_periods=1,center=True).mean()            
sZDR = specH.rolling(vel=10,min_periods=1,center=True).mean() - specV.rolling(vel=10,min_periods=1,center=True).mean()
p2 = ax[1].pcolormesh(outputAll['vel']*np.sin(30/180.*np.pi),outputAll['temperature [°C]'],sZDR,vmin=-0.5, vmax=3,cmap=getNewNipySpectral())
#fig.subplots_adjust(bottom=0.8)
ax[1].grid(True,ls='-.')
ax[1].set_xlim(-1,0.05)
#ax[1].set_ylim(0,-25)
ax[1].set_title('Simulation',fontsize=18)
ax[1].set_xlabel(r'Doppler vel [ms$^{-1}$]',fontsize=18)
ax[1].tick_params(axis='x',labelsize=16)
#plt.draw()
fig.subplots_adjust(bottom=0.35, top=0.9, left=0.15, right=0.95, wspace=0.12, hspace=0.12)
ax_cbar = fig.add_axes([0.15,0.15,0.8,0.05])
cb = fig.colorbar(p2,cax=ax_cbar,orientation='horizontal')
#pos1 = ax[0].get_position().get_points().flatten()
#pos2 = ax[1].get_position().get_points().flatten()
#ax_cbar = fig.add_axes([pos1[0], 0, pos2[2]-pos1[0], 0.05])
#cb = plt.colorbar(p2,ax=ax_cbar,orientation='horizontal')
cb.ax.tick_params(labelsize=16)
cb.set_label('sZDR [dB]',fontsize=18)

#plt.tight_layout()
plt.savefig(date2proc.strftime('%Y%m%d_%H%M%S')+'_sZDR_obs_vs_McSnow_cbar_horizontal.png')
plt.close()
