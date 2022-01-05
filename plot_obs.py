'''
this script is meant to plot the observations in a specific time slot in order to compare the McSnow simulations to that. This plots (pol-) moments and (pol-)spectra
'''


import numpy as np
import xarray as xr
import glob 
import pandas as pd
import matplotlib.pyplot as plt
import plotRoutines as plotRout
date = '20190103'
dateStart = pd.to_datetime(date+' 22:00:00'); dateEnd = pd.to_datetime(date+' 23:00:00')
dataLV2Path = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_2/'
dataPolPath = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_1/'
dataLV0Path = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_0/'
dataPolLV0Path = '/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_0/'
# read in files 
dataLV2 = xr.open_dataset(dataLV2Path+date+'_tripex_pol_3fr_L2_mom.nc') # LV2 moments
dataPol = xr.open_dataset(dataPolPath+date+'_tripex_pol_poldata_L1_mom.nc') # polarimetric moments
dataLV0List = glob.glob(dataLV0Path+dateStart.strftime('%Y')+'/'+dateStart.strftime('%m')+'/'+dateStart.strftime('%d')+'/'+date+'*_tripex_pol_3fr_spec_filtered_regridded.nc') # LV0, so spectra
dataLV0 = xr.open_mfdataset(dataLV0List)#[int(dateStart.strftime('%H')):int(dateEnd.strftime('%H'))+1])
dataLV0.time.attrs['units']='seconds since 1970-01-01 00:00:00 UTC' # time attribute is weird, and python does not recognize it
dataLV0 = xr.decode_cf(dataLV0)
dataPolLV0List = glob.glob(dataLV0Path+dateStart.strftime('%Y')+'/'+dateStart.strftime('%m')+'/'+dateStart.strftime('%d')+'/'+date+'*_tripex_pol_poldata_L0_spec_regridded_dealized.nc') #polarimetric spectra
dataPolLV0 = xr.open_mfdataset(dataPolLV0List)#[int(dateStart.strftime('%H')):int(dateEnd.strftime('%H'))+1])
#print(dataPolLV0)
#quit()
####################
# plot moments timeseries
####################
#plotRout.plotMomentsObs(dataLV2.sel(time=slice(dateStart,dateEnd)),dataPol.sel(time=slice(dateStart,dateEnd)),'20190122/',dateStart.strftime('%Y%m%d_%H%M%S')+'_'+dateEnd.strftime('%Y%m%d_%H%M%S')+'_moments')


####################
# calculate and plot 5 min mean
###################
polMean = dataPol.sel(time=slice(dateStart,dateEnd)).resample(time='5min').mean()
LV2mean = dataLV2.sel(time=slice(dateStart,dateEnd)).resample(time='5min').mean()
LV0mean = dataLV0.sel(time=slice(dateStart,dateEnd)).resample(time='5min').mean()
LV0Polmean = dataPolLV0[['sZDR','Vel2ZeroH']].sel(time=slice(dateStart,dateEnd)).resample(time='5min').mean()
range_in_cloud = LV2mean.range.where(~np.isnan(LV2mean.Ka_DBZ_H)) # calculate all ranges that have a cloud
cloudtop = range_in_cloud.max(dim='range') # now find the maximum range of the cloud, alas the cloudtop range
#cloudtop.to_netcdf(date+'/cloudtop_range.nc')
#range_in_cloud = LV2mean.range.where(~np.isnan(LV2mean.Ka_DBZ_H.sel(range=slice(0,2000)))) # calculate all ranges that have a cloud
#cloudtop = range_in_cloud.max(dim='range')
#print(LV2mean.Ka_DBZ_H.sel(range=cloudtop,method='nearest').sel(time='20190103 22:50:00',method='nearest'))#.sel(range=1300,method='nearest')
#LV2mean.Ka_DBZ_H.plot(y='range')
#cloudtop.plot(x='time',c='k')
#plt.grid()
#plt.ylim([250,1500])
#plt.show()
#quit()
'''
LV2mean.Ka_DBZ_H.plot(y='range')
cloudtop.plot(x='time',c='k')
plt.tight_layout()
plt.savefig(date+'/test_cloudtop.png')
plt.close()
#print(cloudtop)
cloudPD = cloudtop.to_dataframe()
cloudPD.to_csv(date+'/cloudtop_range.csv')
#cloudtop.to_dataset(date+'/cloudtop_range.nc')
#quit()

#print(LV2mean.Ka_DBZ_H.sel(range=4000,method='nearest'))
LV2mean.Ka_DBZ_H.sel(range=cloudtop,method='nearest').plot()
plt.grid()
plt.tight_layout()
plt.savefig(date+'/Ka_DBZ_H_CT_5min_mean.png')
plt.show()
quit()
'''
#plotRout.plotMomentsObs(LV2mean,polMean,date+'/',dateStart.strftime('%Y%m%d_%H%M%S')+'_'+dateEnd.strftime('%Y%m%d_%H%M%S')+'_moments')
#plotRout.plotDWRsObs(LV2mean.sel(time=slice(dateStart,dateEnd)),date+'/',dateStart.strftime('%Y%m%d_%H%M%S')+'_'+dateEnd.strftime('%Y%m%d_%H%M%S')+'_DWRs')
#quit()
#############################
# plot spectra and profiles:
# plot a spectra every 5 minutes, otherwise it takes too long 
# plot profiles of the observations, but 5 minute mean to get rid of KDP uncertainty
############################
plotRout.plotSpectraObs(LV0Polmean,LV0mean,date+'/',ylim=([0,1500]))#.resample(time='5min').nearest(tolerance='4s')#.resample(time='5min').nearest(tolerance='6min')
plotRout.plotProfilesObs(LV2mean,polMean,date+'/',ylim=([0,1500]))
plotRout.plotProfilesObsDWR(LV2mean,date+'/',ylim=([0,1500]))
    
