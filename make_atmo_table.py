import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
date = pd.to_datetime('20181204 02:55:00')
year = date.strftime('%Y'); month=date.strftime('%m'); day=date.strftime('%d'); hour=date.strftime('%H');minu=date.strftime('%M');sec=date.strftime('%S')
cnPath = '/data/data_hatpro/jue/cloudnet/juelich/processed/categorize/'+year+'/'
data = xr.open_dataset(cnPath+year+month+day+'_juelich_categorize.nc')
datasel = data.sel(time=date,method='nearest')
print(datasel)
datasel = datasel[['temperature','pressure','specific_humidity']]
#dataInt = datasel.interp(model_height=dataLV2.range.values[::-1])
dataPD = datasel.to_dataframe()
dataPD['model_height'] = dataPD.index
dataPD = dataPD[['model_height','temperature','pressure','specific_humidity']]
#save as txt because that is what is read in in McSnow
np.savetxt('mcsnow/input/'+year+month+day+'_'+hour+minu+sec+'_CN_atmo_profiles.txt', dataPD.values[::-1], fmt='%f')
print(dataPD)
