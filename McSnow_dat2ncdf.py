'''
this script creates ncdf data from McSnows mass2fr..dat file with the original and some postprocessed variables (like temporal averages)
'''
import numpy as np
import pandas as pd
import xarray as xr
import os


experimentID = os.environ['experiment']
print(experimentID)
inputPath = os.environ['MCexp']+'experiments/'+experimentID+'/'
print(inputPath)
#inputPath = '/work/lvonterz/McSnow_habit/experiments/Jan_Niklas_frag_Leonie_setup/'

data = pd.read_csv(inputPath+'mass2fr.dat',header=0)#,skiprows=1, # header=0
                   #names=['time','mTot','sHeight','vel','dia','area','sMice','sVice','sPhi','sRhoIce','sNmono','sMrime','sVrime','sMult', 'sRho_tot'])
data.columns = data.columns.str.replace(' ','')

# converting pandas dataframe into xarray to add attrs and save into netcdf
dataXR = data.to_xarray()
dataXR['time'].attrs={'units':'min','long_name':'minutes since start of simulation'}
dataXR['mTot'].attrs={'units':'kg','long_name':'total mass of superparticle (ice+rime+liquid)'}
dataXR['mTot'].attrs={'units':'kg','long_name':'total mass of superparticle (ice+rime+liquid)'}
dataXR['sHeight'].attrs={'units':'m','long_name':'height of superparticle'}
dataXR['vel'].attrs={'units':'m/s','long_name':'terminal velocity of superparticle'}
dataXR['dia'].attrs={'units':'m','long_name':'maximum dimension of superparticle'}
dataXR['area'].attrs={'units':'m^2','long_name':'projected area of superparticle'}
dataXR['sMice'].attrs={'units':'kg','long_name':'ice mass of superparticle'}
dataXR['sVice'].attrs={'units':'m^3','long_name':'ice volume of superparticle'}
dataXR['sPhi'].attrs={'units':'','long_name':'aspect ratio of superparticle'}
dataXR['sRhoIce'].attrs={'units':'kg/m^3','long_name':'density of ice of superparticle'}
#dataXR['igf'].attrs={'units':'','long_name':'inherent growth function'}
dataXR['sNmono'].attrs={'units':'','long_name':'number of monomers of superparticle'}
dataXR['sMrime'].attrs={'units':'kg','long_name':'rime mass of superparticle'}
dataXR['sVrime'].attrs={'units':'m^3','long_name':'rime volume of superparticle'}
dataXR['sMult'].attrs={'units':'','long_name':'multiplicity of superparticle'}
#dataXR['sMliqu'].attrs={'units':'kg','long_name':'liquid water mass of superparticle'}
#dataXR['sMmelt'].attrs={'units':'kg','long_name':'melted or frozen water mass of superparticle'}
dataXR['sRho_tot'].attrs={'units':'kg/m^3','long_name':'total density of superparticle'}
dataXR = dataXR.fillna(0) 

#save as nc

dataXR.to_netcdf(inputPath+'mass2fr.nc')
print('saved at ', inputPath+'mass2fr.nc')
#remove .dat file
os.remove(inputPath+'mass2fr.dat')
