import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


modelFile = '../20140201_hyytiala_ecmwf.nc' # file to be read in
data = xr.open_dataset(modelFile)

datasel = data.sel(time=data.time[0]) # select correct time here
datasel = datasel[['height','temperature','pressure','q']].sortby('level').to_dataframe()
datasel = datasel[['height','temperature','pressure','q']]
np.savetxt('mcsnow/input/20140201_hyytiala_ecmwf.txt', datasel.values[1:-1], fmt='%f') # needs to be saved to mcsnow/input/

