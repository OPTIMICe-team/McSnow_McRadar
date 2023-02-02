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
def estPercentil(medFreq, cumFreq, centDwrBins, perc):
    
    func = intp.interp1d(cumFreq, centDwrBins)
    intMedian = func(medFreq)
    
    return intMedian

def getPercentil(cumFreq, centDwrBins, perc):
    
    medFreq = cumFreq[:,-1]*perc # this is the number of profiles at median
    medDwr = np.ones_like(medFreq)*np.nan
    
    for i, tmpFreq in enumerate(medFreq):
        try:
            #plt.plot(centDwrBins,cumFreq[i])
            #print(tmpFreq)
            medDwr[i] = estPercentil(tmpFreq, cumFreq[i], centDwrBins, perc)
            #plt.axhline(y=tmpFreq)
            #plt.axvline(x=medDwr[i],c='r')
            #plt.savefig('prof_'+str(i)+'.png')
            #plt.close()
        except:
            medDwr[i] = np.nan #print(medDwr[i])
        
    return medDwr,medFreq
def getNewNipySpectral(r=False):

    numEnt = 15
    if r == True:
      viridis = cm.get_cmap('nipy_spectral_r', 256)
      newcolors = viridis(np.linspace(0, 1, 256))
      colorSpace = np.linspace(144, 198, numEnt)/256
      colorTest=np.zeros((numEnt,4))
      colorTest[:,3] = 1
      colorTest[:,0]=colorSpace
      newcolors[0:15, :] = colorTest
    else:
      viridis = cm.get_cmap('nipy_spectral', 256)
      colorSpace = np.linspace(198, 144, numEnt)/(256)
      newcolors = viridis(np.linspace(0, 1, 256))
      colorTest=np.zeros((numEnt,4))
      colorTest[:,3] = 1
      colorTest[:,0]=colorSpace
      newcolors[- numEnt:, :] = colorTest
    newcmp = ListedColormap(newcolors)

    return newcmp 
dataPath = '/project/meteo/work/L.Terzi/McSnowoutput/habit/trajectories/'

# plot depogrowth in dependence of Tini
RHi_array = ['10','50','100','150','200']
var2plot = ['dia','mTot','sPhi','sRho_tot']
unitvar = ['[m]','[kg]','',r'[kgm$^{-3}$]']
varNamevec = [r'D$_{\mathrm{max}}$','m',r'$\phi$',r'$\rho$']
IGF = 1
expID = '1d_habit_trajectories_habit1_IGF{0}_xi1_nz200_iwc1_nugam10_mugam10_dtc5_nrp1_vt3_coll_kern0_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo1_ssat'.format(IGF)
fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(10,7))
ax_vec = ([0,0],[0,1],[1,0],[1,1])
colors=['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677']
for i_var,i_ax,var,unit,varName in zip(range(len(unitvar)),ax_vec,var2plot,unitvar,varNamevec):
	print(var)
	a=ax[i_ax[0],i_ax[1]]
	
	for i,RHi in enumerate(RHi_array):
		experiment = '{dataPath}{expID}{ssat}'.format(dataPath=dataPath,expID=expID,ssat=RHi)
		print(RHi)
		if os.path.exists(experiment+'/mass2fr.nc'):
			mcTable = xr.open_dataset(experiment+'/mass2fr.nc')
			mcTable = mcTable.to_dataframe()
		else:
			print(experiment+'does not exist, check path!')
			break
		
		
		atmoFile = np.loadtxt(experiment+'/atmo.dat')
		height = atmoFile[:,0]
		Temp = atmoFile[:,2] -273.15
		atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
		atmoPD.index.name='sHeight'
		mcTableNew=mcTable.set_index('sHeight',drop=False)
		mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
		mcTableXR = mcTableNew.to_xarray()
		atmoXR = atmoPD.to_xarray()
		atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
		mcTableTemp = xr.merge([atmoReindex,mcTableXR])
		#print(mcTableTemp)
		mcTableTemp = mcTableTemp.to_dataframe()
		
		mcTableTmp = mcTableTemp
		if var == 'sRho_tot':
			ylog=False
		else:
			ylog=True
		
		
		a=plot.plotInitempVar(mcTableTmp,a,experiment,var,'Temp',unit,varName,ylog=ylog,color=colors[i])#,zoom=[-20,-10])#,zoom=[0,0.2])
		a.axvline(x=-10,c='r',lw=2,ls='--')
		a.axvline(x=-20,c='r',lw=2,ls='--')
		#a.axvline(x=-21.25,c='k',lw=2,ls='--')
	if var=='sRho_tot':
		a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),a.get_ylim()[1]-0.1*(a.get_ylim()[1]-a.get_ylim()[0]),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		a.set_yticks([400,500,600,700,800,900])
	elif var=='dia':
		if IGF == 1:
		    a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),10**(-2),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		else:
		    a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),8*10**(-3),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
	elif var=='mTot':
		a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),2*10**(-6),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		a.set_yticks([10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5])
		y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
		a.yaxis.set_minor_locator(y_minor)
		a.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		if IGF == 2:
		    #a.legend(fontsize=12,loc='upper right',markerscale=1.5)
		    a.legend(bbox_to_anchor=(0.1,0.8),loc='lower left',fontsize=12,ncol=2,markerscale=1.5)
	elif var=='sPhi':
		#a.set_ylim([3*10**-3,10**4])
		if IGF == 1:
		    a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),3*10**(1),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		    a.legend(fontsize=12,loc='upper right',markerscale=1.5)
		    #a.legend(bbox_to_anchor=(0.1,0.8),loc='lower left',fontsize=12,markerscale=1.5)
		else:
		    a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),3*10**(0),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		
		a.tick_params(axis='y', which='minor',direction='out')
		#a.set_yticks([10**-2,10**-1,10**0,10**1,10**2,10**3,10**4])
		y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
		a.yaxis.set_minor_locator(y_minor)
		a.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
figName = 'Ini_Temp_RHi_dependent_IGF{0}'.format(IGF)
plt.tight_layout()
plt.savefig(dataPath+figName+'.png',bbox_inches='tight')
plt.savefig(dataPath+figName+'.pdf',bbox_inches='tight')
plt.close()
quit()

# plot for Defense
RHi_array = ['101','105','110','115','120']
var2plot = ['sPhi']
unitvar = ['']
varNamevec = [r'$\phi$']
expID = '1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo2'
fig,ax_vec = plt.subplots(ncols=2,nrows=1,figsize=(10,4))
#ax_vec = ([0,0],[0,1],[1,0],[1,1])
colors=['#4477AA','#66CCEE','#228833','#CCBB44','#EE6677']
for i_var,a,var,unit,varName in zip(range(len(unitvar)),ax_vec,var2plot,unitvar,varNamevec):
	print(i_var)
	#a=ax[i_ax[0],i_ax[1]]
	a = ax_vec[1]
	
	
	for i,RHi in enumerate(RHi_array):
		experiment = '{dataPath}{expID}_RHi{RHi}_nucleating_10i_interpPPTrue'.format(dataPath=dataPath,
		                                                                         expID=expID,
		                                                                         RHi=RHi)
		print(RHi)
		if os.path.exists(experiment+'/mass2fr.nc'):
			mcTable = xr.open_dataset(experiment+'/mass2fr.nc')
			mcTable = mcTable.to_dataframe()
		else:
			print(experiment+'does not exist, check path!')
			break
		
		try:
 			minmax = os.environ['minmax']
 			vmin=int(minmax.split('_')[0]); vmax=int(minmax.split('_')[1])
		except:
			print('no minmax')
			minmax=False
		atmoFile = np.loadtxt(experiment+'/atmo.dat')
		height = atmoFile[:,0]
		Temp = atmoFile[:,2] -273.15
		atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
		atmoPD.index.name='sHeight'
		mcTableNew=mcTable.set_index('sHeight',drop=False)
		mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
		mcTableXR = mcTableNew.to_xarray()
		atmoXR = atmoPD.to_xarray()
		atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
		mcTableTemp = xr.merge([atmoReindex,mcTableXR])
		#print(mcTableTemp)
		mcTableTemp = mcTableTemp.to_dataframe()
		
		mcTableTmp = mcTableTemp
		if var == 'sRho_tot':
			ylog=False
		else:
			ylog=True
		a=plot.plotInitempVar(mcTableTmp,a,experiment,var,'Temp',unit,varName,ylog=ylog,color=colors[i])#,zoom=[-20,-10])#,zoom=[0,0.2])
		#a.axvline(x=-10,c='r',lw=2,ls='--')
		#a.axvline(x=-20,c='r',lw=2,ls='--')
		a.axhline(y=-21.25,c='k',lw=2,ls='--')
		a.axvline(x=1,ls=':',lw=2,c='k')
	if var=='sRho_tot':
		a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),a.get_ylim()[1]-0.1*(a.get_ylim()[1]-a.get_ylim()[0]),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		a.set_yticks([400,500,600,700,800,900])
	elif var=='dia':
		a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),3*10**(-1),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
	elif var=='mTot':
		a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),8*10**(-6),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		a.set_yticks([10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5])
		y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
		a.yaxis.set_minor_locator(y_minor)
		a.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	elif var=='sPhi':
		a.set_xlim([3*10**-3,10**4])
		a.set_ylim([-9,-31.5])
		#a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),2*10**(3),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
		a.legend(fontsize=12,loc='lower right',markerscale=1.5)
		a.tick_params(axis='x', which='minor',direction='out')
		a.set_xticks([10**-2,10**-1,10**0,10**1,10**2,10**3,10**4])
		y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
		a.xaxis.set_minor_locator(y_minor)
		a.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
#figName = 'Ini_Temp_RHi_dependent'
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/'+figName+'_fordefense.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/'+figName+'new.pdf',bbox_inches='tight')
#plt.close()
#quit()
'''
'''
# plot particle evolution
RHi = '105'
expID = '1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo2'
experiment = '{dataPath}{expID}_RHi{RHi}_nucleating_10i_interpPPTrue'.format(dataPath=dataPath,
		                                                                         expID=expID,
		                                                                         RHi=RHi)

mcTable = xr.open_dataset(experiment+'/mass2fr.nc')
mcTable = mcTable.to_dataframe()
atmoFile = np.loadtxt(experiment+'/atmo.dat')
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
mcTableNew=mcTable.set_index('sHeight',drop=False)
mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
mcTableXR = mcTableNew.to_xarray()
atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
mcTableTemp = xr.merge([atmoReindex,mcTableXR])
#print(mcTableTemp)
mcTableTemp = mcTableTemp.to_dataframe()
'''
fig = plt.figure(figsize=(15,4))
gs = mpl.gridspec.GridSpec(1, 3)        
gs.update(wspace=0.1,hspace=0.35)
ax00 = fig.add_subplot(gs[0]) 
ax01 = fig.add_subplot(gs[1]) 
ax10 = fig.add_subplot(gs[2]) 
#ax11 = fig.add_subplot(gs[1, 1]) 
ax = np.array([ax00,ax01,ax10])
'''
mcTableTmp = mcTableTemp
#vmin1 = 0; vmax1 = 180; vmin2=180; vmax2=350
#fig,ax = plt.subplots(ncols=2,nrows=1,figsize=(10,4))

#ax[0],s_m = plot.plot_var_particle_ID_temp(mcTableTmp,'dia','sMult',r'D$_{\rm max}$','[m]',log=True,ax=ax[0],cbar=False)#,zoom=[-20,-10])

#print('dia')
#ax[0].text(6*10**-3,-28,'(a)',fontsize=18)

#plot.plot_var_particle_ID_temp(mcTableTemp,inputPath,'dia_mum','sMult')
#ax[0],s_m =plot.plot_var_particle_ID_temp(mcTableTmp,'mTot','sMult','mass','[kg]',log=True,ax=ax[0],cbar=False,xlabel=True)#,zoom=[-20,-10])
#print('mTot')
#ax[0].text(1.5*10**-7,-28,'(a)',fontsize=18)
#ax[0].set_xticks([10**-12,10**-11,10**-10,10**-9,10**-8,10**-7,10**-6])
#y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
#ax[0].xaxis.set_minor_locator(y_minor)
#ax[0].xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

ax_vec[0],s_m = plot.plot_var_particle_ID_temp(mcTableTmp,'sPhi','sMult',r'$\phi$','',log=True,ax=ax_vec[0],cbar=False,xlabel=True)#,zoom=[-30,-10])#,zoom=[0,0.2])
#ax_vec[0].text(19,-28,'(b)',fontsize=18)
ax_vec[0].set_xticks([10**-2,10**-1,10**0,10**1])
y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
ax_vec[0].xaxis.set_minor_locator(y_minor)
ax_vec[0].xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

#ax[2],s_m =plot.plot_var_particle_ID_temp(mcTableTmp,'sRho_tot','sMult',r'$\rho$',r'[kgm$^{-3}$]',ax=ax[2],cbar=False,xlabel=False)#,zoom=[-20,-10])
#ax[2].text(860,-28,'(c)',fontsize=18)
#ax[2].set_xticks([900,700,500,300])
#print('sRhoice')

#for a in ax:
ax_vec[0].set_ylim([-9,-31.5])
ax_vec[0].axhline(y=-10,ls='--',lw=2,c='k')
ax_vec[0].axvline(x=1,ls=':',lw=2,c='k')
#ax_vec[0].axhline(y=-20,ls='--',lw=2,c='r')
ax_vec[0].set_yticks([-10,-15,-20,-25,-30])
	
cbar = fig.colorbar(s_m,ax=ax_vec[0],pad=0.025,aspect=30)
cbar.set_label('particle ID',fontsize=18)
cbar.set_ticks(np.arange(0,mcTableTmp['sMult'].max()+1,50))
cbar.ax.tick_params(labelsize=16)
# now plot two particles with black color for better illustration

'''
mcTableTmp = mcTableTmp.set_index('sMult',drop=False)
outputSel=mcTableTmp.loc[180]
outputSel=mcTableTmp.loc[10]
for var2plot,axes in zip(['dia','sPhi','sRho_tot'],[ax00,ax01,ax10]):
	prop2plot = outputSel[var2plot]
	axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,ls=':',label='Particle A')
	
outputSel=mcTableTmp.loc[100]
for var2plot,axes in zip(['dia','sPhi','sRho_tot'],[ax00,ax01,ax10]):
	prop2plot = outputSel[var2plot]
	axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,ls='-.',label='Particle B')

outputSel=mcTableTmp.loc[180]
for var2plot,axes in zip(['dia','sPhi','sRho_tot'],[ax00,ax01,ax10]):
	prop2plot = outputSel[var2plot]
	axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,label='Particle C')
	

outputSel=mcTableTmp.loc[250]
for var2plot,axes in zip(['dia','sPhi','sRho_tot'],[ax00,ax01,ax10]):
	prop2plot = outputSel[var2plot]
	axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,ls='--',label='Particle D')



ax00.legend(fontsize=12,ncol=4,bbox_to_anchor=(0,1,1.5,0.05))

'''
#ax = np.array([[ax00,ax01],[ax10,ax11]])
#plt.show()
plt.tight_layout()
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/particleEvolution_RHi105_newColor_fordefense2.png',bbox_inches='tight')

plt.close()
quit()
'''
'''
# plot integrated moments
McRadar_Outname = experiment+'/9.6GHz_output_full_0_350_singleParticle.nc'
output = xr.open_dataset(McRadar_Outname)
atmoFile = np.loadtxt(experiment+'/atmo.dat')
height = atmoFile[:,0]
Temp = atmoFile[:,2] -273.15
atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
atmoPD.index.name='sHeight'
print(atmoPD)

atmoXR = atmoPD.to_xarray()
atmoReindex = atmoXR.reindex_like(output,method='nearest')
output = xr.merge([atmoReindex,output])
output['sZDR'] = 10*np.log10(output['sZeH_3.12e+01']) - 10*np.log10(output['sZeV_3.12e+01'])
print(output)
#quit()

datasZDRmax = pd.read_csv('sZDRmax_median_DWRclass_2.txt',delimiter=' ',header=0)
print(datasZDRmax)
#quit()

fig,ax = plt.subplots(ncols=1,figsize=(6,5),sharey=True)

vmin=0;vmax=90
outputSel = output.sel(sMult=slice(vmin,vmax))
outputSel['ZeH'] = outputSel['sZeH_3.12e+01'].sum(dim='sMult')
outputSel['ZeV'] = outputSel['sZeV_3.12e+01'].sum(dim='sMult')
outputSel['ZDR'] = 10*np.log10(outputSel.ZeH)-10*np.log10(outputSel.ZeV)
sZDRmax = outputSel.sZDR.max(dim='sMult')
sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
ZDR = outputSel.ZDR.rolling(sHeight=20,center=True).mean()
#ls1=ax[0].plot(ZDR,outputSel.Temp,c='#66CCEE',label='cold',lw=2)#c='#4477AA'
ls2=ax.plot(sZDRmax,outputSel.Temp,ls='-',c='#66CCEE',label=r'cold',lw=3)#c='#66CCEE',

output['ZeH'] = output['sZeH_3.12e+01'].sum(dim='sMult')
output['ZeV'] = output['sZeV_3.12e+01'].sum(dim='sMult')
output['ZDR'] = 10*np.log10(output.ZeH)-10*np.log10(output.ZeV)
sZDRmax = output.sZDR.max(dim='sMult')
sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
ZDR = output.ZDR.rolling(sHeight=20,center=True).mean()
#ls3=ax[0].plot(ZDR,output.Temp,lw=2,c='#228833',label='all')
ls4=ax.plot(sZDRmax,output.Temp,lw=3,c='#228833',label=r'all')#c='#CCBB44'


vmin=180;vmax=350
outputSel = output.sel(sMult=slice(vmin,vmax))
outputSel['ZeH'] = outputSel['sZeH_3.12e+01'].sum(dim='sMult')
outputSel['ZeV'] = outputSel['sZeV_3.12e+01'].sum(dim='sMult')
outputSel['ZDR'] = 10*np.log10(outputSel.ZeH)-10*np.log10(outputSel.ZeV)
sZDRmax = outputSel.sZDR.max(dim='sMult')
sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
ZDR = outputSel.ZDR.rolling(sHeight=20,center=True).mean()
#ls5=ax[0].plot(ZDR,outputSel.Temp,c='#EE6677',label='DGL',lw=2)
ls6=ax.plot(sZDRmax,outputSel.Temp,ls='-',c='#EE6677',label=r'DGL',lw=3)#c='#AA3377'

ax1 = ax.twiny()
ls7=ax1.plot(datasZDRmax.med,datasZDRmax.Temp,c='k',lw=3,label=r'obs')
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
ls = ls2+ls4+ls6+ls7
labs = [l.get_label() for l in ls]
ax.legend(ls,labs,ncol=2,fontsize=12,loc='upper right')
ax.set_ylabel('T [°C]',fontsize=18)
#ax[0].set_xlabel('simulated ZDR [dB]',fontsize=18)
ax.set_xlabel(r'simulated sZDR$_{\rm max}$ [dB]',fontsize=18)

ax.set_xlim([0,7])
ax.set_ylim([0,-30])
ax.axhline(y=-10,ls='--',lw=2,c='r')
ax.axhline(y=-20,ls='--',lw=2,c='r')
ax.tick_params(axis='both',labelsize=16)
ax.set_yticks([0,-5,-10,-15,-20,-25,-30])
ax.set_xticks([0,1,2,3,4,5,6,7])
#ax.text(0.1,-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
ax.grid()
#ax[0].set_zorder(1)
#ax[0].legend(ls,labs,bbox_to_anchor=(-0,0.8),loc='lower left',fontsize=12,ncol=3)
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/sZDRmax_int.png',bbox_inches='tight')
#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/ZDR_int3.pdf',bbox_inches='tight')
plt.close()
quit()
'''
domTop = [3000, 2400, 1900]
expID1 ='1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop'
expID2 = '._atmo2_RHi105_nucleating_sPhiDmax_allDomTop_interpPPTrue'
#fig,ax = plt.subplots(nrows=3,ncols=3,figsize=(12,9))
fig = plt.figure(figsize=(15,9))
gs = mpl.gridspec.GridSpec(3, 3)
gs.update(wspace=0.45,hspace=0.1)        
ax00 = fig.add_subplot(gs[0, 0]) 
ax01 = fig.add_subplot(gs[0, 1]) 
ax02 = fig.add_subplot(gs[0, 2]) 
ax10 = fig.add_subplot(gs[1, 0]) 
ax11 = fig.add_subplot(gs[1, 1]) 
ax12 = fig.add_subplot(gs[1, 2]) 
ax20 = fig.add_subplot(gs[2, 0]) 
ax21 = fig.add_subplot(gs[2, 1]) 
ax22 = fig.add_subplot(gs[2, 2]) 

ax = np.array([[ax00,ax01,ax02],[ax10,ax11,ax12],[ax20,ax21,ax22]])

i_var=0
for i,dT in enumerate(domTop):
	experiment = '{dataPath}{expID1}{dT}{expID2}'.format(dataPath=dataPath,
		                                                   expID1=expID1,
		                                                   dT=dT,
		                                                   expID2=expID2)
	
	mcTable = xr.open_dataset(experiment+'/mass2fr.nc')
	mcTable = mcTable.to_dataframe()
	atmoFile = np.loadtxt(experiment+'/atmo.dat')
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
	#print(mcTableTemp)
	mcTableTmp = mcTableTmp.to_dataframe()
	mcTableTmp = mcTableTmp.set_index('sMult',drop=False)
	print(mcTableTmp['sMult'].max())
	print('dia max sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmax()]['dia'].max())
	print('mass max sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmax()]['mTot'].max())
	print('Phi max sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmax()]['sPhi'].min())
	print('Rho max sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmax()]['sRho_tot'].min())
	print('dia min sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmin()]['dia'].max())
	print('mass min sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmin()]['mTot'].max())
	print('Phi min sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmin()]['sPhi'].min())
	print('Rho min sMult',mcTableTmp.loc[mcTableTmp['sMult'].idxmin()]['sRho_tot'].min())
	#print(mcTableTmp.sel[1]['dia'])
	
	ax[i,0],s_m = plot.plot_var_particle_ID_temp(mcTableTmp,experiment,'sPhi','sMult',r'$\phi$','',log=True,Dini=True,ax=ax[i,0],cbar=False)#,zoom=[-30,-10])#,zoom=[0,0.2])
	ax[i,0].set_xticks([-15,-10,-5,0])
	ax[i,0].set_xlim([-19,0])
	ax[i,0].set_ylim([3*10**-3,7*10**-1])
	ax[i,0].text(ax[i,0].get_xlim()[0]+0.04*(ax[i,0].get_xlim()[1]-ax[i,0].get_xlim()[0]),ax[i,0].get_ylim()[1]-0.45*(ax[i,0].get_ylim()[1]-ax[i,0].get_ylim()[0]),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
	print('sMult of sPhiMin',mcTableTmp.loc[mcTableTmp['sPhi'].idxmin()]['sMult'].max())
	print('sPhimin',mcTableTmp['sPhi'].min())
	print('sMult of sPhiMax',mcTableTmp.loc[mcTableTmp['sPhi'].idxmax()]['sMult'].max())
	print('sPhimax',mcTableTmp['sPhi'].max())
	i_var +=1
	ax[i,1],s_m = plot.plot_var_particle_ID_temp(mcTableTmp,experiment,'mTot','sMult',r'mass','[kg]',log=True,Dini=True,ax=ax[i,1],cbar=False)
	ax[i,1].set_xlim([-19,0])
	ax[i,1].set_ylim([10**-11,1.5*10**-6])
	ax[i,1].set_yticks([10**-11,10**-10,10**-9,10**-8,10**-7,10**-6])
	y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
	ax[i,1].yaxis.set_minor_locator(y_minor)
	ax[i,1].yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
	ax[i,1].set_xticks([-15,-10,-5,0])
	ax[i,1].text(ax[i,1].get_xlim()[0]+0.04*(ax[i,1].get_xlim()[1]-ax[i,1].get_xlim()[0]),ax[i,1].get_ylim()[1]-0.7*(ax[i,1].get_ylim()[1]-ax[i,1].get_ylim()[0]),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
	i_var +=1
	print('massMax',mcTableTmp['mTot'].max())
	print('sMult of massMax',mcTableTmp.loc[mcTableTmp['mTot'].idxmax()]['sMult'].max())
	print('massMin',mcTableTmp['mTot'].min())
	print('sMult of massMin',mcTableTmp.loc[mcTableTmp['mTot'].idxmin()]['sMult'].max())
	ax[i,2],s_m = plot.plot_var_particle_ID_temp(mcTableTmp,experiment,'sRho_tot','sMult',r'$\rho$',r'[kgm$^{-3}$]',Dini=True,ax=ax[i,2],cbar=False)
	ax[i,2].set_xlim([-19,0])
	ax[i,2].set_ylim([250,950])
	ax[i,2].set_yticks([900,800,700,600,500,400,300])
	ax[i,2].set_xticks([-15,-10,-5,0])
	ax[i,2].text(ax[i,2].get_xlim()[0]+0.04*(ax[i,2].get_xlim()[1]-ax[i,2].get_xlim()[0]),ax[i,2].get_ylim()[1]-0.1*(ax[i,2].get_ylim()[1]-ax[i,2].get_ylim()[0]),'('+string.ascii_lowercase[i_var]+')',fontsize=18)
	i_var +=1
	print('sMult of sRhototMin',mcTableTmp.loc[mcTableTmp['sRho_tot'].idxmin()]['sMult'].max())
	print('sRho_totMin',mcTableTmp['sRho_tot'].min())
	print('sMult of sRhototMax',mcTableTmp.loc[mcTableTmp['sRho_tot'].idxmax()]['sMult'].max())
	print('sRho_totMax',mcTableTmp['sRho_tot'].max())
	
	if i != 2:
	  for j in range(3):
	    ax[i,j].set_xlabel('')
	    ax[i,j].set_xticklabels('')
cbar = fig.colorbar(s_m,ax=ax,pad=0.025,aspect=30,shrink=0.75)
cbar.set_label('particle ID',fontsize=18)
cbar.set_ticks(np.arange(0,mcTableTmp['sMult'].max()+1,2))
cbar.ax.tick_params(labelsize=16)
#plt.tight_layout()
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/initialized_Dmax_Rho_mass_-18_-15_-12.png',bbox_inches='tight')
plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/initialized_Dmax_Rho_mass_-18_-15_-12.pdf',bbox_inches='tight')
plt.close()
quit()

vmin_array=[180];vmax=350
for vmin in vmin_array:
	#- read in obs data to compare with zdr 
	dataPathObs = '/work/lvonterz/tripex_pol/classification/classification_output/DWR_classes/'
	dwrFile = 'histograms_trpol_masked_mean_total_cont_profile_2010_3classes_classDWR_4.0_9.5.nc'
	dwrData = xr.open_dataset(dataPathObs+dwrFile)

	cumFreq = np.nancumsum(dwrData['sZDRmaxHistProf'].values.T, axis=1)
	profMedFast,medFreq = getPercentil(cumFreq, dwrData['sZDRmax'].values, 0.5)

	expID = '1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo2_RHi105_nucleating_10i_interpPPTrue'

	dataXR = xr.open_dataset(dataPath+expID+'/9.6GHz_output_full_180_350_singleParticle.nc')
	print(dataXR) # this contains already single scattering properties
	
	dataPD = dataXR.to_dataframe()# TODO: not very nice to switch between xarray and pandas all the time, but for now the easiest
	dataPD = dataPD[(dataPD['sMult']>vmin) & (dataPD['sMult']<=vmax)] 
	dataXR = dataPD.to_xarray()
	
	mcTable = dataXR.to_dataframe()
	mcTable['sHeight'] = mcTable.index


	#- input for calculating spectra from single scattering results
	freq = 9.6*1e9
	maxVel=3; minVel=-3; nfft=512
	velRes = (maxVel - minVel)/nfft
	velBins = np.arange(minVel, maxVel, velRes)
	wl=[(constants.c / freq) * 1e3]
	velCenterBin = velBins[0:-1]+np.diff(velBins)/2.
	heightRange = np.arange(0, 5000, 50)

	#-- calculate spectra and moments
	specXR = xr.Dataset()
	#specXR_turb = xr.Dataset()
	counts = np.ones_like(heightRange)*np.nan
	vol = 5.0 * 50
	for i, heightEdge0 in enumerate(heightRange): 
		heightEdge1 = heightEdge0 +50
		mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &
		                         (mcTable['sHeight']<=heightEdge1)].copy()
		tmpSpecXR = mcr.getMultFrecSpec(wl, mcTableTmp, velBins,
		                                velCenterBin, heightEdge1,False,19,10**(-40/10),
		                                1e-6, 10.0, 2.0, 0.6/2./180.*np.pi, 
		                                scatSet={'mode':'full'})
		tmpSpecXR = tmpSpecXR/vol
		specXR = xr.merge([specXR, tmpSpecXR])
	print(specXR)
	wlStr = '{:.2e}'.format(wl[0])
	specXR['Ze_H_{0}'.format(wlStr)] = specXR['spec_H_{0}'.format(wlStr)].sum(dim='vel')
	specXR['Ze_V_{0}'.format(wlStr)] = specXR['spec_V_{0}'.format(wlStr)].sum(dim='vel')

	#- now make Temp also variable in the two dataframes (they can not be the same one because of different height coordinate)
	atmoFile = np.loadtxt(dataPath+expID+'/atmo.dat') 
	height = atmoFile[:,0] 
	Temp = atmoFile[:,2] -273.15 
	atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp']) 
	atmoPD.index.name='sHeight' 
	atmoXR = atmoPD.to_xarray() 
	atmoReindex = atmoXR.reindex_like(dataXR,method='nearest') 
	dataXR  = xr.merge([atmoReindex,dataXR])                  

	atmoPD.index.name='range'   
	atmoXR = atmoPD.to_xarray() 
	atmoReindex = atmoXR.reindex_like(specXR,method='nearest') 
	specXR  = xr.merge([atmoReindex,specXR])  
	
	sZDR = mcr.lin2db(specXR['spec_H_3.12e+01']) - mcr.lin2db(specXR['spec_V_3.12e+01'])
	sZDRmax = sZDR.max(dim='vel')
'''	
'''
	dataPol = xr.open_dataset('/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_0/2019/01/30/20190130_13_tripex_pol_poldata_L0_spec_regridded_dealized.nc') 
	timesel = pd.to_datetime('20190130 13:30:20')
	dataPol = dataPol.sel(time=timesel)
	dataPol = dataPol[['sZDR','HSpec']].where(dataPol['sSNR_H'] > 10.0)
    # move sZDR velocity to 0 for visualisation
	maxVelH,minVelH = post.calcOffset(dataPol,'HSpec')
	dataPol['maxVelH'] = maxVelH
	dataPol['minVelH'] = minVelH
	dataPol = post.removeOffset(dataPol)
	
	dataLV2 = xr.open_dataset('/data/obs/campaigns/tripex-pol/processed/tripex_pol_level_2/20190130_tripex_pol_3fr_L2_mom.nc')    
	data2 = dataLV2.sel(time=timesel)
	
	fig,ax = plt.subplots(ncols=2,figsize=(10,7),sharey=True)
	p1 = ax[0].pcolormesh(dataPol.Vel2ZeroH.fillna(0).T,
                     data2['ta'].values,
                     dataPol.sZDR,
                     vmin=0,vmax=4, cmap=getNewNipySpectral())
	ax[0].set_title('Observation',fontsize=22)
	ax[0].set_ylabel('T [°C]',fontsize=20)
	p2 = ax[1].pcolormesh(sZDR.vel,specXR.Temp,sZDR,
	                 vmin=0,vmax=4, cmap=getNewNipySpectral())
	ax[1].set_title('Simulation',fontsize=24)
	cb = plt.colorbar(p2,ax=ax[1])#,pad=0.02,aspect=20)#,ticks=v1)
	cb.set_label('sZDR [dB]',fontsize=20)
	cb.ax.tick_params(labelsize=18)
	for a in ax:
		a.set_xlim([-2,0.1])
		a.set_ylim([0,-20])
		a.set_xticks([-2,-1.5,-1,-0.5,0])
		a.set_xlabel(r'Doppler velocity [ms$^{-1}]$',fontsize=20)
		a.tick_params(axis='both', which='major', labelsize=18)
		a.grid()
	plt.tight_layout()
	plt.savefig('sZDR_obs_mod_20190130.png')
	plt.show()
	print(specXR)
	quit()
'''
'''	
	
	dataPD = dataXR.to_dataframe()
	dataPD = dataPD.set_index('sMult',drop=False)

	norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
	c_m = mpl.cm.gist_rainbow 
	# create a ScalarMappable and initialize a data structure 
	s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm) 
	s_m.set_array([]) 

	fig,ax=plt.subplots(ncols=2,figsize=(10,5),sharey=True) 
	for sID in sorted(dataPD['sMult'].unique()): 
		dataSel = dataPD.loc[sID] 
		ax[0].plot(dataSel['sKDP_3.12e+01'],dataSel['Temp'],color=s_m.to_rgba(sID-1),alpha=0.5)       
		ax[1].plot(mcr.lin2db(dataSel['sZeH_3.12e+01'])-mcr.lin2db(dataSel['sZeV_3.12e+01']),dataSel['Temp'],color=s_m.to_rgba(sID-1),alpha=0.5)
		#plt.show()
	#quit()
	ax[1].plot(mcr.lin2db(specXR['Ze_H_3.12e+01'])-mcr.lin2db(specXR['Ze_V_3.12e+01']),specXR.Temp,lw=3,c='k',ls='--',label='int ZDR') 
	ax[1].plot(sZDRmax,specXR.Temp,lw=3,c='k',label='sZDRmax')
	ax[1].plot(profMedFast*3,dwrData.ta,lw=3,c='k',ls=':',label='median sZDRmax')
	cbar = plt.colorbar(s_m,ax=ax[1]) 
	cbar.set_label('particle ID',fontsize=16) 
	cbar.ax.tick_params(labelsize=14)     
	ax[0].set_ylim([0,-20])        
	ax[0].set_xlim([0,0.0012])     
	ax[0].set_xticks([0,0.0005,0.001])
	ax[0].grid(True,ls='-.')  
	ax[0].set_ylabel('T [°C]',fontsize=18)                  
	ax[0].set_xlabel('sKDP [°/km]',fontsize=18)                   
	ax[1].set_xlabel('ZDR [dB]',fontsize=18)  
	 
	ax[1].grid(True,ls='-.')    
	ax[1].set_xlim([0,6])                                                                                                                                 
	ax[1].legend(loc='lower left',fontsize=16)  
	for a in ax:
		a.tick_params(axis='both',labelsize=16)     
	#fig.suptitle('vmin = {0}, vmax = {1}'.format(vmin,vmax))
	plt.tight_layout()     
	plt.savefig(dataPath+expID+'/single_particle_moments_int_{0}_{1}_color.png'.format(vmin,vmax))  
	plt.show() 	
'''	
