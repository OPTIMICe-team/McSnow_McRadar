'''
This is the beginning of plotting routines for McSnow output
'''
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import mcradar as mcr
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
import matplotlib as mpl
from scipy import constants
import math

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

def wavg(group):
    d = group['sPhi']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgD(group):
    d = group['dia_cm']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgDm(group):
    d = group['dia']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgM(group):
    d = group['mTot_g']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgMkg(group):
    d = group['mTot']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgRho(group):
    d = group['sRho_tot']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgNmono(group):
    d = group['sNmono']
    w = group['sMult']
    return (d * w).sum() / w.sum()



def plotOverview(output,dicSettings,inputPath,wl1,wl2):
	freq1 = (constants.c / wl1)  *1e3 / 1e9
	freq2 = (constants.c / wl2)  *1e3 / 1e9
	freq1 = '{:.1f}'.format(freq1)
	freq2 = '{:.1f}'.format(freq2)
	
	specH90 = mcr.lin2db(output['spec_H'].sel(wavelength=wl1,elevation=90)) # Ka-Band reflectivity
	specH30 = mcr.lin2db(output['spec_H'].sel(wavelength=wl2,elevation=30))
	specV30 = mcr.lin2db(output['spec_V'].sel(wavelength=wl2,elevation=30))
	specH90 = specH90.where(specH90 > -40)
	specH30 = specH30.where(specH30 > -40)
	specV30 = specV30.where(specV30 > -40)
	ZDR = specH30 - specV30
	fig,ax = plt.subplots(ncols=4,figsize=(20,5),sharey=True)
	DWR = mcr.lin2db(output['Ze_H'].sel(wavelength=wl1,elevation=90)) - mcr.lin2db(output['Ze_H'].sel(wavelength=wl2,elevation=90))
	ax[0].plot(DWR,output['Temp'],lw=2)
	ax[0].set_ylim([0,np.min(output.Temp)-1])	
	ax[0].set_ylabel('T [°C]',fontsize=24)
	ax[0].set_xlabel('DWR$_{{{freq1},{freq2}}}$ [dB]'.format(freq1=freq1,freq2=freq2),fontsize=24)
	
	ax[1].plot(output['KDP'].sel(wavelength=wl2,elevation=30),output['Temp'],lw=2)
	ax[1].set_xlabel(r'KDP [°km$^{-1}$]',fontsize=24)
	
	p1=ax[2].pcolormesh(output.vel,output.Temp,specH90,vmin=-30,vmax=10,cmap=getNewNipySpectral(),shading='auto')
	cb = plt.colorbar(p1,ax=ax[2])
	cb.set_label('sZeH [dBz]',fontsize=24)
	cb.set_ticks([-30,-25,-20,-15,-10,-5,0,5,10])
	cb.ax.tick_params(labelsize=20)
	ax[2].set_xlabel(r'Doppler velocity [ms$^{-1}$]',fontsize=24)
	ax[2].set_xlim([-2,0])
	p1=ax[3].pcolormesh(output.vel,output.Temp,ZDR,vmin=-1,vmax=5,cmap=getNewNipySpectral(),shading='auto')
	cb = plt.colorbar(p1,ax=ax[3])
	cb.set_label('sZDR [dB]',fontsize=24)
	cb.ax.tick_params(labelsize=20)
	ax[3].set_xlabel(r'Doppler velocity [ms$^{-1}$]',fontsize=24)
	ax[3].set_xlim([-2,0])
	for a in ax:
		a.tick_params(which='both',labelsize=20)
		a.grid()
	plt.tight_layout()
	plt.savefig(inputPath+'Profiles_overview.png')
	plt.close()
	
def plotPropSpecThesis(ax,heightRange,heightRes,mcTable,velBins,prop,savefig=False,cbar=True,gridBaseArea=5.0,diff=False,mcTable1=None):
# plots the aspect ratios depending on the velocity (so sort of like a Doppler Spectrum)
    
    for i, heightEdge0 in enumerate(heightRange[0:-1]):
        #print(mcTable)
        heightEdge1 = heightEdge0 + heightRes
        height = heightEdge0+heightRes/2
        
        mcTableTmp = mcTable.where((mcTable['sHeight']>heightEdge0) &
                             (mcTable['sHeight']<=heightEdge1))
        mcTableTmp = mcTableTmp.to_dataframe()
        if diff == True:
            mcTableTmp1 = mcTable1.where((mcTable1['sHeight']>heightEdge0.values) &
                     (mcTable1['sHeight']<=heightEdge1.values))

            binVel,sVel = pd.cut(mcTableTmp1['vel'],bins=velBins,retbins=True)
            #group according to velocity bins
            grouped1 = mcTableTmp1.groupby(binVel)
            #mcTableTmp[var] = mcTableTmp[var]-mcTableTmp1[var]
		    
        Temp = np.asarray(mcTableTmp['Temp'].mean())
        # cut into velocity bins
        
        binVel,sVel = pd.cut(mcTableTmp['vel'],bins=velBins,retbins=True)
        #group according to velocity bins
        grouped = mcTableTmp.groupby(binVel)
        
        # apply sMult as weight
        if prop == 'dia':
            if diff:
                weighted = grouped.apply(wavgDm)
                weighted1 = grouped1.apply(wavgDm)
                weighted = weighted1.sub(weighted,fill_value=0)*1e3
                #print(weighted.max())
                
                cmap='coolwarm'
                vmin=-7;vmax=7
                cbar_label = r'$\Delta$D$_{\rm max}}$ [mm]'
            else:
                weighted = grouped.apply(wavgDm)*1e3
                vmin = -2;vmax = 1.5
                cmap = getNewNipySpectral()
                cbar_label = r'D$_{\rm max}}$ [mm]'
        
        elif prop == 'sPhi':
            weighted = grouped.apply(wavg)
            #print(height)
            #print(weighted.min())
            vmin = 0.;vmax = 3.
            bottom = cm.get_cmap('binary',21)
            top = cm.get_cmap('copper',41)
            bottom =  bottom(np.linspace(0.2, 1, 21))
            top = top(np.linspace(0.3, 1, 41))
            bottom[16] = colors.to_rgba('cyan')
            bottom[17] = colors.to_rgba('mediumspringgreen')
            bottom[18] = colors.to_rgba('lime')
            bottom[19] = colors.to_rgba('green')
            bottom[20] = colors.to_rgba('magenta')

            top[0] = colors.to_rgba('red')
            top[1] = colors.to_rgba('darkred')
            top[2] = colors.to_rgba('darkorange')
            top[3] = colors.to_rgba('gold')
      
            newcolors = np.vstack((bottom,
                                   top))
            cmap = ListedColormap(newcolors)
            cbar_label = 'aspect ratio'
        elif prop == 'sNmono':
            if diff:
                weighted = grouped.apply(wavgNmono)
                weighted1 = grouped1.apply(wavgNmono)
                weighted = weighted1.sub(weighted,fill_value=0)
                #print(weighted.max())
                vmin=-50; vmax= 50
                cmap = 'coolwarm'
                cbar_label = r'$\Delta$ Number of monomers'
            else:
                weighted = grouped.apply(wavgNmono)
                vmin = 1;vmax = 50
            #viridis = cm.get_cmap('nipy_spectral', vmax)
            #newcolors = viridis(np.linspace(0, 1, vmax))
                newNipy = getNewNipySpectral()
                newcolors = newNipy(np.linspace(0, 1, vmax))
                red = ([1,0/vmax,0/vmax,1])
                pink = ([248/vmax, 24/vmax, 148/vmax, 1])
                newcolors[:1, :] = colors.to_rgba('cyan')
                cmap = ListedColormap(newcolors)
                cbar_label = 'Number of monomers'
        
        elif prop == 'number_conc':
            volume = heightRes*gridBaseArea
            if diff:
                weighted = grouped['sMult'].sum()/volume
                weighted1 = grouped1['sMult'].sum()/volume
                weighted = weighted1.sub(weighted,fill_value=0)
                weighted = weighted.replace(0,np.nan)
                #print(weighted.max())
                vmin = -100;vmax = 100
                cmap='coolwarm'
                cbar_label = r'$\Delta$ Number conc. [m$^{-3}$]'
            else:
                weighted =grouped['sMult'].sum()/volume #grouped['sMult'].count()# grouped['sMult'].sum() #np.log10(grouped['sMult'].sum())
                vmin = 0;vmax = 4
            # we need concentration per m3, not per volume which is spanned by gridbase area and height res. So we need to divide by gridbaseArea*heightRes
                cmap = getNewNipySpectral()
                cbar_label = r'Number conc. [m$^{-3}$]'
        elif prop == 'mTot':
            if diff:
                weighted = grouped.apply(wavgMkg)*1e3
                weighted1 = grouped1.apply(wavgMkg)*1e3
                weighted = weighted1.sub(weighted,fill_value=0)
                #print(weighted.max())
                
                cmap='coolwarm'
                vmin=-1e-5;vmax=1e-5
                cbar_label = r'$\Delta$ mass [g]'
            else:
                weighted = grouped.apply(wavgMkg)*1e3
                vmin = -8;vmax = -2
                cmap = getNewNipySpectral()
                cbar_label = r'mass [g]'
        elif prop == 'sRho_tot':
            if diff:
                weighted = grouped.apply(wavgRho)*1e3
                weighted1 = grouped1.apply(wavgRho)*1e3
                weighted = weighted1.sub(weighted,fill_value=0)
                #print(weighted.max())
                
                cmap='coolwarm'
                vmin=0;vmax=1000
                cbar_label = r'$\Delta$ mass [g]'
            else:
                weighted = grouped.apply(wavgRho)
                vmin = 0;vmax = 1000
                cmap = getNewNipySpectral()
                cbar_label = r'$\rho$ [kgm$^{-3}$]'
        else:
            print(prop+' not yet defined')
        if i == 0:
            binnedXR = xr.DataArray(weighted.values.reshape(len(weighted),1),
                                    dims=('vel','T'),
                                    coords={'T':Temp.reshape(1),'vel':velBins[0:-1]})
        else:
            tmpXR = xr.DataArray(weighted.values.reshape(len(weighted),1),
                                    dims=('vel','T'),
                                    coords={'T':Temp.reshape(1),'vel':velBins[0:-1]})
            
            binnedXR = xr.concat([binnedXR,tmpXR],dim='T')
    binnedXR = binnedXR.where(~np.isnan(binnedXR['T']), drop=True)
    if diff==False: 
        if (prop == 'dia' or prop == 'number_conc' or prop == 'mTot'):
            
            plot=ax.pcolormesh(binnedXR.vel,binnedXR['T'],binnedXR.T,cmap=cmap,norm=colors.LogNorm(vmin=10**vmin,vmax=10**vmax))#,cmap=cmap,add_colorbar=False)
        else:
            
            plot=ax.pcolormesh(binnedXR.vel,binnedXR['T'],binnedXR.T,cmap=cmap,vmin=vmin,vmax=vmax)
            #plot=binnedXR.plot(ax=ax,x='vel',vmin=vmin,vmax=vmax,cmap=cmap,add_colorbar=False)
    else:
        
        plot=binnedXR.plot(ax=ax,x='vel',vmin=vmin,vmax=vmax,cmap=cmap,add_colorbar=False)
   #if (prop == 'dia' or prop == 'number_conc' or prop == 'mTot'):
   #     plot=ax.pcolormesh(binnedXR.vel,binnedXR['T'],binnedXR.T,cmap=cmap,norm=colors.LogNorm(vmin=10**vmin,vmax=10**vmax))#,cmap=cmap,add_colorbar=False)
   # else:
   #     plot=binnedXR.plot(ax=ax,x='vel',vmin=vmin,vmax=vmax,cmap=cmap,add_colorbar=False)
    if cbar==True:
        cb = plt.colorbar(plot,ax=ax,pad=0.02,aspect=20)#,ticks=v1)
        cb.set_label(cbar_label,fontsize=18)
        cb.ax.yaxis.offsetText.set_fontsize(16)
        cb.ax.tick_params(labelsize=16)
    ax.set_xlabel(r'Velocity [ms$^{-1}$]',fontsize=18)
    ax.set_xticks([-2,-1.5,-1,-0.5,0])
    #ax.pcolormesh(binnedXR.vel,binnedXR['T'],binnedXR,vmin=vmin,vmax=vmax,cmap=cmap)
    #plt.show()
    #ax.set_xlabel('vel [m/s]',fontsize=18)
    
    if savefig==True:
        plt.grid(True,ls='-.')
        plt.ylim([0,dicSettings['maxHeight']])
        plt.tight_layout()
        plt.savefig(inputPath+'1d_habit_spec_'+prop+'_weighted_alltimes.png')
        plt.close()
    else:
        return ax

def plotMoments(dicSettings,output,inputPath,plotTemp=False,mult_conc=False):
	for wl in dicSettings['wl']:
		for elv in dicSettings['elv']:
			#wlStr = '{:.2e}'.format(wl)
			freq = (constants.c / wl)  *1e3 / 1e9
			#print(freq)
			freq = '{:.1f}'.format(freq)

			if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):

				if plotTemp == True:
					saveName = '1d_habit_moments_freq{0}_elv{1}_{2}_Temp.png'.format(freq,elv,dicSettings['scatSet']['mode'])
					vary='T'; varUnit = '[°C]'
					ylim = [0,-30]
				else:
					saveName = '1d_habit_moments_freq{0}_elv{1}_{2}.png'.format(freq,elv,dicSettings['scatSet']['mode'])
					vary='range'; varUnit = '[m]'
					ylim = [0,dicSettings['maxHeight']]

				#print(output['Ze_H_{0}'.format(wlStr)])

				fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
				axes[0].plot(output['MDV_H'].sel(wavelength=wl,elevation=elv),output[vary],lw=2)
				axes[0].set_title('Freq: {0} elv: {1}, MDV'.format(freq, elv),fontsize=16)
				axes[0].grid(True,ls='-.')
				axes[0].set_xlabel('MDV [m/s]',fontsize=16)
				axes[0].set_ylim(ylim)
				axes[0].set_ylabel(vary+' '+varUnit,fontsize=16)
				axes[0].tick_params(axis='both',labelsize=14)

				axes[1].plot(mcr.lin2db(output['Ze_H'].sel(wavelength=wl,elevation=elv)),output[vary],linewidth=2) # TODO: change back to ZeH
				axes[1].set_xlabel('Ze [dB]',fontsize=16)
				axes[1].grid(True,ls='-.')
				axes[1].set_ylim(ylim)
				axes[1].tick_params(axis='both',labelsize=14)
				plt.tight_layout()
			else:
				if plotTemp == True:
					saveName = '1d_habit_moments_freq{freq}_elv{elv}_{mode}_Temp.png'.format(freq=freq,elv=elv,mode=dicSettings['scatSet']['mode'])
					vary = 'Temp';label = ' T [°C]'
					ylim=([0,np.min(output.Temp)-1])

				else:
					saveName = '1d_habit_moments_freq{freq}_elv{elv}_{mode}.png'.format(freq=freq,elv=elv,mode=dicSettings['scatSet']['mode'])
					vary = 'range';varUnit = '[m]'
					ylim= ([0,dicSettings['maxHeight']])

				if 'LDR' in output:
					fig,axes = plt.subplots(ncols=5,figsize=(22,5),sharey=True)
				else:
					fig,axes = plt.subplots(ncols=4,figsize=(20,5),sharey=True)

				if 'KDPAgg' in output:
					axes[0].plot(output['KDPAgg'].sel(wavelength=wl,elevation=elv),output[vary],lw=2,label='Agg')
					axes[0].plot(output['KDPMono'].sel(wavelength=wl,elevation=elv),output[vary],lw=2,label='Mono')
					axes[0].plot(output['KDP'].sel(wavelength=wl,elevation=elv),output[vary],lw=2,label='total')
					axes[0].legend()
				else:
					axes[0].plot(output['KDP'].sel(wavelength=wl,elevation=elv),output[vary],lw=2)
				axes[0].set_title('freq: {0} elv: {1}'.format(freq, elv))
				axes[0].set_ylabel(label, fontsize=16)
				axes[0].grid(True,ls='-.')
				axes[0].set_xlabel(r'KDP [°km$^{-1}$]',fontsize=16)
				axes[0].set_ylim(ylim)
				axes[0].tick_params(axis='both',labelsize=14)

				# plot ZDR
				axes[1].plot(mcr.lin2db(output['Ze_H'].sel(wavelength=wl,elevation=elv))-mcr.lin2db(output['Ze_V'].sel(wavelength=wl,elevation=elv)),output[vary],linewidth=2)#
				axes[1].set_xlabel('ZDR [dB]',fontsize=16)
				axes[1].grid(True,ls='-.')
				axes[1].set_ylim(ylim)
				axes[1].tick_params(axis='both',labelsize=14)
				# plot Ze
				axes[2].plot(mcr.lin2db(output['Ze_H'].sel(wavelength=wl,elevation=elv)),output[vary],linewidth=2) #
				axes[2].set_xlabel('Ze_H [dB]',fontsize=16)
				axes[2].grid(True,ls='-.')
				axes[2].set_ylim(ylim)
				axes[2].tick_params(axis='both',labelsize=14)
				# plot MDV
				axes[3].plot(output['MDV_H'].sel(wavelength=wl,elevation=elv),output[vary],linewidth=2) #
				axes[3].set_xlabel(r'MDV [ms$^{-1}$]',fontsize=16)
				axes[3].grid(True,ls='-.')
				axes[3].set_ylim(ylim)
				axes[3].tick_params(axis='both',labelsize=14)
				
				if 'LDR' in output:
					axes[4].plot(output['LDR'].sel(wavelength=wl,elevation=elv),output[vary],linewidth=2) #
					axes[4].set_xlabel(r'LDR [dB]',fontsize=16)
					axes[4].grid(True,ls='-.')
					axes[4].set_ylim(ylim)
					axes[4].tick_params(axis='both',labelsize=14)
				
				plt.tight_layout()

			plt.savefig(inputPath+saveName)
			plt.close()

def plotDWR(dicSettings,wl1,wl2,output,inputPath,plotTemp=False):
	# input: dicSettings, the two wavelengts, the McRadar output, the axes for the plot
	freq1 = (constants.c / wl1)  *1e3 / 1e9
	freq2 = (constants.c / wl2)  *1e3 / 1e9
	freq1 = '{:.1f}'.format(freq1)
	freq2 = '{:.1f}'.format(freq2)
	for elv in dicSettings['elv']:
		if plotTemp == True:
			saveName = '1d_habit_DWR_{freq1}_{freq2}_elv{elv}_{mode}_Temp.png'.format(freq1=freq1,freq2=freq2,elv=elv,mode=dicSettings['scatSet']['mode'])
			vary = 'Temp'; varUnit = '[°C]'
			ylim = [0,np.min(output.Temp)-1]
		else:
			saveName = '1d_habit_DWR_{freq1}_{freq2}_elv{elv}_{mode}.png'.format(freq1=freq1,freq2=freq2,elv=elv,mode=dicSettings['scatSet']['mode'])
			vary = 'range'; varUnit = '[m]'
			ylim = [0,dicSettings['maxHeight']]

		fig,axes = plt.subplots(figsize=(5,4))
		DWR = mcr.lin2db(output['Ze_H'].sel(wavelength=wl1,elevation=elv)) - mcr.lin2db(output['Ze_H'].sel(wavelength=wl2,elevation=elv))
		axes.plot(DWR,output[vary],lw=2)
		axes.grid(True,ls='-.')
		axes.set_xlabel(r'DWR$_{{{freq1},{freq2}}}$ [dB]'.format(freq1=freq1,freq2=freq2),fontsize=16)
		axes.set_ylim(ylim)
		axes.set_ylabel(vary+' '+varUnit,fontsize=16)
		axes.tick_params(axis='both',labelsize=14)	   
		plt.tight_layout()        
		plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
		plt.close()

def plotDWRspectra(dicSettings,wl1,wl2,output,inputPath,plotTemp=False):
	#for wl in dicSettings['wl']:
	#wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
	#wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
	freq1 = (constants.c / wl1)  *1e3 / 1e9
	freq2 = (constants.c / wl2)  *1e3 / 1e9
	freq1 = '{:.1f}'.format(freq1)
	freq2 = '{:.1f}'.format(freq2)
	for elv in dicSettings['elv']:
		if plotTemp == True:
			saveName = '1d_habit_sDWR_{freq1}_{freq2}_elv{elv}_{mode}_Temp.png'.format(freq1=freq1,freq2=freq2,elv=elv,mode=dicSettings['scatSet']['mode'])
			vary = 'Temp';varUnit = '[°C]'
			ylim = [0,np.min(output.Temp)-1]
		else:
			saveName = '1d_habit_sDWR_{freq1}_{freq2}_elv{elv}_{mode}.png'.format(freq1=freq1,freq2=freq2,elv=elv,mode=dicSettings['scatSet']['mode'])
			vary = 'range'
			varUnit = '[m]'
			ylim = [0,dicSettings['maxHeight']]

		fig,axes = plt.subplots(figsize=(5,4))
		specH1 = mcr.lin2db(output['spec_H'].sel(wavelength=wl1,elevation=elv))
		specH1 = specH1.where(specH1 > -40)
		specH2 = mcr.lin2db(output['spec_H'].sel(wavelength=wl2,elevation=elv))
		specH2 = specH2.where(specH2 > -40)

		DWR =specH1 - specH2
		plot = axes.pcolormesh(output.vel,output[vary],DWR,vmin=0, vmax=15,shading='auto', cmap=getNewNipySpectral())
		cb = plt.colorbar(plot,ax=axes,pad=0.02,aspect=20)#,ticks=v1)
		cb.set_label(r'sDWR$_{{{freq1},{freq2}}}$'.format(freq1=freq1,freq2=freq2),fontsize=16)
		cb.ax.tick_params(labelsize=14)
		axes.set_xlabel('Doppler velocity [m/s]',fontsize=16)
		axes.set_ylabel(vary+' '+varUnit,fontsize=16)
		axes.tick_params(axis='both',labelsize=14)
		axes.grid(True,ls='-.')
		axes.set_ylim(ylim)
		axes.set_xlim(-2, 0)                     
		
		plt.tight_layout()        
		plt.savefig(inputPath+saveName, format='png', dpi=200)#, bbox_inches='tight')
		plt.close()


def plotSpectra(dicSettings,output,inputPath,minmax=None,plotTemp=False):
	for wl in dicSettings['wl']:
		for elv in dicSettings['elv']:
			wlStr = '{:.2e}'.format(wl)
			freq = (constants.c / float(wlStr))  *1e3 / 1e9
			freq = '{:.1f}'.format(freq)

			if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
				
				if plotTemp == True:
					saveName = '1d_habit_spectra_freq{wl}_elv{elv}_{mode}_{part}_Temp.png'.format(wl=freq,elv=elv,mode=dicSettings['scatSet']['mode'],
																							part=dicSettings['scatSet']['particle_name'])
					specH = mcr.lin2db(output['spec_H_{0}_elv{1}'.format(wlStr,elv)])
					vary='Temp';label = 'T [°C]'
					ylim = [0,np.min(output.Temp)-1]
				else:
					saveName = '1d_habit_spectra_freq{wl}_elv{elv}_{mode}_{part}.png'.format(wl=freq,elv=elv,mode=dicSettings['scatSet']['mode'],
																						part=dicSettings['scatSet']['particle_name'])
																						
					specH = mcr.lin2db(output['spec_H'].sel(elevation=elv,wavelength=wl))
					vary = 'range'; label = 'range [m]'
					ylim = [0,dicSettings['maxHeight']]
				fig,ax = plt.subplots(figsize=(5,4))
				specH = specH.where(specH > -40)
				plot = ax.pcolormesh(output.vel,output[vary],specH,cmap=getNewNipySpectral(),vmin=-30,vmax=10,shading='auto')
				cb = plt.colorbar(plot,ax=ax,pad=0.02,aspect=20)#,ticks=v1)
				cb.set_label(r'sZe$_{\rm Ka}$',fontsize=16)
				cb.ax.tick_params(labelsize=14)
				ax.set_xlabel('Doppler velocity [m/s]',fontsize=16)
				ax.set_ylabel(label,fontsize=16)
				ax.tick_params(axis='both',labelsize=14)
				ax.set_ylim(ylim)
				ax.set_xlim(-2, 0)
				ax.grid(True,ls='-.')      
				plt.tight_layout()
			else:
				
				if plotTemp == True:
					if minmax:
						saveName = '1d_habit_spectra_freq{wl}_elv{elv}_{mode}_{minmax}_Temp.png'.format(wl=freq,elv=elv,mode=dicSettings['scatSet']['mode'],minmax=minmax)
					else:
						saveName = '1d_habit_spectra_freq{wl}_elv{elv}_{mode}_Temp.png'.format(wl=freq,elv=elv,mode=dicSettings['scatSet']['mode'])
				else:
					if minmax:
						saveName = '1d_habit_spectra_freq{wl}_elv{elv}_{mode}_{minmax}.png'.format(wl=freq,elv=elv,mode=dicSettings['scatSet']['mode'],minmax=minmax)
					else:
						saveName = '1d_habit_spectra_freq{wl}_elv{elv}_{mode}.png'.format(wl=freq,elv=elv,mode=dicSettings['scatSet']['mode'])


				specH = mcr.lin2db(output['spec_H'].sel(wavelength=wl,elevation=elv))
				specV = mcr.lin2db(output['spec_V'].sel(wavelength=wl,elevation=elv))
				
				specH = specH.where(specH > -40)
				specV = specV.where(specV > -40)
				if 'spec_HV' in output:
					sLDR = mcr.lin2db(output['spec_HV'].sel(wavelength=wl,elevation=elv)) - specH
				ZDR = specH - specV
				if plotTemp == True:
					vary = 'Temp'; label='T [°C]'
					ylim=([0,np.min(output.Temp)-1])
				else:
					vary = 'range'; label='range m'
					ylim= ([0,dicSettings['maxHeight']])
				fig,axes = plt.subplots(ncols=3,figsize=(12,5),sharey=True)
				#specH.plot(ax=axes[0],vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})

				p1 = axes[0].pcolormesh(output.vel,
								output[vary],
								specH,vmin=-30,vmax=10,cmap=getNewNipySpectral(),shading='auto')
				cb = plt.colorbar(p1,ax=axes[0])
				cb.set_label('sZeH [dBz]')
				axes[0].set_xlabel('vel [m/s]')
				
				axes[0].set_title('Ze_H_spec freq: {0} GHz, elv: {1}'.format(freq, elv))
				#axes[0].set_ylim([0,dicSettings['maxHeight']])
				axes[0].set_ylim(ylim)
				axes[0].set_ylabel(label)
				axes[0].set_xlim(-2, 0)
				axes[0].grid(True,ls='-.')

				#specV.plot(ax=axes[1],vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
				if 'spec_HV' in output:
					p2 = axes[1].pcolormesh(output.vel,
									output[vary],
									sLDR,vmin=-35,vmax=-20,cmap=getNewNipySpectral(),shading='auto')
					cb = plt.colorbar(p2,ax=axes[1])
					cb.set_label('sLDR [dB]')
					axes[1].set_title('LDR, freq: {0} GHz, elv: {1}'.format(freq, elv))
				else:
					p2 = axes[1].pcolormesh(output.vel,
									output[vary],
									specV,vmin=-30,vmax=10,cmap=getNewNipySpectral(),shading='auto')
					cb = plt.colorbar(p2,ax=axes[1])
					cb.set_label('sZeH [dBz]')
					axes[1].set_title('Ze_V_spec, freq: {0} GHz, elv: {1}'.format(freq, elv))
				#axes[1].set_ylim([0,dicSettings['maxHeight']])
				axes[1].set_xlim(-2, 0)
				axes[1].grid(True,ls='-.')
				axes[1].set_ylabel('')
				axes[1].set_xlabel('vel [m/s]')
				axes[1].set_ylim(ylim)

				#ZDR.plot(ax=axes[2],vmin=-0.5, vmax=3,cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZDR [dB]'}) 
				p3 = axes[2].pcolormesh(output.vel,
								output[vary],
								ZDR,cmap=getNewNipySpectral(),vmin=-1,vmax=5,shading='auto')
				cb = plt.colorbar(p3,ax=axes[2])
				cb.set_label('sZDR [dB]')
				axes[2].set_title('ZDR freq: {0} GHz, elv: {1}'.format(freq, elv))
				axes[2].set_xlim(-2, 0)
				#axes[2].set_ylim([0,dicSettings['maxHeight']])
				axes[2].grid(True,ls='-.')
				axes[2].set_xlabel('vel [m/s]')
				axes[2].set_ylabel('')
				axes[2].set_ylim(ylim)
				plt.tight_layout()

			plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
			plt.close()

def plotAtmo(atmoFile,inputPath):
    height = atmoFile[:,0]
    rho = atmoFile[:,1]
    Temp = atmoFile[:,2] -273.15
    p = atmoFile[:,3]
    eta = atmoFile[:,4]
    ssat = atmoFile[:,5]
    rh = atmoFile[:,6]
    psatw = atmoFile[:,7]
    psati = atmoFile[:,8]
    lwc = atmoFile[:,9]
    #qv = atmoFile[:,10]
    #qc = atmoFile[:,11]
    fig,axes=plt.subplots(ncols=2,figsize=(10,5),sharey=True)
    axes[0].plot(Temp,height,lw=2)
    axes[0].grid(True,ls='-.')
    axes[0].set_xlabel('Temp [°C]',fontsize=18)
    axes[0].set_ylabel('height [m]',fontsize=18)
    axes[0].tick_params(labelsize=16)
    axes[1].plot(rh,height,label='RHw',lw=2)
    axes[1].plot(100+100*ssat,height,label='RHi',lw=2)
    axes[1].legend(fontsize=18)
    axes[1].grid(True,ls='-.')
    axes[1].set_xlabel('rh [%]',fontsize=18)
    axes[1].tick_params(labelsize=16)
    #axes[2].plot(qv,height)
    #axes[2].grid(True,ls='-.')
    #axes[2].set_xlabel('qv [%]')
    
    #axes[2].plot(ssat,height)
    #axes[2].grid(True,ls='-.')
    #axes[2].set_xlabel('ssat')
    
    plt.tight_layout()
    plt.savefig(inputPath+'atmo_test.png')
    plt.close()
    
    
def plotMomentsObs(dataLV2,dataPol,outPath,outName):
    fig,axes = plt.subplots(nrows=3,figsize=(18,18),sharey=True)
    plot = dataLV2.Ka_DBZ_H.plot(ax=axes[0],x='time',vmin=-35,vmax=25,cmap='jet',add_colorbar=False)
    cb = plt.colorbar(plot,ax=axes[0]) 
    cb.set_label('Ze [dB]',fontsize=18)
    cb.ax.tick_params(labelsize=16)
    CS = dataLV2.ta.plot.contour(ax=axes[0],y='range',levels=[-15,0],colors='k',linewidths=[3,3],linestyles=['--','-'])
    plt.clabel(CS, CS.levels, fontsize=22, fmt='%1.f '+'C')
    axes[0].set_title('Ka-Band Ze',fontsize=20)
    plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
    axes[0].grid()
    axes[0].set_xlabel('')
    axes[0].tick_params(axis='y',labelsize=16)
    axes[0].set_ylabel('range [m]',fontsize=20) 

    plot = dataPol.KDP.plot(ax=axes[1],x='time',vmin=-1,vmax=4,cmap='jet',add_colorbar=False)
    cb = plt.colorbar(plot,ax=axes[1]) 
    cb.set_label('KDP [°/km]',fontsize=18)
    cb.ax.tick_params(labelsize=16)
    CS = dataLV2.ta.plot.contour(ax=axes[1],y='range',levels=[-15,0],colors='k',linewidths=[3,3],linestyles=['--','-'])
    plt.clabel(CS, CS.levels, fontsize=22, fmt='%1.f '+'C')
    axes[1].set_title('W-Band KDP',fontsize=20)
    plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
    axes[1].grid()
    axes[1].set_xlabel('')
    axes[1].tick_params(axis='both',labelsize=16)
    axes[1].set_ylabel('range [m]',fontsize=20) 

    plot = dataPol.ZDR.plot(ax=axes[2],x='time',vmin=-1,vmax=3,cmap='jet',add_colorbar=False)
    cb = plt.colorbar(plot,ax=axes[2]) 
    cb.set_label('ZDR [dB]',fontsize=18)
    cb.ax.tick_params(labelsize=16)
    

    CS = dataLV2.ta.plot.contour(ax=axes[2],y='range',levels=[-15,0],colors='k',linewidths=[3,3],linestyles=['--','-'])
    plt.clabel(CS, CS.levels, fontsize=22, fmt='%1.f '+'C')
    axes[2].set_title('W-Band KDP',fontsize=20)
    plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
    axes[2].grid()
    axes[2].set_xlabel('')
    axes[2].tick_params(axis='both',labelsize=16)
    axes[2].set_ylabel('range [m]',fontsize=20) 
    
    plt.savefig(outPath+outName+'.png',dpi=200,bbox_inches='tight')
    plt.close()
    

def plot_agg_kernel_xr(ax,data,kernel,height,cmap='viridis',noN=False):
    
    #if noN==True:
    #    im = data.plot(ax=ax,cmap=cmap,vmax=0.72,vmin=0.0,add_colorbar=False)
    #else:
    im = data.plot(ax=ax,cmap=cmap,add_colorbar=False)
    cbar = plt.colorbar(im,ax=ax)
    if kernel == 'Ak':
        ax.set_xlabel(r'D$_j$ [m]')
        ax.set_ylabel(r'D$_i$ [m]')
        if noN == False:
            cbar.set_label(r'$(A_i^{0.5}+A_j^{0.5})^2*sqrt(N(D_i)*N(D_j))$')
        else:
            cbar.set_label(r'$(A_i^{0.5}+A_j^{0.5})^2)$')
        plt.tight_layout()
    elif kernel =='Vk':
        ax.set_xlabel(r'D$_j$ [m]')
        ax.set_ylabel(r'D$_i$ [m]')
        if noN == False:
            cbar.set_label(r'$abs(v_i-v_j)*sqrt(N(D_i)*N(D_j))$')
        else: 
            cbar.set_label(r'$abs(v_i-v_j)$')
        plt.tight_layout()
    elif kernel == 'Aggk':
        ax.set_xlabel(r'D$_j$ [m]')
        ax.set_ylabel(r'D$_i$ [m]')
        if noN == False:
            cbar.set_label(r'$Ak*Vk*Ec*Es*sqrt(N(D_i)*N(D_j))$')
        else:
            cbar.set_label(r'$Ak*Vk*Ec*Es)$')
        plt.tight_layout()
    elif kernel == 'Ec':
        ax.set_xlabel(r'D$_j$ [m]')
        ax.set_ylabel(r'D$_i$ [m]')
        if noN == False:
            cbar.set_label(r'$Ec*N(D_i)*N(D_j)$')
        else:
            cbar.set_label(r'$Ec$')
        plt.tight_layout()
    else:
        print('no kernel has been specified, so no axis or colorbar labels will be applied')
    ax.set_title('height='+str(height)+'m')
    
    return ax        
 

  
def plotHeightProf(nz,mcTable,inputPath,dicSettings):
  '''
  plot height profiles of superparticle number concentration. 
  '''
  
  dz = dicSettings['maxHeight']/nz
  Heightrange = np.arange(0,dicSettings['maxHeight'],dz)
  heightCenterBin = Heightrange[0:-1] + dz/2
  #binheight,sheight = pd.cut(mcTable['sHeight'],bins=Heightrange,retbins=True)
  #group according to velocity bins
  grouped = mcTable.groupby_bins('sHeight', Heightrange)
  #grouped = mcTable.groupby(binheight)
  Nsuper = grouped.count()['sMult']#.assign_coords({'sHeight':heightCenterBin})
  
  #height = sheight+dz/2
  fig,ax = plt.subplots(figsize=(5,5))
  ax.plot(Nsuper.values,heightCenterBin)
  ax.set_xlabel('# superparticle')
  ax.set_ylabel('height [m]')
  ax.grid()
  plt.tight_layout()
  plt.savefig(inputPath+'number_sParticle.png')
  plt.close()
  #print(Heightrange)
  #quit()  
  
def plotPSD(mcTable,dicSettings,inputPath,bins,var,gam=None,fac=1,xlim=None,ticks=None,heightEdge0=2900,unit='',sepMono=False,yscale=None,xscale=None):
  #for i, heightEdge0 in enumerate(dicSettings['heightRange'][::-1]):
    heightEdge1 = heightEdge0 + dicSettings['heightRes']
    height = heightEdge0+dicSettings['heightRes']/2 
    mcTableTmp = mcTable.where((mcTable['sHeight']>heightEdge0) &(mcTable['sHeight']<=heightEdge1),drop=True)
    
    fig,ax = plt.subplots()
    # make it per m3, so I need to multiply dh with boxwidth
    volume = dicSettings['gridBaseArea']*dicSettings['heightRes']
    if sepMono:
        mcTableMono = mcTableTmp.where(mcTableTmp.sNmono == 1,drop=True)
        mcTableAgg = mcTableTmp.where(mcTableTmp.sNmono > 1,drop=True)
         # this is the volume that the particles are in
        concMono = mcTableMono.sMult.sum()/volume
        concAgg = mcTableAgg.sMult.sum()/volume
        ax.hist(mcTableAgg[var]*fac,bins=100,weights=mcTableAgg.sMult/volume,label='Agg, tot conc: {0}'.format(int(concAgg)))
        ax.hist(mcTableMono[var]*fac,bins=100,weights=mcTableMono.sMult/volume,label='Mono, tot conc: {0}'.format(int(concMono)))
        ax.legend()
    else:
        ax.hist(mcTableTmp[var]*fac,bins=100,weights=mcTableTmp.sMult/volume)
      
    if yscale:
        ax.set_yscale(yscale)
    if xscale:
        ax.set_xscale(xscale)
    if xlim:
      ax.set_xlim(xlim)
    if ticks:
      ax.set_xticks(ticks)
    #ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel(var+' '+unit)
    ax.set_ylabel(r'#m$^{-3}$')
    ax.set_title('height: {}'.format(height))
    #ax.set_xlim(0,0.1)#[10**-3*10,0.5*10**0])
    plt.tight_layout()
    plt.savefig(inputPath+'{height}_{var}_hist_weigth.png'.format(height=height,var=var))
    plt.close()
    
    print(height,' done')
def calc_gamma_distr(nu_gam,mu_gam,IWC,nrp0,Mbins):
    #calculate gamma distribution how it is written in McSnow
    lam, gam_n1, gam_n2, mean_mass = calc_lam(nu_gam, mu_gam, IWC, nrp0)
  # calculate A
    A = calc_A(nu_gam, mu_gam, lam, nrp0, gam_n1)
 #calculate gamma distribution:
    return calc_gam(nu_gam, mu_gam, A, lam, Mbins)
    
def calc_lam(nu_gam, mu_gam, IWC, nrp0):
  # calculate lambda of gamma distribution: lambda = (gam_n1/(gam_n2*mean_mass))
  mean_mass = IWC/nrp0
  print((nu_gam+1)/mu_gam)
  print((nu_gam+2)/mu_gam)
  print(nu_gam,mu_gam)
  gam_n1 = math.gamma((nu_gam+1)/mu_gam)
  gam_n2 = math.gamma((nu_gam+2)/mu_gam)
  lam = (gam_n1/(gam_n2*mean_mass))**(-mu_gam)
  return lam, gam_n1, gam_n2, mean_mass
  
def calc_A(nu_gam, mu_gam, lam, nrp0, gam_n1):
  # calculate A of gamma distribution
  return (mu_gam*nrp0/gam_n1)*lam**((nu_gam+1)/mu_gam)
  
def calc_gam(nu_gam, mu_gam, A, lam, x):
    return A*x**nu_gam*np.exp(-lam*x**mu_gam)
    
    

		
