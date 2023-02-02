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



def plotPropSpecThesis(ax,heightRange,heightRes,mcTable,velBins,prop,savefig=False,cbar=True,gridBaseArea=5.0,diff=False,mcTable1=None):
# plots the properties depending on the fall velocity (so sort of like a Doppler Spectrum), so the data gets binned into velocity bins and weighted by the multiplicity and then plotted with height/temp on y-axis and velocity on x-axis
    
    for i, heightEdge0 in enumerate(heightRange[0:-1]):
        
        heightEdge1 = heightEdge0 + heightRes
        height = heightEdge0+heightRes/2
        
        mcTableTmp = mcTable[(mcTable['height']>heightEdge0) &
                             (mcTable['height']<=heightEdge1)].copy()
        if diff == True:
            mcTableTmp1 = mcTable1[(mcTable1['height']>heightEdge0.values) &
                     (mcTable1['height']<=heightEdge1.values)].copy()

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
                vmin=200;vmax=1000
                cbar_label = r'$\Delta$ mass [g]'
            else:
                weighted = grouped.apply(wavgRho)
                vmin = 200;vmax = 1000
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
    if diff==False: 
        if (prop == 'dia' or prop == 'number_conc' or prop == 'mTot'):
            plot=ax.pcolormesh(binnedXR.vel,binnedXR['T'],binnedXR.T,cmap=cmap,norm=colors.LogNorm(vmin=10**vmin,vmax=10**vmax))#,cmap=cmap,add_colorbar=False)
        else:
            plot=binnedXR.plot(ax=ax,x='vel',vmin=vmin,vmax=vmax,cmap=cmap,add_colorbar=False)
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

def plotMoments(dicSettings,output,inputPath,convoluted=False,minmax=None,plotTemp=False,mult_conc=False):
    for wl in dicSettings['wl']:
        wlStr = '{:.2e}'.format(wl)
        freq = (constants.c / float(wlStr))  *1e3 / 1e9
        print(freq)
        freq = '{:.1f}'.format(freq)
     
        if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
          
          if plotTemp == True:
              saveName = '1d_habit_moments_{0}_{1}_Temp.png'.format(freq,dicSettings['scatSet']['mode'])
              vary='Temp'; varUnit = '[°C]'
              ylim = [0,-30]
          else:
              saveName = '1d_habit_moments_{0}_{1}.png'.format(freq,dicSettings['scatSet']['mode'])
              vary='range'; varUnit = '[m]'
              ylim = [0,dicSettings['maxHeight']]
          
          print(output['Ze_H_{0}'.format(wlStr)])
              
          fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
          #output['MDV_H_{0}'.format(wlStr)].plot(ax=axes[0],y='range', lw=2)
          axes[0].plot(output['MDV_H_{0}'.format(wlStr)],output[vary],lw=2)
          axes[0].set_title('Freq: {0} elv: {1}, MDV'.format(freq, dicSettings['elv']),fontsize=16)
          #axes[0].set_ylim(0, 5000)
          axes[0].grid(True,ls='-.')
          axes[0].set_xlabel('MDV [m/s]',fontsize=16)
          axes[0].set_ylim(ylim)
          axes[0].set_ylabel(vary+' '+varUnit,fontsize=16)
          axes[0].tick_params(axis='both',labelsize=14)
          
          axes[1].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)]),output[vary],linewidth=2) # TODO: change back to ZeH
          axes[1].set_xlabel('Ze [dB]',fontsize=16)
          #axes[1].set_title('Ze_H',fontsize=16)
          #plt.ylim(0, 5000)
          #plt.xlim(-3, 0)
          axes[1].grid(True,ls='-.')
          axes[1].set_ylim(ylim)
          axes[1].tick_params(axis='both',labelsize=14)
          plt.tight_layout()
          
              
        else:
            if convoluted == True:
              saveName = '1d_habit_moments_{0}_convoluted_nace_40.png'.format(freq)
              fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
              # plot ZDR
              axes[0].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)])-mcr.lin2db(output['Ze_V_{0}'.format(wlStr)]),output['range'],linewidth=2)# TODO: change bach to ZeH-ZeV
              axes[0].set_xlabel('ZDR [dB]')
              axes[0].set_title('ZDR')
              #plt.xlim(-3, 0)
              axes[0].grid(True,ls='-.')
              axes[0].set_ylim([0,dicSettings['maxHeight']])
        #plt.savefig(inputPath+'1d_habit_ZDR_{0}.png'.format(wlStr), format='png', dpi=200, bbox_inches='tight')
        #plt.close()

              axes[1].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)]),output['range'],linewidth=2) # TODO: change back to ZeH
              axes[1].set_xlabel('Z_H [dB]')
              axes[1].set_title('Ze_H')
              #plt.ylim(0, 5000)
              #plt.xlim(-3, 0)
              axes[1].grid(True,ls='-.')
              axes[1].set_ylim([0,dicSettings['maxHeight']])
              plt.tight_layout()
              plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
              plt.close()
            else: 
              
              if plotTemp == True:
                  if minmax:
                      if mult_conc == True:
                          saveName = '1d_habit_moments_{wl}_{mode}_{minmax}_Temp_mult_conc.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],minmax=minmax)
                      else:
                          saveName = '1d_habit_moments_{wl}_{mode}_{minmax}_Temp.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],minmax=minmax)
                  else:
                      if mult_conc == True:
                          saveName = '1d_habit_moments_{wl}_{mode}_Temp_mult_conc.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'])
                      else:
                          saveName = '1d_habit_moments_{wl}_{mode}_Temp.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'])
              else:
                  if minmax:
                      saveName = '1d_habit_moments_{wl}_{mode}_{minmax}.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],minmax=minmax)
                  else:
                      saveName = '1d_habit_moments_{wl}_{mode}.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'])
        
        #plot KDP 
              if plotTemp == True:
                  vary = 'Temp'
                  ylim=([0,-30])
              else:
                  vary = 'range'
                  ylim= ([0,dicSettings['maxHeight']])
              fig,axes = plt.subplots(ncols=3,figsize=(15,5),sharey=True)
              #output['kdpInt_{0}'.format(wlStr)].plot(ax=axes[0],y='range', lw=2)
              axes[0].plot(output['kdpInt_{0}'.format(wlStr)],output[vary],lw=2)
              axes[0].set_title('rad: {0} elv: {1}, KDP'.format(freq, dicSettings['elv']))
              #axes[0].set_ylim(0, 5000)
              axes[0].set_ylabel(vary)
              axes[0].grid(True,ls='-.')
              axes[0].set_xlabel('KDP [°/km]')
              axes[0].set_ylim(ylim)
        
        # plot ZDR
              axes[1].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)])-mcr.lin2db(output['Ze_V_{0}'.format(wlStr)]),output[vary],linewidth=2)# TODO: change bach to ZeH-ZeV
              axes[1].set_xlabel('ZDR [dB]')
              axes[1].set_title('ZDR')
              #plt.xlim(-3, 0)
              axes[1].grid(True,ls='-.')
              #axes[1].set_ylim([0,dicSettings['maxHeight']])
              axes[1].set_ylim(ylim)
              
              axes[2].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)]),output[vary],linewidth=2) # TODO: change back to ZeH
              axes[2].set_xlabel('Z_H [dB]')
              axes[2].set_title('Ze_H')
              #plt.ylim(0, 5000)
              #plt.xlim(-3, 0)
              axes[2].grid(True,ls='-.')
              axes[2].set_ylim(ylim)
              #axes[2].set_ylim([0,dicSettings['maxHeight']])
              plt.tight_layout()
            
        plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
        plt.close()

def plotDWR(dicSettings,wlStr1,wlStr2,output,inputPath,convoluted=False,plotTemp=False):
    #for wl in dicSettings['wl']:
    #wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
    #wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
	freq1 = (constants.c / float(wlStr1))  *1e3 / 1e9
	freq2 = (constants.c / float(wlStr2))  *1e3 / 1e9
	freq1 = '{:.1f}'.format(freq1)
	freq2 = '{:.1f}'.format(freq2)
	if convoluted == True:
		if plotTemp == True:
			saveName = '1d_habit_DWR_{freq1}_{freq2}_convoluted_{mode}_Temp.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'Temp'; varUnit = '[°C]'
			ylim = [0,-30]
		else:
			saveName = '1d_habit_DWR_{freq1}_{freq2}_convoluted_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'range'; varUnit = '[m]'
			ylim = [0,dicSettings['maxHeight']]

		fig,axes = plt.subplots(figsize=(5,4))
		DWR = mcr.lin2db(output['Ze_H_{0}'.format(wlStr1)]) - mcr.lin2db(output['Ze_H_{0}'.format(wlStr2)])
		#DWR.plot(ax=axes,y='range', lw=2)
		axes.plot(DWR,output[vary],lw=2)
		#axes.set_title('DWR_{freq1}_{freq2}'.format(freq1=freq1,freq2=freq2))
		#axes[0].set_ylim(0, 5000)
		axes.grid(True,ls='-.')
		axes.set_xlabel(r'DWR$_{{{freq1},{freq2}}}$ [dB]'.format(freq1=freq1,freq2=freq2),fontsize=16)
		axes.set_ylim(ylim)
		axes.set_ylabel(vary+' '+varUnit,fontsize=16)
		axes.tick_params(axis='both',labelsize=14)
	  
	else: 
		if plotTemp == True:
			saveName = '1d_habit_DWR_{freq1}_{freq2}_{mode}_Temp.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'Temp'; varUnit = '[°C]'
			ylim = [0,-30]
		else:
			saveName = '1d_habit_DWR_{freq1}_{freq2}_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'range'; varUnit = '[m]'
			ylim = [0,dicSettings['maxHeight']]
		

		#saveName = '1d_habit_DWR_{freq1}_{freq2}_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
		   
		fig,axes = plt.subplots(figsize=(5,4))
		Ze1 = mcr.lin2db(output['spec_H_{0}'.format(wlStr1)].sum(dim='vel'))
		Ze2 = mcr.lin2db(output['spec_H_{0}'.format(wlStr2)].sum(dim='vel'))
		DWR = Ze1 - Ze2
		#DWR.plot(ax=axes,y='range', lw=2)
		axes.plot(DWR,output[vary],lw=2)
		axes.set_ylabel(vary+' '+varUnit,fontsize=16)
		axes.tick_params(axis='both',labelsize=14)
		#axes.set_title('DWR_{freq1}_{freq2}'.format(freq1=freq1,freq2=freq2))
		#axes[0].set_ylim(0, 5000)
		axes.grid(True,ls='-.')
		axes.set_xlabel(r'DWR$_{{{freq1},{freq2}}}$ [dB]'.format(freq1=freq1,freq2=freq2),fontsize=16)
		axes.set_ylim(ylim)
		                  
	plt.tight_layout()        
	plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
	plt.close()
def plotDWRspectra(dicSettings,wlStr1,wlStr2,output,inputPath,convoluted=False,plotTemp=False):
    #for wl in dicSettings['wl']:
    #wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
    #wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
	freq1 = (constants.c / float(wlStr1))  *1e3 / 1e9
	freq2 = (constants.c / float(wlStr2))  *1e3 / 1e9
	freq1 = '{:.1f}'.format(freq1)
	freq2 = '{:.1f}'.format(freq2)
	if convoluted == True:
		if plotTemp == True:
			saveName = '1d_habit_sDWR_{freq1}_{freq2}_convoluted_{mode}_Temp.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'Temp';varUnit = '[°C]'
			ylim = [0,-30]
		else:
			saveName = '1d_habit_sDWR_{freq1}_{freq2}_convoluted_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'range'
			varUnit = '[m]'
			ylim = [0,dicSettings['maxHeight']]

		fig,axes = plt.subplots(figsize=(5,4))
		specH1 = mcr.lin2db(output['spec_H_{0}'.format(wlStr1)])
		specH1 = specH1.where(specH1 > -40)
		specH2 = mcr.lin2db(output['spec_H_{0}'.format(wlStr2)])
		specH2 = specH2.where(specH2 > -40)

		DWR =specH1 - specH2
		#DWR.plot(ax=axes,y='range', vmin=0, vmax=15, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sDWR [dB]'})
		plot = axes.pcolormesh(output.vel,output[vary],DWR,vmin=0, vmax=15, cmap=getNewNipySpectral())
		cb = plt.colorbar(plot,ax=axes,pad=0.02,aspect=20)#,ticks=v1)
		cb.set_label(r'sDWR$_{{{freq1},{freq2}}}$'.format(freq1=freq1,freq2=freq2),fontsize=16)
		cb.ax.tick_params(labelsize=14)
		axes.set_xlabel('Doppler velocity [m/s]',fontsize=16)
		axes.set_ylabel(vary+' '+varUnit,fontsize=16)
		axes.tick_params(axis='both',labelsize=14)
		#axes[0].set_ylim(0, 5000)
		axes.grid(True,ls='-.')
		#axes.set_xlabel('DWR [dB]')
		axes.set_ylim(ylim)
		axes.set_xlim(-2, 0)
	else: 
		if plotTemp == True:
			saveName = '1d_habit_sDWR_{freq1}_{freq2}_{mode}_Temp.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'Temp';varUnit = '[°C]'
			ylim = [0,-30]
		else:
			saveName = '1d_habit_sDWR_{freq1}_{freq2}_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
			vary = 'range';varUnit = '[m]'
			ylim = [0,dicSettings['maxHeight']]

		#saveName = '1d_habit_sDWR_{freq1}_{freq2}_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])

		fig,axes = plt.subplots(figsize=(5,4))
		specH1 = mcr.lin2db(output['spec_H_{0}'.format(wlStr1)])
		specH1 = specH1.where(specH1 > -40)
		specH2 = mcr.lin2db(output['spec_H_{0}'.format(wlStr2)])
		specH2 = specH2.where(specH2 > -40)

		DWR =specH1 - specH2
		plot = axes.pcolormesh(output.vel,output[vary],DWR,vmin=0, vmax=15, cmap=getNewNipySpectral())
		cb = plt.colorbar(plot,ax=axes,pad=0.02,aspect=20)#,ticks=v1)
		cb.set_label(r'sDWR$_{{{freq1},{freq2}}}$'.format(freq1=freq1,freq2=freq2),fontsize=16)
		cb.ax.tick_params(labelsize=14)
		axes.set_xlabel('Doppler velocity [m/s]',fontsize=16)
		axes.set_ylabel(vary+' '+varUnit,fontsize=16)
		axes.tick_params(axis='both',labelsize=14)
		#axes[0].set_ylim(0, 5000)
		axes.grid(True,ls='-.')
		#axes.set_xlabel('DWR [dB]')
		axes.set_ylim(ylim)
		axes.set_xlim(-2, 0)                      
	plt.tight_layout()        
	plt.savefig(inputPath+saveName, format='png', dpi=200)#, bbox_inches='tight')
	plt.close()

def plotSpectra(dicSettings,output,inputPath,convoluted=False,minmax=None,plotTemp=False):
    for wl in dicSettings['wl']:
        wlStr = '{:.2e}'.format(wl)
        freq = (constants.c / float(wlStr))  *1e3 / 1e9
        freq = '{:.1f}'.format(freq)
        
        if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
            
            if plotTemp == True:
                saveName = '1d_habit_spectra_{wl}_{mode}_{part}_{elv}_Temp.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],
                																			part=dicSettings['scatSet']['particle_name'],elv=dicSettings['elv'])
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                vary='Temp';varUnit = '[°C]'
                ylim = [0,-30]
            else:
                saveName = '1d_habit_spectra_{wl}_{mode}_{part}_{elv}.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],
                																		part=dicSettings['scatSet']['particle_name'],
                																		elv=dicSettings['elv'])
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                vary = 'range'; varUnit = '[m]'
                ylim = [0,dicSettings['maxHeight']]
            fig,ax = plt.subplots(figsize=(5,4))
            #specH.plot(ax=ax,vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
            
            plot = ax.pcolormesh(output.vel,output[vary],specH,cmap=getNewNipySpectral(),vmin=-30,vmax=5)
            cb = plt.colorbar(plot,ax=ax,pad=0.02,aspect=20)#,ticks=v1)
            cb.set_label(r'sZe$_{\rm Ka}$',fontsize=16)
            cb.ax.tick_params(labelsize=14)
            ax.set_xlabel('Doppler velocity [m/s]',fontsize=16)
            ax.set_ylabel(vary+' '+varUnit,fontsize=16)
            ax.tick_params(axis='both',labelsize=14)
            #ax.set_title('Ze_H')
            #ax.set_title('Ze_H_spec rad: {0} GHz, elv: {1}'.format(freq, dicSettings['elv']))# TODO: change back to ZeH or to ZeV if you use first LUT of dendrites
            #axes[0].set_ylim([0,dicSettings['maxHeight']])
            ax.set_ylim(ylim)
            ax.set_xlim(-2, 0)
            ax.grid(True,ls='-.')      
            plt.tight_layout()
        else:
            if convoluted == True:
                saveName = '1d_habit_spectra_{0}_convoluted.png'.format(freq)
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                specV = mcr.lin2db(output['spec_V_{0}'.format(wlStr)])
                specH = specH.where(specH > -40)
                specV =  specV.where(specV > -40)
                #dataSmooth = specTable.rolling(vel=10,min_periods=1,center=True).mean()            
                ZDR = specH.rolling(vel=10,min_periods=1,center=True).mean() - specV.rolling(vel=10,min_periods=1,center=True).mean()
            else: 
                if plotTemp == True:
                    if minmax:
                        saveName = '1d_habit_spectra_{wl}_{mode}_{minmax}_Temp.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],minmax=minmax)
                    else:
                        saveName = '1d_habit_spectra_{wl}_{mode}_Temp.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'])
                else:
                    if minmax:
                        saveName = '1d_habit_spectra_{wl}_{mode}_{minmax}.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'],minmax=minmax)
                    else:
                        saveName = '1d_habit_spectra_{wl}_{mode}.png'.format(wl=freq,mode=dicSettings['scatSet']['mode'])
        
                
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                specV = mcr.lin2db(output['spec_V_{0}'.format(wlStr)])
                ZDR = specH - specV
            if plotTemp == True:
                vary = 'Temp'
                ylim=([0,-30])
            else:
                vary = 'range'
                ylim= ([0,dicSettings['maxHeight']])
            fig,axes = plt.subplots(ncols=3,figsize=(12,5),sharey=True)
            #specH.plot(ax=axes[0],vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
            
            p1 = axes[0].pcolormesh(output.vel,
                                    output[vary],
                                    specH,vmin=-30,vmax=5,cmap=getNewNipySpectral())
            cb = plt.colorbar(p1,ax=axes[0])
            cb.set_label('sZeH [dBz]')
            axes[0].set_xlabel('vel [m/s]')
            axes[0].set_title('Ze_H')
            axes[0].set_title('Ze_H_spec rad: {0} GHz, elv: {1}'.format(freq, dicSettings['elv']))
            #axes[0].set_ylim([0,dicSettings['maxHeight']])
            axes[0].set_ylim(ylim)
            axes[0].set_ylabel(vary)
            axes[0].set_xlim(-2, 0)
            axes[0].grid(True,ls='-.')
        
            #specV.plot(ax=axes[1],vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
            p2 = axes[1].pcolormesh(output.vel,
                                    output[vary],
                                    specV,vmin=-30,vmax=5,cmap=getNewNipySpectral())
            cb = plt.colorbar(p2,ax=axes[1])
            cb.set_label('sZeV [dBz]')
            axes[1].set_title('Ze_V_spec rad: {0} GHz, elv: {1}'.format(freq, dicSettings['elv']))
            #axes[1].set_ylim([0,dicSettings['maxHeight']])
            axes[1].set_xlim(-2, 0)
            axes[1].grid(True,ls='-.')
            axes[1].set_ylabel('')
            axes[1].set_xlabel('vel [m/s]')
            axes[1].set_ylim(ylim)
            
            #ZDR.plot(ax=axes[2],vmin=-0.5, vmax=3,cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZDR [dB]'}) 
            p3 = axes[2].pcolormesh(output.vel,
                                    output[vary],
                                    ZDR,vmin=-0.5,vmax=3,cmap=getNewNipySpectral())
            cb = plt.colorbar(p3,ax=axes[2])
            cb.set_label('sZDR [dB]')
            axes[2].set_title('ZDR rad: {0} GHz, elv: {1}'.format(freq, dicSettings['elv']))
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
    
def plotDWRsObs(dataLV2,outPath,outName):
    DWR_KaW = dataLV2.Ka_DBZ_H - dataLV2.W_DBZ_H
    DWR_XKa = dataLV2.X_DBZ_H - dataLV2.Ka_DBZ_H
    fig,axes = plt.subplots(nrows=2,figsize=(18,12),sharey=True)
    plot = DWR_KaW.plot(ax=axes[0],x='time',vmin=-1,vmax=10,cmap='jet',add_colorbar=False)
    cb = plt.colorbar(plot,ax=axes[0]) 
    cb.set_label('DWR-KaW [dB]',fontsize=18)
    cb.ax.tick_params(labelsize=16)
    CS = dataLV2.ta.plot.contour(ax=axes[0],y='range',levels=[-15,0],colors='k',linewidths=[3,3],linestyles=['--','-'])
    plt.clabel(CS, CS.levels, fontsize=22, fmt='%1.f '+'C')
    axes[0].set_title('DWR-KaW',fontsize=20)
    plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
    axes[0].grid()
    axes[0].set_xlabel('')
    axes[0].tick_params(axis='y',labelsize=16)
    axes[0].set_ylabel('range [m]',fontsize=20) 

    plot = DWR_XKa.plot(ax=axes[1],x='time',vmin=-1,vmax=10,cmap='jet',add_colorbar=False)
    cb = plt.colorbar(plot,ax=axes[1]) 
    cb.set_label('DWR-XKa [dB]',fontsize=18)
    cb.ax.tick_params(labelsize=16)
    CS = dataLV2.ta.plot.contour(ax=axes[1],y='range',levels=[-15,0],colors='k',linewidths=[3,3],linestyles=['--','-'])
    plt.clabel(CS, CS.levels, fontsize=22, fmt='%1.f '+'C')
    axes[1].set_title('DWR-XKa',fontsize=20)
    plt.setp(plot.axes.xaxis.get_majorticklabels(), rotation=0)
    axes[1].grid()
    axes[1].set_xlabel('')
    axes[1].tick_params(axis='both',labelsize=16)
    axes[1].set_ylabel('range [m]',fontsize=20) 

    plt.savefig(outPath+outName+'.png',dpi=200,bbox_inches='tight')
    plt.close()
        

def plotSpectraObs(dataPol,dataLV0,dataLV2,outPath,ylim=False):
    for t in dataLV0.time:
        print(t)
        ti = pd.to_datetime(str(t.values)).strftime('%Y%m%d_%H%M%S')
        spec = 10*np.log10(dataLV0.KaSpecH.sel(time=t)) 
        #specV = 10*np.log10(dataLV0.KaSpecV.sel(time=t))
        #LDR = specV - specH
        #print(spec)
        #quit()
        dataPolSel = dataPol.sel(time=t)
        dataLV2Sel = dataLV2.sel(time=t)
        fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
        #spec.plot(ax = axes[0],y='range',vmin=-30,vmax=5,cmap=getNewNipySpectral())
        plot1=axes[0].pcolormesh(spec.dopplerKa, # this is needed because we measure at an elevation angle, therefore if we would use the Doppler Velocity as x-axis, the spectra looks like a snake, so everything is moved to zero. Does not always work because of uncertainty of getting the spectral edge velocity correctly
                       #dataPolSel.height,
                       dataLV2Sel.ta.values,
                       spec.values,
                       vmin=-30,vmax=10,cmap=getNewNipySpectral())
        cb1 = plt.colorbar(plot1,ax=axes[0]) 
        cb1.set_label('Ze [dB]')
        axes[0].set_ylabel('height [m]')
        axes[0].set_xlabel('Doppler vel [m/s]')
        axes[0].grid(True,ls='-.')
        axes[0].set_xlim([-2.5,1.5])
        if ylim:
            axes[0].set_ylim(ylim)
        plot=axes[1].pcolormesh(dataPolSel.Vel2ZeroH.fillna(0), # this is needed because we measure at an elevation angle, therefore if we would use the Doppler Velocity as x-axis, the spectra looks like a snake, so everything is moved to zero. Does not always work because of uncertainty of getting the spectral edge velocity correctly
                       #dataPolSel.height,
                       dataLV2Sel.ta.values,
                       dataPolSel.sZDR.values,
                       vmin=-1,vmax=3,cmap=getNewNipySpectral())
        cb = plt.colorbar(plot,ax=axes[1]) 
        cb.set_label('ZDR [dB]')
    #cb.ax.tick_params(labelsize=16)
        axes[1].set_xlabel('corrected Doppler vel [m/s]')
        axes[1].grid(True,ls='-.')
        axes[1].set_xlim([-2.5,0.5])
        if ylim:
            axes[1].set_ylim(ylim)
        plt.tight_layout()
        plt.savefig(outPath+ti+'_spectra_5min_mean.png')
        plt.close()
        print(ti,' finished')    
        
def plotProfilesObs(dataLV2,dataPol,outPath,ylim=False):
    for t in dataPol.time:
        KDP = dataPol.KDP.sel(time=t)
        ZDR = dataPol.ZDR.sel(time=t)
        Ze = dataLV2.Ka_DBZ_H.sel(time=t)
        ti = pd.to_datetime(str(t.values)).strftime('%Y%m%d_%H%M%S')
        fig,axes = plt.subplots(ncols=3,figsize=(15,5),sharey=True)
        Ze.plot(ax=axes[0],y='range',LineWidth=3)
        axes[0].set_ylabel('height [m]')
        axes[0].set_xlabel('Ze [dB]')
        axes[0].grid(True,ls='-.')
        ZDR.plot(ax=axes[1],y='height',LineWidth=3)
        axes[1].set_xlabel('ZDR [dB]')
        axes[1].set_ylabel('')
        axes[1].grid(True,ls='-.')
        KDP.plot(ax=axes[2],y='height',LineWidth=3)
        axes[2].set_xlabel('KDP [°/km]')
        axes[2].set_ylabel('')
        axes[2].grid(True,ls='-.')
        axes[2].set_xlim([0,4])
        if ylim:
            axes[0].set_ylim(ylim)
            axes[1].set_ylim(ylim)
            axes[2].set_ylim(ylim)
        plt.tight_layout()
        plt.savefig(outPath+ti+'_moments_profile_5min_mean.png')
        plt.close()
        print(ti,' finished')    
    
def plotProfilesObsDWR(dataLV2,outPath,ylim=False):
    for t in dataLV2.time:
        DWR_KaW = dataLV2.Ka_DBZ_H.sel(time=t) - dataLV2.W_DBZ_H.sel(time=t)
        DWR_XKa = dataLV2.X_DBZ_H.sel(time=t) - dataLV2.Ka_DBZ_H.sel(time=t)
        ti = pd.to_datetime(str(t.values)).strftime('%Y%m%d_%H%M%S')
        fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
        DWR_KaW.plot(ax=axes[0],y='range',LineWidth=3)
        axes[0].set_ylabel('height [m]')
        axes[0].set_xlabel('DWR-KaW [dB]')
        axes[0].grid(True,ls='-.')
        axes[0].set_xlim([0,9])
        DWR_XKa.plot(ax=axes[1],y='range',LineWidth=3)
        axes[1].set_xlabel('DWR-XKa [dB]')
        axes[1].set_ylabel('')
        axes[1].grid(True,ls='-.')
        axes[1].set_xlim([0,9])
        if ylim:
            axes[0].set_ylim(ylim)
            axes[1].set_ylim(ylim)
        plt.tight_layout()
        plt.savefig(outPath+ti+'_DWR_profile_5min_mean.png')
        plt.close()
        print(ti,' finished')    
        

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
 


def plot_var_particle_ID_temp(mcTable,prop,var2col,varName,unit,zoom=False,onlySpecPart=False,log=False,Dini=False,ax=None,cbar=True,xlabel=True,save=False,inputPath=None):
  # this plots the particles according to their ID and height. This is to test the evolution of different parameters with respect to height. Input: mcTable, path to save the plot and the variable to plot (prop)
  # zoom needs to  be array of ylim
  # first lets make temp like height
    
  if onlySpecPart:
    n = onlySpecPart[1]-onlySpecPart[0]
    norm = mpl.colors.Normalize(vmin=onlySpecPart[0], vmax=onlySpecPart[1])
    figName = '_Temp_index_selSpecPart{0}_{1}.png'.format(onlySpecPart[0],onlySpecPart[1])
  else:
    n = max(mcTable[var2col])
    norm = mpl.colors.Normalize(vmin=0, vmax=n)
    figName = '_Temp_index.png'
  mcTable = mcTable.set_index(var2col,drop=False)
  if log==True:
    figName = '_log'+figName
  #c_m = mpl.cm.gist_rainbow
  top = cm.get_cmap('cool',170)
  bottom = cm.get_cmap('autumn',180)
  top =  top(np.linspace(0, 1, 170))
  bottom = bottom(np.linspace(0, 1, 180))
  

  newcolors = np.vstack((bottom,
                           top))
  cmap = ListedColormap(newcolors)
  # create a ScalarMappable and initialize a data structure
  s_m = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
  s_m.set_array([])
  step = np.floor(n/5)
  step = step 
  
  tickmarks = np.arange(0,n,step)
  
  if not ax:
    fig, ax = plt.subplots(figsize=(9,5))
    save=True
  else:
    save=False
  for sID in np.unique(mcTable[var2col]):
    mcTableTmp = mcTable.loc[sID]
    plot=ax.plot(mcTableTmp[prop],mcTableTmp['Temp'],color = s_m.to_rgba(sID-1),alpha=0.5)
  if ('dia' in prop) or (log==True):
    ax.set_xscale('log')
  #cax = fig.add_axes([0.95, 0.1, 0.02, 0.75])
  if var2col == 'time':
    var2col='time'
  if var2col == 'sMult':
    var2col = 'particle_ID'
  if cbar == True:
    cbar = plt.colorbar(s_m,ax=ax)#, pad=0.15)
    if Dini:
      DmaxIni = (10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,300.,400.,500.,600.,700.,800.,900.,1000.)
      cbar.set_ticks([t for t in tickmarks])
      cbar.set_ticklabels([DmaxIni[int(tickmarks[t])] for t in range(len(tickmarks))])
      var2col='Dini'
      cbar.set_label(r'D$_{\rm ini}$ [$\mu$m]',fontsize=16)
    else:
      cbar.set_label(var2col,fontsize=16)
    cbar.ax.tick_params(labelsize=16)
  if xlabel:
    ax.set_ylabel('T [°C]',fontsize=18)
  else:
    ax.set_yticklabels('')
 
  ax.set_xlabel(varName+' '+unit,fontsize=18)
  #ax1.set_xlim([-20,0])
  if zoom:
    plt.xlim([zoom[0],zoom[1]])
    prop=prop+'_zoom'
  #plt.xlim([-20,0])
  ax.grid()
  #ax2.grid(ls='-.')
  ax.tick_params(axis='both',labelsize=16)
  
  if save==True:
    plt.tight_layout()
    plt.savefig(inputPath+prop+'_'+var2col+figName)
    plt.close()
  else:
    return ax,s_m
    
def plotInitempVar(mcTable,ax,inputPath,prop,var2col,unit,varName,zoom=False,ylog=False,xlog=False,color='C0',xlim=None,ylim=None):
  try:
    RHi=np.float(inputPath.split('ssat')[1].split('/')[0])*0.1+100
  except:
    RHi=np.float(inputPath.split('ssat')[1].split('_')[0])*0.1+100
  #RHi=np.float(inputPath.split('RHi')[1].split('_')[0])
  #print(RHi)
  #quit()
  mcTableSmult = mcTable.set_index('sMult',drop=False)
  iniDia = np.array([10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,160.,180.,200.,300.,400.,500.,600.,700.,800.,900.,1000.])*1e-6 
  
  for i,sID in enumerate(np.sort(mcTableSmult['sMult'].unique())):
    
      
    mcTableTmp = mcTableSmult.loc[sID]
    
    #mcTableTmp = mcTableTmp.sort_values('Temp')
    mcTableTmp = mcTableTmp.where(mcTableTmp['Temp']<=-10)
    mcTableTmp = mcTableTmp.sort_values('time')
    mcTableTmp = mcTableTmp.set_index('time')
    
    #iniTemp = mcTableTmp.Temp.min()
    iniTemp = mcTableTmp.Temp.iloc[0]
    #print(mcTableTmp.iloc[0])
    
    iniDF = mcTableTmp.loc[mcTableTmp.Temp.idxmin()]
    #print(iniDF)
    #quit()
    if var2col == 'Temp':
      ini = iniTemp
    elif var2col == 'dia':
      ini = iniDia[i]      
    else:
      ini = iniDF[var2col]
    
    if i == 2:
      label=r'RH$_{\rm i}$'+'={0}%'.format(int(RHi))
    else:
      label = '_nolegend_'
    #label='_nolegend_'
    
    dfAt10 = mcTableTmp.loc[np.abs(mcTableTmp.Temp - (-10)).idxmin()]
    
    if ~np.isnan(ini):
      ax.plot(ini,dfAt10[prop],'.',ms=6,color=color,label=label)
      #ax.plot(dfAt10[prop],ini,'.',ms=6,color=color,label=label)
      
     
    
  if ('dia' in prop) or (ylog==True):
    ax.set_yscale('log')
    #ax1.set_yscale('log')
  if (ylog==True):
    ax.set_yscale('log')
  if xlim:
    ax.set_xlim(xlim)
  if ylim:
    #ax1.set_ylim(ylim)
    ax.set_ylim(ylim)
  
  
  ax.set_ylabel(varName+' at -10°C '+unit,fontsize=18)
  #ax1.set_ylabel(r'$|\Delta$('+varName+')| '+unit,fontsize=18)
  ax.set_xlabel(r'T$_{\rm ini}$ [°C]',fontsize=18)
  ax.tick_params(axis='both',labelsize=16)
  #ax1.tick_params(axis='both',labelsize=16)
  ax.grid()
  
  
  return ax

  
def plotHeightProf(nz,mcTable,inputPath,dicSettings):
  '''
  plot height profiles of superparticle number concentration. 
  '''
  
  dz = dicSettings['maxHeight']/nz
  Heightrange = np.arange(0,dicSettings['maxHeight'],dz)
  binheight,sheight = pd.cut(mcTable['sHeight'],bins=Heightrange,retbins=True)
  #group according to velocity bins
  grouped = mcTable.groupby(binheight)
  Nsuper = grouped['sMult'].count()
  
  height = sheight+dz/2
  fig,ax = plt.subplots(figsize=(5,5))
  ax.plot(Nsuper.values,height[0:-1])
  ax.set_xlabel('# superparticle')
  ax.set_ylabel('height [m]')
  ax.grid()
  plt.tight_layout()
  plt.savefig(inputPath+'number_sParticle.png')
  plt.close()
  #print(Heightrange)
  #quit()  
  
def plotPSD(mcTable,dicSettings,inputPath,bins,var,gam=None,fac=1,xlim=None,ticks=None):
  #for i, heightEdge0 in enumerate(dicSettings['heightRange'][::-1]):
    heightEdge0 = 2900
    heightEdge1 = heightEdge0 + dicSettings['heightRes']
    height = heightEdge0+dicSettings['heightRes']/2 
    mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &(mcTable['sHeight']<=heightEdge1)].copy()
    
    fig,ax = plt.subplots()
    ax.hist(mcTableTmp[var]*fac,bins=200,weights=mcTableTmp.sMult)
    #if gam.any():
      # you need to plot it in mass
    #  ax.plot(bins,gam)
      
    ax.set_yscale('log')
    if xlim:
      ax.set_xlim(xlim)
    if ticks:
      ax.set_xticks(ticks)
    #ax.set_yscale('log')
    ax.grid()
    ax.set_xlabel(var)
    ax.set_ylabel('#')
    #ax.set_xlim(0,0.1)#[10**-3*10,0.5*10**0])
    plt.tight_layout()
    plt.savefig(inputPath+'{height}_{var}_hist_weigth_ylog.png'.format(height=height,var=var))
    plt.close()
    
    print(height,' done')
    
def plot_scat_particle_ID_temp(mcTable,inputPath,prop,var2col,varName,unit,zoom=False,onlySpecPart=False,log=False,xlog=False,ax=None,cbar=True,y_label=True,separate=None):
  # this plots the particles according to their ID and height. This is to test the evolution of different parameters with respect to height. Input: mcTable, path to save the plot and the variable to plot (prop)
  # zoom needs to  be array of ylim
  # first lets make temp like height
    
  if onlySpecPart:
    n = onlySpecPart[1]-onlySpecPart[0]
    norm = mpl.colors.Normalize(vmin=onlySpecPart[0], vmax=onlySpecPart[1])
    figName = '_Temp_index_selSpecPart{0}_{1}.png'.format(onlySpecPart[0],onlySpecPart[1])
  else:
    n = 350
    norm = mpl.colors.Normalize(vmin=0, vmax=n)
    figName = '_Temp_index.png'
  #mcTable = mcTable.set_index(var2col,drop=False)
  if log==True:
    figName = '_log'+figName
  #c_m = mpl.cm.nipy_spectral
  #c_m = getNewNipySpectral()
  top = cm.get_cmap('cool',170)
  bottom = cm.get_cmap('autumn',180)
  top =  top(np.linspace(0, 1, 170))
  bottom = bottom(np.linspace(0, 1, 180))
  

  newcolors = np.vstack((bottom,
                           top))
  c_m = ListedColormap(newcolors)
  # create a ScalarMappable and initialize a data structure
  s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
  s_m.set_array([])
  step = np.floor(n/5)
  step = step 
  
  tickmarks = np.arange(0,n,step)
  
  if not ax:
    fig, ax = plt.subplots(figsize=(9,5))
    save=True
  else:
    save=False
  if separate:
    ax1 = ax.twiny()
    #mcTableTmp1 = mcTable.sel(sMult=slice(0,separate))
    #mcTableTmp2 = mcTable.sel(sMult=slice(separate,350))
    for sID in mcTable[var2col]:
      mcTableTmp = mcTable.sel(sMult=sID)
      prop2plot = mcTableTmp[prop]
      TempNaN = mcTableTmp.Temp.where(~np.isnan(prop2plot))
      if sID <= separate:
        plot=ax.plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),color = s_m.to_rgba(sID-1),alpha=0.5)
      else:
        plot1=ax1.plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),color = s_m.to_rgba(sID-1),alpha=0.5)
    plot = plot + plot1
  else:
    for sID in mcTable[var2col]:
      mcTableTmp = mcTable.sel(sMult=sID)
      prop2plot = mcTableTmp[prop]
      TempNaN = mcTableTmp.Temp.where(~np.isnan(prop2plot))
    
      plot=ax.plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),color = s_m.to_rgba(sID-1),alpha=0.5)
    
  if ('dia' in prop) or (log==True):
    ax.set_yscale('log')
  if (xlog==True):
    ax.set_xscale('log')
  #cax = fig.add_axes([0.95, 0.1, 0.02, 0.75])
  if var2col == 'time':
    var2col='time'
  if var2col == 'sMult':
    var2col = 'particle_ID'
  if cbar == True:
    cbar = plt.colorbar(s_m,ax=ax)#, pad=0.15)
    cbar.set_label(var2col,fontsize=16)
    cbar.ax.tick_params(labelsize=16)
  ax.set_ylim([0,-31.5])
  if y_label:
    ax.set_ylabel('T [°C]',fontsize=18)
  else:
    ax.set_yticklabels('')
  
  ax.set_xlabel(varName+' '+unit,fontsize=18)
  if prop=='sZDR':
    ax.set_xlim([-0.2,7])
  ax.grid()
  ax.tick_params(axis='both',labelsize=16)
  #plt.tight_layout()
  if save==True:
    plt.savefig(inputPath+prop+'_'+var2col+figName)
    plt.close()
  else:
    return ax,s_m
    
def plotRadar(data,var,ax,diff=False,data1=None):
	if var =='sZe':
		wl = (constants.c / (35.5*1e9)) * 1e3
		wlStr = '{:.2e}'.format(wl)
		if diff:
			specH = data['spec_H_{0}'.format(wlStr)]
			specH = specH.where(specH > 10**(-40/10))
			specH1 = data1['spec_H_{0}'.format(wlStr)]
			specH1 = specH1.where(specH1 > 10**(-40/10))
			specH = specH1.fillna(0) - specH.fillna(0)
			specH = specH.where(specH!=0,np.nan) #.sub(specH1, fill_value=0)
			specHpos = mcr.lin2db(specH.where(specH>0))
			specHneg = mcr.lin2db(-1*specH.where(specH<0))
			clabel=r'$\Delta$sZe'
			vmin=-20;vmax=0
			cmap = 'coolwarm'
			p1 = ax.pcolormesh(data.vel,data.Temp,specHpos,vmin=vmin,vmax=vmax,cmap='Reds')
			p2 = ax.pcolormesh(data.vel,data.Temp,specHneg,vmin=vmin,vmax=vmax,cmap='Blues')
			cb = plt.colorbar(p1,ax=ax,pad=0.001,aspect=40)
			cb.set_label('$\Delta$sZe [dB]',fontsize=16)
			cb.ax.tick_params(labelsize=14)
			cb.set_ticks([-20,-15,-10,-5,0])
			cb = plt.colorbar(p2,ax=ax,pad=0.02,aspect=40)
			#cb.set_label('$\Delta$sZe, negative',fontsize=16)
			cb.ax.set_yticklabels('')#tick_params(labelsize=14)
			cb.set_ticks([-20,-15,-10,-5,0])
		else:
			specH = data['spec_H_{0}'.format(wlStr)]
			specH = mcr.lin2db(specH)
			specH = specH.where(specH > -40)
			vmin=-30;vmax=5
			clabel='sZe [dBz]'
			cmap=getNewNipySpectral()
			
			p1 = ax.pcolormesh(data.vel,data.Temp,specH,vmin=vmin,vmax=vmax,cmap=cmap)
		
			cb = plt.colorbar(p1,ax=ax,pad=0.02,aspect=20)
			cb.set_label(clabel,fontsize=16)
			cb.ax.tick_params(labelsize=14)
		ax.set_xlabel(r'$v_D$ [ms$^{-1}$]',fontsize=18)
		ax.set_xlim([-2,0])
	elif var=='DWR':
		wl1 = (constants.c / (35.5*1e9)) * 1e3
		wlStr1 = '{:.2e}'.format(wl1)
		wl2 = (constants.c / (94.0*1e9)) * 1e3
		wlStr2 = '{:.2e}'.format(wl2)
		wl3 = (constants.c / (9.6*1e9)) * 1e3
		wlStr3 = '{:.2e}'.format(wl3)
		if diff:
			ZeH11 = data1['Ze_H_{0}'.format(wlStr1)]
			ZeH1 = data['Ze_H_{0}'.format(wlStr1)]
			#ZeH1 = ZeH1.fillna(0) - ZeH11.fillna(0) #.sub(ZeH11,fill_value=0)
			#ZeH1 = ZeH1.where(ZeH1!=0,np.nan)
			
			ZeH21 = data1['Ze_H_{0}'.format(wlStr2)]
			ZeH2 = data['Ze_H_{0}'.format(wlStr2)]
			#ZeH2 = ZeH2.fillna(0) - ZeH21.fillna(0) #.sub(ZeH11,fill_value=0)
			#ZeH2 = ZeH2.where(ZeH2!=0,np.nan)
			DWR1 = mcr.lin2db(ZeH1)-mcr.lin2db(ZeH2)
			DWR2 = mcr.lin2db(ZeH11)-mcr.lin2db(ZeH21)
			DWR = DWR2-DWR1
			xlim=[-7,7]
			xlabel=r'$\Delta$DWR$_{\rm KaW}$ [dB]'
			xticks = [-5,-2.5,0,2.5,5]
			ax.plot(DWR,data['Temp'],lw=2)
		else:
			ZeH1 = data['Ze_H_{0}'.format(wlStr1)]
			ZeH2 = data['Ze_H_{0}'.format(wlStr2)]
			ZeH3 = data['Ze_H_{0}'.format(wlStr3)]
			DWR1 = mcr.lin2db(ZeH1) - mcr.lin2db(ZeH2)
			DWR2 = mcr.lin2db(ZeH3) - mcr.lin2db(ZeH1)
			xlim=[0,10]
			xticks = [0,2.5,5,7.5,10]
			xlabel=r'DWR$_{\rm KaW}$ [dB]'
			ax.plot(DWR1,data['Temp'],lw=2,label=r'DWR$_{\rm KaW}$')
			ax.plot(DWR2,data['Temp'],lw=2,label=r'DWR$_{\rm XKa}$')
			ax.legend(fontsize=14)
		#DWR.plot(ax=axes,y='range', lw=2)
		
		ax.set_xlabel(xlabel,fontsize=18)
		ax.set_xlim(xlim)
	elif var=='ZeMDV':
		wl = (constants.c / (35.5*1e9)) * 1e3
		wlStr = '{:.2e}'.format(wl)
		if diff:
			ZeH1 = mcr.lin2db(data1['Ze_H_{0}'.format(wlStr)])
			ZeH = mcr.lin2db(data['Ze_H_{0}'.format(wlStr)])
			ZeH = ZeH1 - ZeH #.sub(ZeH11,fill_value=0)
			
			xlabel=r'$\Delta$Ze [dB]'
			xlim=[-15,10]
			MDV = data['MDV_H_{0}'.format(wlStr)]
			MDV1 = data1['MDV_H_{0}'.format(wlStr)]
			MDV = MDV1-MDV#.sub(ZeH11,fill_value=0)
			xlabel1 = r'$\Delta$MDV [ms$^{-1}$]'
			x1lim=[-0.3,0.3]
			loc='upper right'
			
		else:
			ZeH = mcr.lin2db(data['Ze_H_{0}'.format(wlStr)])
			MDV = data['MDV_H_{0}'.format(wlStr)]
			xlim=[-20,10]
			x1lim=[-1.5,-0.5]
			xlabel='Ze [dB]'
			xlabel1 = r'MDV [ms$^{-1}$]'
			loc='upper center'
			#ls1 = ax.plot(ZeH,data.Temp,lw=2,c='C0',label='Ze')
		ls1 = ax.plot(ZeH,data.Temp,lw=2,c='C0',label='Ze')
		ax.set_xlabel(xlabel,fontsize=18)
		ax.set_xlim(xlim)
		ax1 = ax.twiny()
		ls2 = ax1.plot(MDV,data.Temp,lw=2,c='C1',label='MDV')
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
		ax1.set_xlabel(xlabel1,fontsize=18)
		ax1.set_xlim(x1lim)
		ls = ls1+ls2
		labs = [l.get_label() for l in ls]
		ax.legend(ls,labs,loc=loc,fontsize=14)
		
	elif var=='Ze':
		wl = (constants.c / (35.5*1e9)) * 1e3
		wlStr = '{:.2e}'.format(wl)
		if diff:
			ZeH1 = mcr.lin2db(data1['Ze_H_{0}'.format(wlStr)])
			ZeH = mcr.lin2db(data['Ze_H_{0}'.format(wlStr)])
			ZeH = ZeH1 - ZeH #.sub(ZeH11,fill_value=0)
			xlabel=r'$\Delta$Ze [dB]'
			xlim=[-10,10]
			xticks = [-10,-5,0,5,10]
		else:
			ZeH = mcr.lin2db(data['Ze_H_{0}'.format(wlStr)])
			xlim=[-15,10]
			xlabel='Ze [dB]'
			xticks = [-15,-10,-5,0,5,10]
		ax.plot(ZeH,data.Temp,lw=2,c='C0')
		ax.set_xlabel(xlabel,fontsize=18)
		ax.set_xlim(xlim)
		ax.set_xticks
		
	return ax
		
		
