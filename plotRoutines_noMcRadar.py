'''
This is the beginning of plotting routines for McSnow output
'''
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as colors
from cmcrameri import cm as colmap
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



def plotPropSpecThesis(ax,heightRange,heightRes,mcTable,velBins,prop,cbar=True,gridBaseArea=5.0,diff=False,mcTable1=None):
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
            mcTableTmp1 = mcTableTmp1.to_dataframe()
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
                
                vmin = -2;vmax = 1.4 #1.7
                #cmap = colmap.batlow.copy()#getNewNipySpectral()
                #cmap.set_extremes(under='blue', over='red')
                cmap = colmap.managua.copy()
                cmap.set_extremes(under='bisque', over='Grey') 
                cbar_label = r'D$_{\rm max}}$ [mm]'
                extend = 'both'
        elif prop == 'sPhi':
            weighted = grouped.apply(wavg)
            #print(height)
            #print(weighted.min())
            vmin=-2;vmax=1.5
            
            top = cm.get_cmap('cool',76)
            bottom = cm.get_cmap('autumn',101)
            top =  top(np.linspace(0, 1, 76))
            bottom = bottom(np.linspace(0, 1, 101))
            top[0] = colors.to_rgba('k')
            bottom[-1] = colors.to_rgba('k')
            newcolors = np.vstack((bottom,
                                   top))
            cmap = ListedColormap(newcolors)
            cmap.set_extremes(under='crimson', over='mediumvioletred')
            cbar_label = 'aspect ratio'
            extend = 'both'
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
                #newNipy = getNewNipySpectral()
                #newcolors = newNipy(np.linspace(0, 1, vmax))
                #newcolors = colmap.batlow(np.linspace(0,1,vmax))
                newcolors = colmap.managua(np.linspace(0,1,vmax))
                newcolors[:1, :] = colors.to_rgba('cyan')
                cmap = ListedColormap(newcolors)
                cmap.set_extremes(under='blue', over='Grey') 
                
                #cmap.set_extremes(under='blue', over='red')
                
                cbar_label = 'Number of monomers'
                extend = 'max'
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
                #cmap = colmap.batlow.copy()#getNewNipySpectral()
                #cmap.set_extremes(under='blue', over='red')
                cmap = colmap.managua.copy()
                cmap.set_extremes(under='bisque', over='Grey') 
                cbar_label = r'Number conc. [m$^{-3}$]'
                extend = 'both'
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
                #cmap = colmap.batlow.copy()#getNewNipySpectral()
                #cmap.set_extremes(under='blue', over='red')
                cmap = colmap.managua.copy()
                cmap.set_extremes(under='bisque', over='Grey') 
                cbar_label = r'mass [g]'
                extend = 'both'
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
                #cmap = colmap.batlow.copy()#getNewNipySpectral()
                cmap = colmap.managua.copy()
                cbar_label = r'$\rho$ [kgm$^{-3}$]'
                extend = 'neither'
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
        if (prop == 'dia' or prop == 'number_conc' or prop == 'mTot' or prop == 'sPhi'):
            
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
        cb = plt.colorbar(plot,ax=ax,pad=0.02,aspect=20,extend=extend)#,ticks=v1)
        cb.set_label(cbar_label,fontsize=18)
        cb.ax.yaxis.offsetText.set_fontsize(16)
        #if prop == 'sPhi':
        #	cb.set_ticks(np.arange(11))
        cb.ax.tick_params(labelsize=16)
        
    ax.set_xlabel(r'Velocity [ms$^{-1}$]',fontsize=18)
    ax.set_xticks([-2,-1.5,-1,-0.5,0])
    #ax.pcolormesh(binnedXR.vel,binnedXR['T'],binnedXR,vmin=vmin,vmax=vmax,cmap=cmap)
    #plt.show()
    #ax.set_xlabel('vel [m/s]',fontsize=18)
    
    
    return ax


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
    
    
    

def plotconcHeightProf(particle,mcTable,inputPath,gridBaseArea,maxHeight,nz=200,vmax=10000):
	'''
	plot height profiles of superparticle number concentration. 
	or of real particles
	'''
	if particle == 'Super':
		dz = maxHeight/nz
		Heightrange = np.arange(0,maxHeight,dz)
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
	else:
		dz = 10
		Heightrange = np.arange(0,maxHeight,dz)
		heightCenterBin = Heightrange[0:-1] + dz/2
		#binheight,sheight = pd.cut(mcTable['sHeight'],bins=Heightrange,retbins=True)
		#group according to velocity bins
		mcTableMono = mcTable.where(mcTable['sNmono']==1)
		mcTableAgg = mcTable.where(mcTable['sNmono']>1)
		groupedMono = mcTableMono.groupby_bins('sHeight', Heightrange)
		groupedAgg = mcTableAgg.groupby_bins('sHeight', Heightrange)
		groupedTot = mcTable.groupby_bins('sHeight', Heightrange)
		concMono = groupedMono.sum()['sMult']
		concAgg = groupedAgg.sum()['sMult']
		concTot = groupedTot.sum()['sMult']
		Temp = groupedTot.mean()['Temp']
		fig,ax = plt.subplots(figsize=(5,5))
		vmax = (concTot/(gridBaseArea*dz)).max()
		ax.plot((concTot/(gridBaseArea*dz)).values,Temp.values,label='All')
		ax.plot((concAgg/(gridBaseArea*dz)).values,Temp.values,label='Aggregates')
		ax.plot((concMono/(gridBaseArea*dz)).values,Temp.values,label='Mono')
		ax.set_xlabel('concentration per m$^{3}$')
		#ax.set_ylabel('height [m]')
		ax.set_ylabel('T [°C]')
		ax.set_ylim([0,np.min(Temp)-1])
		ax.legend()
		ax.grid()
		ax.set_xlim([0,vmax])
		plt.tight_layout()
		plt.savefig(inputPath+'concentration_realParticle_Temp.png')
		plt.close()
  #print(Heightrange)
  #quit()  
  
def plotPSD(mcTable,gridBaseArea,heightRes,inputPath,bins,var,gam=None,fac=1,xlim=None,ticks=None,heightEdge0=2900,unit='',sepMono=False,yscale=None,xscale=None):
  #for i, heightEdge0 in enumerate(dicSettings['heightRange'][::-1]):
    heightEdge1 = heightEdge0 + heightRes
    height = heightEdge0+heightRes/2 
    mcTableTmp = mcTable.where((mcTable['sHeight']>heightEdge0) &(mcTable['sHeight']<=heightEdge1),drop=True)
    
    fig,ax = plt.subplots()
    # make it per m3, so I need to multiply dh with boxwidth
    volume = gridBaseArea*heightRes
    if sepMono:
        mcTableMono = mcTableTmp.where(mcTableTmp.sNmono == 1,drop=True)
        mcTableAgg = mcTableTmp.where(mcTableTmp.sNmono > 1,drop=True)
         # this is the volume that the particles are in
        concMono = mcTableMono.sMult.sum()/volume
        concAgg = mcTableAgg.sMult.sum()/volume
        ax.hist(mcTableAgg[var]*fac,bins=bins,weights=mcTableAgg.sMult/volume,label='Agg, tot conc: {0}'.format(int(concAgg)))
        ax.hist(mcTableMono[var]*fac,bins=bins,weights=mcTableMono.sMult/volume,label='Mono, tot conc: {0}'.format(int(concMono)))
        ax.legend()
    else:
        ax.hist(mcTableTmp[var]*fac,bins=bins,weights=mcTableTmp.sMult/volume)
      
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
    plt.savefig(inputPath+'{height}_{var}_hist_weigth_perm3.png'.format(height=height,var=var))
    plt.close()
    
    print(height,' done')

