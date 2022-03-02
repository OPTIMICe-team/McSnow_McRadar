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
def wavgM(group):
    d = group['mTot_g']
    w = group['sMult']
    return (d * w).sum() / w.sum()
def wavgNmono(group):
    d = group['sNmono']
    w = group['sMult']
    return (d * w).sum() / w.sum()


def plotArSpec(dicSettings,mcTable,velBins,inputPath):
# plots the aspect ratios depending on the velocity (so sort of like a Doppler Spectrum)
    for i, heightEdge0 in enumerate(dicSettings['heightRange']):
        heightEdge1 = heightEdge0 + dicSettings['heightRes']
        height = heightEdge0+dicSettings['heightRes']/2
        mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &
                             (mcTable['sHeight']<=heightEdge1)].copy()
        
        # cut into velocity bins
        binVel,sVel = pd.cut(mcTableTmp['vel'],bins=velBins,retbins=True)
        #group according to velocity bins
        grouped = mcTableTmp.groupby(binVel)
        # apply sMult as weight
        weighted = grouped.apply(wavg)
        if i == 0:
            binnedXR = xr.DataArray(weighted.values.reshape(len(weighted),1),
                                    dims=('vel','height'),
                                    coords={'height':height.reshape(1),'vel':velBins[0:-1]})
        else:
            tmpXR = xr.DataArray(weighted.values.reshape(len(weighted),1),
                                    dims=('vel','height'),
                                    coords={'height':height.reshape(1),'vel':velBins[0:-1]})
            binnedXR = xr.concat([binnedXR,tmpXR],dim='height')
    top = cm.get_cmap('autumn', 200)
    bottom = cm.get_cmap('winter', 100)

    newcolors = np.vstack((bottom(np.linspace(0, 1, 100)),
                       top(np.linspace(0, 1, 200))))
    newcmp = ListedColormap(newcolors, name='OrangeBlue')
    binnedXR.plot(x='vel',y='height',vmin=0,vmax=3,cmap=newcmp,cbar_kwargs={'label':'ar'})
    plt.grid(True,ls='-.')
    plt.xlabel('vel [m/s]')
    plt.tight_layout()
    plt.savefig(inputPath+'1d_habit_spec_ar_weighted.png')
    plt.close()
def plotPropSpec(dicSettings,mcTable,velBins,inputPath,prop):
# plots the aspect ratios depending on the velocity (so sort of like a Doppler Spectrum)
    for i, heightEdge0 in enumerate(dicSettings['heightRange']):
        heightEdge1 = heightEdge0 + dicSettings['heightRes']
        height = heightEdge0+dicSettings['heightRes']/2
        mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &
                             (mcTable['sHeight']<=heightEdge1)].copy()
        
        # cut into velocity bins
        #mcTableTmp = mcTableTmp[(mcTableTmp['sNmono'] > 1.0)]
        #print(mcTableTmp['dia_cm'].max())
        #quit()
        #mcTableTmp = mcTableTmp[(mcTableTmp['sPhi']>=0.01)]
        binVel,sVel = pd.cut(mcTableTmp['vel'],bins=velBins,retbins=True)
        #group according to velocity bins
        grouped = mcTableTmp.groupby(binVel)
        # apply sMult as weight
        if prop == 'dia_cm':
            weighted = np.log10(grouped.apply(wavgD))
            vmin = -2;vmax = 1
            cmap = getNewNipySpectral()
            cbar_label = 'log(dia) [cm]'
        elif prop == 'sPhi':
            weighted = grouped.apply(wavg)
            #print(height)
            #print(weighted.min())
            vmin = 0.;vmax = 3.
            #pink = np.array([248/300, 24/300, 148/300, 1])
            #top = cm.get_cmap('autumn', 200)
            #bottom = cm.get_cmap('winter', 99)
            #bottom = bottom(np.linspace(0, 1, 100))
            #bottom[99] = pink           
            #top = top(np.linspace(0, 1, 200))
            #top[0] = pink
            #newcolors = np.vstack((bottom,
            #                       top))
            #cmap = ListedColormap(newcolors, name='OrangeBlue')
            
            #cmap = colors.ListedColormap(['sandybrown','lightsalmon','lightcoral','salmon','darksalmon','coral','darkorange','orangered','darkred','red',
            #                              'magenta',
            #                              'green','lime','aquamarine','turquoise','lightseagreen','mediumturquoise','teal','darkcyan','darkturquoise',
            #                              'cadetblue','deepskyblue','steelblue','dodgerblue','cornflowerblue','royalblue','blue','mediumblue','darkblue','navy','midnightblue'])
            
            #boundaries = np.linspace(0,3,31)
            #norm = colors.BoundaryNorm(boundaries, cmap.N, clip=True)
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
        elif prop == 'mTot_g':
            weighted = np.log10(grouped.apply(wavgM))
            vmin = -5;vmax = 1.5
            cmap = getNewNipySpectral()
            cbar_label = 'log(mTot) [g]'
        elif prop == 'sNmono':
            weighted = grouped.apply(wavgNmono)
            vmin = 1;vmax = 50
            viridis = cm.get_cmap('viridis', vmax)
            newcolors = viridis(np.linspace(0, 1, vmax))
            red = ([1,0/vmax,0/vmax,1])
            pink = np.array([248/vmax, 24/vmax, 148/vmax, 1])
            newcolors[:1, :] = red
            cmap = ListedColormap(newcolors)
            cbar_label = 'Number of monomers'
        elif prop == 'sNmono_min':
            weighted = grouped['sNmono'].min()
            vmin = 1;vmax = 50
            viridis = cm.get_cmap('viridis', vmax)
            newcolors = viridis(np.linspace(0, 1, vmax))
            red = ([1,0/vmax,0/vmax,1])
            pink = np.array([248/vmax, 24/vmax, 148/vmax, 1])
            newcolors[:1, :] = red
            cmap = ListedColormap(newcolors)
            cbar_label = 'minimum number of monomers'
        elif prop == 'sNmono_max':
            weighted = grouped['sNmono'].max()
            vmin = 1;vmax = 50
            viridis = cm.get_cmap('viridis', vmax)
            newcolors = viridis(np.linspace(0, 1, vmax))
            red = ([1,0/vmax,0/vmax,1])
            pink = np.array([248/vmax, 24/vmax, 148/vmax, 1])
            newcolors[:1, :] = red
            cmap = ListedColormap(newcolors)
            cbar_label = 'maximum number of monomers'
        elif prop == 'number_conc':
            weighted = np.log10(grouped['sMult'].sum())
            vmin = 0;vmax = 6
            cmap = getNewNipySpectral()
            cbar_label = 'Number concentration [log]'
        else:
            print(prop+' not yet defined')
        if i == 0:
            binnedXR = xr.DataArray(weighted.values.reshape(len(weighted),1),
                                    dims=('vel','height'),
                                    coords={'height':height.reshape(1),'vel':velBins[0:-1]})
        else:
            tmpXR = xr.DataArray(weighted.values.reshape(len(weighted),1),
                                    dims=('vel','height'),
                                    coords={'height':height.reshape(1),'vel':velBins[0:-1]})
            binnedXR = xr.concat([binnedXR,tmpXR],dim='height')
    binnedXR.plot(x='vel',y='height',vmin=vmin,vmax=vmax,cmap=cmap,cbar_kwargs={'label':cbar_label})
    plt.grid(True,ls='-.')
    plt.xlabel('vel [m/s]')
    plt.ylim([0,dicSettings['maxHeight']])
    plt.tight_layout()
    plt.savefig(inputPath+'1d_habit_spec_'+prop+'_weighted_alltimes.png')
    plt.close()
def plotMoments(dicSettings,output,inputPath,convoluted=False):
    for wl in dicSettings['wl']:
        wlStr = '{:.2e}'.format(wl)
        if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
          if convoluted == True:
              saveName = '1d_habit_moments_{0}_convoluted.png'.format(wlStr)
              fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
              # plot ZDR
              output['MDV_H_{0}'.format(wlStr)].plot(ax=axes[0],y='range', lw=2)
              axes[0].set_title('rad: {0} elv: {1}, MDV'.format(wlStr, dicSettings['elv']))
              #axes[0].set_ylim(0, 5000)
              axes[0].grid(True,ls='-.')
              axes[0].set_xlabel('MDV [m/s]')
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
              saveName = '1d_habit_moments_{wl}_{mode}.png'.format(wl=wlStr,mode=dicSettings['scatSet']['mode'])
           
              fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
              output['MDV_H_{0}'.format(wlStr)].plot(ax=axes[0],y='range', lw=2)
              axes[0].set_title('rad: {0} elv: {1}, MDV'.format(wlStr, dicSettings['elv']))
              #axes[0].set_ylim(0, 5000)
              axes[0].grid(True,ls='-.')
              axes[0].set_xlabel('MDV [m/s]')
              axes[0].set_ylim([0,dicSettings['maxHeight']])
        
              axes[1].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)]),output['range'],linewidth=2) # TODO: change back to ZeH
              axes[1].set_xlabel('Z_H [dB]')
              axes[1].set_title('Ze_H')
              #plt.ylim(0, 5000)
              #plt.xlim(-3, 0)
              axes[1].grid(True,ls='-.')
              axes[1].set_ylim([0,dicSettings['maxHeight']])
              plt.tight_layout()
              
              
        else:
            if convoluted == True:
              saveName = '1d_habit_moments_{0}_convoluted_nace_40.png'.format(wlStr)
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
              saveName = '1d_habit_moments_{0}_ar_0.01.png'.format(wlStr)
        
        #plot KDP    
              fig,axes = plt.subplots(ncols=3,figsize=(15,5),sharey=True)
              output['kdpInt_{0}'.format(wlStr)].plot(ax=axes[0],y='range', lw=2)
              axes[0].set_title('rad: {0} elv: {1}, KDP'.format(wlStr, dicSettings['elv']))
              #axes[0].set_ylim(0, 5000)
              axes[0].grid(True,ls='-.')
              axes[0].set_xlabel('KDP [°/km]')
              axes[0].set_ylim([0,dicSettings['maxHeight']])
        
        # plot ZDR
              axes[1].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)])-mcr.lin2db(output['Ze_V_{0}'.format(wlStr)]),output['range'],linewidth=2)# TODO: change bach to ZeH-ZeV
              axes[1].set_xlabel('ZDR [dB]')
              axes[1].set_title('ZDR')
              #plt.xlim(-3, 0)
              axes[1].grid(True,ls='-.')
              axes[1].set_ylim([0,dicSettings['maxHeight']])
      
              axes[2].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)]),output['range'],linewidth=2) # TODO: change back to ZeH
              axes[2].set_xlabel('Z_H [dB]')
              axes[2].set_title('Ze_H')
              #plt.ylim(0, 5000)
              #plt.xlim(-3, 0)
              axes[2].grid(True,ls='-.')
              axes[2].set_ylim([0,dicSettings['maxHeight']])
              plt.tight_layout()
            
        plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
        plt.close()
def plotDWR(dicSettings,wlStr1,wlStr2,output,inputPath,convoluted=False):
    #for wl in dicSettings['wl']:
    #wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
    #wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
    if convoluted == True:
      saveName = '1d_habit_DWR_{freq1}_{freq2_convoluted_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
          
      fig,axes = plt.subplots(figsize=(5,5))
      DWR = mcr.lin2db(output['Ze_H_{0}'.format(wlStr1)]) - mcr.lin2db(output['Ze_H_{0}'.format(wlStr2)])
      DWR.plot(ax=axes,y='range', lw=2)
      axes.set_title('DWR_{freq1}_{freq2}'.format(freq1=freq1,freq2=freq2))
      #axes[0].set_ylim(0, 5000)
      axes.grid(True,ls='-.')
      axes.set_xlabel('DWR {freq1}_{freq2} [dB]'.format(freq1=freq1,freq2=freq2))
      axes.set_ylim([0,dicSettings['maxHeight']])
          
    else: 
      freq1 = (constants.c / float(wlStr1))  *1e3 / 1e9
      freq2 = (constants.c / float(wlStr2))  *1e3 / 1e9
      freq1 = '{:.1f}'.format(freq1)
      freq2 = '{:.1f}'.format(freq2)
      
      saveName = '1d_habit_DWR_{freq1}_{freq2}_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
           
      fig,axes = plt.subplots(figsize=(5,5))
      Ze1 = mcr.lin2db(output['spec_H_{0}'.format(wlStr1)].sum(dim='vel'))
      Ze2 = mcr.lin2db(output['spec_H_{0}'.format(wlStr2)].sum(dim='vel'))
      DWR = Ze1 - Ze2
      DWR.plot(ax=axes,y='range', lw=2)
      axes.set_title('DWR_{freq1}_{freq2}'.format(freq1=freq1,freq2=freq2))
      #axes[0].set_ylim(0, 5000)
      axes.grid(True,ls='-.')
      axes.set_xlabel('DWR {freq1}_{freq2} [dB]'.format(freq1=freq1,freq2=freq2))
      axes.set_ylim([0,dicSettings['maxHeight']])
                              
    plt.tight_layout()        
    plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
    plt.close()
def plotDWRspectra(dicSettings,wlStr1,wlStr2,output,inputPath,convoluted=False):
    #for wl in dicSettings['wl']:
    #wlStr1 = '{:.2e}'.format(dicSettings['wl'][0])
    #wlStr2 = '{:.2e}'.format(dicSettings['wl'][1])
    if convoluted == True:
      saveName = '1d_habit_sDWR_{freq1}_{freq2_convoluted_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
          
      fig,axes = plt.subplots(figsize=(5,5))
      specH1 = mcr.lin2db(output['spec_H_{0}'.format(wlStr1)])
      specH1 = specH1.where(specH1 > -40)
      specH2 = mcr.lin2db(output['spec_H_{0}'.format(wlStr2)])
      specH2 = specH2.where(specH2 > -40)
      
      DWR =specH1 - specH2
      DWR.plot(ax=axes,y='range', vmin=0, vmax=15, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sDWR [dB]'})
      axes.set_title('sDWR_{freq1}_{freq2}'.format(freq1=freq1,freq2=freq2))
      #axes[0].set_ylim(0, 5000)
      axes.grid(True,ls='-.')
      axes.set_xlabel('DWR {freq1}_{freq2} [dB]'.format(freq1=freq1,freq2=freq2))
      axes.set_ylim([0,dicSettings['maxHeight']])
      axes.set_xlim(-2, 0)
    else: 
      freq1 = (constants.c / float(wlStr1))  *1e3 / 1e9
      freq2 = (constants.c / float(wlStr2))  *1e3 / 1e9
      freq1 = '{:.1f}'.format(freq1)
      freq2 = '{:.1f}'.format(freq2)
      
      saveName = '1d_habit_sDWR_{freq1}_{freq2}_{mode}.png'.format(freq1=freq1,freq2=freq2,mode=dicSettings['scatSet']['mode'])
           
      fig,axes = plt.subplots(figsize=(5,5))
      specH1 = mcr.lin2db(output['spec_H_{0}'.format(wlStr1)])
      specH1 = specH1.where(specH1 > -40)
      specH2 = mcr.lin2db(output['spec_H_{0}'.format(wlStr2)])
      specH2 = specH2.where(specH2 > -40)
      
      DWR =specH1 - specH2
      DWR.plot(ax=axes,y='range', vmin=0, vmax=15, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sDWR [dB]'})
      axes.set_title('sDWR_{freq1}_{freq2}'.format(freq1=freq1,freq2=freq2))
      #axes[0].set_ylim(0, 5000)
      axes.grid(True,ls='-.')
      axes.set_xlabel('DWR {freq1}_{freq2} [dB]'.format(freq1=freq1,freq2=freq2))
      axes.set_ylim([0,dicSettings['maxHeight']])
      axes.set_xlim(-2, 0)                        
    plt.tight_layout()        
    plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
    plt.close()

def plotSpectra(dicSettings,output,inputPath,convoluted=False):
    for wl in dicSettings['wl']:
        wlStr = '{:.2e}'.format(wl)
        print(wlStr)
        if (dicSettings['scatSet']['mode'] == 'SSRGA') or (dicSettings['scatSet']['mode'] == 'Rayleigh') or (dicSettings['scatSet']['mode'] == 'SSRGA-Rayleigh'):
            if convoluted == True:
                saveName = '1d_habit_spectra_{wl}_convoluted_{mode}.png'.format(wl=wlStr,mode=dicSettings['scatSet']['mode'])
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                specH = specH.where(specH > -40)
            else: 
                saveName = '1d_habit_spectra_{wl}_{mode}.png'.format(wl=wlStr,mode=dicSettings['scatSet']['mode'])
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
            fig,ax = plt.subplots(figsize=(5,4))
            specH.plot(ax=ax,vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
            ax.set_xlabel('Z_H [dB]')
            ax.set_title('Ze_H')
            ax.set_title('Ze_H_spec rad: {0} elv: {1}'.format(wlStr, dicSettings['elv']))# TODO: change back to ZeH or to ZeV if you use first LUT of dendrites
            #axes[0].set_ylim([0,dicSettings['maxHeight']])
            ax.set_xlim(-2, 0)
            ax.grid(True,ls='-.')      
            plt.tight_layout()
        else:
            if convoluted == True:
                saveName = '1d_habit_spectra_{0}_convoluted.png'.format(wlStr)
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                specV = mcr.lin2db(output['spec_V_{0}'.format(wlStr)])
                specH = specH.where(specH > -40)
                specV =  specV.where(specV > -40)
                #dataSmooth = specTable.rolling(vel=10,min_periods=1,center=True).mean()            
                ZDR = specH.rolling(vel=10,min_periods=1,center=True).mean() - specV.rolling(vel=10,min_periods=1,center=True).mean()
            else: 
                saveName = '1d_habit_spectra_{0}.png'.format(wlStr)
                specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
                specV = mcr.lin2db(output['spec_V_{0}'.format(wlStr)])
                ZDR = specH - specV
            fig,axes = plt.subplots(ncols=3,figsize=(12,5),sharey=True)
            specH.plot(ax=axes[0],vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
            axes[0].set_xlabel('Z_H [dB]')
            axes[0].set_title('Ze_H')
            axes[0].set_title('Ze_H_spec rad: {0} elv: {1}'.format(wlStr, dicSettings['elv']))# TODO: change back to ZeH or to ZeV if you use first LUT of dendrites
            #axes[0].set_ylim([0,dicSettings['maxHeight']])
            axes[0].set_xlim(-2, 0)
            axes[0].grid(True,ls='-.')
        
            specV.plot(ax=axes[1],vmin=-30, vmax=5, cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZe [dB]'})
            axes[1].set_title('Ze_V_spec rad: {0} elv: {1}'.format(wlStr, dicSettings['elv']))# TODO: change back to ZeV or to ZeV if you use first LUT of dendrites
            axes[1].set_ylim([0,dicSettings['maxHeight']])
            axes[1].set_xlim(-2, 0)
            axes[1].grid(True,ls='-.')
            axes[1].set_ylabel('')
            ZDR.plot(ax=axes[2],vmin=-0.5, vmax=3,cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZDR [dB]'}) # TODO: change back to ZeH-ZeV or to ZeV if you use first LUT of dendrites
            axes[2].set_title('ZDR rad: {0} elv: {1}'.format(wlStr, dicSettings['elv']))
            axes[2].set_xlim(-2, 0)
            axes[2].set_ylim([0,dicSettings['maxHeight']])
            axes[2].grid(True,ls='-.')
            axes[2].set_ylabel('')
            plt.tight_layout()
        
        plt.savefig(inputPath+saveName, format='png', dpi=200, bbox_inches='tight')
        plt.close()

def plotProposal(dicSettings,output,inputPath):
    for wl in dicSettings['wl']:
        
        wlStr = '{:.2e}'.format(wl)
        specH = mcr.lin2db(output['spec_H_{0}'.format(wlStr)])
        specV = mcr.lin2db(output['spec_V_{0}'.format(wlStr)])
        specH = specH.where(specH > -40)
        specV =  specV.where(specV > -40)
        #dataSmooth = specTable.rolling(vel=10,min_periods=1,center=True).mean()            
        sZDR = specH.rolling(vel=10,min_periods=1,center=True).mean() - specV.rolling(vel=10,min_periods=1,center=True).mean()
        fig,axes = plt.subplots(ncols=4,figsize=(15,5),sharey=True)
        axes[0].plot(output['conc_agg']/250,output['temperature [°C]'], lw=2,label='Aggregate')
        axes[0].plot(output['conc_mono']/250,output['temperature [°C]'], lw=2,label='Crystal')
        axes[0].grid(True,ls='-.')
        axes[0].set_ylim([-2.5,-20])
        #axes[0].set_xticks([-5,-10,-15,-20])
        axes[0].set_xlabel(r'conc [m$^{-3}$]',fontsize=18)
        axes[0].set_ylabel(r'Temp [°C]',fontsize=18)
        axes[0].tick_params(labelsize=16)
        axes[0].legend(fontsize=16)
        axes[0].text(1,-18.6,'a)',fontsize=22)
        plt.tight_layout()        
        axes[1].plot(output['kdpInt_{0}'.format(wlStr)],output['temperature [°C]'], lw=2)
        axes[1].grid(True,ls='-.')
        #axes[1].set_ylim([-20,0])
        axes[1].set_xlabel(r'KDP [°km$^{-1}$]',fontsize=18)
        axes[1].set_ylabel('')
        axes[1].tick_params(labelsize=16)
        axes[1].text(0.002,-18.6,'b)',fontsize=22)
        plt.tight_layout()
        axes[2].plot(mcr.lin2db(output['Ze_H_{0}'.format(wlStr)])-mcr.lin2db(output['Ze_V_{0}'.format(wlStr)]),output['temperature [°C]'], lw=2)
        axes[2].set_xlabel('ZDR [dB]',fontsize=18)
        #axes[2].set_ylim([-20,0])
        axes[2].grid(True,ls='-.')
        axes[2].set_ylabel('')
        axes[2].tick_params(labelsize=16)
        axes[2].text(0.48,-18.6,'c)',fontsize=22)
        plt.tight_layout()
        #sZDR.plot(ax=axes[3],x='vel',vmin=-0.5, vmax=3,cmap=getNewNipySpectral(),cbar_kwargs={'label':'sZDR [dB]'}) # TODO: change back to ZeH-ZeV or to ZeV if you use first LUT of dendrites
        p=axes[3].pcolormesh(output['vel']*np.sin(30/180.*np.pi),output['temperature [°C]'],sZDR,vmin=-0.5, vmax=3,cmap=getNewNipySpectral())
        cb=fig.colorbar(p,ax=axes[3])
        cb.set_label('sZDR [dB]',fontsize=16)#,fontsize=22)
        cb.ax.tick_params(labelsize=16)
        axes[3].set_xlim(-1, 0)
        #axes[3].set_ylim([-20,0])
        axes[3].set_ylabel('')
        axes[3].set_xlabel(r'Doppler velocity [ms$^{-1}$]',fontsize=18)
        axes[3].grid(True,ls='-.')
        axes[3].tick_params(labelsize=16)
        axes[3].text(-0.96,-18.6,'d)',fontsize=22)
        plt.tight_layout()
        
        plt.savefig(inputPath+'KDP_ZDR_sZDR_proposal_sin30.png', format='png', dpi=200, bbox_inches='tight')
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
    qv = atmoFile[:,10]
    qc = atmoFile[:,11]
    fig,axes=plt.subplots(ncols=4,figsize=(15,5),sharey=True)
    axes[0].plot(Temp,height)
    axes[0].grid(True,ls='-.')
    axes[0].set_xlabel('Temp [°C]')
    axes[0].set_ylabel('height [m]')
    
    axes[1].plot(rh,height)
    axes[1].grid(True,ls='-.')
    axes[1].set_xlabel('rh [%]')
    
    axes[2].plot(qv,height)
    axes[2].grid(True,ls='-.')
    axes[2].set_xlabel('qv [%]')
    
    axes[3].plot(ssat,height)
    axes[3].grid(True,ls='-.')
    axes[3].set_xlabel('ssat [%]')
    
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
        

def plotSpectraObs(dataPol,dataLV0,outPath,ylim=False):
    for t in dataLV0.time:
        print(t)
        ti = pd.to_datetime(str(t.values)).strftime('%Y%m%d_%H%M%S')
        spec = 10*np.log10(dataLV0.KaSpecH.sel(time=t)) 
        #specV = 10*np.log10(dataLV0.KaSpecV.sel(time=t))
        #LDR = specV - specH
        
        dataPolSel = dataPol.sel(time=t)
        fig,axes = plt.subplots(ncols=2,figsize=(10,5),sharey=True)
        spec.plot(ax = axes[0],y='range',vmin=-30,vmax=5,cmap=getNewNipySpectral())
        axes[0].set_ylabel('height [m]')
        axes[0].set_xlabel('Doppler vel [m/s]')
        axes[0].grid(True,ls='-.')
        axes[0].set_xlim([-2.5,1.5])
        if ylim:
            axes[0].set_ylim(ylim)
        plot=axes[1].pcolormesh(dataPolSel.Vel2ZeroH.fillna(0), # this is needed because we measure at an elevation angle, therefore if we would use the Doppler Velocity as x-axis, the spectra looks like a snake, so everything is moved to zero. Does not always work because of uncertainty of getting the spectral edge velocity correctly
                       dataPolSel.height,
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
        
        
def plotPSD(dicSettings,mcTable,dBins,inputPath):
# plots the PSD binned in D bins, dBins in cm
    for i, heightEdge0 in enumerate(dicSettings['heightRange'][::-1]):
        heightEdge1 = heightEdge0 + dicSettings['heightRes']
        height = heightEdge0+dicSettings['heightRes']/2
        mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &
                             (mcTable['sHeight']<=heightEdge1)].copy()
        
        
        mcTableTmpMono = mcTableTmp[mcTableTmp['sNmono'] == 1].copy()
        mcTableTmpAgg = mcTableTmp[mcTableTmp['sNmono'] > 1].copy()
        plt.xscale('log')
        #plt.yscale('log')
        plt.hist(mcTableTmpMono.dia_cm,bins=dBins,label='Crystal')
        plt.hist(mcTableTmpAgg.dia_cm,bins=dBins,label='Agg')
        plt.ylabel('#(SP)')
        plt.xlabel(r'D [cm]')
        plt.title('Height: ' + str(height))
        plt.grid(which='both')
        plt.legend()
        outfileName = 'histD_height_' + str(height)
        plt.savefig(inputPath + '/' + outfileName + '.png')
        plt.show()
        quit()
        
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
 
def plot_var_particle_ID(mcTable,inputPath,prop,var2col,atmoFile,zoom=False):
  # this plots the particles according to their ID (which in this first setup is sNmono) and height. This is to test the evolution of different parameters with respect to height. Input: mcTable, path to save the plot and the variable to plot (prop)
  # zoom needs to  be array of ylim
  # first lets make temp like height
  height = atmoFile[:,0]
  Temp = atmoFile[:,2] -273.15
  idx15 = (np.abs(Temp - (-15))).argmin()
  idx20 = (np.abs(Temp - (-16))).argmin()
  idx10 = (np.abs(Temp - (-10))).argmin()
  
  n = max(mcTable[var2col])
  mcTable = mcTable.set_index(var2col,drop=False)
  
  norm = mpl.colors.Normalize(vmin=0, vmax=n)
  c_m = mpl.cm.gist_rainbow
  # create a ScalarMappable and initialize a data structure
  s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
  s_m.set_array([])
  step = n
  step = step-1
  tickmarks = np.arange(0,n)
  fig, ax1 = plt.subplots(figsize=(9,5))

  ax2 = ax1.twinx()
  for sID in mcTable[var2col]:
    mcTableTmp = mcTable.loc[sID]
    if 'dia' in prop:
      ax1.semilogy(mcTableTmp['sHeight'],mcTableTmp[prop],color = s_m.to_rgba(sID-1))
      
    else:
      ax1.plot(mcTableTmp['sHeight'],mcTableTmp[prop],color = s_m.to_rgba(sID-1))
  
  ax1.axvline(x=height[idx20],c='k',linestyle=':',label='-16°C')    
  ax1.axvline(x=height[idx15],c='k',linestyle='-',label='-15°C')
  ax1.axvline(x=height[idx10],c='k',linestyle='--',label='-10°C')
  #mcTable.plot.scatter(x='sHeight',y=prop,c='dia_cm',cmap='gist_rainbow') 
  #print(mcTableTmp)
  #cax = fig.add_axes([0.95, 0.1, 0.02, 0.75])
  if var2col == 'time':
    var2col='time [min]'
  if var2col == 'sMult':
    var2col = 'particle_ID'
  cbar = plt.colorbar(s_m, pad=0.15)
  cbar.ax.tick_params(labelsize=12)
  cbar.set_label(var2col,fontsize=13)
  ax2.plot(height,Temp,c='b')
  ax1.set_xlabel('sHeight [m]',fontsize=13)
  ax2.set_ylim([0,-20])
  ax1.set_ylabel(prop,fontsize=13)
  ax2.set_ylabel('Temp [°C]',color='b',fontsize=13)
  ax2.set_xlim([5000,2000])
  ax1.set_xlim([5000,2000])
  if zoom:
    plt.ylim([zoom[0],zoom[1]])
    prop=prop+'_zoom'
  ax1.grid()
  #ax2.grid(ls='-.')
  ax1.legend()
  plt.tight_layout()
  plt.savefig(inputPath+prop+'_'+var2col+'_Temp.png')
  plt.close()
  

def plot_var_particle_ID_temp(dicSettings,atmoPD,mcTable,inputPath,prop,var2col,zoom=False):
  # this plots the particles according to their ID (which in this first setup is sNmono) and height. This is to test the evolution of different parameters with respect to height. Input: mcTable, path to save the plot and the variable to plot (prop)
  # zoom needs to  be array of ylim
  # first lets make temp like height
  Temp = np.empty(len(dicSettings['heightRange']))

  propvec = Temp.copy()
  var2colvec = Temp.copy()
  #dia_mum = Temp.copy()
  #sPhi = Temp.copy()
  #dia = Temp.copy()
  for i, heightEdge0 in enumerate(dicSettings['heightRange']):
    heightEdge1 = heightEdge0 + dicSettings['heightRes']
    height = heightEdge0+dicSettings['heightRes']/2
    atmoTmp = atmoPD[(atmoPD['height']>heightEdge0) & (atmoPD['height']<=heightEdge1)].copy()
    Temp[i] = atmoTmp.temp.mean()
    mcTableTmp = mcTable[(mcTable['sHeight']>heightEdge0) &(mcTable['sHeight']<=heightEdge1)].copy()
    propvec[i] = mcTableTmp[prop].mean()
    var2colvec[i] = mcTableTmp[var2col].mean()
   # dia[i] = mcTableTmp['dia'].mean()
   # dia_mum[i] = mcTableTmp['dia_mum'].mean()
   # mTot[i] = mcTableTmp['mTot'].mean()
  
    
          
