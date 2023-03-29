import mcradar as mcr
from scipy import constants
import pandas as pd
import xarray as xr
import numpy as np
import plotRoutines as plot
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
def str2bool(v):
  return v.lower() in ("yes", "True", "t", "1","true")
  
freqEnv = os.environ['freq'].split('_')
print(freqEnv)
elvEnv = os.environ['elv'].split('_')
#elv = float(os.environ['elv'])
elv = np.array([float(e) for e in elvEnv])
particle_name = os.environ['particle']
freq = np.array([float(f)*1e9 for f in freqEnv])
experimentID = os.environ['experiment']
print(experimentID)
#quit()
inputPath = os.environ['MCexp']+'/'+experimentID+'/'
scatMode = os.environ['scatMode']
print(inputPath)
#convolute=os.environ['convolute']
print('loading the settings')
splitPath = inputPath.split('domtop')[1]
domTop = splitPath[0:4]
lutPath = os.environ['LUT_dir']

# decide what you want to plot
plot_initemp = False
plot_inidia = False
plot_thesis_particle_evolution = False
mult_conc = False
single_particle=str2bool(os.environ['singleParticle'])
convoluted=str2bool(os.environ['convolute'])

#-- load the settings of McSnow domain, as well as elevation you want to plot:
#In order to avoid volume sampling problems, you have to insert the gridBaseArea as it was defined in the McSnow simulation
if ('trajectories' not in experimentID) and ('trajectories' not in inputPath):
	heightRes = 36
else:
	heightRes = 2
dicSettings = mcr.loadSettings(dataPath=inputPath+'mass2fr.nc',
                               elv=elv, freq=freq,gridBaseArea=5.0,maxHeight=int(domTop),
                               ndgsVal=50,heightRes=heightRes,scatSet={'mode':scatMode, 'lutPath':lutPath,'particle_name':particle_name,'safeTmatrix':False})

print('loading the McSnow output')
# now generate a table from the McSnow output. You can specify xi0, if it is not stored in the table (like in my cases)
mcTable = mcr.getMcSnowTable(dicSettings['dataPath'])


try:
  minmax = os.environ['minmax']
  vmin=int(minmax.split('_')[0]); vmax=int(minmax.split('_')[1])
except:
  minmax=False
#minmax='0_350'; vmin=0; vmax=350
plotTemp=True
#- now reading in McSnow output. If we did not run McSnow, need to comment that out TODO: make that automatic
McRadar_Outname = os.environ['McRadarfileName']
print(McRadar_Outname)
# plot property spectra from McSnow output. THis is done for all simulations that are NOT a trajectory
if ('trajectories' not in experimentID) and ('trajectories' not in inputPath):
	#- read in atmo file to get Temperature information if you want to plot it with that
	selTime = mcTable['time'].max()
	times = mcTable['time']
	#mcTable = mcTable.sort_values('sHeight')
	mcTable = mcTable.where(times==selTime,drop=True)
	mcTable = mcTable.sortby('sHeight').set_index(index='sHeight').rename({'index':'sHeight'})
	
	atmoFile = np.loadtxt(inputPath+'atmo.dat')
	plot.plotAtmo(atmoFile,inputPath)
	height = atmoFile[:,0]
	Temp = atmoFile[:,2] -273.15
	atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
	atmoPD.index.name='sHeight'
	atmoXR = atmoPD.to_xarray()
	atmoReindex = atmoXR.reindex_like(mcTable,method='nearest')
	mcTableTmp = xr.merge([atmoReindex,mcTable])
	mcTableTmp = mcTableTmp.to_dataframe()

	mcTableMono = mcTableTmp[mcTableTmp.sNmono==1]
	#quit()
	# now plotting stuff directly from McSnow output but in the shape of a velocity spectrum:
	#print('plotting aspect ratios')
	velBins = np.linspace(-3,0,100)
	dBins = 10**(np.linspace(-3,0,100))
	#print(dBins)
	'''
	print('plotting particle properties')
	fig,ax=plt.subplots(ncols=3,nrows=2,figsize=(15,10))
	varVec = ['dia','mTot','sRho_tot']
	for i,var in enumerate(varVec):	
		print(var)
		ax[0,i]=plot.plotPropSpecThesis(ax[0,i],dicSettings['heightRange'],dicSettings['heightRes'],mcTableTmp,velBins,var)
		ax[0,i].set_ylim([mcTableTmp['Temp'].max()+0.1,mcTableTmp['Temp'].min()-0.1])
		ax[0,i].set_xlim([-2,0])
		ax[0,i].tick_params(axis='both',labelsize=16)
		ax[0,i].text(ax[0,i].get_xlim()[0]+0.04*(ax[0,i].get_xlim()[1]-ax[0,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
		ax[0,i].grid(ls='-.')
		if mcTableTmp['Temp'].min() < -20:
			ax[0,i].axhline(y=-20,ls='--',color='r',linewidth=2)
			ax[0,i].axhline(y=-10,ls='--',color='r',linewidth=2)
		if i == 0:
			ax[0,i].set_ylabel('T [°C]',fontsize=18)
		else:
			ax[0,i].set_ylabel('')
			ax[0,i].set_yticklabels('')
	#plot.plotPropSpecThesis(dicSettings,mcTable,velBins,inputPath,'dia_cm')
	varVec = ['sNmono','sPhi','number_conc']
	for i,var in enumerate(varVec):	
		print(var)
		ax[1,i]=plot.plotPropSpecThesis(ax[1,i],dicSettings['heightRange'],dicSettings['heightRes'],mcTableTmp,velBins,var)
		ax[1,i].set_xlim([-2,0])
		ax[1,i].set_ylim([mcTableTmp['Temp'].max()+.1,mcTableTmp['Temp'].min()-.1])
		ax[1,i].tick_params(axis='both',labelsize=16)
		ax[1,i].text(ax[1,i].get_xlim()[0]+0.04*(ax[1,i].get_xlim()[1]-ax[1,i].get_xlim()[0]),-27,'('+string.ascii_lowercase[i+3]+')',fontsize=18)
		ax[1,i].grid(ls='-.')
		if mcTableTmp['Temp'].min() < -20:
			ax[1,i].axhline(y=-20,ls='--',color='r',linewidth=2)
			ax[1,i].axhline(y=-10,ls='--',color='r',linewidth=2)
		if i == 0:
			ax[1,i].set_ylabel('T [°C]',fontsize=18)
		else:
			ax[1,i].set_ylabel('')
			ax[1,i].set_yticklabels('')
	plt.tight_layout()
	plt.savefig(inputPath+'properties.png')
	plt.close()
	#quit()
	'''
	#-- plot number of superparticles per grid cell
	
	nz = float(inputPath.split('nz')[1].split('_')[0])
	#heightProf = pd.read_csv(inputPath+'hei2massdens.dat',header=0)#,skiprows=1, # header=0
	#print(heightProf)
	print('plotting number of sp')
	plot.plotHeightProf(nz,mcTable,inputPath,dicSettings)
	print('plotting PSD')
	mBins = 10**(np.linspace(-12,-9,100))
	plot.plotPSD(mcTable,dicSettings,inputPath,mBins,'mTot',heightEdge0=1900,unit='[kg]',sepMono=True,yscale='log',xscale='log')
	dBins = 10**(np.linspace(-6,-2.5,100))
	plot.plotPSD(mcTable,dicSettings,inputPath,dBins,'dia',heightEdge0=1900,unit='[m]',sepMono=True,yscale='log',xscale='log')

	if os.path.exists(inputPath+McRadar_Outname):
		output = xr.open_dataset(inputPath+McRadar_Outname)
		#outPut1 = xr.open_dataset(inputPath+os.environ['freq']+'GHz_output_{mode}_180_350_singleParticle.nc'.format(mode=dicSettings['scatSet']['mode']))
		#outPut = xr.merge([outPut,outPut1])
		print('now plotting McRadar')
		#else:			
		if plotTemp == True:
			atmoFile = np.loadtxt(inputPath+'atmo.dat')
			height = atmoFile[:,0]
			Temp = atmoFile[:,2] -273.15
			atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
			atmoPD.index.name='range'
			atmoXR = atmoPD.to_xarray()
			atmoReindex = atmoXR.reindex_like(output,method='nearest')
			output = xr.merge([atmoReindex,output])
		print('plotting spectra')
		if minmax:
			plot.plotSpectra(dicSettings,output,inputPath,minmax=minmax,plotTemp=plotTemp)#,mult_conc=mult_conc)#,convoluted=True)
		else:
			plot.plotSpectra(dicSettings,output,inputPath,plotTemp=plotTemp)#,mult_conc=mult_conc)#,convoluted=True)
		print('plotting moments')

		if minmax:
			plot.plotMoments(dicSettings,output,inputPath,minmax=minmax,plotTemp=plotTemp)
		else:
			plot.plotMoments(dicSettings,output,inputPath,plotTemp=plotTemp)
		print(freq)
		if len(freq)==2:
			print('plotting DWR')
			plot.plotDWR(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
			plot.plotDWRspectra(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
		elif len(freq)==3:
			print('plotting DWR')
			plot.plotDWR(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
			plot.plotDWR(dicSettings,dicSettings['wl'][0],dicSettings['wl'][2],output,inputPath,plotTemp=plotTemp)
			plot.plotDWRspectra(dicSettings,dicSettings['wl'][0],dicSettings['wl'][1],output,inputPath,plotTemp=plotTemp)
			plot.plotDWRspectra(dicSettings,dicSettings['wl'][0],dicSettings['wl'][2],output,inputPath,plotTemp=plotTemp)


#- plot ar test setup (with bnd_type==3)
if ('trajectories' in experimentID) or ('trajectories' in inputPath):
	atmoFile = np.loadtxt(inputPath+'atmo.dat')
	plot.plotAtmo(atmoFile,inputPath)
	height = atmoFile[:,0]
	Temp = atmoFile[:,2] -273.15
	atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
	atmoPD.index.name='sHeight'
	mcTableNew=mcTable.set_index('sHeight',drop=False)
	mcTableNew = mcTableNew.rename(columns={'sHeight':'height'})
	mcTableXR = mcTableNew.to_xarray()
	atmoXR = atmoPD.to_xarray()
	atmoReindex = atmoXR.reindex_like(mcTableXR,method='nearest')
	mcTableTempXR = xr.merge([atmoReindex,mcTableXR])
	#print(mcTableTemp)
	mcTableTemp = mcTableTempXR.to_dataframe()
	print(mcTableTemp)

	###############################################################
	# plot particle properties dependend on ini temp or ini diam (Figure 5.1 in my thesis)
	############################################################## 
	if plot_initemp:
		fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(10,9))
		if minmax:
			mcTableTmp = mcTableTemp[(mcTableTemp['sMult']>vmin) & (mcTableTemp['sMult']<=vmax)]
		else:
			mcTableTmp = mcTableTemp
		
		a=plot.plotInitempVar(mcTableTmp,ax[0,0],inputPath,'dia','Temp','[m]',r'D$_{\rm max}$',ylog=True)#,zoom=[-20,-10])#,zoom=[0,0.2])
		#a.set_xlim([-9,-31.5])
		#a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),2.5*10**(-3),'(a)',fontsize=18)
		print('dia')
		
		a=plot.plotInitempVar(mcTableTmp,ax[0,1],inputPath,'mTot','Temp','[kg]','mass',ylog=True)#,zoom=[-20,-10])#,zoom=[0,0.2])
		#a.set_xlim([-9,-31.5])
		#a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),7*10**(-8),'(b)',fontsize=18)
		#a.set_yticks([10**-11,10**-10,10**-9,10**-8,10**-7,10**-6,10**-5])
		print('mTot')
		
		a=plot.plotInitempVar(mcTableTmp,ax[1,0],inputPath,'sPhi','Temp','',r'$\phi$',ylog=True)#,zoom=[-20,-10])#,zoom=[0,0.2])
		#a.set_xlim([-9,-31.5])
		#a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),1.2*10**0,'(c)',fontsize=18)
		
		print('sPhi')
		
		a=plot.plotInitempVar(mcTableTmp,ax[1,1],inputPath,'sRho_tot','Temp',r'[kgm$^{-3}$]',r'$\rho$')#,zoom=[-20,-10])#,zoom=[0,0.2])
		#a.set_xlim([-9,-31.5])
		#a.text(a.get_xlim()[0]+0.04*(a.get_xlim()[1]-a.get_xlim()[0]),a.get_ylim()[1]-0.1*(a.get_ylim()[1]-a.get_ylim()[0]),'(d)',fontsize=18)
		a.set_yticks([400,500,600,700,800,900])
		print('sRho_tot')
		
		plt.tight_layout()
		plt.savefig(inputPath+'ini_temp.png')
		plt.savefig(inputPath+'ini_temp.pdf')
		plt.close()
		#quit()
	if plot_inidia:
		fig,ax = plt.subplots(ncols=2,nrows=2,figsize=(10,7))
		if minmax:
			mcTableTmp = mcTableTemp[(mcTableTemp['sMult']>vmin) & (mcTableTemp['sMult']<=vmax)]
		else:
			mcTableTmp = mcTableTemp
		a=plot.plotInitempVar(mcTableTmp,ax[0,0],inputPath,'dia','dia','[m]',r'D$_{\rm max}$',xlog=True)#,zoom=[-20,-10])#,zoom=[0,0.2])
		a.text(10**-5,1.75*10**(-3),'(a)',fontsize=18)
		print('dia')
		a=plot.plotInitempVar(mcTableTmp,ax[0,1],inputPath,'mTot','dia','[kg]','mass',ylog=True,xlog=True)#,zoom=[-20,-10])#,zoom=[0,0.2])
		a.text(10**-5,3.4*10**(-7),'(b)',fontsize=18)
		print('mTot')
		#quit()
		a=plot.plotInitempVar(mcTableTmp,ax[1,0],inputPath,'sPhi','dia','',r'$\phi$',ylog=True,xlog=True,ylim=[6*10**-2,10**0])#,zoom=[-20,-10])#,zoom=[0,0.2])
		a.text(10**-5,8*10**(-1),'(c)',fontsize=18)
		print('sPhi')
		a=plot.plotInitempVar(mcTableTmp,ax[1,1],inputPath,'sRho_tot','dia',r'[kgm$^{-3}$]',r'$\rho$',xlog=True)#,zoom=[-20,-10])#,zoom=[0,0.2])
		a.text(10**-5,a.get_ylim()[1]-0.1*(a.get_ylim()[1]-a.get_ylim()[0]),'(d)',fontsize=18)
		print('sRho_tot')

		plt.tight_layout()
		plt.savefig(inputPath+'ini_dia.png')
		plt.savefig(inputPath+'ini_dia.pdf')
		plt.close()
	
	
	###################################################################
	# plot particle evolution of one case (Figure 5.2 in my thesis)
	###################################################################
	
	fig = plt.figure(figsize=(11,8))
	gs = mpl.gridspec.GridSpec(2, 2)        
	gs.update(wspace=0.1,hspace=0.35)
	ax00 = fig.add_subplot(gs[0, 0]) 
	ax01 = fig.add_subplot(gs[0, 1]) 
	ax10 = fig.add_subplot(gs[1, 0]) 
	ax11 = fig.add_subplot(gs[1, 1]) 
	ax = np.array([[ax00,ax01],[ax10,ax11]])
	mcTableTmp = mcTableTemp
	vmin1 = 0; vmax1 = 180; vmin2=180; vmax2=350
	#mcTableTmp1 = mcTableTemp[(mcTableTemp['sMult']>vmin1) & (mcTableTemp['sMult']<=vmax1)]
	#mcTableTmp2 = mcTableTemp[(mcTableTemp['sMult']>vmin2) & (mcTableTemp['sMult']<=vmax2)]
	ax00,s_m = plot.plot_var_particle_ID_temp(mcTableTmp,'dia','sMult',r'D$_{\rm max}$','[m]',log=True,ax=ax00,cbar=False)#,zoom=[-20,-10])
	#ax00,s_m = plot.plot_var_particle_ID_temp(mcTableTmp2,inputPath,'dia','sMult',r'D$_{\rm max}$','[m]',onlySpecPart=[vmin2,vmax2],log=True,ax=ax00,cbar=False)#,zoom=[-20,-10])

	print('dia')
	#ax00.text(3*10**-3,-28,'(a)',fontsize=18)

	#plot.plot_var_particle_ID_temp(mcTableTemp,inputPath,'dia_mum','sMult')
	ax01,s_m =plot.plot_var_particle_ID_temp(mcTableTmp,'mTot','sMult','mass','[kg]',log=True,ax=ax01,cbar=False,xlabel=False)#,zoom=[-20,-10])
	print('mTot')
	#ax01.text(2*10**-7,-28,'(b)',fontsize=18)
	ax01.set_xticks([10**-12,10**-10,10**-8,10**-6])
	y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
	ax01.xaxis.set_minor_locator(y_minor)
	ax01.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

	ax10,s_m = plot.plot_var_particle_ID_temp(mcTableTmp,'sPhi','sMult',r'$\phi$','',log=True,ax=ax10,cbar=False)#,zoom=[-30,-10])#,zoom=[0,0.2])
	#ax10.text(4,-28,'(c)',fontsize=18)
	ax10.set_xticks([10**-2,10**-1,10**0,10**1])
	y_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
	ax10.xaxis.set_minor_locator(y_minor)
	ax10.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())

	ax11,s_m =plot.plot_var_particle_ID_temp(mcTableTmp,'sRho_tot','sMult',r'$\rho$',r'[kgm$^{-3}$]',ax=ax11,cbar=False,xlabel=False)#,zoom=[-20,-10])
	#ax11.text(875,-28,'(d)',fontsize=18)
	ax11.set_xticks([900,700,500,300])
	print('sRhoice')
	for a in [ax00,ax01,ax10,ax11]:
		a.set_ylim([0,-31.5])
		a.axhline(y=-10,ls='--',lw=2,c='r')
		a.axhline(y=-20,ls='--',lw=2,c='r')
	for a in [ax00,ax01,ax10,ax11]:
		a.set_yticks([0,-5,-10,-15,-20,-25,-30])
		a.axhline(y=-10,ls='--',lw=2,c='r')
		a.axhline(y=-20,ls='--',lw=2,c='r')
	cbar = fig.colorbar(s_m,ax=ax,pad=0.05,aspect=30,shrink=0.75)
	cbar.set_label('particle ID',fontsize=18)
	cbar.set_ticks(np.arange(0,mcTableTmp['sMult'].max()+1,50))
	cbar.ax.tick_params(labelsize=16)
    # now plot two particles with black color for better illustration
	if plot_thesis_particle_evolution:
		print(mcTableTempXR)
		mcTableTmp = mcTableTmp.set_index('sMult',drop=False)
		outputSel=mcTableTmp.loc[180]
		outputSel=mcTableTmp.loc[10]
		for var2plot,axes in zip(['dia','mTot','sPhi','sRho_tot'],[ax00,ax01,ax10,ax11]):
			prop2plot = outputSel[var2plot]
			axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,ls=':',label='Particle A')
			
		outputSel=mcTableTmp.loc[100]
		for var2plot,axes in zip(['dia','mTot','sPhi','sRho_tot'],[ax00,ax01,ax10,ax11]):
			prop2plot = outputSel[var2plot]
			axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,ls='-.',label='Particle B')

		outputSel=mcTableTmp.loc[180]
		for var2plot,axes in zip(['dia','mTot','sPhi','sRho_tot'],[ax00,ax01,ax10,ax11]):
			prop2plot = outputSel[var2plot]
			axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,label='Particle C')
			

		outputSel=mcTableTmp.loc[250]
		for var2plot,axes in zip(['dia','mTot','sPhi','sRho_tot'],[ax00,ax01,ax10,ax11]):
			prop2plot = outputSel[var2plot]
			axes.plot(prop2plot,outputSel.Temp,c='k',lw=2,ls='--',label='Particle D')



		ax00.legend(fontsize=12,ncol=4,bbox_to_anchor=(0,1,1.5,0.05))
    #ax = np.array([[ax00,ax01],[ax10,ax11]])
    #plt.show()
		plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/particleEvolution_RHi105_particleK_Ty_newColor2.png',bbox_inches='tight')
		plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/particleEvolution_RHi105_particleK_Ty_newColor2.pdf',bbox_inches='tight')
	else:
		plt.savefig(inputPath+'particleEvolution.png',bbox_inches='tight')
	plt.close()
	
	if os.path.exists(inputPath+McRadar_Outname):
		vmin=0;vmax=350
		print(inputPath+McRadar_Outname)
		#quit()mcTableTmp = mcTableTemp[(mcTableTemp['sMult']>vmin) & (mcTableTemp['sMult']<=vmax)]
		output = xr.open_dataset(inputPath+McRadar_Outname)
		
		#quit()
		#outPut1 = xr.open_dataset(inputPath+os.environ['freq']+'GHz_output_{mode}_180_350_singleParticle.nc'.format(mode=dicSettings['scatSet']['mode']))
		#outPut = xr.merge([outPut,outPut1])
		
		print('now plotting McRadar')

		if plotTemp == True:
			atmoFile = np.loadtxt(inputPath+'atmo.dat')
			height = atmoFile[:,0]
			Temp = atmoFile[:,2] -273.15
			atmoPD = pd.DataFrame(data=Temp,index=height,columns=['Temp'])
			atmoPD.index.name='sHeight'
			

			atmoXR = atmoPD.to_xarray()
			atmoReindex = atmoXR.reindex_like(output,method='nearest')
			output = xr.merge([atmoReindex,output])
		#print(dicSettings['wl'])
		wlStr = '{:.2e}'.format(dicSettings['wl'][0])
		output['sZDR'] = 10*np.log10(output['sZeH_{0}'.format(wlStr)]) - 10*np.log10(output['sZeV_{0}'.format(wlStr)])
		#output = output.to_dataframe()
		
		##################################################################################################
		# plot single particle ZDR and single particle KDP with potentially four particles for better visibility (when setting plot_four_particles to true) (Figure 5.3 in my thesis)
		##################################################################################################
		fig = plt.figure(figsize=(10,4))
		gs = mpl.gridspec.GridSpec(1, 2)        
		axes = [fig.add_subplot(gs[0, col]) for col in range(2)]

		axes[0],s_m = plot.plot_scat_particle_ID_temp(output,inputPath,'sZDR','sMult',r'ZDR$_{\rm particle}$','[dB]',ax=axes[0],cbar=False)
		axes[0].text(0,-29,'(a)',fontsize=18)
		axes[0].set_yticks([0,-5,-10,-15,-20,-25,-30])

		outputSel = output.sel(sMult=slice(0,180))
		axes[1],s_m1=plot.plot_scat_particle_ID_temp(output,inputPath,'sKDP_{0}'.format(wlStr),'sMult',r'KDP$_{\rm particle}$',r'[°km$^{-3}$]',ax=axes[1],cbar=False,y_label=False)
		outputSel = output.sel(sMult=slice(180,350))
		#ax1 = axes[1].twiny()
		#ax1,s_m2=plot.plot_scat_particle_ID_temp(outputSel,inputPath,'sKDP_3.12e+01','sMult',r'KDP$_{\rm particle}$',r'[°km$^{-3}$]',ax=ax1,cbar=False,y_label=False)
		axes[1].set_yticks([0,-5,-10,-15,-20,-25,-30])
		#ax1.grid(False)
		#x_minor = mpl.ticker.LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 10)
		#axes[1].xaxis.set_minor_locator(x_minor)
		#axes[1].xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
		#ax1.set_xlim([-0.0011,0.0015])

		cbar = fig.colorbar(s_m,ax=axes,pad=0.025)
		cbar.set_label('particle ID',fontsize=16)
		cbar.ax.tick_params(labelsize=14)
		cbar.set_ticks(np.arange(vmin,vmax+1,50))
		axes[1].text(0,-29,'(b)',fontsize=18)

		if plot_thesis_particle_evolution:
			outputSel=output.sel(sMult=10)
			prop2plot = outputSel['sZDR']
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			axes[0].plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',ls=':',lw=2,label='Particle A')
			prop2plot = outputSel['sKDP_3.12e+01']
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			axes[1].plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',ls=':',lw=2,label='Particle A')
			
			outputSel=output.sel(sMult=100)
			prop2plot = outputSel['sZDR']
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			axes[0].plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',ls='-.',lw=2,label='Particle B')
			prop2plot = outputSel['sKDP_3.12e+01']
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			axes[1].plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',ls='-.',lw=2,label='Particle B')
			
			outputSel=output.sel(sMult=180)
			prop2plot = outputSel['sZDR']
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			axes[0].plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',lw=2,label='Particle C')
			
			prop2plot = outputSel['sKDP_{0}'.format(wlStr)]
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			ax1.plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',lw=2,label='Particle C')
			
			outputSel=output.sel(sMult=250)
			prop2plot = outputSel['sZDR']
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			axes[0].plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',ls='--',lw=2,label='Particle D')
			prop2plot = outputSel['sKDP_{0}'.format(wlStr)]
			TempNaN = outputSel.Temp.where(~np.isnan(prop2plot))
			ax1.plot(prop2plot.dropna(dim='sHeight'),TempNaN.dropna(dim='sHeight'),c='k',ls='--',lw=2,label='Particle D')		
			ax1.axhline(y=-10,ls='--',lw=2,c='r')
			ax1.axhline(y=-20,ls='--',lw=2,c='r')
			axes[0].legend(bbox_to_anchor=(-0.1,1),loc='lower left',fontsize=12,ncol=2)

		# now calculate integrated ZDR and sZDRmax
		for i,a in enumerate(axes):
			a.axhline(y=-10,ls='--',lw=2,c='r')
			a.axhline(y=-20,ls='--',lw=2,c='r')
		#plt.tight_layout()
		plt.savefig(inputPath+'sZDR_sKDP_forward_sim_{mode}_200bins.png'.format(mode=scatMode),bbox_inches='tight')
		#plt.savefig(inputPath+'sZDR_sKDP_forward_sim.pdf',bbox_inches='tight')
		plt.close()
		quit()
		##################################################
		#- plot integrated ZDR and compare to observations (Figure 5.4 in my thesis)
		##################################################
		datasZDRmax = pd.read_csv('sZDRmax_median_DWRclass_2.txt',delimiter=' ',header=0)
		print(datasZDRmax)
		#quit()


		fig,ax = plt.subplots(ncols=2,figsize=(10,4),sharey=True)
		
		vmin=0;vmax=90
		outputSel = output.sel(sMult=slice(vmin,vmax))
		outputSel['ZeH'] = outputSel['sZeH_{0}'.format(wlStr)].sum(dim='sMult')
		outputSel['ZeV'] = outputSel['sZeV_{0}'.format(wlStr)].sum(dim='sMult')
		outputSel['ZDR'] = 10*np.log10(outputSel.ZeH)-10*np.log10(outputSel.ZeV)
		sZDRmax = outputSel.sZDR.max(dim='sMult')
		sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
		ZDR = outputSel.ZDR.rolling(sHeight=20,center=True).mean()
		ls1=ax[0].plot(ZDR,outputSel.Temp,c='#66CCEE',label='cold',lw=2)#c='#4477AA'
		ls2=ax[1].plot(sZDRmax,outputSel.Temp,ls='-',c='#66CCEE',label=r'sZDR$_{\rm max}$ cold',lw=2)#c='#66CCEE',
		
		output['ZeH'] = output['sZeH_{0}'.format(wlStr)].sum(dim='sMult')
		output['ZeV'] = output['sZeV_{0}'.format(wlStr)].sum(dim='sMult')
		output['ZDR'] = 10*np.log10(output.ZeH)-10*np.log10(output.ZeV)
		sZDRmax = output.sZDR.max(dim='sMult')
		sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
		ZDR = output.ZDR.rolling(sHeight=20,center=True).mean()
		ls3=ax[0].plot(ZDR,output.Temp,lw=2,c='#228833',label='all')
		ls4=ax[1].plot(sZDRmax,output.Temp,lw=2,c='#66CCEE',label=r'sZDR$_{\rm max}$ all')#c='#CCBB44'


		vmin=180;vmax=350
		outputSel = output.sel(sMult=slice(vmin,vmax))
		outputSel['ZeH'] = outputSel['sZeH_{0}'.format(wlStr)].sum(dim='sMult')
		outputSel['ZeV'] = outputSel['sZeV_{0}'.format(wlStr)].sum(dim='sMult')
		outputSel['ZDR'] = 10*np.log10(outputSel.ZeH)-10*np.log10(outputSel.ZeV)
		sZDRmax = outputSel.sZDR.max(dim='sMult')
		sZDRmax = sZDRmax.rolling(sHeight=20,center=True).mean()
		ZDR = outputSel.ZDR.rolling(sHeight=20,center=True).mean()
		ls5=ax[0].plot(ZDR,outputSel.Temp,c='#EE6677',label='DGL',lw=2)
		ls6=ax[1].plot(sZDRmax,outputSel.Temp,ls='-',c='#EE6677',label=r'sZDR$_{\rm max}$ DGL',lw=2)#c='#AA3377'

		ax1 = ax[1].twiny()
		ls7=ax1.plot(datasZDRmax.med,datasZDRmax.Temp,c='k',lw=2,label=r'obs')
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
		ls = ls1+ls3+ls5+ls7
		labs = [l.get_label() for l in ls]
		ax[1].legend(ls,labs,ncol=2,fontsize=12)
		ax[0].set_ylabel('T [°C]',fontsize=18)
		ax[0].set_xlabel('simulated ZDR [dB]',fontsize=18)
		ax[1].set_xlabel(r'simulated sZDR$_{\rm max}$ [dB]',fontsize=18)
		for i,a in enumerate(ax):
			a.set_xlim([0,7])
			a.set_ylim([0,-30])
			a.axhline(y=-10,ls='--',lw=2,c='r')
			a.axhline(y=-20,ls='--',lw=2,c='r')
			a.tick_params(axis='both',labelsize=16)
			a.set_yticks([0,-5,-10,-15,-20,-25,-30])
			a.set_xticks([0,1,2,3,4,5,6,7])
			a.text(0.1,-27,'('+string.ascii_lowercase[i]+')',fontsize=18)
			a.grid()
		#ax[0].set_zorder(1)
		#ax[0].legend(ls,labs,bbox_to_anchor=(-0,0.8),loc='lower left',fontsize=12,ncol=3)
		plt.savefig(inputPath+'ZDR_sZDRmax_obs_{mode}.png'.format(mode=scatMode),bbox_inches='tight')
		#plt.savefig('/home/lvonterz/thesis/plots/McSnow_chapter/ZDR_int3.pdf',bbox_inches='tight')
		plt.close()
		#quit()
		#ax.text(6,-29,'(c)',fontsize=18)

	
	
		
	








