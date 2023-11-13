#!/bin/bash



# THIS SCRIPT DOES ALL THE WORK: 
#- compiling and running McSnow 
#- run McRadar on the output
#- plot the output of McSNow and PAMTRA
# 

#############
# you might not always need to do all together (McSnow+McRadar+ plotting)
#, so you can give some input arguments to deactivate some of the procedures: e.g. "bash McSnow_McRadar.sh 1 0 0 0 0" will only compile McSnow
#INPUT ARGUMENTS: Boolean (0) dont do (1) do
#$1 RECOMPILE
#$2 RUN MCSNOW
#$3 CREATE NCDF FROM MCSNOW RUN
#$4 RUN McRadar
#$5 DO PLOTTING
#$6 run ncl postskript
#

#############
source /home/l/L.Terzi/.bashrc

set -e #exit if error occurs

#module load python/3.9-2021.11 # needed to run xarray without mistakes. TODO check regularly if the standard python changes to python 3.9!! 

#if cheops is shut down use old folder for tests
export MC=/project/meteo/work/L.Terzi/McSnow_habit/mcsnow/ #foldder of mcsnow
export McR=/project/meteo/work/L.Terzi/McRadar/notebooks/ # folder were my McRadar notebooks are
export MCexp=/project/meteo/work/L.Terzi/McSnowoutput/habit/case_studies/20220206/  # folder where to store output 
export postMcSnowDir=/project/meteo/work/L.Terzi/McSnow_habit/mcsnow/postprocessing/
export cur_dir=$(pwd -P) #current directory
export LUT_dir=/project/meteo/work/L.Terzi/McRadar/LUT/
export freq=9.6_35.5_94.0 #35.5 #9.6
export scatMode=DDA
export particle=vonTerzi_dendrite
export elv=30_90
export convolute=True
export singleParticle=False
export attenuation=False
export McRadarfileName="${freq}GHz_output_${scatMode}_${particle}_${elv}_convolute${convolute}_singleParticle${singleParticle}_attenuation${attenuation}_testKDP.nc"
#define what this script should do
if [ -z "$1" ]; then
    execwhat="r0Mc0dat2nc0McRad0plot1plotNCL0genTS0"  #recompile, McSnow, create ncdf, McRadar, plot Output, plotNCL output #set 1 to activate and 0 to deactivate one of these steps
else #the following allows the user to define what the script should do by (currently 5) booleans
    execwhat="r"$1"Mc"$2"dat2nc"$3"McRad"$4"plot"$5"plotNCL"$6"genTS"$7
    echo $execwhat
fi
######
# choose setup (testcase)

export McSnow_testcase=("habit_case") #"habit_secondMode" #habit_trajectories
# if second mode and agg from above: in mo_bnd_cond.f90 line 120, 121 set monomer number to two and Vi=m/rho. Otherwise comment these lines out!!!

for testcase in "${McSnow_testcase[@]}"
do
    ###########################################################
    #loop over different namelist settings 
    ##########################################################
    domTop_array=(9000) #(2700 2500 2300 2100) # change domain top
    atmo="2" # 1: idealized profiles; 2: own profiles
    ssat_array=(50) # 50 100 200) #(0 10 20 30 40 50)  #supersaturation over ice: [1/1000] -> 1 is 0.001
    atmoFileName="radiosondes_juelich_20220206_042141" #"radiosondes_juelich_20220206_042141" #"radiosondes_juelich_20220215_194355"
    xi0="50" 
    nz="200"
    bndtype="2" #0: zero flux, 1: const xi=xi0, 2: xi ~ flux cdf, 3: discrete
    timeend="10800" #"6000" #"10800"
    dt_1dprof="3600" #"10" #"60" #"3600"
    
    #aggregation parameters:
    agg_array=(2) # 0: no aggregation. 1: Golovin kernel 2: hydrodynamic kernel
    stick="2" #sticking efficiency: 0: E_s=1, 1: PK97, 2: Connolly12, 3: 0.5*Connolly12, 10: Markus sticking efficiency, 11: 0.5*Markus sticking efficiency
    agggeo="2" # 2: mitchell, 3: monodependent, 4: dendrite aggregates
  	am="0.03" # TODO: adapt to follow mD of our dendrites
    bm="2.2"
    aa="0.08"
    ba="1.9"
    sp_kern_sig_array=(10) # 20 30) # 50)
    McSnow_geom="2" # 1.constant 2.binary 3.monodep #we need binary here
    coll_kern="0" #0:Area, 1:Dmax
    
    #depo growth
    dep="1" # 0: no diffusion, 1: vapour diffusion
    habit="1"
    IGF="2"
   
    #- setup of nucleation
    iwc_array=(2) # 5 7 10) #2 3) # [1/100000 kg/m3]
    ncl_array=(100) # 50 100 500) #(10 20 50) #nucleation rate [10^-4 SP/sm3] #setting a high numer gets expensive (CPU-time and memory!); compensate this by a high multiplicity (xi0 in runscript: McSnow_runscripts/1d_bimodal2mode)
    nclmass_array=(4.8) # 10 100 1000) # *10**-13 if using mass distr. change nclmass in 1d_habit # 
    nugam_array=(0) # 3.5 5 10) #"3.5" #"10" #"3.5" #"3.545"
    mugam_array=(0) # 3.5 5 10) #"0.5" #"10" #"0.5" #"0.455" 3.5 here makes it much more narrow!!
    nrp0_array=(0) # 2 5 7) # 2 3) # number of real particles
    #nh1_array=(1800 2000 2200 2300)
    #nh2_array=(3000 3500 4000) # 3500 4000) # 3000
    nhwidth_array=(1500) #250 500 750 1000) # half width of layer
    nucl_type="5" #nucleation: 0) off 1) constant rate 2) Meyers 5) simple nucleation with 2 layers
    
    #- fragmentation
    iceicebreak_fpm2_array=(0) # 12 24 48) # default value is 6
    iceicebreak_fpm1=0 # default value is 1
    Dmode="150"
    
    for (( i_ssat = 0 ; i_ssat < ${#ssat_array[@]} ; i_ssat++ )); do
    	ssat=${ssat_array[$i_ssat]}
		for (( i_ncl = 0 ; i_ncl < ${#ncl_array[@]} ; i_ncl++ )); do
			ncl=${ncl_array[$i_ncl]}
			for (( i_nclmass = 0 ; i_nclmass < ${#nclmass_array[@]} ; i_nclmass++ )); do
				nclmass=${nclmass_array[$i_nclmass]}
				for (( i_domtop = 0 ; i_domtop < ${#domTop_array[@]} ; i_domtop++ )); do
					domTop=${domTop_array[$i_domtop]} 
					for (( i_ice = 0 ; i_ice < ${#iceicebreak_fpm2_array[@]} ; i_ice++ )); do
						iceicebreak_fpm2=${iceicebreak_fpm2_array[$i_ice]}  
						#for (( i_nh = 0 ; i_nh < ${#nhwidth_array[@]} ; i_nh++ )); do
						#	nhwidth=${nhwidth_array[$i_nh]} 
						for (( i_agg = 0 ; i_agg < ${#agg_array[@]} ; i_agg++ )); do
							agg=${agg_array[$i_agg]} 
							#nh1=$(echo "2420-${nhwidth}" | bc -l) # define nh1 and nh2, which is centered around 2420 meters, which is the -15°C isotherm
							#nh2=$(echo "2420+${nhwidth}" | bc -l)
							#nh1=$(echo "3300-${nhwidth}" | bc -l) # define nh1 and nh2, which is centered around 2420 meters, which is the -15°C isotherm
							#nh2=$(echo "3300+${nhwidth}" | bc -l) # 1650=-10°,3200=-20°,4950=-30°
							#nh1=1650 
							#nh2=4950
							# looking at 20220206 042141, -10°C=3000, -20=4900; -30=6250; -40=7500; -50: 8500
							nh1=4900
							nh2=8500
							#looking at 20220215 19:43:55: -10=3100 -20=4600 -30=6000 -40=7500 -50=9000
							#nh1=3100
							#nh2=4600
							for (( i_nu = 0 ; i_nu < ${#nugam_array[@]} ; i_nu++ )); do
								nugam=${nugam_array[$i_nu]} 
								for (( i_mu = 0 ; i_mu < ${#mugam_array[@]} ; i_mu++ )); do
									mugam=${mugam_array[$i_mu]} 
									for (( i_iwc = 0 ; i_iwc < ${#iwc_array[@]} ; i_iwc++ )); do
										iwc=${iwc_array[$i_iwc]} 
										for (( i_nrp = 0 ; i_nrp < ${#nrp0_array[@]} ; i_nrp++ )); do
											nrp0=${nrp0_array[$i_nrp]} 
							 				for (( i_kern = 0 ; i_kern < ${#sp_kern_sig_array[@]} ; i_kern++ )); do
												sp_kern_sig=${sp_kern_sig_array[$i_kern]}
											
												echo "#################"
												echo "analyzing testcase: "$testcase "with ssat=" $ssat ", stick=" $stick ",ncl=," $ncl ",nclmass," $nclmass ",model_setup," $model_setup
												echo "#################"
												export testcase; export McSnow_geom; export model_setup; export switch_off_processes; export nclmass #make them visible for called routines
														
												if [[ "$execwhat" == *r1* ]] ; then #recompile for a change in the model
													echo "############"
													echo "recompiling"
													echo "############"

													cd $MC

													#alternative if cheops is shut down
													make release_omp #alternative if cheops is shut down
												fi
												#get foldername from current runscript (even when it is not executed
												cd $cur_dir
												cp McSnow_runscripts/1d_${testcase} $MC/run/1d #McSnow_runscripts/1d_${testcase%_*} $MC/run/1d
												#if [[ $atmo == 2 ]]; then
												#	export atmoFileName="atmoLeonie_nh1${nh1}_nh2${nh2}_ssat${ssat}"
												#	python3 make_atmo_table.py $nh1 $nh2 $ssat $atmoFileName
												#fi
												# here the name of the folder is written into the environment for the other skripts to find
												export experiment=$($MC/run/1d "onlyname" $testcase $domTop $atmo $ssat $atmoFileName $xi0 $nz $bndtype $timeend $dt_1dprof $agg $stick $agggeo $am $bm $aa $ba $sp_kern_sig $coll_kern $dep $habit $IGF $iwc $ncl $nclmass $nugam $mugam $nrp0 $nh1 $nh2 $nucl_type $iceicebreak_fpm2 $iceicebreak_fpm1 $Dmode)
												echo "analyze experiment:" $experiment
												
												if [ -f "${MCexp}${experiment}.tar.gz" ]; then
													echo "experiment already exists, not rerunning"
												else
																										
												if [[ "$execwhat" == *Mc1* ]] ; then #run McSnow (f.e 1d-model)
													echo "############"
													echo "run McSnow"
													echo "############"
													# if you are running with ssat enhanced in second nucleation layer:
													#atmoFileName="atmoLeonie_nh1${nh1}_nh2${nh2}_ssat${ssat}.txt"
												
													
													cd $MC/run #alternative if cheops is shut down
													
													# here McSnow is run with the setup you specify here.
													./1d "fullrun" $testcase $domTop $atmo $ssat $atmoFileName $xi0 $nz $bndtype $timeend $dt_1dprof $agg $stick $agggeo $am $bm $aa $ba $sp_kern_sig $coll_kern $dep $habit $IGF $iwc $ncl $nclmass $nugam $mugam $nrp0 $nh1 $nh2 $nucl_type $iceicebreak_fpm2 $iceicebreak_fpm1 $Dmode $MCexp 
													fi

													if [[ "$execwhat" == *dat2nc1* ]] ; then #run McSnow (f.e 1d-model)
														echo "############"
														echo "compiling netcdf file from mass2fr.dat file"
														echo "############"

														cd $cur_dir 
														python3 McSnow_dat2ncdf.py        
													fi


													if [[ "$execwhat" == *McRad1* ]] ; then #let McRadar run on output
														echo "############"
														echo "run McRadar"
														echo "############"
														cd $cur_dir
														
														python3 calc_McRadar_output.py
														
													fi
													 if [[ "$execwhat" == *plot1* ]] ; then #produce quicklook of McSnow and the corresponding McRadar run
														echo "############"
														echo "plotting"
														echo "############"

														cd $cur_dir
														#python3 plot_PSD.py
														#python3 plotModObs.py
														#python3 plot_agg_kernel.py
														#for (( i_minmax = 0 ; i_minmax < ${#minmax_array[@]} ; i_minmax++ )); do
														#export minmax=${minmax_array[$i_minmax]}
														python3 plot_output.py # sofar this produces plots of the McRadar moments and spectra, aswell as the aspect ratios in spectral form
														#done # for minmax
													fi
													if [[ "$execwhat" == *plotNCL1* ]] ; then #produce quicklook of McSnow and the corresponding McRadar run
														echo "############"
														echo "plotting nclscript"
														echo "############"

														#cd $postMcSnowDir
														#echo $(pwd -P)
														#./ncl_post.sh $MCexp'/'$experiment
														cd $cur_dir
														python3 calc_PSD.py
														#cd $McR
														#python3 generate_ts.py
														#rm -v !("filename1"|"filename2")
												 		#cd ${MCexp}${experiment}
												 		#tar -czvf ${MCexp}${experiment}.tar.gz ${MCexp}${experiment} && rm -r ${MCexp}${experiment}
												 		#shopt -s extglob
												 		#rm -v !("mass2fr.nc"|"atmo.dat"|"gamma_psd_parameters.txt"|"parameters.feather")
													fi
													if [[ "$execwhat" == *genTS1* ]] ; then #produce quicklook of McSnow and the corresponding McRadar run
														echo "############"
														echo "generating ts"
														echo "############"
														cd $McR
														python3 generate_ts.py
												 		tar -czvf ${MCexp}${experiment}.tar.gz ${MCexp}${experiment} && rm -r ${MCexp}${experiment}
												fi
							 					fi # for the check if tar.gz already exists!
							 				done
							 			done # for nrp0
							 		done # for iwc
		             			done # for mu_gam
							done # for nu_gam
						done # for nhwidth
    				done #for icebreak
				done #for domTop
			done #for nclmass
		done #for ncl
    done #for ssat
    
done #for testcase


cd $MCexp #go back to start directory
#1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo2_RHi105_nucleating_10i_interpPPTrue

