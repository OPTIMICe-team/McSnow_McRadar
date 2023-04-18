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
source /home/L.Terzi/.bashrc

if [ "$2" == "1" ]; then #McSnow can somehow not be run with open evince instances in certain folders
    killall evince #otherwise there is an "directory not empty" error
fi

set -e #exit if error occurs

#if cheops is shut down use old folder for tests
export MC=/project/meteo/work/L.Terzi/McSnow_habit/mcsnow/ #foldder of mcsnow
export McR=/project/meteo/work/L.Terzi/McRadar/notebooks/ # folder were my McRadar notebooks are
export MCexp=/project/meteo/work/L.Terzi/McSnowoutput/habit/fragmentation #second_nucleation  # folder where to store output 
export postMcSnowDir=/project/meteo/work/L.Terzi/McSnow_habit/mcsnow/postprocessing/
export cur_dir=$(pwd -P) #current directory
export LUT_dir=/project/meteo/work/L.Terzi/McRadar/LUT/
export freq=94.0 #35.5 #5.6 #9.6 #35.5 #9.6
export scatMode=DDA #full #SSRGA-Rayleigh
export particle=vonTerzi_dendrite
export elv=30
export convolute=False
export singleParticle=False
export McRadarfileName="${freq}GHz_output_${scatMode}_${particle}_${elv}_convolute${convolute}_singleParticle${singleParticle}.nc"
#define what this script should do
if [ -z "$1" ]; then
    execwhat="r0Mc0dat2nc0McRad0plot1plotNCL0"  #recompile, McSnow, create ncdf, McRadar, plot Output, plotNCL output #set 1 to activate and 0 to deactivate one of these steps
else #the following allows the user to define what the script should do by (currently 5) booleans
    execwhat="r"$1"Mc"$2"dat2nc"$3"McRad"$4"plot"$5"plotNCL"$6
    echo $execwhat
fi
######
# choose setup (testcase)

export McSnow_testcase=("habit_fragmentation") #"habit_secondMode" #habit_trajectories
# if second mode and agg from above: in mo_bnd_cond.f90 line 120, 121 set monomer number to two and Vi=m/rho. Otherwise comment these lines out!!!
export adapt_version=3
#switch of processes (for McSnow this is so far not used)
export switch_off_processes_list=("_") # "_noicesnow_nosnowself" "_nosnowself" "_noicesnow") # _noiceself, _noicesnow, _nosnowself

for testcase in "${McSnow_testcase[@]}"
do
    ###########################################################
    #loop over different namelist settings 
    ##########################################################
    ssat_array=(50) #(0 10 20 30 40 50)  #supersaturation over ice: [1/1000] -> 1 is 0.001
    stick="2" #sticking efficiency: 0: E_s=1, 1: PK97, 2: Connolly12, 3: 0.5*Connolly12
    ncl_array=(0) #(10 20 50) #nucleation rate [10^-4 SP/sm3] #setting a high numer gets expensive (CPU-time and memory!); compensate this by a high multiplicity (xi0 in runscript: McSnow_runscripts/1d_bimodal2mode)
    nclmass_array=(1) #(100 1000 5000) # *10**-11 if using mass distr. change nclmass in 1d_habit # TODO: right now the setup calculated mtot from distribution of Dmax and aspect ratio!!
    McSnow_geom="2" # 1.constant 2.binary 3.monodep #we need binary here
    
    domTop_array=(5000) #(2700 2500 2300 2100) # change domain top
    #minmax_array=(180_350) # 0_180 180_350) # 190_200 210_220 250_260 280_290 310_320 340_350) #210_220 220_230 230_240 240_250 250_260 260_270 270_280 280_290 290_300)
    agg="2" # 0: no aggregation. 1: Golovin kernel 2: hydrodynamic kernel
    dep="1" # 0: no diffusion, 1: vapour diffusion
    xi0="50" 
    nz="200"
    iwc="3"
    coll_kern="0" #0:Area, 1:Dmax
    nugam="3.5" #"10" #"3.5" #"3.545"
    mugam="3.5" #"10" #"0.5" #"0.455"
    bndtype="2" #0: zero flux, 1: const xi=xi0, 2: xi ~ flux cdf, 3: discrete
    atmo="1" # 1: idealized profiles; 2: own profiles
    nrp0="2" # number of real particles
    nh1="0" #2300
    nh2="0" # 3000
    habit="1"
    IGF="2"
    iceicebreak_fpm2_array=(6) # 12 24 48) # default value is 6
    iceicebreak_fpm1=0 # default value is 1
    sp_kern_sig=10
    timeend="10800" #"6000" #"10800"
    dt_1dprof="3600" #"60" #"3600"
    
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
    
							
						echo "#################"
						echo "analyzing testcase: "$testcase "with ssat=" $ssat ", stick=" $stick ",ncl=," $ncl ",nclmass," $nclmass ",model_setup," $model_setup
						echo "#################"
						export testcase; export McSnow_geom; export model_setup; export switch_off_processes; export nclmass #make them visible for called routines
								
						export tstep="300" #start of timestep analyzed by McRadar and plotting routine in minutes # not implemented in MCRadar
						export tstep_end="600" #end timestep analyzed by McRadar and plotting routine in minutes (defines end of averaging period)# not implemented in MCRadar
						export av_tstep="5" #average the McSnow output in the PAMTRA adaption by ... s #do not set this lower than dtc (f.e. "5"s)# not implemented in MCRadar
						
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
						cp McSnow_runscripts/1d_${testcase%_*} $MC/run/1d
						# here the name of the folder is written into the environment for the other skripts to find
						export experiment=$($MC/run/1d "onlyname" $ssat $stick $testcase  $ncl $nclmass $domTop $agg $xi0 $nz $iwc $coll_kern $nugam $mugam $bndtype $atmo $nrp0 $nh1 $nh2 $habit $IGF $iceicebreak_fpm2 $iceicebreak_fpm1 $sp_kern_sig $dep $timeend $dt_1dprof)
						echo "analyze experiment:" $experiment
						
						if [[ "$execwhat" == *Mc1* ]] ; then #run McSnow (f.e 1d-model)
							echo "############"
							echo "run McSnow"
							echo "############"

							cd $MC/run #alternative if cheops is shut down
							
							# here McSnow is run with the setup you specify here.
							./1d "fullrun" $ssat $stick $testcase  $ncl $nclmass $domTop $agg $xi0 $nz $iwc $coll_kern $nugam $mugam $bndtype $atmo $nrp0 $nh1 $nh2 $habit $IGF $iceicebreak_fpm2 $iceicebreak_fpm1 $sp_kern_sig $dep $timeend $dt_1dprof $MCexp 
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

							cd $postMcSnowDir
							echo $(pwd -P)
							./ncl_post.sh $MCexp'/'$experiment
							
					 
						fi
                 


    				done #for icebreak
				done #for domTop
			done #for nclmass
		done #for ncl
    done #for ssat
    
done #for testcase


cd $MCexp #go back to start directory
#1d_habit_habit1_xi1_nz300_iwc0.01_nugam3.545_mugam0.455_dtc5_nrp1_vt3_coll_kern1_at0_stick2_colleffi1_dt1_bndtype3_ba500_domtop5000._atmo2_RHi105_nucleating_10i_interpPPTrue

