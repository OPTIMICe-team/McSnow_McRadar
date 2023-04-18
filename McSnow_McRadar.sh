#!/bin/bash



# THIS SCRIPT DOES ALL THE WORK: 
#- compiling and running McSnow 
#- run McRadar on the output
#- plot the output of McSNow and PAMTRA
# 

#############
# you might not always need to do all together (McSnow+McRadar+ plotting)
#, so you can give some input arguments to deactivate some of the procedures: e.g. "bash McSnow_McRadar.sh 1 0 0 0 0 0" will only compile McSnow
#INPUT ARGUMENTS: Boolean (0) dont do (1) do
#$1 RECOMPILE
#$2 RUN MCSNOW
#$3 CREATE NCDF FROM MCSNOW RUN
#$4 RUN McRadar
#$5 DO PLOTTING
#$6 run ncl postskript
#

#############

#TODO change all paths to your setup

#if cheops is shut down use old folder for tests
export MC=/project/meteo/work/L.Terzi/McSnow_habit/mcsnow/ #foldder of mcsnow {Path2McSnow} #
export McR=/project/meteo/work/L.Terzi/McRadar/notebooks/ # folder were my McRadar notebooks are {Path2McRadar} #
export MCexp=/project/meteo/work/L.Terzi/McSnowoutput/habit/second_nucleation/ #second_nucleation  # folder where to store output  {Path2McSnowoutput} #
export postMcSnowDir=${MC}postprocessing/
export cur_dir=$(pwd -P) #current directory
export LUT_dir=/project/meteo/work/L.Terzi/McRadar/LUT/ #{Path2McRadarLUT} 
export freq=94.0 #This controls the frequency that is input into McRadar. Separate the frequencies with "_". Also, only a single frequency is possible.
export scatMode=DDA #This controls which scattering mode to use. For you SSRGA-Rayleigh makes the most sense This uses SSRGA LUTs generated from snowScatt for aggregates and Rayleigh for monomers. 
					#Other possibilities: full (full Tmatrix calculation), DDA (DDA LUT) among others
export particle=dendrite # this is the name of the particle you want to use for scattering. In case of SSRGA, currently available are vonTerzi_dendrite, vonTerzi_column,.. (see List in mcradar/settings.py)
export elv=90 # elevation to use. Here you can also choose multiple elevations, again separated by "_" (e.g. 30_90)
export convolute=True # convolute noise and broadening terms onto Doppler spectrum. Suggested: True (the variables such as noise power, eddy dissipation rate etc can be adjusted in McRadar)

export McRadarfileName="${freq}GHz_output_${scatMode}_${particle}_${elv}.nc" # this is the name of your McRadar output file. It will be stored in MCexp
#define what this script should do
if [ -z "$1" ]; then
    execwhat="r0Mc0dat2nc0McRad0plot1plotNCL0"  #recompile, McSnow, create ncdf, McRadar, plot Output, plotNCL output #set 1 to activate and 0 to deactivate one of these steps
else #the following allows the user to define what the script should do by (currently 5) booleans
    execwhat="r"$1"Mc"$2"dat2nc"$3"McRad"$4"plot"$5"plotNCL"$6
    echo $execwhat
fi
######
# choose setup (testcase)

export McSnow_testcase=("habit") #"habit_secondMode" #habit_trajectories

for testcase in "${McSnow_testcase[@]}"
do
  ###########################################################
    #loop over different namelist settings 
    ##########################################################
    ssat_array=(50) #(0 10 20 30 40 50)  #supersaturation over ice: [1/1000] -> 1 is 0.001
    stick="2" #sticking efficiency: 0: E_s=1, 1: PK97, 2: Connolly12, 3: 0.5*Connolly12
    McSnow_geom="2" # 1.constant 2.binary 3.monodep #we need binary here
    
    domTop_array=(5000) #(2700 2500 2300 2100) # change domain top
    
    agg="2" # 0: no aggregation. 1: Golovin kernel 2: hydrodynamic kernel
    dep="1" # 0: no diffusion, 1: vapour diffusion
    xi0="50" # multiplicity
    nz="200" # number of vertical cells
    iwc="3" # IWC at top [1/100000 kg/m3]
    coll_kern="0" #0:Area, 1:Dmax
    bndtype="2" #0: zero flux, 1: const xi=xi0, 2: xi ~ flux cdf, 3: discrete
    atmo="1" # 1: idealized profiles; 2: own profiles
    nrp0="2" # number of real particles
    habit="1" # 0: m-D-relation, 1: habit prediciton TODO: check if you have with or without habit prediction
    rt="0" #0: no riming, 1: continuous, 2: stochastic
    h1="2000" # bottom of liquid water layer
    h2="3000" # top of liquid water layer
    lwc="2" # lwc [1/10000 kg/m3]
    timeend="10800" # end time of McSnow simulations in s
    dt_1dprof="3600" # output step of McSnow to mass2fr.dat
    
    
    for (( i_ssat = 0 ; i_ssat < ${#ssat_array[@]} ; i_ssat++ )); do # loop over relative humidities
	    ssat=${ssat_array[$i_ssat]}
			for (( i_domtop = 0 ; i_domtop < ${#domTop_array[@]} ; i_domtop++ )); do # loop over domain tops
				domTop=${domTop_array[$i_domtop]} 
						
						echo "#################"
						echo "analyzing testcase: "$testcase "with ssat=" $ssat ", stick=" $stick ",model_setup," $model_setup
						echo "#################"
						export testcase; export McSnow_geom; export model_setup;  #make them visible for called routines
								
						
						if [[ "$execwhat" == *r1* ]] ; then #recompile for a change in the model
							echo "############"
							echo "recompiling"
							echo "############"

							cd $MC

							make release_omp 
						fi
						#get foldername from current runscript 
						cd $cur_dir
						cp McSnow_runscripts/1d_${testcase%_*} $MC/run/1d
						
						# here the name of the folder is written into the environment for the other skripts to find
						export experiment=$($MC/run/1d "onlyname" $ssat $stick $testcase  $domTop $agg $xi0 $nz $iwc $coll_kern $bndtype $atmo $nrp0 $habit $dep $rt $h1 $h2 $lwc $timeend $dt_1dprof )
						echo "analyze experiment:" $experiment
						
						if [[ "$execwhat" == *Mc1* ]] ; then #run McSnow (f.e 1d-model)
							echo "############"
							echo "run McSnow"
							echo "############"

							cd $MC/run 
							
							# here McSnow is run with the setup you specify here.
							./1d "fullrun" $ssat $stick $testcase  $domTop $agg $xi0 $nz $iwc $coll_kern $bndtype $atmo $nrp0 $habit $dep $rt $h1 $h2 $lwc $timeend $dt_1dprof $MCexp 
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
							
							python3 plot_output.py # this produces plots of the McRadar moments and spectra, aswell as microphysical spectra (i.e. aspect ratio in dependency of fall velocity)
							
						fi
						if [[ "$execwhat" == *plotNCL1* ]] ; then #produce quicklook of McSnow and the corresponding McRadar run
							echo "############"
							echo "plotting nclscript"
							echo "############"

							cd $postMcSnowDir
							echo $(pwd -P)
							./ncl_post.sh $MCexp'/'$experiment
							
					 
						fi
                 


    				
				done #for domTop
    done #for ssat
    
done #for testcase


cd $MCexp #go back to start directory


