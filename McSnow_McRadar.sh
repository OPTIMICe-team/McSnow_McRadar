#!/bin/bash



# THIS SCRIPT DOES ALL THE WORK: 
#- compiling and running McSnow 
#- run McRadar on the output
#- plot the output of McSNow and McRadar
# 

#############
# you might not always need to do all together (McSnow+McRadar+ plotting)
#, so you can give some input arguments to deactivate some of the procedures: e.g. "bash McSnow_Pamtra_quicklooks_bimodalityexp.sh 1 0 0 0 0" will only compile McSnow
#INPUT ARGUMENTS: Boolean (0) dont do (1) do
#$1 RECOMPILE
#$2 RUN MCSNOW
#$3 CREATE NCDF FROM MCSNOW RUN
#$4 RUN McRadar
#$5 DO PLOTTING
#
# 
#############
if [ "$2" == "1" ]; then #McSnow can somehow not be run with open evince instances in certain folders
    killall evince #otherwise there is an "directory not empty" error
fi

set -e #exit if error occurs

export MC=/work/lvonterz/McSnow_habit/mcsnow/ #foldder of mcsnow
export McR=/work/lvonterz/McRadar/notebooks/ # folder were my McRadar notebooks are
export MCexp=/work/lvonterz/McSnow_habit/  # folder where to store output 
export cur_dir=$(pwd -P) #current directory

#define what this script should do
if [ -z "$1" ]; then
    execwhat="r0Mc0dat2nc0McRad0plot1" #recompile, McSnow, create ncdf, McRadar, plot Output #set 1 to activate and 0 to deactivate one of these steps
else #the following allows the user to define what the script should do by (currently 5) booleans
    execwhat="r"$1"Mc"$2"dat2nc"$3"McRad"$4"plot"$5
    echo $execwhat
fi
######
# choose setup of McRadar
######
export freq=35.5 #5.6 #9.6 #35.5 #9.6
export scatMode=SSRGA
###
#use next line to plot McSnow runs with different geometries
###
export McSnow_testcase=("habit") 
export adapt_version=3

for testcase in "${McSnow_testcase[@]}"
do
    ###########################################################
    #loop over different namelist settings (relevant for running McSnow not for plotting)
    ##########################################################
    ssat_array=(50) #(0 10 20 30 40 50)  #supersaturation over ice: [1/1000] -> 1 is 0.001
    stick_array=(2) #sticking efficiency: 0: E_s=1, 1: PK97, 2: Connolly12, 3: 0.5*Connolly12
    ncl_array=(10) #(10 20 50) #nucleation rate [3/100 SP/sm3] #setting a high numer gets expensive (CPU-time and memory!); compensate this by a high multiplicity (xi0 in runscript: McSnow_runscripts/1d_bimodal2mode)
    nclmass_array=(10) #(100 1000 5000) #
    McSnow_geom=2 # 1.constant 2.binary 3.monodep #we need binary here
    len_namelistcomb=${#ssat_array[@]}
    for (( i_ssat = 0 ; i_ssat < ${#ssat_array[@]} ; i_ssat++ )); do
        ssat=${ssat_array[$i_ssat]}
        for (( i_stick = 0 ; i_stick < ${#stick_array[@]} ; i_stick++ )); do
            stick=${stick_array[$i_stick]}
            for (( i_ncl = 0 ; i_ncl < ${#ncl_array[@]} ; i_ncl++ )); do
                ncl=${ncl_array[$i_ncl]}
                for (( i_nclmass = 0 ; i_nclmass < ${#nclmass_array[@]} ; i_nclmass++ )); do
                    nclmass=${nclmass_array[$i_nclmass]}
                #for (( i_lwc = 0 ; i_lwc < ${#lwc_array[@]} ; i_lwc++ )); do
                #    lwc=${lwc_array[$i_lwc]}

                #bash doesnt care about indentation; for looping over other variables just yank the for loop and the following line

        
    echo "#################"
    echo "analyzing testcase: "$testcase "with ssat=" $ssat ", stick=" $stick ",ncl=," $ncl ",nclmass," $nclmass ",model_setup," $model_setup
    echo "#################"
    export testcase; export McSnow_geom; export model_setup; export switch_off_processes; export nclmass #make them visible for called routines
            
    export tstep="300" #start of timestep analyzed by McRadar and plotting routine in minutes
    export tstep_end="600" #end timestep analyzed by McRadar and plotting routine in minutes (defines end of averaging period)
     export av_tstep="5" #average the McSnow output in the PAMTRA adaption by ... s #do not set this lower than dtc (f.e. "5"s)
    
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
    export experiment=$($MC/run/1d "onlyname" $ssat $stick $model_setup $McSnow_geom $depos_from0mto $riming_fromdeposto  $ncl $nclmass )
    echo "analyze experiment:" $experiment
    
    if [[ "$execwhat" == *Mc1* ]] ; then #run McSnow (f.e 1d-model)
        echo "############"
        echo "run McSnow"
        echo "############"

        cd $MC/run #alternative if cheops is shut down
        echo "SETTINGS: " $ssat $stick $model_setup $McSnow_geom $depos_from0mto $riming_fromdeposto $testcase $habit $ncl $nclmass
        ./1d "fullrun" $ssat $stick $model_setup $McSnow_geom $depos_from0mto $riming_fromdeposto  $ncl $nclmass
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
        python3 plot_output.py # sofar this produces plots of the McRadar moments and spectra, aswell as the aspect ratios in spectral form
        
    fi
    
    
    

    done #for ssat
    done #for stick
    done #for ncl
    done #for nclmass
done #for testcase


