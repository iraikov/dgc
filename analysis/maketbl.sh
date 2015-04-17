#!/bin/bash

##
## This script collects the results from all DGC test result files and
## assembles them in a single table.
##
##

exec_path=`dirname $0`

tblpath=$1

nworkers=20
ncells=1000

path_prefix="$2"

echo "%% gid,DC_input_resistance,dendarea,vmin,vtau0,tau0,maximum_Vsoma,minimum_Vsoma,threshold,AP_amplitude_rel_threshold,AHP,Rel_AP_amplitude_dendrite_0,Rel_AP_amplitude_dendrite_1,Rel_AP_amplitude_dendrite_2,Rel_AP_amplitude_dendrite_3,Rel_AP_amplitude_dendrite_4,number_of_spikes,mean_FR,mean_ISI,stdev_ISI,adaptation_ISI1,adaptation_ISI2,adaptation_ISI3,adaptation_ISI3,adaptation_ISI4,HC_amp,HC_rise,HC_decay,BC_amp,BC_rise,BC_decay,HCC_amp,HCC_rise,HCC_decay,AA_amp,AA_rise,AA_decay,NGFC_GABAA_amp,NGFC_GABAA_rise,NGFC_GABAA_decay,NGFC_GABAB_amp,NGFC_GABAB_rise,NGFC_GABAB_decay" > $tblpath

for (( i = 0; i < ncells; i++ )) do

    id=`printf '%06d' $i`
    path="$path_prefix"
    passive_results_path="$path/DGC_passive_results_$id.dat"
    single_ap_results_path="$path/DGC_single_ap_results_$id.dat"
    threshold_results_path="$path/DGC_threshold_results_$id.dat"
    spikes_path="$path/DGC_spikes_$id.dat"
    synapse_path="$path/DGC_synapse_$id.dat"

    awk -v "ID=$i" -f $exec_path/csv-filter.awk $passive_results_path $single_ap_results_path $threshold_results_path $spikes_path $synapse_path >> $tblpath

done

