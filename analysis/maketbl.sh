#!/bin/bash

##
## This script collects the results from all DGC test result files and
## assembles them in a single table.
##
##

tblpath=$1

nworkers=20
ncells=1000

path_prefix="$2"

echo "## gid,DC_input_resistance,dendarea,vmin,vtau0,tau0,maximum_Vsoma,minimum_Vsoma,threshold,AP_amplitude_rel_threshold,AHP,Rel_AP_amplitude_dendrite_0,Rel_AP_amplitude_dendrite_1,Rel_AP_amplitude_dendrite_2,Rel_AP_amplitude_dendrite_3,Rel_AP_amplitude_dendrite_4,number_of_spikes,mean_FR,mean_ISI,stdev_ISI,adaptation_ISI1,adaptation_ISI2,adaptation_ISI3,adaptation_ISI3,adaptation_ISI4,HC_amp,HC_rise,HC_decay,BC_amp,BC_rise,BC_decay,HCC_amp,HCC_rise,HCC_decay,AA_amp,AA_rise,AA_decay,NGFC_GABAA_amp,NGFC_GABAA_rise,NGFC_GABAA_decay" > $tblpath

for (( i = 0; i < ncells; i++ )) do

    id=`printf '%06d' $i`
    path="$path_prefix"
    passive_results_path="$path/DGC_passive_results_$id.dat"
    single_ap_results_path="$path/DGC_single_ap_results_$id.dat"
    threshold_results_path="$path/DGC_threshold_results_$id.dat"
    spikes_path="$path/DGC_spikes_$id.dat"
    synapse_path="$path/DGC_synapse_$id.dat"

    DC_input_resistance=`grep "DC input resistance" $passive_results_path | cut -f4 -d' '`
    vmin=`grep "vmin" $passive_results_path | cut -f2 -d' '`
    vtau0=`grep "vtau0" $passive_results_path | cut -f2 -d' '`
    tau0=`grep "^tau0" $passive_results_path | cut -f2 -d' '`
    dendarea=`grep "dendritic surface area" $passive_results_path | cut -f5 -d' '`

    maximum_Vsoma=`grep "maximum Vsoma" $single_ap_results_path | cut -f3 -d' '`
    minimum_Vsoma=`grep "minimum Vsoma" $single_ap_results_path | cut -f3 -d' '`
    threshold=`grep "^threshold:" $threshold_results_path | cut -f2 -d' '`
    AP_amplitude_rel_threshold=`grep "AP amplitude relative" $threshold_results_path | cut -f6 -d' '`
    AHP=`grep "AHP relative" $threshold_results_path | cut -f5 -d' '`
    Rel_AP_amplitude_dendrite_0=`grep "Relative amplitude of AP in dendrite 0" $single_ap_results_path | cut -f8 -d' '`
    Rel_AP_amplitude_dendrite_1=`grep "Relative amplitude of AP in dendrite 1" $single_ap_results_path | cut -f8 -d' '`
    Rel_AP_amplitude_dendrite_2=`grep "Relative amplitude of AP in dendrite 2" $single_ap_results_path | cut -f8 -d' '`
    Rel_AP_amplitude_dendrite_3=`grep "Relative amplitude of AP in dendrite 3" $single_ap_results_path | cut -f8 -d' '`
    Rel_AP_amplitude_dendrite_4=`grep "Relative amplitude of AP in dendrite 4" $single_ap_results_path | cut -f8 -d' '`

    number_of_spikes=`grep "number of spikes" $spikes_path | cut -f5 -d' '`
    mean_FR=`grep "FR mean" $spikes_path | cut -f4 -d' '`
    mean_ISI=`grep "ISI mean" $spikes_path | cut -f4 -d' '`
    stdev_ISI=`grep "ISI stdev" $spikes_path | cut -f4 -d' '`
    adaptation_ISI1=`grep "ISI adaptation 1" $spikes_path | cut -f5 -d' '`
    adaptation_ISI2=`grep "ISI adaptation 2" $spikes_path | cut -f5 -d' '`
    adaptation_ISI3=`grep "ISI adaptation 3" $spikes_path | cut -f5 -d' '`
    adaptation_ISI4=`grep "ISI adaptation 4" $spikes_path | cut -f5 -d' '`
    
    HC_amp=`grep "HC-GC synapses" $synapse_path | cut -f6 -d' '`
    HC_rise=`grep "HC-GC synapses" $synapse_path | cut -f11 -d' '`
    HC_decay=`grep "HC-GC synapses" $synapse_path | cut -f16 -d' '`
    BC_amp=`grep "BC-GC synapses" $synapse_path | cut -f6 -d' '`
    BC_rise=`grep "BC-GC synapses" $synapse_path | cut -f11 -d' '`
    BC_decay=`grep "BC-GC synapses" $synapse_path | cut -f16 -d' '`
    HCC_amp=`grep "HCC-GC synapses" $synapse_path | cut -f6 -d' '`
    HCC_rise=`grep "HCC-GC synapses" $synapse_path | cut -f11 -d' '`
    HCC_decay=`grep "HCC-GC synapses" $synapse_path | cut -f16 -d' '`
    AA_amp=`grep "AA-GC synapses" $synapse_path | cut -f6 -d' '`
    AA_rise=`grep "AA-GC synapses" $synapse_path | cut -f11 -d' '`
    AA_decay=`grep "AA-GC synapses" $synapse_path | cut -f16 -d' '`
    NGFC_GABAA_amp=`grep "NGFC-GC GABAA synapses" $synapse_path | cut -f7 -d' '`
    NGFC_GABAA_rise=`grep "NGFC-GC GABAA synapses" $synapse_path | cut -f12 -d' '`
    NGFC_GABAA_decay=`grep "NGFC-GC GABAA synapses" $synapse_path | cut -f17 -d' '`


    echo $i,$DC_input_resistance,$dendarea,$vmin,$vtau0,$tau0,$maximum_Vsoma,$minimum_Vsoma,$threshold,$AP_amplitude_rel_threshold,$AHP,$Rel_AP_amplitude_dendrite_0,$Rel_AP_amplitude_dendrite_1,$Rel_AP_amplitude_dendrite_2,$Rel_AP_amplitude_dendrite_3,$Rel_AP_amplitude_dendrite_4,$number_of_spikes,$mean_FR,$mean_ISI,$stdev_ISI,$adaptation_ISI1,$adaptation_ISI2,$adaptation_ISI3,$adaptation_ISI4,$HC_amp,$HC_rise,$HC_decay,$BC_amp,$BC_rise,$BC_decay,$HCC_amp,$HCC_rise,$HCC_decay,$AA_amp,$AA_rise,$AA_decay,$NGFC_GABAA_amp,$NGFC_GABAA_rise,$NGFC_GABAA_decay >> $tblpath
    
done

