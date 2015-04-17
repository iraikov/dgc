#
#

END { 
   if (PREVIOUS != 0)
   {
        printf ("%d,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g\n",ID,DC_input_resistance,dendarea,vmin,vtau0,tau0,maximum_Vsoma,minimum_Vsoma,threshold,AP_amplitude_rel_threshold,AHP,Rel_AP_amplitude_dendrite_0,Rel_AP_amplitude_dendrite_1,Rel_AP_amplitude_dendrite_2,Rel_AP_amplitude_dendrite_3,Rel_AP_amplitude_dendrite_4,number_of_spikes,mean_FR,mean_ISI,stdev_ISI,adaptation_ISI1,adaptation_ISI2,adaptation_ISI3,adaptation_ISI4,HC_amp,HC_rise,HC_decay,BC_amp,BC_rise,BC_decay,HCC_amp,HCC_rise,HCC_decay,AA_amp,AA_rise,AA_decay,NGFC_GABAA_amp,NGFC_GABAA_rise,NGFC_GABAA_decay,NGFC_GABAB_amp,NGFC_GABAB_rise,NGFC_GABAB_decay)
    
   }
}

BEGIN { 
     IGNORECASE = 1 
     PREVIOUS = 0
}

/^DC input resistance: / \
{
   DC_input_resistance=$4
   PREVIOUS=1
}

/^vmin: / \
{
    vmin=$2
    PREVIOUS=1
}

/^vtau0: / \
{
    vtau0=$2
    PREVIOUS=1
}

/^tau0: / \
{
    tau0=$2
    PREVIOUS=1
}

/^total dendritic surface area: / \
{
    dendarea=$5
    PREVIOUS=1
}

/^maximum Vsoma: / \
{
    maximum_Vsoma=$3
    PREVIOUS=1
}

/^minimum Vsoma: / \
{
    minimum_Vsoma=$3
    PREVIOUS=1
}

/^threshold: / \
{
    threshold=$2
    PREVIOUS=1
}

/^AP amplitude relative / \
{
    AP_amplitude_rel_threshold=$6
    PREVIOUS=1
}

/^AHP relative / \
{
    AHP=$5
    PREVIOUS=1
}

/^Relative amplitude of AP in dendrite 0/ \
{
    Rel_AP_amplitude_dendrite_0=$8
    PREVIOUS=1
}

/^Relative amplitude of AP in dendrite 1/ \
{
    Rel_AP_amplitude_dendrite_1=$8
    PREVIOUS=1
}

/^Relative amplitude of AP in dendrite 2/ \
{
    Rel_AP_amplitude_dendrite_2=$8
    PREVIOUS=1
}

/^Relative amplitude of AP in dendrite 3/ \
{
    Rel_AP_amplitude_dendrite_3=$8
    PREVIOUS=1
}

/^Relative amplitude of AP in dendrite 4/ \
{
    Rel_AP_amplitude_dendrite_4=$8
    PREVIOUS=1
}

/number of spikes/ \
{
    number_of_spikes=$5
    PREVIOUS=1
}

/number of spikes/ \
{
    number_of_spikes=$5
    PREVIOUS=1
}

/FR mean/ \
{
    mean_FR=$4
    PREVIOUS=1
}

/ISI mean/ \
{
    mean_ISI=$4
    PREVIOUS=1
}

/ISI stdev/ \
{
    stdev_ISI=$4
    PREVIOUS=1
}

/ISI adaptation 1/ \
{
    adaptation_ISI1=$5
    PREVIOUS=1
}

/ISI adaptation 1/ \
{
    adaptation_ISI1=$5
    PREVIOUS=1
}

/ISI adaptation 2/ \
{
    adaptation_ISI2=$5
    PREVIOUS=1
}

/ISI adaptation 3/ \
{
    adaptation_ISI3=$5
    PREVIOUS=1
}

/ISI adaptation 4/ \
{
    adaptation_ISI4=$5
    PREVIOUS=1
}

/HC-GC synapses/ \
{
    HC_amp=$4
    HC_rise=$8
    HC_decay=$12
    PREVIOUS=1
}

/BC-GC synapses/ \
{
    BC_amp=$4
    BC_rise=$8
    BC_decay=$12
    PREVIOUS=1
}

/HCC-GC synapses/ \
{
    HCC_amp=$4
    HCC_rise=$8
    HCC_decay=$12
    PREVIOUS=1
}

/AA-GC synapses/ \
{
    AA_amp=$4
    AA_rise=$8
    AA_decay=$12
    PREVIOUS=1
}

/NGFC-GC GABAA synapses/ \
{
    NGFC_GABAA_amp=$5
    NGFC_GABAA_rise=$9
    NGFC_GABAA_decay=$13
    PREVIOUS=1
}

/NGFC-GC GABAB synapses/ \
{
    NGFC_GABAB_amp=$5
    NGFC_GABAB_rise=$9
    NGFC_GABAB_decay=$13
    PREVIOUS=1
}
