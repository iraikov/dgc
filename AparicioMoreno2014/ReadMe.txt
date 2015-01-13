This NEURON model of a dentate granule cell is taken from the 
following paper:

Pedro Mateos-Aparicio, Ricardo Murphy and Johan F. Storm (2014). 
Complementary functions of SK and Kv7/M potassium channels in 
excitability control and synaptic integration in rat hippocampal 
dentate granule cells. Journal of Physiology 592, 669-693.

It is based on the model of Aradi & Holmes (1999; Journal of 
Computational Neuroscience 6, 215-235) which uses an idealized 
morphology (DGC_Morphology.hoc). The model was used to help 
understand the contribution of M and SK channels to the medium 
afterhyperpolarization (mAHP) following one or seven spikes, as 
well as the contribution of M channels to the slow 
afterhyperpolarization (sAHP). We found that SK channels are the 
main determinants of the mAHP, in contrast to CA1 pyramidal 
cells where the mAHP is primarily caused by the opening of M 
channels. The model reproduced these experimental results, but 
we were unable to reproduce the effects of the M-channel blocker 
XE991 on the sAHP. It is suggested that either the XE991-
sensitive component of the sAHP is not due to M channels, or that 
when contributing to the sAHP, these channels operate in a mode 
different from that associated with the mAHP.

After compiling the .mod files, run the program mosinit.hoc. Two 
panels will appear. Click on "Figure 6A" or "Figure 6B" in the 
smaller panel to reproduce the corresponding figures in the paper 
(details in DGC_Figures.hoc). The larger panel allows you to 
inject a test pulse into the soma and change various model 
parameters (see DGC_Biophysics.hoc, DGC_Parameters.hoc, 
DGC_SetUp.hoc, the .mod files and our paper's supplementary 
information for details).

