#! /usr/bin/octave -qf
args = argv ();

DGC_results = load(args{1});

input_resistance=DGC_results(:,2);
membrane_tau=DGC_results(:,6);
spike_threshold=DGC_results(:,9);
spike_amplitude=DGC_results(:,10);
spike_ahp=DGC_results(:,11);
spike_amplitude_dend1=DGC_results(:,12);
spike_amplitude_dend2=DGC_results(:,13);
spike_amplitude_dend3=DGC_results(:,14);
number_of_spikes=DGC_results(:,15);
mean_firing_rate=DGC_results(:,16);
mean_isi=DGC_results(:,17);
isi_adaptation2=DGC_results(:,20);

h = figure(1);

subplot(3, 2, 1)
hist(input_resistance,50);
title(sprintf("Input resistance [MOhm]; mean = %g",mean(input_resistance)));

subplot(3, 2, 2)
hist(membrane_tau,50);
title(sprintf("Membrane time constant [ms]; mean = %g",mean(membrane_tau)));

subplot(3, 2, 3)
hist(spike_amplitude,50);
title(sprintf("Rel. AP amplitude [mV]; mean = %g",mean(spike_amplitude)));

subplot(3, 2, 4)
hist(spike_threshold,50);
title(sprintf("AP threshold [mV]; mean = %g",mean(spike_threshold)));

subplot(3, 2, 5)
hist(spike_ahp,50);
title(sprintf("Fast AHP [mV]; mean = %g",mean(spike_ahp)));

subplot(3, 2, 6)
hist(isi_adaptation2,50);
title(sprintf("ISI adaptation 2; mean = %g",mean(isi_adaptation2)));

print (h, "DGC_results1.pdf", "-dpdf")

h = figure(2);

s = subplot(2, 2, 1)
spike_amplitude_dend1_mean = mean(spike_amplitude_dend1);
spike_amplitude_dend1_stdev = std(spike_amplitude_dend1);
spike_amplitude_dend2_mean = mean(spike_amplitude_dend2);
spike_amplitude_dend2_stdev = std(spike_amplitude_dend2);
spike_amplitude_dend3_mean = mean(spike_amplitude_dend3);
spike_amplitude_dend3_stdev = std(spike_amplitude_dend3);
spike_amplitude_dend_means = [spike_amplitude_dend1_mean,
                              spike_amplitude_dend2_mean,
                              spike_amplitude_dend3_mean];
spike_amplitude_dend_stds = [spike_amplitude_dend1_stdev,
                             spike_amplitude_dend2_stdev,
                             spike_amplitude_dend3_stdev];
errorbar([80,160,240],spike_amplitude_dend_means, ...
         spike_amplitude_dend_stds);
axis([0 300 -0.25 1.25]);
set(s,'XTick',0:80:240);
set(s,'YTick',0:0.2:1);
title("Amplitude of dendritic AP rel. to soma [mV]");

subplot(2, 2, 2)
hist(number_of_spikes,20);
title(sprintf("Number of spikes; mean = %g",mean(number_of_spikes)));

subplot(2, 2, 3)
hist(mean_firing_rate,50);
title(sprintf("Mean firing rate [Hz]"));

subplot(2, 2, 4)
hist(mean_isi,50);
title(sprintf("Mean ISI [ms]"));


print (h, "DGC_results2.pdf", "-dpdf")
