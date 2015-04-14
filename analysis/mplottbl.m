DGC_results = load('/home/igr/src/model/DGC/Mateos-Aparicio2014/results/forest/DGC_forest_110_test_results_passive_201504141201.dat');

input_resistance=DGC_results(:,2);
membrane_tau=DGC_results(:,6);
spike_threshold=DGC_results(:,9);
spike_amplitude=DGC_results(:,10);
spike_ahp=DGC_results(:,11);
rel_amplitude_dend1=DGC_results(:,12);
rel_amplitude_dend2=DGC_results(:,13);
rel_amplitude_dend3=DGC_results(:,14);
rel_amplitude_dend4=DGC_results(:,15);
rel_amplitude_dend5=DGC_results(:,16);
number_of_spikes=DGC_results(:,17);
mean_firing_rate=DGC_results(:,18);
mean_isi=DGC_results(:,19);
isi_adaptation3=DGC_results(:,23);
isi_adaptation4=DGC_results(:,24);

h = figure(1);

subplot(3, 2, 1)
hist(input_resistance,50);
title(['Input resistance [MOhm]', sprintf('mean = %.2f std = %.2f', mean(input_resistance),std(input_resistance))]);

subplot(3, 2, 2)
hist(membrane_tau,50);
title(['Membrane time constant [ms]', sprintf('mean = %.2f std = %.2f', mean(membrane_tau),std(membrane_tau))]);

subplot(3, 2, 3)
hist(spike_amplitude,50);
title(['Rel. AP amplitude [mV]', sprintf('mean = %.2f std = %.2f', mean(spike_amplitude),std(spike_amplitude))]);

subplot(3, 2, 4)
hist(spike_threshold,50);
title(sprintf('AP threshold [mV]; mean = %.2f std = %.2f', mean(spike_threshold),std(spike_threshold)));

subplot(3, 2, 5)
hist(spike_ahp,50);
title(sprintf('Fast AHP [mV]; mean = %.2f std = %.2f',mean(spike_ahp),std(spike_ahp)));
 
print (h, 'DGC_results1.pdf', '-dpdf')

h = figure(2);

subplot(2, 2, 1)
hist(isi_adaptation4,50);
title(sprintf('ISI adaptation 4; mean = %.2f std = %.2f', mean(isi_adaptation4),std(isi_adaptation4)));

subplot(2, 2, 2)
hist(number_of_spikes,20);
title(sprintf('Number of spikes; mean = %.2f',mean(number_of_spikes)));

s = subplot(2, 2, 3)
rel_amplitude_dend1_mean = mean(rel_amplitude_dend1);
rel_amplitude_dend1_stdev = std(rel_amplitude_dend1);
rel_amplitude_dend2_mean = mean(rel_amplitude_dend2);
rel_amplitude_dend2_stdev = std(rel_amplitude_dend2);
rel_amplitude_dend3_mean = mean(rel_amplitude_dend3);
rel_amplitude_dend3_stdev = std(rel_amplitude_dend3);
rel_amplitude_dend4_mean = mean(rel_amplitude_dend4);
rel_amplitude_dend4_stdev = std(rel_amplitude_dend4);
rel_amplitude_dend5_mean = mean(rel_amplitude_dend5);
rel_amplitude_dend5_stdev = std(rel_amplitude_dend5);


rel_amplitude_dend_means = [rel_amplitude_dend1_mean, rel_amplitude_dend2_mean, rel_amplitude_dend3_mean, rel_amplitude_dend4_mean, rel_amplitude_dend5_mean,];
rel_amplitude_dend_stds = [rel_amplitude_dend1_stdev, rel_amplitude_dend2_stdev, rel_amplitude_dend3_stdev, rel_amplitude_dend4_stdev, rel_amplitude_dend5_stdev];
errorbar([50,100,150,200,250],rel_amplitude_dend_means, rel_amplitude_dend_stds);
axis([0 300 0 1]);
set(s,'XTick',0:50:250);
set(s,'YTick',0:0.2:1);
title('Amplitude of dendritic AP rel. to soma');

print (h, 'DGC_results2.pdf', '-dpdf')

% subplot(2, 3, 1)
% hist(mean_firing_rate,50);
% title(sprintf('Mean firing rate [Hz]'));

% subplot(2, 2, 4)
% hist(mean_isi,50);
% title(sprintf('Mean ISI [ms]'));

h = figure(3);

HC_amp=DGC_results(:,25);
HC_rise=DGC_results(:,26);
HC_decay=DGC_results(:,27);
BC_amp=DGC_results(:,28);
BC_rise=DGC_results(:,29);
BC_decay=DGC_results(:,30);
HCC_amp=DGC_results(:,31);
HCC_rise=DGC_results(:,32);
HCC_decay=DGC_results(:,33);
AA_amp=DGC_results(:,34);
AA_rise=DGC_results(:,35);
AA_decay=DGC_results(:,36);
NGFC_GABAA_amp=DGC_results(:,36);
NGFC_GABAA_rise=DGC_results(:,37);
NGFC_GABAA_decay=DGC_results(:,38);

subplot(3, 3, 1)
hist(HC_amp,50);
title(['HIPP syn amp. ', sprintf('mean = %.2f std = %.2f', mean(HC_amp),std(HC_amp))]);

subplot(3, 3, 2)
hist(HC_rise,50);
title(['HIPP syn rise ', sprintf('mean = %.2f std = %.2f', mean(HC_rise),std(HC_rise))]);

subplot(3, 3, 3)
hist(HC_decay,50);
title(['HIPP syn decay ', sprintf('mean = %.2f std = %.2f', mean(HC_decay),std(HC_decay))]);


subplot(3, 3, 4)
hist(BC_amp,50);
title(['PVBC syn amp. ', sprintf('mean = %.2f std = %.2f', mean(BC_amp),std(BC_amp))]);

subplot(3, 3, 5)
hist(BC_rise,50);
title(['PVBC syn rise ', sprintf('mean = %.2f std = %.2f', mean(BC_rise),std(BC_rise))]);

subplot(3, 3, 6)
hist(BC_decay,50);
title(['PVBC syn decay ', sprintf('mean = %.2f std = %.2f', mean(BC_decay),std(BC_decay))]);


subplot(3, 3, 7)
hist(AA_amp,50);
title(['AA syn amp. ', sprintf('mean = %.2f std = %.2f', mean(AA_amp),std(AA_amp))]);

subplot(3, 3, 8)
hist(AA_rise,50);
title(['AA syn rise ', sprintf('mean = %.2f std = %.2f', mean(AA_rise),std(AA_rise))]);

subplot(3, 3, 9)
hist(AA_decay,50);
title(['AA syn decay ', sprintf('mean = %.2f std = %.2f', mean(AA_decay),std(AA_decay))]);


print (h, 'DGC_results3.pdf', '-dpdf')

system('pdftk DGC_results1.pdf DGC_results2.pdf DGC_results3.pdf cat output DGC_results.pdf');
system('rm DGC_results1.pdf DGC_results2.pdf DGC_results3.pdf');
