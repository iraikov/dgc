#! /usr/bin/octave -qf
args = argv ();

DGC_results = load(args{1});

h = figure(1);

subplot(3, 3, 1)
hist(DGC_results(:,2),50);
title("Input resistance [MOhm]");


subplot(3, 3, 2)
hist(DGC_results(:,6),50);
title("Membrane time constant [ms]");

subplot(3, 3, 3)
hist(DGC_results(:,8),50);
title("Threshold [ms]");

subplot(3, 3, 4)
hist(DGC_results(:,20),50);
title("ISI adaptation 2");

print (h, "DGC_results.pdf", "-dpdf")

