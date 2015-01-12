
DGC_results = load("DGC_results_table.dat" );

figure(1);
hist(DGC_results(:,2),50);
title("Input resistance [MOhm]");
figure(2);
hist(DGC_results(:,6),50);
title("Membrane time constant [ms]");
figure(3);
hist(DGC_results(:,8),50);
title("Threshold [ms]");
figure(4);
hist(DGC_results(:,20),50);
title("ISI adaptation 3");
