
DGC_results = load("DGC_results_table.dat" );

figure(1);
hist(DGC_results(:,2));
title("Input resistance [MOhm]");
figure(2);
hist(DGC_results(:,6));
title("Membrane time constant [ms]");
figure(3);
hist(DGC_results(:,9));
