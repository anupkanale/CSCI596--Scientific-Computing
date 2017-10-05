clear; close all

set(figure(), 'position', [0 0 866 600], 'color', 'w');

dat1 = dlmread('pv_1.dat');
nbins1 = dat1(:,1);
vel1 = dat1(:,2);
plot(nbins1, vel1, 'linewidth', 2);
hold on

dat2 = dlmread('pv_2.dat');
nbins2 = dat2(:,1);
vel2 = dat2(:,2);
plot(nbins2, vel2, 'linewidth', 2);
hold on

dat3 = dlmread('pv_3.dat');
nbins3 = dat3(:,1);
vel3 = dat3(:,2);
plot(nbins3, vel3, 'linewidth', 2);

xlim([0 3]);
% ylim();
xlabel('v', 'fontsize', 18, 'interpreter', 'latex'); 
ylabel('P(v)', 'fontsize', 18, 'interpreter', 'latex');
set(gca, 'fontsize', 18)