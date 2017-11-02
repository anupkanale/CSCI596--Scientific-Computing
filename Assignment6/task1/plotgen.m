clear; close all

set(figure(), 'position', [2500 500 866 600], 'color', 'w');

dat1 = dlmread('pdf_cpu.d');
nbins_cpu = dat1(:,1);
pdf_cpu = dat1(:,2);

dat2 = dlmread('pdf_gpu.d');
nbins_gpu = dat2(:,1);
pdf_gpu = dat2(:,2);

plot(nbins_cpu, pdf_cpu, 'b-', 'linewidth', 3);
hold on
plot(nbins_gpu, pdf_gpu, 'r--', 'linewidth', 2);

% ylim();
xlabel('r', 'fontsize', 18, 'interpreter', 'latex'); 
ylabel('g(r)', 'fontsize', 18, 'interpreter', 'latex');
set(gca, 'fontsize', 18)

legend('cpu', 'gpu')