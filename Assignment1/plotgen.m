% Script to generate plot and fit to a curve
clear
close all;

nAtoms = 4* (4:10).^3;
timelmd = [0.017377, 0.022644, 0.031774, 0.06952,...
    0.084479, 0.108944, 0.180908];

p = polyfit(nAtoms, timelmd, 1);
timeFit = polyval(p, nAtoms);

set(figure, 'position', [1000 300 800 600])
plot(nAtoms, timelmd, 'ro', 'linewidth', 2)
hold on;
plot(nAtoms, timeFit, 'b-', 'linewidth', 2)
txt = strcat(num2str(p(1)), 'x+', num2str(p(2)));
text(nAtoms(5)+100, timelmd(5), txt);
xlabel('nAtoms', 'fontsize', 14)
ylabel('Time', 'fontsize', 14)

% MD plot
timemd = [0.012428, 0.037953, 0.105062, 0.251538,...
    0.530956, 1.068831, 1.9434];
p = polyfit(nAtoms, timemd, 2);
timeFit = polyval(p, nAtoms);

set(figure, 'position', [1000 00 800 600])
plot(nAtoms, timemd, 'ro', 'linewidth', 2)
hold on;
plot(nAtoms, timeFit, 'b-', 'linewidth', 2)
txt = strcat(num2str(p(1)), 'x^2+', num2str(p(2)), 'x+', num2str(p(3)));
text(nAtoms(5)+100, timemd(5), txt);
xlabel('nAtoms', 'fontsize', 14)
ylabel('Time', 'fontsize', 14)