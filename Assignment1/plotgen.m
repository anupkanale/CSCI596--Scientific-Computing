% Script to generate plot and fit to a curve
clear
close all;

nAtoms = log(4* (4:10).^3);
timelmd = log([0.017377, 0.022644, 0.031774, 0.06952,...
    0.084479, 0.108944, 0.180908]);

plmd = polyfit(nAtoms, timelmd, 1);
timeFit = polyval(plmd, nAtoms);

set(figure, 'position', [1000 300 800 600], 'color', 'w')
h1 = plot(nAtoms, timelmd, 'ko', 'linewidth', 2);
hold on;
h2 = plot(nAtoms, timeFit, 'b-', 'linewidth', 2);

% MD plot
timemd = log([0.012428, 0.037953, 0.105062, 0.251538,...
    0.530956, 1.068831, 1.9434]);
pmd = polyfit(nAtoms, timemd, 1);
timeFit = polyval(pmd, nAtoms);

hold on;
h3 = plot(nAtoms, timemd, 'k^', 'linewidth', 2);
hold on;
h4 = plot(nAtoms, timeFit, 'r-', 'linewidth', 2);

txt = strcat(num2str(plmd(1)), 'ln(nAtom)+', num2str(plmd(2)));
text(nAtoms(4)+0.1, timelmd(4)-0.2, txt);
txt = strcat(num2str(pmd(1)), 'ln(nAtom)+', num2str(pmd(2)));
text(nAtoms(2)+0.3, timemd(4)+0.5, txt);

xlabel('ln(nAtoms)', 'fontsize', 14)
ylabel('ln(time)', 'fontsize', 14)
legn = legend([h2,h4], 'Linked-list MD', 'MD', 'Location', 'best');

set(gca, 'fontsize', 14)