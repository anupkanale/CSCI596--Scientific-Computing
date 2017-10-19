clear;
% close all;

nbin = 1e7; % 10 mil points used in the c-program
nProcs = [1,2,4,8];

%% Fixed problem-size scaling
% execTime1 = [2.32e-1, 9.99e-2, 5.02e-2, 3.06e-2];%>1
% execTime1 = [2.23e-1, 1.30e-1, 6.68e-2, 3.56e-2]; % good
% execTime1 = [2.23e-1, 1.17e-1, 6.56e-2, 3.26e-2]; % correct but not good
% execTime1 = [2.23e-1, 1.06e-1, 5.09e-2, 3.25e-2]; %>1
execTime1 = [3.24e-1, 1.95e-1, 1.25e-1, 7.94e-2]; % best
% execTime1 = [2.386839e-01, 1.158540e-01, 6.487799e-02, 3.309393e-02]; %>1
% execTime1 = [1.582389e-01, 8.969712e-02, 5.591893e-02, 4.283094e-02]; @
% hpc1118 node

eff1 = zeros(length(nProcs),1);
for ii=1:length(nProcs)
    eff1(ii) = execTime1(1)/(nProcs(ii)*execTime1(ii));
end

set(figure, 'position', [0 0 1000 800], 'color', 'w')
plot(nProcs,eff1, 'ro--', 'linewidth', 2);
hold on;
plot(nProcs,ones(length(eff1),1), 'k-', 'linewidth', 1.5);
ylim([0 1.1]);
xlabel('number of processors', 'Interpreter','latex');
ylabel('parallel efficiency', 'Interpreter','latex');
set(gca, 'fontsize', 30);

%% Isogranular scaling
execTime2 = [1.632669e-01, 1.760671e-01, 1.816571e-01, 1.894770e-01];
% execTime2 = [1.633959e-01, 1.671591e-01, 2.110820e-01, 1.981881e-01]; %
% correct, not great though

eff2 = zeros(length(nProcs),1);
for ii=1:length(nProcs)
    eff2(ii) = execTime2(1)/execTime2(ii);
end

set(figure, 'position', [0 0 800 600], 'color', 'w')
plot(nProcs,eff2, 'ro--', 'linewidth', 2);
hold on;
plot(nProcs,ones(length(eff2),1), 'k-', 'linewidth', 1.5);
ylim([0 1.1]);
xlabel('number of processors', 'Interpreter','latex');
ylabel('parallel efficiency', 'Interpreter','latex');
set(gca, 'fontsize', 30)