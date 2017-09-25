clear;
close all;

nProcs = [1,2,4,8];

%% Fixed problem-size scaling
execTime1 = [2.62855, 1.3267, 6.9056e-1, 4.0104e-1];

eff1 = zeros(length(nProcs),1);
for ii=1:length(nProcs)
    eff1(ii) = execTime1(1)/nProcs(ii)/execTime1(ii);
end

plot(nProcs,eff1, 'ro--', 'linewidth', 2);
hold on;
plot(nProcs,ones(length(eff1),1), 'k-', 'linewidth', 1.5);
ylim([0 1.1]);
xlabel('number of processors', 'Interpreter','latex');
ylabel('parallel efficiency', 'Interpreter','latex');
set(gca, 'fontsize', 20)
set(gcf,'color','w');

%% Isogranular scaling
execTime2 = [2.62855,2.57204,2.59,1];

eff2 = zeros(length(nProcs),1);
for ii=1:length(nProcs)
    eff2(ii) = execTime2(1)/execTime2(ii);
end

hold off;
figure()
plot(nProcs,eff2, 'ro--', 'linewidth', 2);
hold on;
plot(nProcs,ones(length(eff2),1), 'k-', 'linewidth', 1.5);
ylim([0 1.1]);
xlabel('number of processors', 'Interpreter','latex');
ylabel('parallel efficiency', 'Interpreter','latex');
set(gca, 'fontsize', 20)
set(gcf,'color','w');