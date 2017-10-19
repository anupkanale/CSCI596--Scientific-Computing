clear
close all;

thrds = [1,2,4,8];
tn = [4.574044e+01, 2.979513e+01, 2.975431e+01, 2.234813e+01 ];

Sp = tn./tn(1);
Eff = Sp./thrds;

set(figure(), 'position', [1000 600 800 600], 'color', 'w');
plot(thrds, Eff, 'bo-', 'linewidth', 2);
hold on;
plot(0:10, ones(11,1), '-k','linewidth',2);

xlabel('# of threads', 'fontsize', 18, 'interpreter', 'latex');
ylabel('Efficiency', 'fontsize', 18, 'interpreter', 'latex');
gca('fontsize',18)