%dir = 'Test4-Galileo-Deimos';
% Betts no shadow and betts shadow representations proble
load(['../',dir,'/output.out']);
% alpha =load([dir,'/alpha.out']);
% beta =load('beta.out');
% load('at.out');
% L =load('L.out');
% Lends =load('Lends.out');
% at =load('at.out');
% T =load('T.out');
signoid = load(['../',dir,'/signoid.out']);

t     = 365*output(:,1);
p     = output(:,2);
f     = output(:,3);
g     = output(:,4);
h     = output(:,5);
k     = output(:,6);


signoid = signoid>0.3;
sum(t(diff(signoid)==1)- t(diff(signoid)==-1))

oe = mee2oe([p,f, g h, k, zeros(size(k))]);

fig = figure(1);
plot(t,oe(:,1))
legend({'Step2'},'FontSize',12,'Interpreter','latex')
title('Semi-major axis')
xlabel('Time (days)','FontSize',12)
savefig(['../',dir,'/Semimajor_axis.fig'])
saveas(fig,['../',dir,'/Semimajor_axis.png']);
% % 
fig = figure(2);
plot(t,oe(:,2))
legend({'Step2'},'FontSize',12,'Interpreter','latex')
title('Eccentricity')
xlabel('Time (days)','FontSize',12)
savefig(['../',dir,'/Eccentricity.fig'])
saveas(fig,['../',dir,'/Eccentricity.png']);
% % 
fig = figure(3);
plot(t,oe(:,3)*180/pi)
legend({'Step2'},'FontSize',12,'Interpreter','latex')
title('Inclination')
xlabel('Time (days)','FontSize',12)
savefig(['../',dir,'/Inclination.fig'])
saveas(fig,['../',dir,'/Inclination.png']);

fig = figure(4);
plot(t,signoid)
legend({'Step2'},'FontSize',12,'Interpreter','latex')
title('Eclipse')
xlabel('Time (days)','FontSize',12)
savefig(['../',dir,'/signoid.fig'])
saveas(fig,['../',dir,'/signoid.png']);
% 
fig = figure(5);
plot(t,oe(:,4)*180/pi)
legend({'Step2'},'FontSize',12,'Interpreter','latex')
title('Eclipse')
xlabel('Time (days)','FontSize',12)
savefig(['../',dir,'/RAAN.fig'])
saveas(fig,['../',dir,'/RAAN.png']);
