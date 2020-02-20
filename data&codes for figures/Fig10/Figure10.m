close all
clear all
load('../csgA/sim_frac_csga');
all=[csga1;csga2;csga3];
mu=mean(all);
load('../WT/sim_frac_wt');
all=[WT1;WT2;WT3];
wt=mean(all);
%% Panel A
figure
xi = (1:length(mu)) / 2 / 60-1.5;
plot(xi,mu,'b');
hold on
xi = (1:length(wt)) / 2 / 60-1.5;
plot(xi,wt,'r');
load('csgA_wtnp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_csganp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'k');
hold off
xlim([0,5]);
ylim([0,0.5]);
%% Panel B
figure
xi = (1:length(mu)) / 2 / 60-1.5;
plot(xi,mu,'b');
hold on
xi = (1:length(wt)) / 2 / 60-1.5;
plot(xi,wt,'r');
load('csgA_wtnpb');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_csganpb');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'k');
hold off
xlim([0,5]);
ylim([0,0.5]);
%% Panel C
figure
xi = (1:length(mu)) / 2 / 60-1.5;
plot(xi,mu,'b');
hold on
xi = (1:length(wt)) / 2 / 60-1.5;
plot(xi,wt,'r');
load('csgA_wtpdur');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_csgapdur');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'k');
hold off
xlim([0,5]);
ylim([0,0.56]);
%% Panel D
figure
xi = (1:length(mu)) / 2 / 60-1.5;
plot(xi,mu,'b');
hold on
xi = (1:length(wt)) / 2 / 60-1.5;
plot(xi,wt,'r');
load('csgA_wtpsp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_csgapsp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'k');
hold off
xlim([0,5]);
ylim([0,0.5]);