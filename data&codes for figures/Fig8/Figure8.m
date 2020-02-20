close all
clear all
load('../pilC/sim_frac_pilc');
all=[pilc1;pilc2;pilc3];
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
load('pilC_wtnp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_pilcnp');
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
load('pilC_wtnpb');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_pilcnpb');
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
load('pilC_wtpdur');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_pilcpdur');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'k');
hold off
xlim([0,5]);
ylim([0,0.5]);
%% Panel D
figure
xi = (1:length(mu)) / 2 / 60-1.5;
plot(xi,mu,'b');
hold on
xi = (1:length(wt)) / 2 / 60-1.5;
plot(xi,wt,'r');
load('pilC_wtpsp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'g');
load('WT_pilcpsp');
xi = (1:length(m)) / 2 / 60-1.5;
plot(xi,m,'k');
hold off
xlim([0,5]);
ylim([0,0.5]);