close all
clear all
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');
org_color=[0 0 1];
rd_color=[0 0 0];
colormap = [org_color; rd_color];

%% Panel A
load('../pilC/sim_frac_pilc');
all=[pilc1;pilc2;pilc3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
hold on 
load('pilC_npbrddis');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.4]);
hold off
%% Panel B
figure
load('../csgA/sim_frac_csga');
all=[csga1;csga2;csga3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('csgA_npbrddis');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.5]);
%% Panel C
figure
load('../WT/sim_frac_wt');
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('WT_npbrddis');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.5]);
%% Panel D
figure
load('../pilC/sim_frac_pilc');
all=[pilc1;pilc2;pilc3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
hold on 
load('pilC_nprobrddis');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.4]);
hold off
%% Panel E
figure
load('../csgA/sim_frac_csga');
all=[csga1;csga2;csga3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('csgA_nprobrddis');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.5]);
%% Panel F
figure
load('../WT/sim_frac_wt');
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('WT_nprobrddis');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.5]);
