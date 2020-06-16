close all
clear all
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');
org_color=[0 0 1];
rd_color=[1 0 0];
blk=[0,0,0];
colormap = [org_color; rd_color;blk];
load('../WT/sim_frac_wt');
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('frequentNon');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
load('longdur');
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(3,colormap),'alpha');
xlim([0,7]);
ylim([0,0.5]);