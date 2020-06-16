close all
clear all
colormap=[[0,0,1];[0 0 0]];
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');
load('../csgA/sim_frac_csga');
all=[csga1;csga2;csga3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('shortdur')
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.5]);