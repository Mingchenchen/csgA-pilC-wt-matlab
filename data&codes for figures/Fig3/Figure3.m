clear all
close all
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');
CSGA_COLOR = [1 0 0];
WT_COLOR = [0 0 1];
PILC_COLOR=[0 0 0];
% Auto calculated based on above settings
colormap = [WT_COLOR; CSGA_COLOR; PILC_COLOR];
%% Panel A
figure
load('../WT/frac_exp1')
load('../WT/frac_exp2')
load('../WT/frac_exp3')
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
hold on
load('../csgA/frac_exp1')
load('../csgA/frac_exp2')
load('../csgA/frac_exp3')
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
load('../pilC/frac_exp1')
load('../pilC/frac_exp2')
load('../pilC/frac_exp3')
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(3,colormap),'alpha');
xlim([0,5]);
hold off
%% Panel B
figure
load('../WT/sim_frac_wt');
all=[WT1;WT2;WT3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(1,colormap),'alpha');
load('../csgA/sim_frac_csga');
all=[csga1;csga2;csga3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(2,colormap),'alpha');
load('../pilC/sim_frac_pilc');
all=[pilc1;pilc2;pilc3];
m=mean(all);
s=std(all,[],1)' ;
xi = (1:length(m)) / 2 / 60-1.5;
boundedline(xi,m,[s,s],'-','cmap',color_chooser(3,colormap),'alpha');
xlim([0,5]);
hold off
