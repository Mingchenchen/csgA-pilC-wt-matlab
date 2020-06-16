clear all
load('data')
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');
colormap=[0 0 1;1 0 0;0 0 0];
xi = (1:length(upilc)) / 2 / 60-1.5;
boundedline(xi,upilc,[spilc,spilc],'-','cmap',color_chooser(3,colormap),'alpha');
xi = (1:length(ucsga)) / 2 / 60-1.5;
boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2,colormap),'alpha');
xlim([0,5]);
ylim([0,0.5]);