close all
clear all

addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');

%% Panel A
%%load pilC only experimental figures
fname='pilCend.tif';
load(fname)
 info = imfinfo(fname);
 %plot and detect cells
detectCellsRegionProps(fname,'MinLevel',30,'check',length(info));

%% Panel B
%%load csgA only experimental figures
fname='csgAend.tif';
load(fname)
 info = imfinfo(fname);
 %plot and detect cells
detectCellsRegionProps(fname,'MinLevel',30,'check',length(info));

%% Panel C
%%load pilC with WT experimental figures, tif image can be found in Dryad
fname='4223w2_1_MMStack_Pos0.ome.tif';
load(fname)
 info = imfinfo(fname);
 %plot and detect cells
detectCellsRegionProps(fname,'MinLevel',30,'check',length(info));

%% Panel D
%%load csgA with WT experimental figures, tif image can be found in Dryad
fname='A_30_MMStack_Pos0.ome.tif';
load(fname)
 info = imfinfo(fname);
 %plot and detect cells
detectCellsRegionProps(fname,'MinLevel',30,'check',length(info));