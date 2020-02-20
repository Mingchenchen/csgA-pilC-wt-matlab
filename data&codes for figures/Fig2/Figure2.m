clear all
close all

load('../pilC/pilCExpData')
load('../csgA/csgAExpData')
load('../WT/WTExpData')
%% Panel A
load('../csgA/aggfrac1')
csgAAggFrac=aggfrac1;
load('../csgA/aggfrac2')
csgAAggFrac=[csgAAggFrac;aggfrac2];
load('../csgA/aggfrac3')
csgAAggFrac=[csgAAggFrac;aggfrac3];
aggfrac=mean(csgAAggFrac);
xi=1:length(aggfrac);
plot(xi/120,aggfrac,'r');
hold on
load('../pilC/aggfrac1')
pilCAggFrac=aggfrac1;
load('../pilC/aggfrac2')
pilCAggFrac=[pilCAggFrac;aggfrac2];
load('../pilC/aggfrac3')
pilCAggFrac=[pilCAggFrac;aggfrac3];
aggfrac=mean(pilCAggFrac);
xi=1:length(aggfrac);
plot(xi/120,aggfrac,'k');

load('../WT/aggfrac1')
WTAggFrac=aggfrac1;
load('../WT/aggfrac2')
WTAggFrac=[WTAggFrac;aggfrac2];
load('../WT/aggfrac3')
WTAggFrac=[WTAggFrac;aggfrac3];
aggfrac=mean(WTAggFrac);
xi=1:length(aggfrac);
plot(xi/120,aggfrac,'b');
hold off
xlim([0,5]);
%% Panel B
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');
xi_agg = floor([2 3 4] * 60 * 2);

props = FigureCreation.Development.AggPlots();

WT =  props.combine_agg_props(AllData);
pilC =  props.combine_agg_props(AllData1);
csgA =  props.combine_agg_props(AllData2);

figure
    props.agg_count(csgA,WT,pilC,xi_agg)

ylabel('Agg Count')
%% Panel C
figure
    props.plot_property(csgA,WT,pilC,xi_agg,'area')

ylabel('Area (\mum^2)')
  
    
   