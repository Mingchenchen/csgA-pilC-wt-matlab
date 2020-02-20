%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

% How many bootstrap samples
% 1000 for publication, set to 10 to keep things moving quickly for testing
NBOOT = 1000;

% Density cutoff defining inside a ripple
DENSITY_CUTOFF = 0;

%
RIPPLE_START = 200;

% % Publication: setting to true supresses some labeling during plotting
% % to make plot publication ready, but less visually readable. 
% PUB_READY = false;
% 
% % Load previous bootstrapped data
% LOAD_PREVIOUS = false;

SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/rippleing/RevInRipple/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

BAR_WIDTH = 0.9;

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [15, 25];

%%%
%% END SETABLE PARAMS %%

%%%
% Real
%%%

%%%% WT %%%%

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataTableWrappedStops.mat']);
load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataWrappedStops.mat']);
WT_SETS = unique(AllDataTable.set(ismember(cellstr(AllDataTable.movie),keys(WT_RIPPLE_FOLDERS))))';

Rwt = FigureCreation.Rippleing.revs_in_ripple(AllDataTable,AllDataWrappedStops,WT_SETS,RIPPLE_START,NBOOT,DENSITY_CUTOFF);

%%%% CsgA %%%%
CSGA_SETS = unique(AllDataTable.set(ismember(cellstr(AllDataTable.movie),keys(CSGA_RIPPLE_FOLDERS))))';

RcsgA = FigureCreation.Rippleing.revs_in_ripple(AllDataTable,AllDataWrappedStops,CSGA_SETS,RIPPLE_START,NBOOT,DENSITY_CUTOFF);

%%%% DifA %%%%
% DIFA_SETS = unique(AllDataTable.set(ismember(cellstr(AllDataTable.movie),keys(DIFA_RIPPLE_FOLDERS))))';
% 
% RdifA = FigureCreation.Rippleing.revs_in_ripple(AllDataTable,AllDataWrappedStops,DIFA_SETS,RIPPLE_START,NBOOT,DENSITY_CUTOFF);

%%
% Random
%%%
load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataRandomWrappedStops.mat']);
load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataTableRandomWrappedStops.mat']);

%%%% WT %%%%
WT_SETS = unique(AllDataTable.set(ismember(cellstr(AllDataTable.movie),keys(WT_RIPPLE_FOLDERS))))';

Rwt_random = FigureCreation.Rippleing.revs_in_ripple(AllDataTable,AllDataRandomWrappedStops,WT_SETS,RIPPLE_START,NBOOT,DENSITY_CUTOFF);

%%%% CsgA %%%%
CSGA_SETS = unique(AllDataTable.set(ismember(cellstr(AllDataTable.movie),keys(CSGA_RIPPLE_FOLDERS))))';

RcsgA_random = FigureCreation.Rippleing.revs_in_ripple(AllDataTable,AllDataRandomWrappedStops,CSGA_SETS,RIPPLE_START,NBOOT,DENSITY_CUTOFF);


%%
figure, hold on
    %%%% Reversals %%%%
    WTmean = mean(Rwt.revs.samples);
    CsgAmean = mean(RcsgA.revs.samples);
    %DifAmean = mean(RcsgA.revs.samples);
    
    bar_values = [WTmean,CsgAmean];
    ci_neg = [WTmean - Rwt.revs.ci(1), ...
              CsgAmean - RcsgA.revs.ci(1)];

    ci_pos = [Rwt.revs.ci(2) - WTmean, ...
              RcsgA.revs.ci(2) - CsgAmean];
    
    points = 1:length(ci_pos);
    xi = [1:length(bar_values)] * BAR_WIDTH;
    
    l1 = plot(points,bar_values,'or');
    e1 = errorbar(points,bar_values,ci_neg,ci_pos,'r');
    e1.LineStyle = 'none';
    

    %%%% All positions %%%%
    WTmean = mean(Rwt.trajectories.samples);
    CsgAmean = mean(RcsgA.trajectories.samples);
    %DifAmean = mean(RdifA.trajectories.samples);
    
    bar_values = [WTmean,CsgAmean];
    ci_neg = [WTmean - Rwt.trajectories.ci(1), ...
              CsgAmean - RcsgA.trajectories.ci(1)];

    ci_pos = [Rwt.trajectories.ci(2) - WTmean, ...
              RcsgA.trajectories.ci(2) - CsgAmean];
    
    l2 = plot(points,bar_values,'ob');
    e2 = errorbar(points,bar_values,ci_neg,ci_pos,'b');
    e2.LineStyle = 'none';
    
    %%% 
    % Randomized 
    %%%
    WTmean = mean(Rwt_random.revs.samples);
    CsgAmean = mean(RcsgA_random.revs.samples);    
    
    bar_values = [WTmean,CsgAmean];
    ci_neg = [WTmean - Rwt_random.revs.ci(1), ...
              CsgAmean - RcsgA_random.revs.ci(1)];

    ci_pos = [Rwt_random.revs.ci(2) - WTmean, ...
              RcsgA_random.revs.ci(2) - CsgAmean];
    
    l3 = plot(points,bar_values,'*r');
    e3 = errorbar(points,bar_values,ci_neg,ci_pos,'r');
    e3.LineStyle = 'none';
    
    WTmean = mean(Rwt_random.trajectories.samples);
    CsgAmean = mean(RcsgA_random.trajectories.samples);
    %DifAmean = mean(RdifA.trajectories.samples);
    
    bar_values = [WTmean,CsgAmean];
    ci_neg = [WTmean - Rwt_random.trajectories.ci(1), ...
              CsgAmean - RcsgA_random.trajectories.ci(1)];

    ci_pos = [Rwt_random.trajectories.ci(2) - WTmean, ...
              RcsgA_random.trajectories.ci(2) - CsgAmean];
    
    l4 = plot(points,bar_values,'*b');
    e4 = errorbar(points,bar_values,ci_neg,ci_pos,'b');
    e4.LineStyle = 'none';
    
    %%% Figure Style %%%
    ax = gca;
    ax.YLim = [0,1];
    ax.XLim = [points(1) - points(1),points(end) + points(1)];
    ax.XTick = points;
    ax.XTickLabel = {'WT','csgA'};
    xlabel('Genotype')
    ylabel('Fraction Inside Ripple Crests')
    %legend([l1 l2, e1 e3],{'Reversals','All Cell Positions','Exp. Data','Randomized Controls'})
       
    LINE_WIDTH = 1.5;
    CAP_SIZE = 25;
    e1.LineWidth = LINE_WIDTH;
    e1.CapSize = CAP_SIZE;
    e2.LineWidth = LINE_WIDTH;
    e2.CapSize = CAP_SIZE;
    e3.LineWidth = LINE_WIDTH;
    e3.CapSize = CAP_SIZE;
    e4.LineWidth = LINE_WIDTH;
    e4.CapSize = CAP_SIZE;

if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER 'RevInRipple'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end
    