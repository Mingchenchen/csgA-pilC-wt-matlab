%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
addpath('Libraries/TrajectoryAnalysis/')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllData.mat']);

%AllData cell indicies for each movie type
WTinWT = [1:3];
CsgAinWT = [4:6];

% PERSISTENT = true caluclates graph for persistent runs
% PERSISTENT = false calculates graph for non-persistent runs
PERSISTENT = false;

% How many bootstrap samples
% 1000 for publication, set to 10 to keep things moving quickly for testing
NBOOT = 1000;

%Number of frames to skip prior to calculating inside aggreagte values
%Set to a frame after aggreagte tracking becomes fairly stable. 
STABLE_AGG_CUTOFF = 120;

% Sliding window size, in frames
WINDOW = 20;

% Inside aggreagte density cutoff (cells/um^2)
DENSITY_CUTOFF = 2.3;

% Only measuremnts taken after this frame are included in graphs
TIME_CUTOFF = 3 * 2 * 60; % Set to 3 hours after start of experiment

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = true;

SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/Instantaious/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

xi = [-50:1:150];

CSGA_COLOR = [1 0 0];

CSGA_IN_CSGA_COLOR = [0 0 1];

WT_COLOR = [0 0 0];

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [14, 13];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

%%
WT.inside_agg = []; %Speeds measured inside aggs;
WT.outside_agg = []; %Speeds measured outside aggs;

CsgA.inside_agg = [];
CsgA.outside_agg = [];

for mov = WTinWT
    tracks = AllData{mov}.m_tracks;
    tracks = subStruct(tracks,tracks.frame > TIME_CUTOFF);
    
    [spd, distance, state, dens] = calculateSpeed(tracks);
    
    WT.inside_agg = [WT.inside_agg; spd(dens >= DENSITY_CUTOFF & state < 3)];
    WT.outside_agg = [WT.outside_agg; spd(dens < DENSITY_CUTOFF & state < 3)];
end

for mov = CsgAinWT
    tracks = AllData{mov}.m_tracks;
    tracks = subStruct(tracks,tracks.frame > TIME_CUTOFF);
    
    [spd, distance, state, dens] = calculateSpeed(tracks);
    
    CsgA.inside_agg = [CsgA.inside_agg; spd(dens >= DENSITY_CUTOFF & state < 3)];
    CsgA.outside_agg = [CsgA.outside_agg; spd(dens < DENSITY_CUTOFF & state < 3)];
end

%%
text = ['   Cells     |   Speed +/- std (um/min)\n', ...
        '-------------|---------------------------\n', ...
        'CsgA inside  |    %0.2f +/- %0.2f       \n',  ...
        'CsgA outside |    %0.3f +/- %0.2f       \n',  ...
        '   WT inside |    %0.3f +/- %0.2f       \n',  ...
        '  WT outside |    %0.3f +/- %0.2f       \n'];
        
sprintf(text,mean(CsgA.inside_agg),std(CsgA.inside_agg), ...
             mean(CsgA.outside_agg),std(CsgA.outside_agg), ...
             mean(WT.inside_agg),std(WT.inside_agg), ...
             mean(WT.outside_agg),std(WT.inside_agg))