%% Global Vars  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = [pwd, '/../'];
run('ENVS.m')

%% Track cells in rippling movies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% WT
run('+data_processing/cell_tracking/Rippleing/WT in WT/Track_RIPPLE_3.m');
run('+data_processing/cell_tracking/Rippleing/WT in WT/Track_RIPPLE_24.m');
run('+data_processing/cell_tracking/Rippleing/WT in WT/Track_RIPPLE_28.m');

%% CsgA 
run('+data_processing/cell_tracking/Rippleing/csgA in WT/Track_RIPPLE_25.m');
run('+data_processing/cell_tracking/Rippleing/csgA in WT/Track_RIPPLE_26.m');
run('+data_processing/cell_tracking/Rippleing/csgA in WT/Track_RIPPLE_27.m');

%% DifA
run('+data_processing/cell_tracking/Rippleing/difA in WT/Track_RIPPLE_6.m');

%% Normalize Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This would normally be performed prior to estimating cell density from
% the images, but for rippling we don't care about the real cell density,
% so in the following steps we just use the normalized images.
%%%

%%% Choose what movies to analyze
%MAP = ALL_RAW_RIPPLE_MOVIES;
%MAP = WT_RAW_RIPPLE_MOVIES;
%MAP = CSGA_RAW_RIPPLE_MOVIES;
MAP = DIFA_RAW_RIPPLE_MOVIES;

addpath('Libraries/DensityEstimation')
for i = keys(MAP)
    i
    k = i{1};
    [Kum] = normalizeImages([PROJECT_BASENAME MAP(k)], ...
                            [PROJECT_BASENAME, '/data/processed/', MAP(k)]);                
end

%% Smooth Images and futher normalization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subtract off a local gaussian window to perform furlocal background subtractioin that was not performed in 
% filter_images Cotter et al 2017 and apply a large gaussin filter to 
% act as a low pass filter to filter out individual cells and make ripples 
% more pronounced. 
%%%

%%% Choose what movies to analyze
% MAP = All_RIPPLE_FOLDRS;
% MAP = WT_RIPPLE_FOLDERS;
% MAP = CSGA_RIPPLE_FOLDERS
MAP = DIFA_RIPPLE_FOLDERS

addpath('Libraries/Rippleing')
for i = keys(MAP)
    i
    k = i{1};
    
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k)];
    load([folder '/Knormalized.mat'])
    [Kripple] = filter_images(K);   
    save([folder '/Kripple.mat'],'Kripple','-v7.3');
end

%% Extract Wavelengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use a 2D fourier transform to extract the magnitude of various 
% spatial wavelengths in each frame
%
% see ripple_workflow_checks for an exmample on how to use this data
%
%      Hand measured statistics (from 10 measurements)
% Movie      Mean (um)   Min (um)    Max (um)
% Ripple_3      99         91          110
% Ripple_24     108        95          124
% Ripple_28     93         80          108          
%%%
addpath('Libraries/Rippleing')

%%% Choose what movies to analyze
% MAP = All_RIPPLE_FOLDRS;
% MAP = WT_RIPPLE_FOLDERS;
MAP = CSGA_RIPPLE_FOLDERS;
% MAP = DIFA_RIPPLE_FOLDERS

p = Progress(length(MAP));
k = keys(MAP);
for i = 1:length(MAP)
    p.d(i)
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];
    load([folder '/Kripple.mat'])
    wavelengths = rippleWavelength(Kripple);
    save([folder '/wavelengths.mat'],'wavelengths');
end


%% Extract Wave Periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Peforms short time 1D fourier transforms to extract the average 
% wavelength periods in the movie
%%%
addpath('Libraries/Rippleing')

FIRST_FRAME = 200;
LAST_FRAME = 450; %Have to use 400 for DIFA..

Fs = 2; %Sample frequency in samples/min
WINDOW_SIZE = 2^5;
STEP_SIZE = 5;

%%% Choose what movies to analyze
% MAP = All_RIPPLE_FOLDRS;
% MAP = WT_RIPPLE_FOLDERS;
%MAP = CSGA_RIPPLE_FOLDERS;
% MAP = DIFA_RIPPLE_FOLDERS
MAP = [WT_RIPPLE_FOLDERS; CSGA_RIPPLE_FOLDERS];

p = Progress(length(MAP));
k = keys(MAP);
fft_bin = containers.Map();
for i = 1:length(MAP)
    p.d(i)
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];
    load([folder '/Kripple.mat'])
    [fz, centers] = ripplePeriod(Kripple(:,:,FIRST_FRAME:LAST_FRAME),WINDOW_SIZE,STEP_SIZE);
    fft_bin(MAP(k{i})) = fz;
end

centers = (FIRST_FRAME + centers) / 2 / 60;
periods = WINDOW_SIZE ./ cell2mat(values(fft_bin));
save([PROJECT_BASENAME, '/data/processed/Rippleing/' 'RipplePeriods'],'fft_bin','periods','centers');

%% Build ripple_tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% Kripple is filtered using a 1D fourier transform along the time axis
% for each pixel position. This helps with ripple crest detection and 
% normalizes the mean cell density at 0. A good cutoff for inside a ripple
% crest is then density > 0
%%%
addpath('Libraries/Utils')
addpath('Libraries/Rippleing')
addpath('Libraries/TrajectoryAnalysis/')

FIRST_FRAME = 200;
LAST_FRAME = 450; %Have to use 400 for DIFA...

%%% Choose what movies to analyze
% MAP = All_RIPPLE_FOLDRS;
% MAP = WT_RIPPLE_FOLDERS;
% MAP = CSGA_RIPPLE_FOLDERS;
% MAP = DIFA_RIPPLE_FOLDERS;
MAP = [WT_RIPPLE_FOLDERS; CSGA_RIPPLE_FOLDERS];

p = Progress(length(MAP));
k = keys(MAP);
for i = 1:length(MAP)
    p.d(i)
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];
    
    load([folder 'Kripple.mat'])
    load([folder 'tracks'])
    
    %%
    [Kfripple,~,~,~] = fftRippleFilter(Kripple(:,:,FIRST_FRAME:LAST_FRAME),11);

    %Padding with NaNs is the simplest way to still use createMTracks,
    % and all the data of tracks, since sliceing tracks would result in 
    % some very short trajectories.
    paddedKfripple = NaN(size(Kripple));
    paddedKfripple(:,:,FIRST_FRAME:LAST_FRAME) = Kfripple;
    ripple_tracks = createMTracks(tracks,paddedKfripple);
    clear paddedKfripple;
    
    save([folder 'ripple_tracks'],'ripple_tracks')
    save([folder 'Kfripple'],'Kfripple','-v7.3')
end

%% Build random ripple_tracks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
% Same as ripple_tracks, except the frames of Kripple are shuffled to uncoorleate
% cell reversals with ripple crests. This is used as a control since we do
% not have data from a cell type that we are sure does nto particiapte in ripples. 
%%%
addpath('Libraries/Utils')
addpath('Libraries/Rippleing')
addpath('Libraries/TrajectoryAnalysis/')

FIRST_FRAME = 200;
LAST_FRAME = 450; %Have to use 400 for DIFA...

%%% Choose what movies to analyze
% MAP = All_RIPPLE_FOLDRS;
% MAP = WT_RIPPLE_FOLDERS;
% MAP = CSGA_RIPPLE_FOLDERS;
% MAP = DIFA_RIPPLE_FOLDERS;
MAP = [WT_RIPPLE_FOLDERS; CSGA_RIPPLE_FOLDERS];

p = Progress(length(MAP));
k = keys(MAP);
for i = 1:length(MAP)
    p.d(i)
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];
    
    load([folder 'Kfripple.mat'])
    load([folder 'tracks'])
    
    %Shuffle Kfripple
    shuffled = randperm(size(Kfripple,3));
    Kfripple = Kfripple(:,:,shuffled);
    
    %Padding with NaNs is the simplest way to still use createMTracks,
    % and all the data of tracks, since sliceing tracks would result in 
    % some very short trajectories.
    paddedKfripple = NaN([size(Kfripple,1),size(Kfripple,2),max(tracks.frame)]);
    paddedKfripple(:,:,FIRST_FRAME:LAST_FRAME) = Kfripple;

    ripple_tracks = createMTracks(tracks,paddedKfripple);
    clear paddedKfripple;
    
    save([folder 'random_ripple_tracks'],'ripple_tracks')
end

%% Generate AllDatas %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a structure that contains all informatoin about each movie.
% I've called these "AllData". 
%
% AllDataWrappedStops: Parts of cell trajectories that were classified as 
%  stops are wrapped into previous persistent states so that runs only 
%  end at movement reversals.

addpath('Libraries/Utils')
addpath('Libraries/TrajectoryAnalysis/')

MAX_FRAME = 400;

%%% Choose what movies to analyze
% MAP = ALL_RIPPLE_FOLDERS;
% MAP = WT_RIPPLE_FOLDERS;
% MAP = CSGA_RIPPLE_FOLDERS;
%MAP = DIFA_RIPPLE_FOLDERS;
MAP = [WT_RIPPLE_FOLDERS; CSGA_RIPPLE_FOLDERS];

%%% 
AllDataWrappedStops = data_processing.Rippleing.generateAllData(MAP,PROJECT_BASENAME,true,MAX_FRAME,'ripple_tracks.mat');
AllDataRandomWrappedStops =  data_processing.Rippleing.generateAllData(MAP,PROJECT_BASENAME,true,MAX_FRAME,'random_ripple_tracks.mat');
AllData = data_processing.Rippleing.generateAllData(MAP,PROJECT_BASENAME,false,MAX_FRAME,'ripple_tracks.mat');
AllDataFull = data_processing.Rippleing.generateAllData(MAP,PROJECT_BASENAME,false,repmat(1000,1,length(MAP)),'ripple_tracks.mat');

save([PROJECT_BASENAME '/data/processed/Rippleing/' 'AllDataWrappedStops'],'AllDataWrappedStops');
save([PROJECT_BASENAME '/data/processed/Rippleing/' 'AllDataRandomWrappedStops'],'AllDataRandomWrappedStops');
save([PROJECT_BASENAME, '/data/processed/Rippleing/' 'AllData'],'AllData');
save([PROJECT_BASENAME, '/data/processed/Rippleing/' 'AllDataFull'],'AllDataFull');

%% Generate AllDataTables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataWrappedStops.mat']);
AllDataTable = CompileVariables(AllDataWrappedStops);
save([PROJECT_BASENAME, '/data/processed/Rippleing/', 'AllDataTableWrappedStops'],'AllDataTable');

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataRandomWrappedStops.mat']);
AllDataTable = CompileVariables(AllDataRandomWrappedStops);
save([PROJECT_BASENAME, '/data/processed/Rippleing/', 'AllDataTableRandomWrappedStops'],'AllDataTable');

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllData.mat']);
AllDataTable = CompileVariables(AllData);
save([PROJECT_BASENAME, '/data/processed/Rippleing/', 'AllDataTable'],'AllDataTable');

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataFull.mat']);
AllDataTable = CompileVariables(AllDataFull);
save([PROJECT_BASENAME, '/data/processed/Rippleing/', 'AllDataTableFull'],'AllDataTable');