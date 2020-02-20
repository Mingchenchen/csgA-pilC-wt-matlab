
%% Global Vars  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
addpath('Libraries/Utils')
run('ENVS.m')

%% Track Cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%

%% WT in WT
run('+data_processing/cell_tracking/Development/WT in WT/Track_10112015.m');
run('+data_processing/cell_tracking/Development/WT in WT/Track_10152015.m');
run('+data_processing/cell_tracking/Development/WT in WT/Track_10252015.m');

%% CsgA in WT
run('+data_processing/cell_tracking/Development/csgA in WT/Track_A_30.m');
run('+data_processing/cell_tracking/Development/csgA in WT/Track_DEV_1.m');
run('+data_processing/cell_tracking/Development/csgA in WT/Track_Exp1_35.m');

%% CsgA in Csga
run('+data_processing/cell_tracking/Development/csgA in csgA/Exp1_34.m');
run('+data_processing/cell_tracking/Development/csgA in csgA/Exp1_32.m');
run('+data_processing/cell_tracking/Development/csgA in csgA/Exp1_31.m.m');

%% Normalize Images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Performed prior to estimating cell density from
% the images.
%%%

%%% Choose what movies to analyze
MAP = ALL_RAW_DEV_MOVIES;
%MAP = WT_RAW_DEV_MOVIES;
%MAP = CSGA_RAW_DEV_MOVIES;

addpath('Libraries/DensityEstimation')
for i = keys(MAP)
    i
    k = i{1};
    [Kum] = normalizeImages([PROJECT_BASENAME MAP(k)], ...
                            [PROJECT_BASENAME, '/data/processed/', MAP(k)]);
end

%% Estimate Cell Density %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate cell density using the normalized images as described in
% Cotter et al. PNAS 2017
%%%

%%% Choose what movies to analyze
MAP = ALL_RAW_DEV_MOVIES;
%MAP = WT_RAW_DEV_MOVIES;
%MAP = CSGA_RAW_DEV_MOVIES;

addpath('Libraries/DensityEstimation')
for i = keys(MAP)
    i
    k = i{1};
    %%% TODO: ADD CODE
    Knorm = load([PROJECT_BASENAME, '/data/processed/', MAP(k), 'Knormalized.mat'])
    Kum = floro2dens(Knorm); % Estimated with PNAS 2018 Coefficients
    save([PROJECT_BASENAME, '/data/processed/', MAP(k), '/Kum'],'Kum')
end

%% Detect Cell State %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detect the cell state (persistent or non-persistent) and add biofilm
% cell density information to the images
%%%


%%% Choose what movies to analyze
MAP = ALL_RAW_DEV_MOVIES;
%MAP = WT_RAW_DEV_MOVIES;
%MAP = CSGA_RAW_DEV_MOVIES;

addpath('Libraries/TrajectoryAnalysis')
for i = keys(MAP)
    i
    k = i{1};
    %%% TODO: ADD CODE
    load([PROJECT_BASENAME, '/data/processed/', MAP(k), '/Kum.mat'])
    load([PROJECT_BASENAME, '/data/processed/', MAP(k), '/tracks.mat'])
    m_tracks = createMTracks(tracks, Kum);
    save([PROJECT_BASENAME, '/data/processed/', MAP(k), '/m_tracks.mat'], 'm_tracks')
end

%% Align Movies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Library/Utils')
%%% Choose what movies to analyze
MAP = All_DEV_FOLDERS;
%MAP = WT_RAW_DEV_MOVIES;
%MAP = CSGA_RAW_DEV_MOVIES;
%MAP = [WT_RAW_DEV_MOVIES CSGA_RAW_
assert(length(unique(keys(MAP))) == length(MAP),'All keys in MAP must be unique')

aggregation_magnitudes = containers.Map();
p = Progress(length(MAP));
c = 0;
for i = keys(MAP)
    p.d(c); c = c + 1;

    k = i{1};

    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k)];
    load([folder '/Kum.mat']);

    aggregation_magnitudes(k) = data_processing.Development.AggAlignment.extract_aggragte_magnitude_dfft(Kum,50,500);
end

save([PROJECT_BASENAME, '/data/interim/Development/aggregation_magnitudes'],'aggregation_magnitudes')

%% Find centerpoint in aggregate growth by finding magnitude > 0.1
load([PROJECT_BASENAME, '/data/interim/Development/aggregation_magnitudes'])

agg_start = containers.Map();

% Detect where aggreagtion begins
for k = keys(aggregation_magnitudes)
    key = k{1};
    c1 = aggregation_magnitudes(key);

    dmin = min(c1);
    dmax = max(c1);

    data = (c1 - dmin) ./ (dmax-dmin);

    agg_start(key) = find(data > 0.2,1,'first');
end

save([PROJECT_BASENAME, '/data/interim/Development/agg_start'],'agg_start')

%% Crop the data so that they all start at the same time
load([PROJECT_BASENAME, '/data/interim/Development/agg_start'])
load([PROJECT_BASENAME, '/data/interim/Development/aggregation_magnitudes'])

movie_start = containers.Map();
min_start =  min(cell2mat(values(agg_start)));
for k = keys(agg_start)
    key = k{1};

    movie_start(key) = agg_start(key) - min_start + 1;
end

movie_length = min(cellfun(@length,values(aggregation_magnitudes)) - cell2mat(values(movie_start)));

save([PROJECT_BASENAME, '/data/interim/Development/movie_length'],'movie_length')
save([PROJECT_BASENAME, '/data/interim/Development/movie_start'],'movie_start')

%% Track Aggreagtes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%
%% TODO: ADD CODE

%% Generate AllData %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%
addpath('Libraries/Utils')
addpath('Libraries/TrajectoryAnalysis/')

load([PROJECT_BASENAME, '/data/processed/Development/AllData_Combined_Init.mat']);
load([PROJECT_BASENAME, '/data/interim/Development/movie_length'])
load([PROJECT_BASENAME, '/data/interim/Development/movie_start'])
WRAP_STOPS = false;
FILTER_UNSTABLE = true;

%%% Choose what movies to analyze
%MAP = All_DEV_FOLDERS;
%MAP = WT_DEV_FOLDERS;
%MAP = CSGA_DEV_FOLDERS;
MAP = [WT_DEV_FOLDERS;  CSGA_DEV_FOLDERS];

aligned_start = cell2mat(values(movie_start,keys(MAP)));
aligned_stop = aligned_start + movie_length;

%%% c
AllData = {};

assert(length(unique(keys(MAP))) == length(keys(MAP)),'All keys in MAP must be unique')

k = keys(MAP);
p = Progress(length(MAP));
for i = 1:length(MAP)
    p.d(i)

    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];

    load([folder '/m_tracks.mat']);
    load([folder '/agg_tracks_with_area.mat']);

    %%
    AllData{i}.movie = k{i};
    AllData{i}.all_runs = createRunVectors(m_tracks,WRAP_STOPS);
    AllData{i}.m_tracks = m_tracks;
    AllData{i}.AlignedStart = aligned_start(i);
    AllData{i}.AlignedStop = aligned_stop(i);

    agg_tracks.stable = ismember(agg_tracks.id,stable_aggs);
    if(FILTER_UNSTABLE)
        % Sort out only the filtered stable aggreagtes
        % This is the version used by Cotter et al. PNAS 2017
        agg_tracks = subStruct(agg_tracks,ismember(agg_tracks.id,stable_aggs),'r');
    else
        %Use both stable and unstable aggreagtes
        agg_tracks = subStruct(agg_tracks,ismember(agg_tracks.id,union(stable_aggs,unstable_aggs)),'r');
    end
    AllData{i}.agg_tracks = agg_tracks;

    AllData{i}.STOPS_WRAPPED = WRAP_STOPS;
    AllData{i}.FILTER_UNSTABLE_AGGS = FILTER_UNSTABLE;
end

save([PROJECT_BASENAME, '/data/processed/Development/' 'AllData'],'AllData');

%% Generate AllData for CsgA in CsgA cells %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%
addpath('Libraries/Utils')
addpath('Libraries/TrajectoryAnalysis/')

load([PROJECT_BASENAME, '/data/processed/Development/AllData_Combined_Init.mat']);
load([PROJECT_BASENAME, '/data/interim/Development/movie_length'])
load([PROJECT_BASENAME, '/data/interim/Development/movie_start'])
WRAP_STOPS = false;
FILTER_UNSTABLE = true;

%%% Choose what movies to analyze
MAP = CSGA_IN_CSGA_DEV_FOLDERS;

%%%
AllData = {};

assert(length(unique(keys(MAP))) == length(keys(MAP)),'All keys in MAP must be unique')

k = keys(MAP);
p = Progress(length(MAP));
for i = 1:length(MAP)
    p.d(i)

    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];

    load([folder '/m_tracks.mat']);

    %%
    AllData{i}.movie = k{i};
    AllData{i}.all_runs = createRunVectors(m_tracks,WRAP_STOPS);
    AllData{i}.m_tracks = m_tracks;
    AllData{i}.AlignedStart = min(m_tracks.frame);
    AllData{i}.AlignedStop = max(m_tracks.frame);

    AllData{i}.STOPS_WRAPPED = WRAP_STOPS;
    AllData{i}.FILTER_UNSTABLE_AGGS = FILTER_UNSTABLE;
end

save([PROJECT_BASENAME, '/data/processed/Development/' 'AllDataCsgAinCsgA'],'AllData');

%% Add Aggreagte Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%
addpath('Libraries/Utils')
load([PROJECT_BASENAME, '/data/processed/Development/' 'AllData']);

for set = 1:length(AllData)
    agg_tracks = AllData{set}.agg_tracks;

    T = table();

    T.majorAxis    = agg_tracks.majorAxis;
    T.minorAxis    = agg_tracks.minorAxis;
    T.orientation  = agg_tracks.orientation;
    T.agg_id       = agg_tracks.id;
    T.area         = agg_tracks.area;
    T.eccentricity = agg_tracks.eccentricity;
    T.mean_density = agg_tracks.meanIntensity;
    T.frame        = agg_tracks.frame;
    T.movie_id     = repmat(set,length(agg_tracks.frame),1);
    T.index        = agg_tracks.index;

    AllData{set}.AggProps = T;
end

save([PROJECT_BASENAME, '/data/processed/Development/' 'AllData'],'AllData');

%% Generate AllDataTable %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Libraries/TrajectoryAnalysis/')

load([PROJECT_BASENAME, '/data/processed/Development/' 'AllData'])

assert(isfield(AllData{1},'agg_tracks'), 'Must run above code block to add aggreagte properties')
assert(isfield(AllData{1},'AggProps'), 'Must run above code block to add aggreagte properties')

AllDataTable = CompileVariables(AllData);

save([PROJECT_BASENAME, '/data/processed/Development/AllDataTable.mat'],'AllDataTable');

%% Generate AllDataTable for CsgA in CsgA cells %%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Libraries/TrajectoryAnalysis/')

load([PROJECT_BASENAME, '/data/processed/Development/' 'AllDataCsgAinCsgA'])

AllDataTable = CompileVariables(AllData);

save([PROJECT_BASENAME, '/data/processed/Development/AllDataTableCsgAinCsgA.mat'],'AllDataTable');

%% Generate Frac Cells in aggreagte curves
%
%%%
load([PROJECT_BASENAME, '/data/interim/Development/movie_length'])
load([PROJECT_BASENAME, '/data/interim/Development/movie_start'])

MAP = All_DEV_FOLDERS;

AGG_DESNITY_CUTOFF = 2.3; %Cutoff density for detecting aggreagtes
                            %in cells/um^2

frac_curves = containers.Map();
for k = keys(MAP)
    key = k{1};

    first_frame = movie_start(key);
    last_frame = first_frame + movie_length;

    load([PROJECT_BASENAME, '/data/processed/' MAP(key), 'm_tracks.mat'])

    tracks = subStruct(m_tracks,m_tracks.frame > first_frame);
    tracks.frame = tracks.frame - first_frame;

    in_agg = tracks.density > AGG_DESNITY_CUTOFF;

    frac_curves(key) = histc(tracks.frame(in_agg),0:movie_length)' ./ histc(tracks.frame,0:movie_length)';
end

save([PROJECT_BASENAME, '/data/processed/Development/frac_curves.mat'],'frac_curves');
