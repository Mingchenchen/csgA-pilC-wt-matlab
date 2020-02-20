% Generates a table of aggreagte properties from agg_tracks_with_area.mat files
%%

% %WT in WT
% load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/AllDataAligned.mat'

% Csga in WT
load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/AllDataAligned.mat'

%%
addpath('Libraries/Utils')
AggProps = table();

for set = 1:length(AllData)
    agg_tracks = AllData{set}.agg_tracks;
    real_aggs = subStruct(agg_tracks,ismember(agg_tracks.id,union(stable_aggs,unstable_aggs)));

    T = table();

    T.majorAxis    = real_aggs.majorAxis;
    T.minorAxis    = real_aggs.minorAxis;
    T.orientation  = real_aggs.orientation;
    T.id           = real_aggs.id;
    T.area         = real_aggs.area;
    T.eccentricity = real_aggs.eccentricity;
    T.mean_density = real_aggs.meanIntensity;
    T.frame        = real_aggs.frame;
    T.movie_id     = repmat(set,length(real_aggs.frame),1);
    
    AggProps = [AggProps; T];
end

%%

% % WT in WT
% save '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/AggregateProperites.mat' AggProps

% Csga in WT
save '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/AggregateProperites.mat' AggProps