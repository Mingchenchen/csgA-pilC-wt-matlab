addpath('Libraries/TrajectoryAnalysis')
addpath('Libraries/Utils')
addpath('Libraries/AggregateTracking')

BASE_FOLDER = '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/' %don't forget trailing slash

load([BASE_FOLDER, 'AllData_Combined_Init.mat'])

for set = 1:length(AllData)
    set
    load([AllData{set}.set 'agg_tracks_with_area.mat']);
    load([AllData{set}.set 'm_tracks.mat']);
    
    %%
    sp = strsplit(AllData{set}.set,'/')
    AllData{set}.movie = sp(9);

    AllData{set} = CalculateData(AllData{set}, m_tracks, agg_tracks, stable_aggs, unstable_aggs,'IncludeUnstable',false)
    
    %%
    AllData{set}.AggProps = generateAggProps(AllData{set}.agg_tracks);
end

save([BASE_FOLDER, 'AllData'],'AllData');

%%
%generate all data table
%%%
AllDataTable = CompileVariables(AllData);
save([BASE_FOLDER, 'AllDataTable'],'AllDataTable');