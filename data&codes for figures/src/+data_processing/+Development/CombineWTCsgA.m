load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/AllData_Combined.mat'
AllDataCombined = AllData

%%%
% CsgA
%%%
folder_base = '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/'
load([folder_base, 'AllData.mat'])

% %%%
% % WT
% %%%
% folder_base = '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/'
% load([folder_base, 'AllData_init.mat'])

for i = 1:length(AllData)
    cur = cell2mat(cellfun(@(x) strcmp(x.set,AllData{i}.set),AllDataCombined,'UniformOutput',false))
    AllData{i}.AlignedStart = AllDataCombined{cur}.AlignedStart
    AllData{i}.AlignedStop  = AllDataCombined{cur}.AlignedStop
end

save([folder_base, 'AllDataAligned_init.mat'],'AllData')

%%
addpath('Libraries/TrajectoryAnalysis')
addpath('Libraries/Utils')
addpath('Libraries/AggregateTracking')

% %%%
% % WT
% %%%
% BASE_FOLDER = '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/' %don't forget trailing slash

%%%
% Csga
%%%
BASE_FOLDER = '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/' %don't forget trailing slash

load([BASE_FOLDER, 'AllDataAligned_init.mat'])

for set = 1:length(AllData)
    set
    load([AllData{set}.set 'agg_tracks_with_area.mat']);
    load([AllData{set}.set 'm_tracks.mat']);

    %%
    AllData{set}.AggProps = generateAggProps(agg_tracks);
    
    %%
    sp = strsplit(AllData{set}.set,'/')
    AllData{set}.movie = sp(9);

    AllData{set} = CalculateData(AllData{set}, m_tracks, agg_tracks, stable_aggs, unstable_aggs)
    
end

save([BASE_FOLDER, 'AllDataAligned'],'AllData');