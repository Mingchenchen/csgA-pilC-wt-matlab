addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME '/data/processed/Development/AllDataTable.mat'])
T = AllDataTable;

load([PROJECT_BASENAME '/data/processed/Development/AllData.mat'])
A = AllData;

%% WT
AllDataTable = T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:);

AllData = {};
for i = 1:length(A)
    if(ismember(A{i}.movie,keys(WT_DEV_FOLDERS)))
       AllData = [AllData; A{i}];
    end
end

save([PROJECT_BASENAME '/data/processed/Development/WT in WT/AllDataTable'],'AllDataTable');
save([PROJECT_BASENAME '/data/processed/Development/WT in WT/AllData'],'AllData');

%% csgA
AllDataTable = T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:);

AllData = {};
for i = 1:length(A)
    if(ismember(A{i}.movie,keys(CSGA_DEV_FOLDERS)))
       AllData = [AllData; A{i}];
    end
end

save([PROJECT_BASENAME '/data/processed/Development/csgA in WT/AllDataTable'],'AllDataTable');
save([PROJECT_BASENAME '/data/processed/Development/csgA in WT/AllData'],'AllData');
