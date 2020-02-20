% %%%
% % CsgA in WT
% %%%
fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/Exp1_35/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/DEV_1/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/A_30/'};
                  
%%% 
% WT in WT
%%%
% fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10282015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10152015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10112015/'};
      
addpath('Libraries/LAPTracker')
addpath('Libraries/Utils')

for set = 1:length(fnames)
    set
    load([fnames{set} 'Kum.mat']);
    
    %%
    [agg_tracks,stable_aggs,unstable_aggs] = trackAggregates(Kum,2.3);
    save([fnames{set} 'agg_tracks_with_area'],'agg_tracks','stable_aggs','unstable_aggs');
    
    %%
    clear Kum
end
