run('./ENVS.m')
addpath('Libraries/LAPTracker')
addpath('Libraries/Utils')
addpath('Libraries/TrajectoryAnalysis')

tiff_file = WT_RAW_DEV_MOVIES('10282015');

processed_folder = ['../data/processed/', WT_DEV_FOLDERS('10282015')];
check_folder = ['../reports/DataChecks/', WT_DEV_FOLDERS('10282015')];

mkdir(check_folder)
mkdir(processed_folder)

SAVE_RESULTS = true;
%OVERWRITE = true; % Overrides setting in ENVS.m

%%
screePlot(tiff_file,2.0:0.2:5);

if(SAVE_RESULTS)
    if(~OVERWRITE && exist([check_folder 'screePlot.svg'],'file') > 0)
        error('File screePlot.svg already Exists, set OVERWRITE = true to overwrite');
    end
    saveFigures('SaveAs',[check_folder 'screePlot'],'Style','none','Formats',{'svg','fig'});
end

%%
if(~HEADLESS)
    figure
    detectCellsRegionProps(tiff_file,'MinLevel',3.4,'check',1,'large_cutoff',20);
end

%%
posList1 = detectCellsRegionProps(tiff_file,'MinLevel',3.4,'large_cutoff',20);

if(SAVE_RESULTS)
    if(~OVERWRITE && exist([save_folder 'posList.mat'],'file') > 0)
        error('File posList.mat already exists, set OVERWRITE = true to overwrite');
    end
    save([save_folder 'posList'],'posList1')
end

%%
if(~HEADLESSS)
    figure
    for i = min(posList1.frame):max(posList1.frame)
        t = subStruct(posList1,posList1.frame == i);

        cla
        plot(t.x,t.y,'o')
        title(i)
        drawnow
        pause(0.1)
    end
end

%%
%For csgA in csgA tracking
[ tracks, NO_LINKING_COST, CUTOFFS, deltaD,deltaO,displacement] = LAPtracking(posList1, ...
                                                   'Debug',true, ...
                                                   'MinTrackLength',10, ...
                                                   'MaxSliceGap',3, ...
                                                   'GapClosingSearchRadiusCutoff',20, ...
                                                   'EndCostMutiplier',1.1, ...
                                                   'IncludePosListIndex',true, ...
                                                   'NoiseCutoff',1, ...
                                                   'EndCostSeed',15, ...
                                                   'NoGapClosing',false);
                
if(SAVE_RESULTS)
    if(~OVERWRITE && exist([save_folder 'tracks.mat'],'file') > 0)
        error('File tracks.mat already exists, set OVERWRITE = true to overwrite');
    end
    save([save_folder 'tracks'],'tracks')
end

%%
if(~HEADLESS)
    figure
    axis

    for i = (max(tracks.frame)-150):max(tracks.frame)
        t = subStruct(tracks,tracks.frame == i);

        cla
        hold on
        %A = double(bpass(imread(fname1,i,'Info',info),1,10));
        %imagesc(A)
        plot(t.x,t.y,'o')
        xlim([0,max(tracks.x) + 10]);
        ylim([0,max(tracks.y) + 10]);
        title(i)
        drawnow
        pause(0.1)
    end
end
%%
visRuns(tracks);

if(SAVE_RESULTS)
    if(~OVERWRITE && exist([check_folder 'trajectories.png'],'file') > 0)
        error('File trajectories.png already exists, set OVERWRITE = true to overwrite');
    end
    %A .fig format is not saved because it takes way to long to export (due to the large number of plotted lines)
    saveFigures('SaveAs',[check_folder 'trajectories.png'],'Style','none','Formats',{'png'});
end