if(~exist('PROJECT_BASENAME','var'))
    %this should be the project root directory if running this script
    %directly. Otherwise it should be set by calling script.
    PROJECT_BASENAME = '../../../../'; 
end
if(~exist('SAVE_RESULTS','var')) SAVE_RESULTS = true; end
if(~exist('HEADLESS','var'))     HEADLESS = false;    end
if(~exist('OVERWRITE','var'))    OVERWRITE = false;   end

run([PROJECT_BASENAME, '/src/ENVS.m'])

addpath([PROJECT_BASENAME, '/src/Libraries/LAPTracker'])
addpath([PROJECT_BASENAME, '/src/Libraries/Utils'])
addpath([PROJECT_BASENAME, '/src/Libraries/TrajectoryAnalysis'])

tiff_file = CSGA_RAW_RIPPLE_MOVIES('RIPPLE_27');

processed_folder = [PROJECT_BASENAME, '/data/processed/', CSGA_RIPPLE_FOLDERS('RIPPLE_27')];
check_folder = [PROJECT_BASENAME '/reports/DataChecks/', CSGA_RIPPLE_FOLDERS('RIPPLE_27')];

mkdir(check_folder)
mkdir(processed_folder)

%%
screePlot(tiff_file,3:.3:6);

if(SAVE_RESULTS)
    if(~OVERWRITE && exist([check_folder 'screePlot.svg'],'file') > 0)
        error('File screePlot.svg already Exists, set OVERWRITE = true to overwrite');
    end
    saveFigures('SaveAs',[check_folder 'screePlot'],'Style','none','Formats',{'svg','fig'});
end

%%
if(~HEADLESS)
    figure
    detectCellsRegionProps(tiff_file,'MinLevel',5.4,'check',1,'large_cutoff',20);
end

%%
posList1 = detectCellsRegionProps(tiff_file,'MinLevel',5.4,'small_cutoff',3);

if(SAVE_RESULTS)
    if(~OVERWRITE && exist([processed_folder 'posList.mat'],'file') > 0)
        error('File posList.mat already exists, set OVERWRITE = true to overwrite');
    end
    save([processed_folder 'posList'],'posList1')
end

%%
if(~HEADLESS)
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
    if(~OVERWRITE && exist([processed_folder 'tracks.mat'],'file') > 0)
        error('File tracks.mat already exists, set OVERWRITE = true to overwrite');
    end
    save([processed_folder 'tracks'],'tracks')
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

%%
% load(density_file);
% load([processed_folder 'tracks'])
% 
% m_tracks = createMTracks(tracks,Kripple); 
% 
% if(SAVE_RESULTS)
%     if(~OVERWRITE && exist([processed_folder 'm_tracks.mat'],'file') > 0)
%         error('File m_tracs.mat already exists, set OVERWRITE = true to overwrite');
%     end
%     save([processed_folder 'm_tracks'],'m_tracks')
% end
