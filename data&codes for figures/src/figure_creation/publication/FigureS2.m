%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTable.mat']);

xi = [-50:1:150];
zi = [(-pi):pi/16:(pi)];

% PERSISTENT = true caluclates graph for persistent runs
% PERSISTENT = false calculates graph for non-persistent runs
PERSISTENT = true;

%Number of frames to skip prior to calculating inside aggreagte values
%Set to a frame after aggreagte tracking becomes fairly stable. 
STABLE_AGG_CUTOFF = 120;

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = true;

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/Bias/PolarPlot/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = false;

% Sliding window size, in um
WINDOW = [10,pi/4];

% Points at which to center the window for calculation
xi = [-50:1:150];


% Color to use for plotting csgA data
CSGA_COLOR = [1 1 1] * 0.5;

% Color to use for plotting WT data
WT_COLOR = [1 1 1] * 0;

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [20, 20];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

%% Generate or load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(LOAD_PREVIOUS)
    if(PERSISTENT)
        load([SAVE_FOLDER '/BOOT_distance_Persistent']);
        load([SAVE_FOLDER '/BOOT_duration_Persistent'])
    else
        load([SAVE_FOLDER '/BOOT_distance_NonPersistent']);
        load([SAVE_FOLDER '/BOOT_duration_NonPersistent'])
    end
else
    %% Duration
    T = AllDataTable;
    T = T(T.TSS1 > STABLE_AGG_CUTOFF,:);
    T = T(T.state1 < 3,:);

    Twt = T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:);
    [wt_ud,~,~] = movingmeanxyz(Twt.Dn1,Twt.beta1,Twt.Rd1,WINDOW,xi,zi);

    Tcsga = T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:);
    [csga_ud,~,~] = movingmeanxyz(Tcsga.Dn1,Tcsga.beta1,Tcsga.Rd1,WINDOW,xi,zi);

    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER '/BOOT_distance_Persistent'],'wt_ud','csga_ud');
        else
            save([SAVE_FOLDER '/BOOT_distance_NonPersistent'],'wt_ud','csga_ud');
        end
    end

    %% Duration
    T = AllDataTable;
    T = T(T.TSS1 > STABLE_AGG_CUTOFF,:);
    T = T(T.state1 < 3,:);

    Twt = T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:);
    [wt_ut,~,~] = movingmeanxyz(Twt.Dn1,Twt.beta1,Twt.Rt1,WINDOW,xi,zi);

    Tcsga = T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:);
    [csga_ut,~,~] = movingmeanxyz(Tcsga.Dn1,Tcsga.beta1,Tcsga.Rt1,WINDOW,xi,zi);

    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER '/BOOT_duration_Persistent'],'wt_ut','csga_ut');
        else
            save([SAVE_FOLDER '/BOOT_duration_NonPersistent'],'wt_ut','csga_ut');
        end
    end
end

%%
%%% WT
figure, hold on
    polarplot3d(wt_ut,'AngularRange',zi,'RadialRange',xi - xi(1));

    ax = gca
    ax.XTick = -(xi(end) - xi(1)):50:(xi(end) - xi(1));
    ax.XTickLabel = [(150:-50:-50) (0:50:150)];

    ax.YTick = -(xi(end) - xi(1)):50:(xi(end) - xi(1));
    ax.YTickLabel = [(150:-50:-50) (0:50:150)];

    h = colorbar;
    %ylabel(h, 'Run Distance (\mum)')
    h.FontSize = 10;

    box on
    view(2)
    if(PUB_READY)
        xlabel(' ')
        ylabel(' ')
    else
        xlabel('Distance to Agg (\mum)')
        ylabel('Distance to Agg (\mum)')
        %clabel('Run Duration (min)')
        title('(A) WT')
    end

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER '/WT_duration_persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER '/WT_duration_nonpersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end

%%% CsgA
figure, hold on
    polarplot3d(csga_ut,'AngularRange',zi,'RadialRange',xi - xi(1));

    ax = gca
    ax.XTick = -(xi(end) - xi(1)):50:(xi(end) - xi(1));
    ax.XTickLabel = [(150:-50:-50) (0:50:150)];

    ax.YTick = -(xi(end) - xi(1)):50:(xi(end) - xi(1));
    ax.YTickLabel = [(150:-50:-50) (0:50:150)];

    h = colorbar;
    %ylabel(h, 'Run Distance (\mum)')
    h.FontSize = 10;
    
    box on
    view(2)
    if(PUB_READY)
        xlabel(' ')
        ylabel(' ')
    else
        xlabel('Distance to Agg (\mum)')
        ylabel('Distance to Agg (\mum)')
        %clabel('Run Duration (min)')
        title('(C) csgA')
    end

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER '/csga_duration_persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER '/csga_duration_nonpersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end

%%% CsgA Duration
% subplot(2,2,4)
%     polarplot3d(csga_ut,'AngularRange',zi,'RadialRange',xi - xi(1));
% 
%     ax = gca
%     ax.XTick = -(xi(end) - xi(1)):50:(xi(end) - xi(1));
%     ax.XTickLabel = [(150:-50:-50) (0:50:150)];
% 
%     ax.YTick = -(xi(end) - xi(1)):50:(xi(end) - xi(1));
%     ax.YTickLabel = [(150:-50:-50) (0:50:150)];
% 
%     h = colorbar;
%     ylabel(h, 'Run Duration (min)')
%     h.FontSize = 10;
% 
%     xlabel('Distance to Agg (\mum)')
%     ylabel('Distance to Agg (\mum)')
%     title('(D) csgA')
%     view(2)
