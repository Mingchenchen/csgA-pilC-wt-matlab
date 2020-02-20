%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTable.mat']);

% PERSISTENT = true caluclates graph for persistent runs
% PERSISTENT = false calculates graph for non-persistent runs
PERSISTENT = true;

% How many bootstrap samples
% 1000 for publication, set to 10 to keep things moving quickly for testing
NBOOT = 1000;

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = true;

%Number of frames to skip prior to calculating inside aggreagte values
%Set to a frame after aggreagte tracking becomes fairly stable. 
STABLE_AGG_CUTOFF = 160;

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/AggAlignment/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

% Sliding window size, um
WINDOW = 10;

xi = [-30:1:100];

CSGA_COLOR = [1 0 0];

WT_COLOR = [0 0 0];

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [20, 20];


%% Data
if(LOAD_PREVIOUS)
    if(PERSISTENT)
        load([SAVE_FOLDER 'BOOTAggAlign_Persistent']);
    else
        load([SAVE_FOLDER 'BOOTAggAlign_NonPersistent']);
    end
else
    T = AllDataTable;
    T = T(T.frame > STABLE_AGG_CUTOFF,:);

    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end

    wt = FigureCreation.Development.bootstrap_angle_to_agg(xi,T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:),WINDOW,NBOOT);
    csgA = FigureCreation.Development.bootstrap_angle_to_agg(xi,T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:),WINDOW,NBOOT);

    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER 'BOOTAggAlign_Persistent'],'wt','csgA');
        else
            save([SAVE_FOLDER 'BOOTAggAlign_NonPersistent'],'wt','csgA');
        end
    end
end

%% Pub Figure: Distance
figure, hold on 
    clear h
    h(1) = FigureCreation.plot_boot(xi,wt,'-','cmap',WT_COLOR,'alpha');
    h(2) = FigureCreation.plot_boot(xi,csgA,'-','cmap',CSGA_COLOR,'alpha');
    h(3) = FigureCreation.plot_boot(xi,wt.rand,'--','cmap',WT_COLOR,'alpha');
    
    if(PUB_READY)
        ylabel(' ')
        xlabel(' ')
    else
        ylabel('Alignment')
        xlabel('Distance to Agg. Boundary (\mum)')
    end
    
    xlim auto
    ax = gca;
    ylim([-1 1])
    xlim([xi(1), xi(end)])
    ax.XTick = xi(1):20:xi(end);
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'AggAlignemnt_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'AggAlignemnt_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end


    