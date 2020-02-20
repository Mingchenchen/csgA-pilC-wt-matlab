%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTableCsgAinCsgA.mat']);
AllDataTableCsgAinCsgA = AllDataTable;

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

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/NeighborAlignment/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

% Sliding window size, in frames
WINDOW = 20;

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [30, 30];

CSGA_COLOR = [1 0 0];

CSGA_IN_CSGA_COLOR = [0 0 1];

WT_COLOR = [0 0 0];

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [20, 20];

%% Data
if(LOAD_PREVIOUS)
    if(PERSISTENT)
        load([SAVE_FOLDER 'BOOTNeighborAlign_Persistent']);
    else
        load([SAVE_FOLDER 'BOOTNeighborAlign_NonPersistent']);
    end
else
    T = AllDataTable;

    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end

    wt = FigureCreation.Development.bootstrap_neighbor_alignment(T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:),WINDOW,NBOOT);
    csgA =  FigureCreation.Development.bootstrap_neighbor_alignment(T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:),WINDOW,NBOOT);
    
    T = AllDataTableCsgAinCsgA;
    
    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end
    
    T = T(T.TSS1 <= max(AllDataTable.TSS1),:);
    
    csgAincsgA =  FigureCreation.Development.bootstrap_neighbor_alignment(T(ismember(cellstr(T.movie),keys(CSGA_IN_CSGA_DEV_FOLDERS)),:),WINDOW,NBOOT);
    
    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER 'BOOTNeighborAlign_Persistent'],'wt','csgA','csgAincsgA');
        else
            save([SAVE_FOLDER 'BOOTNeighborAlign_NonPersistent'],'wt','csgA','csgAincsgA');
        end
    end
end

%% Pub Figure: Distance
figure, hold on 
    clear h
    h(1) = FigureCreation.plot_boot(wt.xi/2/60,wt,'-','cmap',WT_COLOR,'alpha');
    h(2) = FigureCreation.plot_boot(wt.xi/2/60,csgA,'-','cmap',CSGA_COLOR,'alpha');
    h(3) = FigureCreation.plot_boot(wt.xi/2/60,wt.rand,'--','cmap',WT_COLOR,'alpha');

    h(4) = FigureCreation.plot_boot(csgAincsgA.xi/2/60,csgAincsgA,'-','cmap',CSGA_IN_CSGA_COLOR,'alpha');
    
    if(PUB_READY)
        xlabel(' ')
        ylabel(' ')
    else
        ylabel('Alignment')
        xlabel('Time (hr)')
    end
    
    xlim auto
    ax = gca;
    ylim([-1 1])
    xlim([0, 5]);
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'NeighborAlignemnt_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'NeighborAlignemnt_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end


    