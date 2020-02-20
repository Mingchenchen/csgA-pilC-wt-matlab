%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
addpath('Libraries/')
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataTable.mat']);

% PERSISTENT = true caluclates graph for persistent runs
% PERSISTENT = false calculates graph for non-persistent runs
PERSISTENT = true;

% How many bootstrap samples
% 1000 for publication, set to 10 to keep things moving quickly for testing
NBOOT = 100;

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = false;

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/rippleing/NeighborAlignment/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end
   
% Load previous bootstrapped data
LOAD_PREVIOUS = true;

% Sliding window size, in frames
WINDOW = 20;

% Color to use for plotting csgA data
CSGA_COLOR = [1 1 1] * 0.5;

% Color to use for plotting WT data
WT_COLOR = [1 1 1] * 0;

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [30, 30];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

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

    wt = FigureCreation.Rippleing.bootstrap_neighbor_alignment(T(ismember(cellstr(AllDataTable.movie),keys(WT_RIPPLE_FOLDERS)),:),WINDOW,NBOOT);
    csgA =  FigureCreation.Rippleing.bootstrap_neighbor_alignment(T(ismember(cellstr(AllDataTable.movie),keys(CSGA_RIPPLE_FOLDERS)),:),WINDOW,NBOOT);
    
    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER 'BOOTNeighborAlign_Persistent'],'wt','csgA');
        else
            save([SAVE_FOLDER 'BOOTNeighborAlign_NonPersistent'],'wt','csgA');
        end
    end
end

%% Pub Figure: Alignment
figure, hold on 
    clear h
    h(1) = FigureCreation.plot_boot(wt.xi/2/60,wt,'-','cmap',color_chooser(1,cmap),'alpha');
    h(2) = FigureCreation.plot_boot(csgA.xi/2/60,csgA,'-','cmap',color_chooser(2,cmap),'alpha');
    h(3) = FigureCreation.plot_boot(wt.xi/2/60,wt.rand,'--','cmap',color_chooser(1,cmap),'alpha');
    
    ylabel('Alignment')
    xlabel('Time (hr)')
    
    xlim auto
    ax = gca;
    ylim([-1 1])
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'NeighborAlignemnt_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'NeighborAlignemnt_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end


    