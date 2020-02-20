%%%
% Creates a Figure
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTableCsgAinCsgA.mat']);
AllDataTableCsgAinCsgA = AllDataTable;

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTable.mat']);
AllDataTable = AllDataTable;

% How many bootstrap samples
% 1000 for publication, set to 10 to keep things moving quickly for testing
NBOOT = 1000;

%Number of frames to skip prior to calculating inside aggreagte values
%Set to a frame after aggreagte tracking becomes fairly stable. 
STABLE_AGG_CUTOFF = 120;

% Sliding window size, in frames
WINDOW = 20;

% Inside aggreagte density cutoff (cells/um^2)
DENSITY_CUTOFF = 2.3;

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = true;

SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/Time/State/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

CSGA_COLOR = [1 0 0];

CSGA_IN_CSGA_COLOR = [0 0 1];

WT_COLOR = [0 0 0];

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [14, 13];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

%% Data
if(LOAD_PREVIOUS)
    load([SAVE_FOLDER 'BOOTinTime']);
else

    %
    % CsgA/WT in WT
    %

    %P(stop)
    T = AllDataTable;
    [wt.outside.P3,wt.inside.P3] = FigureCreation.Development.bootstrap_state_transitions_in_time(T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:),NBOOT,STABLE_AGG_CUTOFF,WINDOW,DENSITY_CUTOFF);    
    [csgA.outside.P3,csgA.inside.P3] = FigureCreation.Development.bootstrap_state_transitions_in_time(T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:),NBOOT,STABLE_AGG_CUTOFF,WINDOW,DENSITY_CUTOFF);    
    
    %
    % CsgA in CsgA
    %
    T = AllDataTableCsgAinCsgA;
    [csgAinCsgA.P3] = FigureCreation.bootstrap_state_transitions_in_time(T(ismember(cellstr(T.movie),keys(CSGA_IN_CSGA_DEV_FOLDERS)),:),NBOOT,WINDOW);    

    if(SAVE)
            save([SAVE_FOLDER 'BOOTinTime'],'wt','csgA','csgAinCsgA');
    end
end

%% Pub Figure: Transition Probability
figure, hold on 
    clear h
    FigureCreation.plot_boot(wt.outside.P3.xi / 2 / 60,wt.outside.P3,'-','cmap',WT_COLOR,'alpha')
    FigureCreation.plot_boot(wt.inside.P3.xi / 2 / 60,wt.inside.P3,'--','cmap',WT_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgA.outside.P3.xi / 2 / 60,csgA.outside.P3,'-','cmap',CSGA_COLOR,'alpha')
    FigureCreation.plot_boot(csgA.inside.P3.xi / 2 / 60,csgA.inside.P3,'--','cmap',CSGA_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgAinCsgA.P3.xi / 2 / 60,csgAinCsgA.P3,'-','cmap',CSGA_IN_CSGA_COLOR,'alpha')


    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('P(stop)')
        xlabel('Time (hr)')
    end

    xlim auto
    ax = gca;
    ylim([0 1])
    ax.YTick = 0:0.2:1;
    xlim([0 5])
    ax.XTick = 0:1:5;
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER 'Pstop'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end
