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

%Number of frames to skip prior to calculating inside aggreagte values
%Set to a frame after aggreagte tracking becomes fairly stable. 
STABLE_AGG_CUTOFF = 120;

% Inside aggreagte density cutoff (cells/um^2)
DENSITY_CUTOFF = 2.3;

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = true;

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/Bias/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

% Sliding window size, in um
WINDOW = 10;

% Points at which to center the window for calculation
xi = [-50:1:150];


% Color to use for plotting csgA data
CSGA_COLOR = [1 0 0];

% Color to use for plotting WT data
WT_COLOR = [0 0 0];

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [14, 13];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

%% Data
if(LOAD_PREVIOUS)
    if(PERSISTENT)
        load([SAVE_FOLDER 'BOOTinBeta_Persistent']);
    else
        load([SAVE_FOLDER 'BOOTinBeta_NonPersistent']);
    end
else
    T = AllDataTable;
    T = T(T.TSS1 > STABLE_AGG_CUTOFF,:);
    
    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end

    wt = FigureCreation.Development.bootstrap_runs_in_beta(xi,T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:),WINDOW,NBOOT);
    csgA =  FigureCreation.Development.bootstrap_runs_in_beta(xi,T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:),WINDOW,NBOOT);

    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER 'BOOTinBeta_Persistent'],'wt','csgA');
        else
            save([SAVE_FOLDER 'BOOTinBeta_NonPersistent'],'wt','csgA');
        end
    end
end

%% Pub Figure: Distance
figure, hold on 
    clear h
    h(1) = FigureCreation.plot_boot(xi,wt.towards.Rd,'-','cmap',color_chooser(1,cmap),'alpha');
    h(2) = FigureCreation.plot_boot(xi,wt.away.Rd,'--','cmap',color_chooser(1,cmap),'alpha');
    
    h(3) = FigureCreation.plot_boot(xi,csgA.towards.Rd,'-','cmap',color_chooser(2,cmap),'alpha');
    h(4) = FigureCreation.plot_boot(xi,csgA.away.Rd,'--','cmap',color_chooser(2,cmap),'alpha');
    
    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('Distance (\mum)')
        xlabel('Distance to Nearest Agg. (\mum)')
    end
    
    xlim auto
    ax = gca;
    ylim([0 42])
    xlim([xi(1),xi(end)])
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'Distance_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'Distance_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end
    
%% Pub Figure: Duration
figure, hold on 
    clear h
    h(1) = FigureCreation.plot_boot(xi,wt.towards.Rt,'-','cmap',color_chooser(1,cmap),'alpha');
    h(2) = FigureCreation.plot_boot(xi,wt.away.Rt,'--','cmap',color_chooser(1,cmap),'alpha');
    
    h(3) = FigureCreation.plot_boot(xi,csgA.towards.Rt,'-','cmap',color_chooser(2,cmap),'alpha');
    h(4) = FigureCreation.plot_boot(xi,csgA.away.Rt,'--','cmap',color_chooser(2,cmap),'alpha');
    
    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('Duration (min)')
        xlabel('Distance to Nearest Agg. (\mum)')
    end
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([xi(1),xi(end)])
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'Duration_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'Duration_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end
    
%% Pub Figure: Speed
figure, hold on 
    clear h
    h(1) = FigureCreation.plot_boot(xi,wt.towards.Rs,'-','cmap',color_chooser(1,cmap),'alpha');
    h(2) = FigureCreation.plot_boot(xi,wt.away.Rs,'--','cmap',color_chooser(1,cmap),'alpha');
    
    h(3) = FigureCreation.plot_boot(xi,csgA.towards.Rs,'-','cmap',color_chooser(2,cmap),'alpha');
    h(4) = FigureCreation.plot_boot(xi,csgA.away.Rs,'--','cmap',color_chooser(2,cmap),'alpha');
    
    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('Speed (\mum/min)')
        xlabel('Distance to Nearest Agg. (\mum)')
    end
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([xi(1),xi(end)])
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'Speed_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'Speed_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end


    