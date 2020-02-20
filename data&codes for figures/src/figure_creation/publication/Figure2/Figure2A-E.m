%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTableCsgAinCsgA.mat']);
AllDataTableCsgAinCsgA = AllDataTable;

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTable.mat']);
AllDataTable = AllDataTable;

% PERSISTENT = true caluclates graph for persistent runs
% PERSISTENT = false calculates graph for non-persistent runs
PERSISTENT = true;

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
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/RunBehaviors/Time/DistanceDurationSpeed/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

xi = [-50:1:150];

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
    if(PERSISTENT)
        load([SAVE_FOLDER 'BOOTinTime_Persistent']);
    else
        load([SAVE_FOLDER 'BOOTinTime_NonPersistent']);
    end
else

    %
    % CsgA/WT in WT
    %
    T = AllDataTable;
    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end

    %Run Speed, Distance, and Time
    wt = FigureCreation.Development.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:),NBOOT,STABLE_AGG_CUTOFF,WINDOW,DENSITY_CUTOFF);
    csgA =  FigureCreation.Development.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:),NBOOT,STABLE_AGG_CUTOFF,WINDOW,DENSITY_CUTOFF);

    %
    % CsgA in CsgA
    %
    T = AllDataTableCsgAinCsgA;
    
    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end
    
    %Run Speed, Distance, and Time
    csgAinCsgA =  FigureCreation.Rippleing.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(CSGA_IN_CSGA_DEV_FOLDERS)),:),NBOOT,WINDOW);
    
    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER '/BOOTinTime_Persistent'],'wt','csgA','csgAinCsgA');
        else
            save([SAVE_FOLDER '/BOOTinTime_NonPersistent'],'wt','csgA','csgAinCsgA');
        end
    end
end

%% Pub Figure: Distance
figure, hold on 
    clear h
    FigureCreation.plot_boot(wt.outside.xi / 2 / 60,wt.outside.Rd,'-','cmap',WT_COLOR,'alpha')
    FigureCreation.plot_boot(wt.inside.xi / 2 / 60,wt.inside.Rd,'--','cmap',WT_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgA.outside.xi / 2 / 60,csgA.outside.Rd,'-','cmap',CSGA_COLOR,'alpha')
    FigureCreation.plot_boot(csgA.inside.xi / 2 / 60,csgA.inside.Rd,'--','cmap',CSGA_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgAinCsgA.xi / 2 / 60,csgAinCsgA.Rd,'-','cmap',CSGA_IN_CSGA_COLOR,'alpha')
    
    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('Distance (\mum)');
        xlabel('Time (hr)');
    end
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([0 5])
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER '/Distance_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER '/Distance_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end
    
%% Pub Figure: Duration
figure, hold on 
    clear h
    FigureCreation.plot_boot(wt.outside.xi / 2 / 60,wt.outside.Rt,'-','cmap',WT_COLOR,'alpha')
    FigureCreation.plot_boot(wt.inside.xi / 2 / 60,wt.inside.Rt,'--','cmap',WT_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgA.outside.xi / 2 / 60,csgA.outside.Rt,'-','cmap',CSGA_COLOR,'alpha')
    FigureCreation.plot_boot(csgA.inside.xi / 2 / 60,csgA.inside.Rt,'--','cmap',CSGA_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgAinCsgA.xi / 2 / 60,csgAinCsgA.Rt,'-','cmap',CSGA_IN_CSGA_COLOR,'alpha')

    
    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('Duration (min)');
        xlabel('Time (hr)');
    end
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2) + 1])
    xlim([0 5])
    
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
    FigureCreation.plot_boot(wt.outside.xi / 2 / 60,wt.outside.Rs,'-','cmap',WT_COLOR,'alpha')
    FigureCreation.plot_boot(wt.inside.xi / 2 / 60,wt.inside.Rs,'--','cmap',WT_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgA.outside.xi / 2 / 60,csgA.outside.Rs,'-','cmap',CSGA_COLOR,'alpha')
    FigureCreation.plot_boot(csgA.inside.xi / 2 / 60,csgA.inside.Rs,'--','cmap',CSGA_COLOR,'alpha')
    
    FigureCreation.plot_boot(csgAinCsgA.xi / 2 / 60,csgAinCsgA.Rs,'-','cmap',CSGA_IN_CSGA_COLOR,'alpha')


    if(PUB_READY)
        xlabel(' ');
        ylabel(' ');
    else
        ylabel('Speed (\mum/min)');
        xlabel('Time (hr)');
    end
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([0 5])
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'Speed_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'Speed_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end

%% Data Table
text = ['   Cells              |    Mean Run Distance  \n', ...
        '----------------------|-----------------------\n', ...
        ' Rescued CsgA outside |         %0.2f      \n',  ...
        '      WT outside      |         %0.3f      \n', ...
        '     CsgA             |         %0.3f      \n', ...
        ' Rescued CsgA inside  |         %0.2f      \n',  ...
        '      WT inside       |         %0.3f       \n'];

T = AllDataTable;  
Trcsga = T(ismember(cellstr(T.movie),keys(CSGA_DEV_FOLDERS)),:);
Twt = T(ismember(cellstr(T.movie),keys(WT_DEV_FOLDERS)),:);

T = AllDataTableCsgAinCsgA;
Tcsga = T(ismember(cellstr(T.movie),keys(CSGA_IN_CSGA_DEV_FOLDERS)),:);

R_csgA_outside = mean(Trcsga.Rd1(Trcsga.rho1 <= DENSITY_CUTOFF));
R_csgA_inside = mean(Trcsga.Rd1(Trcsga.rho1 > DENSITY_CUTOFF));
csgA = mean(Tcsga.Rd1);
wt_outside = mean(Twt.Rd1(Twt.rho1 <= DENSITY_CUTOFF));
wt_inside = mean(Twt.Rd1(Twt.rho1 > DENSITY_CUTOFF));
    
sprintf(text,R_csgA_outside,wt_outside, ...
            csgA, ...
            R_csgA_inside,wt_inside)


    