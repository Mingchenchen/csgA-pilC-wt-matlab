%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataTableWrappedStops.mat']);

% How many bootstrap samples
% 1000 for publication, set to 10 to keep things moving quickly for testing
NBOOT = 100;

% Sliding window size, in frames
WINDOW = 20;

% Inside aggreagte density cutoff (cells/um^2)
DENSITY_CUTOFF = 2.3;

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = false;

SAVE = false;

SAVE_FOLDER = [PROJECT_BASENAME '/reports/rippleing/RunBehaviors/Wavelength/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = false;

xi = [-50:1:150];

CSGA_SETS = [4:6];
CSGA_COLOR = [1 1 1] * 0.5;
CSGA_MARKER = '--';

WT_SETS = [1:3];
WT_COLOR = [1 1 1] * 0;
WT_MARKER = '-';

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [20, 20];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

%% Data
if(LOAD_PREVIOUS)
    load([SAVE_FOLDER 'BOOTinTime']);
else
    T = AllDataTable;

    wt = FigureCreation.Rippleing.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(WT_RIPPLE_FOLDERS)),:),NBOOT,WINDOW);
    csgA = FigureCreation.Rippleing.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(CSGA_RIPPLE_FOLDERS)),:),NBOOT,WINDOW);

    if(SAVE)
        save([SAVE_FOLDER 'BOOTinTime'],'wt','csgA');
    end
end

%%
load([PROJECT_BASENAME, '/data/processed/Rippleing/' 'RipplePeriods']) % loads 'periods','centers'
data = cell2mat(values(periods)) ./ 2;
period_u = mean(data,2);
period_s = std(data,[],2);


%% Pub Figure: Distance
figure, hold on 
    clear h
    FigureCreation.plot_boot(wt.xi / 2 / 60,wt.Rd,'-','cmap',color_chooser(1,cmap),'alpha')

    FigureCreation.plot_boot(csgA.xi / 2 / 60,csgA.Rd,'-','cmap',color_chooser(2,cmap),'alpha')
    
    ylabel('Distance (\mum)')
    xlabel('Time (hr)')
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([0 3.5])
    
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
    FigureCreation.plot_boot(wt.xi / 2 / 60,wt.Rt,'-','cmap',color_chooser(1,cmap),'alpha')
    
    FigureCreation.plot_boot(csgA.xi / 2 / 60,csgA.Rt,'-','cmap',color_chooser(2,cmap),'alpha')
    
    hold on;
    plot(centers,period_u)
	plot(centers,period_u + period_s,'--k')
    plot(centers,period_u - period_s,'--k')
    
    ylabel('Duration (min)')
    xlabel('Time (hr)')
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([0 3.5])
    
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
    FigureCreation.plot_boot(wt.xi / 2 / 60,wt.Rs,'-','cmap',color_chooser(1,cmap),'alpha')

    FigureCreation.plot_boot(csgA.xi / 2 / 60,csgA.Rs,'--','cmap',color_chooser(2,cmap),'alpha')
    
    ylabel('Speed (\mum/min)')
    xlabel('Time (hr)')
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    xlim([0 3.5])
    
    ax.FontSize = FONT_SIZE;
    
    box on;

if (SAVE)
    if(PERSISTENT)
        saveFigures('SaveAs',[SAVE_FOLDER 'Speed_Persistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    else
        saveFigures('SaveAs',[SAVE_FOLDER 'Speed_NonPersistent'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
end


    