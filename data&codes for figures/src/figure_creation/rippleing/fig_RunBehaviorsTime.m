%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Rippleing/AllDataTable.mat']);

% PERSISTENT = true caluclates graph for persistent runs
% PERSISTENT = false calculates graph for non-persistent runs
PERSISTENT = true;

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
SAVE = true;

SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/rippleing/RunBehaviors/Time/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Load previous bootstrapped data
LOAD_PREVIOUS = true;

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
    if(PERSISTENT)
        load([SAVE_FOLDER 'BOOTinTime_Persistent']);
    else
        load([SAVE_FOLDER 'BOOTinTime_NonPersistent']);
    end
else
    T = AllDataTable;

    if(PERSISTENT)
        T = T(T.state1 < 3,:);
    else
        T = T(T.state1 == 3,:);
    end

    wt = FigureCreation.Rippleing.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(WT_RIPPLE_FOLDERS)),:),NBOOT,WINDOW);
    csgA = FigureCreation.Rippleing.bootstrap_runs_in_time(T(ismember(cellstr(T.movie),keys(CSGA_RIPPLE_FOLDERS)),:),NBOOT,WINDOW);

    if(SAVE)
        if(PERSISTENT)
            save([SAVE_FOLDER 'BOOTinTime_Persistent'],'wt','csgA');
        else
            save([SAVE_FOLDER 'BOOTinTime_NonPersistent'],'wt','csgA');
        end
    end
end

%%
load([PROJECT_BASENAME, '/data/processed/Rippleing/' 'RipplePeriods']) % loads 'periods','centers'
data = periods / 2;
period_u = mean(data,2);
period_s = std(data,[],2);

%%
NPoints = 2^10;
MAP = [WT_RIPPLE_FOLDERS; CSGA_RIPPLE_FOLDERS];
all_wavelengths = [];
for k = keys(MAP)
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{1})];
    load([folder 'wavelengths.mat'])
    
    Nframes = size(wavelengths,2);
    
    wl = NPoints ./ (0:(NPoints / 2 + 1));
    
    total_magnitude = sum(wavelengths,1);
    D = wavelengths ./ total_magnitude;
    D = wl' .* D;
    major_wavelength = sum(D(2:end,:));
    
    total_magnitude  = (total_magnitude - min(total_magnitude)) ./ range(total_magnitude);

    %Crude conversion to um
    major_wavelength = major_wavelength .* (986 / 2^10 + 739 / 2^10) / 2;
    
    all_wavelengths = [all_wavelengths;  major_wavelength(200:400)];
end
wavelength_u = mean(all_wavelengths) / 2;
wavelength_s = std(all_wavelengths);
wavelength_xi = (200:400) / 2 / 60;

%% Pub Figure: Distance
figure, hold on 
    clear h
    FigureCreation.plot_boot(wt.xi / 2 / 60,wt.Rd,'-','cmap',color_chooser(1,cmap),'alpha')

    FigureCreation.plot_boot(csgA.xi / 2 / 60,csgA.Rd,'-','cmap',color_chooser(2,cmap),'alpha')
    
    ylabel('Distance (\mum)')
    xlabel('Time (hr)')
    
    xlim auto
    ax = gca;
    ylim([0 60])
    xlim([0 3.5])
    
    ax.FontSize = FONT_SIZE;
    
    hold on;
    yyaxis right
    h = boundedline(wavelength_xi,wavelength_u,[wavelength_s',wavelength_s'],'-b')
    
    ax2 = gca;
    ylim(ax2,[0,60]);
    ylabel('1/2 Wavelength (\mum)')
    
    
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
    
    ax1 = gca;
    xlim auto
    ylim(ax1,[0 ax1.YLim(2)])
    ax1.FontSize = FONT_SIZE;
    ylims = ax1.YLim;
    
    ylabel('Duration (min)')
    xlabel('Time (hr)')
    
    hold on;
    yyaxis right
    h = boundedline(centers,period_u,[period_s,period_s],'-b')
    
    ax2 = gca;
    ylim(ax2,ylims);
    ylabel('1/2 Wave Period (min)')
    

    ylim(ax2,[0 ax1.YLim(2)])
    xlim([0 3.5])
    
    
    
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


    