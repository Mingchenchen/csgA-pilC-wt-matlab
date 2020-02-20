%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
addpath('Libraries/')
run('ENVS.m')

xi_agg = floor([1 2.75 4.5] * 60 * 2);

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = true;

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/simulation/development/SupFig_WTDev/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Figure Properties
FONT_SIZE = 12;


%% Panel A
%
SIM_NAME = 'PUB_NEW_1';     %Used to name output files
prerun_length =180;
FIG_SIZE = [13, 13];

load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/simFracCurves.mat']);
load([PROJECT_BASENAME '/data/processed/Development/frac_curves.mat'])

CsgA = cell2mat(values(frac_curves,keys(WT_DEV_FOLDERS))')
ucsga = mean(CsgA,1)' .* 100;
scsga = std(CsgA,[],1)' .* 100;

sim_frac_mean = mean(simFracCurves(prerun_length:end,:),2) .* 100;
sim_frac_std = std(simFracCurves(prerun_length:end,:),[],2) .* 100;

figure, clf, hold on;
    xi = (1:length(sim_frac_mean)) / 2 / 60
    boundedline(xi,sim_frac_mean,[sim_frac_std,sim_frac_std],'--','cmap',color_chooser(2),'alpha')
    
    xi = (1:length(ucsga)) / 2 / 60;
%    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1),'alpha') %Compare with WT
    boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2),'alpha') %Compare with CsgA
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    box on
    ax.XTick = 0:1:5;
    xlim([0, 5]);
    if(PUB_READY)
        ylabel(' ')
        xlabel(' ')
    else
        ylabel('% Cells Inside Aggregates');
        xlabel('Time (hr)')
    end
   
if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER 'frac_cells'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%% Panels B-D
%
AllData = {};
FIG_SIZE = [10, 13];
SIM_NAME = 'PUB_NEW_2';
load([PROJECT_BASENAME, '/data/processed/Development/AllData.mat']);
props = FigureCreation.Development.AggPlots();

REAL =  props.combine_agg_props(AllData(1:3));

load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/AggProps.mat'])
SIM = AggProps;

% Panel B
figure
    props.plot_property(REAL,SIM,xi_agg,'area','jitter')
    
    if(PUB_READY)
        xlabel(' ')
        ylabel(' ')
        legend off
        ax = gca;
        ax.YTick = 0:5e3:2e4;
        ax.YTickLabel = ax.YTick/1e4;
    else
        ylabel('Area (\mum/min)')
        
        l = legend({'Real','Sim'},'Location','northwest')
        entries = l.EntryContainer.Children;
        for i=1:numel(entries)
             entries(i).Icon.Transform.Children.Children(1).Size = 25;
        end
    end
    box on
    
    if SAVE
        saveFigures('SaveAs',[SAVE_FOLDER 'Area'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end

% Panel C
figure
    props.plot_property(REAL,SIM,xi_agg,'mean_density','jitter')
    
    if(PUB_READY)
        xlabel(' ')
        ylabel(' ')
        legend off
    else
        ylabel('Cell Density (cells/\mum)')
        
        l = legend({'Real','Sim'},'Location','northwest')
        entries = l.EntryContainer.Children;
        for i=1:numel(entries)
             entries(i).Icon.Transform.Children.Children(1).Size = 25;
        end
    end
    box on
    
    if SAVE
        saveFigures('SaveAs',[SAVE_FOLDER 'Density'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end

% Panel D
figure,
    props.agg_count(REAL,SIM,xi_agg,'jitter')

    if(PUB_READY)
        xlabel(' ')
        ylabel(' ')
        legend off
    else
        ylabel('Agg Count')
        
        l = legend({'Real','Sim'},'Location','southeast')
        entries = l.EntryContainer.Children;
        for i=1:numel(entries)
             entries(i).Icon.Transform.Children.Children(1).Size = 25;
        end
    end
    ylim([0, 15])
    box on
    
    if SAVE
        saveFigures('SaveAs',[SAVE_FOLDER 'AggCounts'],'Style','none','Formats',{'png','fig'},'Size',FIG_SIZE);
    end
    
%% Panels E-H
prerun_length =180;
FIG_SIZE = [8, 8];

% Panel E
SIM_NAME = 'PUB_NEW_10';
load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/simFracCurves.mat']);

sim_frac_mean = mean(simFracCurves(prerun_length:end,:),2) .* 100;
sim_frac_std = std(simFracCurves(prerun_length:end,:),[],2) .* 100;

figure, clf, hold on;
    xi = (1:length(sim_frac_mean)) / 2 / 60
    boundedline(xi,sim_frac_mean,[sim_frac_std,sim_frac_std],'--','cmap',color_chooser(2),'alpha')
    
    xi = (1:length(ucsga)) / 2 / 60;
%    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1),'alpha') %Compare with WT
    boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2),'alpha') %Compare with CsgA
    
    xlim auto
    ax = gca;
    ylim([0 65])
    box on
    grid on
    ax.XTick = 0:1:5;
    ax.YTick = 0:10:65;
    xlim([0, 5]);
    if(~PUB_READY)
        ylabel('% Cells Inside Aggregates');
        xlabel('Time (hr)')
    end
   
if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '/' SIM_NAME '-frac_cells'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%Panel F
SIM_NAME = 'PUB_NEW_11';
load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/simFracCurves.mat']);

sim_frac_mean = mean(simFracCurves(prerun_length:end,:),2) .* 100;
sim_frac_std = std(simFracCurves(prerun_length:end,:),[],2) .* 100;

figure, clf, hold on;
    xi = (1:length(sim_frac_mean)) / 2 / 60
    boundedline(xi,sim_frac_mean,[sim_frac_std,sim_frac_std],'--','cmap',color_chooser(2),'alpha')
    
    xi = (1:length(ucsga)) / 2 / 60;
%    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1),'alpha') %Compare with WT
    boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2),'alpha') %Compare with CsgA
    
    xlim auto
    ax = gca;
    ylim([0 65])
    box on
    grid on
    ax.XTick = 0:1:5;
    ax.YTick = 0:10:65;
    ax.YTickLabel = '  ';
    xlim([0, 5]);
    if(~PUB_READY)
        ylabel('% Cells Inside Aggregates');
        xlabel('Time (hr)')
    end
   
if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '/' SIM_NAME '-frac_cells'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%Panel G
SIM_NAME = 'PUB_NEW_4';
load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/simFracCurves.mat']);

sim_frac_mean = mean(simFracCurves(prerun_length:end,:),2) .* 100;
sim_frac_std = std(simFracCurves(prerun_length:end,:),[],2) .* 100;

figure, clf, hold on;
    xi = (1:length(sim_frac_mean)) / 2 / 60
    boundedline(xi,sim_frac_mean,[sim_frac_std,sim_frac_std],'--','cmap',color_chooser(2),'alpha')
    
    xi = (1:length(ucsga)) / 2 / 60;
%    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1),'alpha') %Compare with WT
    boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2),'alpha') %Compare with CsgA
    
    xlim auto
    ax = gca;
    ylim([0 65])
    box on
    grid on
    ax.XTick = 0:1:5;
    ax.YTick = 0:10:65;
    ax.YTickLabel = '  ';
    xlim([0, 5]);
    if(~PUB_READY)
        ylabel('% Cells Inside Aggregates');
        xlabel('Time (hr)')
    end
   
if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '/' SIM_NAME '-frac_cells'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%Panel H
SIM_NAME = 'PUB_NEW_3';
load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/simFracCurves.mat']);

sim_frac_mean = mean(simFracCurves(prerun_length:end,:),2) .* 100;
sim_frac_std = std(simFracCurves(prerun_length:end,:),[],2) .* 100;

figure, clf, hold on;
    xi = (1:length(sim_frac_mean)) / 2 / 60
    boundedline(xi,sim_frac_mean,[sim_frac_std,sim_frac_std],'--','cmap',color_chooser(2),'alpha')
    
    xi = (1:length(ucsga)) / 2 / 60;
%    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1),'alpha') %Compare with WT
    boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2),'alpha') %Compare with CsgA
    
    xlim auto
    ax = gca;
    ylim([0 65])
    box on
    grid on
    ax.XTick = 0:1:5;
    ax.YTick = 0:10:65;
    ax.YTickLabel = '  ';
    xlim([0, 5]);
    if(~PUB_READY)
        ylabel('% Cells Inside Aggregates');
        xlabel('Time (hr)')
    end
   
if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '/' SIM_NAME '-frac_cells'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end
