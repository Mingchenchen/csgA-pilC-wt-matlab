%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
addpath('Libraries/')
run('ENVS.m')

SIM_NAME = 'PUB_A1';     %Used to name output files
prerun_length =180;
COMAPRE_WITH = 'WT';

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [30, 30];

xi_agg = floor([1 2.75 4.5] * 60 * 2);

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/simulation/development/' SIM_NAME '/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end
   
%% Plot Sim Frac Curve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/simFracCurves.mat']);
load([PROJECT_BASENAME '/data/processed/Development/frac_curves.mat'])

CsgA = cell2mat(values(frac_curves,keys(CSGA_DEV_FOLDERS))')
ucsga = mean(CsgA,1)' .* 100;
scsga = std(CsgA,[],1)' .* 100;

WT = cell2mat(values(frac_curves,keys(WT_DEV_FOLDERS))')
uwt = mean(WT,1)' .* 100;
swt = std(WT,[],1)' .* 100;

sim_frac_mean = mean(simFracCurves(prerun_length:end,:),2) .* 100;
sim_frac_std = std(simFracCurves(prerun_length:end,:),[],2) .* 100;

figure, clf, hold on;
    xi = (1:length(sim_frac_mean)) / 2 / 60
    boundedline(xi,sim_frac_mean,[sim_frac_std,sim_frac_std],'--','cmap',color_chooser(2),'alpha')
    
    xi = (1:length(ucsga)) / 2 / 60;
    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1),'alpha') %Compare with WT
%    boundedline(xi,ucsga,[scsga,scsga],'-','cmap',color_chooser(2),'alpha') %Compare with CsgA
    
    xlim auto
    ax = gca;
    ylim([0 ax.YLim(2)])
    ylabel('% Cells Inside Aggregates');
    xlabel('Time (hr)')
    box on;


if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER 'frac_cells'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end
    
%% Plot Aggreagte Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([PROJECT_BASENAME, '/data/processed/Development/AllData.mat']);
props = FigureCreation.Development.AggPlots();
WTinWT =  props.combine_agg_props(AllData(1:3));
CsgAinWT =  props.combine_agg_props(AllData(4:end));

load([PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/AggProps.mat'])
SIM = AggProps;

COMPARE_WITH = 'CSGA';

if(strcmp(COMPARE_WITH,'WT'))
    REAL = WTinWT;
elseif(strcmp(COMPARE_WITH,'CSGA'))
    REAL = CsgAinWT;
end

%% Area
figure
    props.plot_property(REAL,SIM,xi_agg,'area','jitter')
    ylabel('Area (\mum/min)')
    
    l = legend({'Real','Sim'},'Location','northwest')
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end

if SAVE
    saveFigures('SaveAs',[SAVE_FOLDER 'Area'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%% Cell Density
figure, 
    props.plot_property(REAL,SIM,xi_agg,'mean_density','jitter')
    ylabel('Cell Density (cells/\mum)')
    grid on
    
    l = legend({'Real','Sim'},'Location','northwest')
    legend boxoff
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end
    
if SAVE
    saveFigures('SaveAs',[SAVE_FOLDER 'CellDensity'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%% Agg Counts
figure,
    props.agg_count(REAL,SIM,xi_agg,'jitter')
    ylabel('Agg Count')
    grid on
    
    l = legend({'Real','Sim'},'Location','southeast')
    legend boxoff
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end

if SAVE
    saveFigures('SaveAs',[SAVE_FOLDER 'AggCounts'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end