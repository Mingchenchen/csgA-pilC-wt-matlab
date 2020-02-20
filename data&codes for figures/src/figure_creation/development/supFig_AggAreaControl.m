%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllData.mat']);

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/AggProps/AggAreaWTCsgAControl/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

xi = floor([1 2.75 4.5] * 60 * 2);

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [20, 20];

CSGA_COLOR = [1 0 0];
WT_COLOR = [0 0 0];

% Auto calculated based on above settings
cmap = [WT_COLOR; CSGA_COLOR];

default_colormap = get(0,'DefaultFigureColormap');
set(0, 'DefaultFigureColormap', cmap)

%%
props = FigureCreation.Development.AggPlots();
WTinWT =  props.combine_agg_props(AllData(1:3));
CsgAinWT =  props.combine_agg_props(AllData(4:end));


%% Area
figure
    props.plot_property(WTinWT,CsgAinWT,xi,'area','jitter')
    ylabel('Area (\mum/min)')
    
    l = legend({'WT','CsgA'},'Location','northwest')
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end

if SAVE
    saveFigures('SaveAs',[SAVE_FOLDER 'Area'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%% Cell Density
figure, 
    props.plot_property(WTinWT,CsgAinWT,xi,'mean_density','jitter')
    ylabel('Cell Density (cells/\mum)')
    grid on
    
    l = legend({'WT','CsgA'},'Location','northwest')
    legend boxoff
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end

if SAVE
    saveFigures('SaveAs',[SAVE_FOLDER 'CellDensity'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
end

%% Agg counts
figure,
    props.agg_count(WTinWT,CsgAinWT,xi,'jitter')
    ylabel('Agg Count')
    grid on
    
    l = legend({'WT','CsgA'},'Location','southeast')
    legend boxoff
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end
    
    if SAVE
        saveFigures('SaveAs',[SAVE_FOLDER 'AggCounts'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
    end
    
figure, 
    props.plot_property(REAL,SIM,xi,'meanIntensity','jitter','colormap',cmap)
    ylabel('Cell Density (cells/\mum)')
    grid on
    
    l = legend({'Real','Sim'},'Location','northwest')
    legend boxoff
    entries = l.EntryContainer.Children;
    for i=1:numel(entries)
         entries(i).Icon.Transform.Children.Children(1).Size = 25;
    end

%%
set(0, 'DefaultFigureColormap', default_colormap)