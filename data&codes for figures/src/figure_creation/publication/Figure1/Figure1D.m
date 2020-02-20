addpath('Libraries/Utils')
addpath('Libraries/')
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/frac_curves.mat']);

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
FOR_PUB = true;

% If to save the data and figures, and the folder to save them in
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME 'reports/figures/development/fig_FracCells/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [13, 13];

CSGA_COLOR = [1 0 0];
WT_COLOR = [0 0 0];

% Auto calculated based on above settings
colormap = [WT_COLOR; CSGA_COLOR];

%%
CsgA = cell2mat(values(frac_curves,keys(CSGA_DEV_FOLDERS))');
ucsga = mean(CsgA,1)' .* 100;
scsga = std(CsgA,[],1)' .* 100;

WT = cell2mat(values(frac_curves,keys(WT_DEV_FOLDERS))');
uwt = mean(WT,1)' .* 100;
swt = std(WT,[],1)' .* 100;

figure; hold on;
    xi = (1:length(ucsga)) / 2 / 60;
    boundedline(xi,uwt,[swt,swt],'-','cmap',color_chooser(1,colormap),'alpha');
    boundedline(xi,ucsga,[scsga,scsga],'--','cmap',color_chooser(2,colormap),'alpha');
    
    xlim auto;
    ax = gca;
    ylim([0 ax.YLim(2)]);
    xlim([0 5]);
    ax.XTick = 0:5;
    ax.YTick = 0:10:70;
    ax.FontSize = FONT_SIZE;
    if(~FOR_PUB)
        ylabel('% Cells Inside Aggs.');
        xlabel('Time (hr)');
    else
        ylabel(' ');
        xlabel(' ');
    end
    box on;

if SAVE
    saveFigures('SaveAs',[SAVE_FOLDER 'FracCellsInAggs'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE),'r',1000;
end

%%
text = ['   Cells      |    %% Cells in Aggs \n', ...
        '--------------|-----------------------\n', ...
        ' Rescued CsgA |    %0.2f +/- %0.2f       \n',  ...
        '      WT      |    %0.3f +/- %0.2f       \n'];
        
sprintf(text,ucsga(end),scsga(end), ...
             uwt(end),swt(end))