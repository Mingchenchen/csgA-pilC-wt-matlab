%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

load([PROJECT_BASENAME, '/data/processed/Development/AllDataTable.mat']);
AllDataTable = AllDataTable;

AGG_FRAME = 540;

DISTCONV = 0.5136;

SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/development/Figure1/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

FIG_SIZE = [13, 13];

SCALE_BAR.LENGTH = 100; %micrometers
SCALE_BAR.X_MARGIN = 10; %micrometers
SCALE_BAR.Y_MARGIN = 15; %micrometers
SCALE_BAR.WIDTH = 3; %micrometers

%% WT
WT_TRAJECTORIES_FNAME = [PROJECT_BASENAME '/data/processed/Development/WT in WT/10282015/'];
WINDOW = [0, 500, 0, 500] % [xmin, xmax, ymin, ymax]

load([WT_TRAJECTORIES_FNAME, '/m_tracks.mat']);
load([WT_TRAJECTORIES_FNAME, '/agg_tracks_with_area.mat']);

tracks = subStruct(m_tracks,m_tracks.x >= WINDOW(1) & m_tracks.x <= WINDOW(2) & ...
                            m_tracks.y >= WINDOW(3) & m_tracks.y <= WINDOW(4));
% 
agg_tracks = subStruct(agg_tracks,agg_tracks.x >= WINDOW(1) & agg_tracks.x <= WINDOW(2) & ...
                            agg_tracks.y >= WINDOW(3) & agg_tracks.y <= WINDOW(4) & ...
                            agg_tracks.frame == AGG_FRAME & ...
                            ismember(agg_tracks.id,stable_aggs));

% Cell Trajectories
figure, hold on
for i = unique(tracks.id)'
    xi = tracks.x(tracks.id == i);
    yi = tracks.y(tracks.id == i);
    plot(xi,yi,'-','Color',color_chooser(i))
end

%Agg locations
hold on
xc = cos(0:0.1:2*pi);
yc = sin(0:0.1:2*pi);
for j = 1:length(agg_tracks.x)
    ph = agg_tracks.orientation(j);
    a = agg_tracks.majorAxis(j)/2;
    b = agg_tracks.minorAxis(j)/2;

    R = [ cos(ph)   sin(ph)
         -sin(ph)   cos(ph)];
    xy = R*[xc .* a; yc .* b];
    xrotated = xy(1,:) + agg_tracks.x(j);
    yrotated = xy(2,:) + agg_tracks.y(j);
    plot(xrotated,yrotated,'--k','LineWidth',1.5);
end

%Scale Bar
plot([WINDOW(2) - SCALE_BAR.LENGTH - SCALE_BAR.X_MARGIN; ...
      WINDOW(2) - SCALE_BAR.X_MARGIN], ...
     [SCALE_BAR.Y_MARGIN; ...
      SCALE_BAR.Y_MARGIN], ...
     'LineWidth',SCALE_BAR.WIDTH, ...
     'Color','k');

%Figure formatting
ax = gca;
xlim(WINDOW(1:2))
ylim(WINDOW(3:4))
ax.XTick = [];
ax.YTick = [];
ax.XColor = 'none';
ax.YColor = 'none';

if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '(A)WT'],'Style','none','Formats',{'png'},'Size',FIG_SIZE,'r',1000);
end

%% pure csgA
CSGA_TRAJECTORIES_FNAME = [PROJECT_BASENAME '/data/processed/Development/csgA in csgA/Exp1_32'];
WINDOW = [0, 500, 0, 500] % [xmin, xmax, ymin, ymax]

load([CSGA_TRAJECTORIES_FNAME, '/m_tracks.mat']);

tracks = subStruct(m_tracks,m_tracks.x >= WINDOW(1) & m_tracks.x <= WINDOW(2) & ...
                            m_tracks.y >= WINDOW(3) & m_tracks.y <= WINDOW(4));
                        
% Cell Trajectories
figure, hold on
for i = unique(tracks.id)'
    xi = tracks.x(tracks.id == i);
    yi = tracks.y(tracks.id == i);
    plot(xi,yi,'-','Color',color_chooser(i))
end

%Figure formatting
ax = gca;
xlim(WINDOW(1:2))
ylim(WINDOW(3:4))
ax.XTick = [];
ax.YTick = [];
ax.XColor = 'none';
ax.YColor = 'none';

if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '(B)csgA'],'Style','none','Formats',{'png'},'Size',FIG_SIZE,'r',1000);
end

%% Rescued Csga
Rcsga_TRAJECTORIES_FNAME = [PROJECT_BASENAME '/data/processed/Development/csgA in WT/Exp1_35'];
WINDOW = [400, 900, 300, 700]; % [xmin, xmax, ymin, ymax]

load([Rcsga_TRAJECTORIES_FNAME, '/m_tracks.mat']);
load([Rcsga_TRAJECTORIES_FNAME, '/agg_tracks_with_area.mat']);

tracks = subStruct(m_tracks,m_tracks.x >= WINDOW(1) & m_tracks.x <= WINDOW(2) & ...
                            m_tracks.y >= WINDOW(3) & m_tracks.y <= WINDOW(4));
% 
agg_tracks = subStruct(agg_tracks,agg_tracks.x >= WINDOW(1) & agg_tracks.x <= WINDOW(2) & ...
                            agg_tracks.y >= WINDOW(3) & agg_tracks.y <= WINDOW(4) & ...
                            agg_tracks.frame == AGG_FRAME & ...
                            ismember(agg_tracks.id,stable_aggs));

% Cell Trajectories
figure, hold on
for i = unique(tracks.id)'
    xi = tracks.x(tracks.id == i);
    yi = tracks.y(tracks.id == i);
    plot(xi,yi,'-','Color',color_chooser(i))
end

%Agg locations
hold on
xc = cos(0:0.1:2*pi);
yc = sin(0:0.1:2*pi);
for j = 1:length(agg_tracks.x)
    ph = agg_tracks.orientation(j);
    a = agg_tracks.majorAxis(j)/2;
    b = agg_tracks.minorAxis(j)/2;

    R = [ cos(ph)   sin(ph)
         -sin(ph)   cos(ph)];
    xy = R*[xc .* a; yc .* b];
    xrotated = xy(1,:) + agg_tracks.x(j);
    yrotated = xy(2,:) + agg_tracks.y(j);
    plot(xrotated,yrotated,'--k','LineWidth',1.5);
end

%Figure formatting
ax = gca;
xlim(WINDOW(1:2))
ylim(WINDOW(3:4))
ax.XTick = [];
ax.YTick = [];
ax.XColor = 'none';
ax.YColor = 'none';

if (SAVE)
    saveFigures('SaveAs',[SAVE_FOLDER '(C)rCsga'],'Style','none','Formats',{'png'},'Size',FIG_SIZE,'r',1000);
end