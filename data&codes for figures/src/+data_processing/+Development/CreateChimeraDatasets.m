addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

%Number of frames to skip prior to calculating inside aggreagte values
%Set to a frame after aggreagte tracking becomes fairly stable. 
STABLE_AGG_CUTOFF = 120;

% Inside aggreagte density cutoff (cells/um^2)
DENSITY_CUTOFF = 2.3;

load([PROJECT_BASENAME '/data/processed/Development/WT in WT/AllDataTable'])
WT_T = AllDataTable;

load([PROJECT_BASENAME '/data/processed/Development/csgA in WT/AllDataTable'])
CSGA_T = AllDataTable;

%% Csga Behaviors outside WT behaivors inside
AllDataTable = WT_T(WT_T.rho1 > DENSITY_CUTOFF,:);
AllDataTable = [AllDataTable; CSGA_T(CSGA_T.rho1 <= DENSITY_CUTOFF,:)];
    
save_folder = [PROJECT_BASENAME '/data/processed/Development/Chimeras/CsgAOutsideWTInside/'];
mkdir(save_folder)
save([save_folder '/AllDataTable'],'AllDataTable')

%Data Check
T = AllDataTable;
T = T(T.state1 < 3,:);

figure,
    subplot(2,1,1)
    Ti = T(T.rho1 <= DENSITY_CUTOFF,:);
    u = movingmeanxy(Ti.TSS1,Ti.Rd1,20,unique(Ti.TSS1));
    plot(unique(Ti.TSS1),u)
    title('csgA Outisde')
    
    subplot(2,1,2)
    Ti = T(T.rho1 > DENSITY_CUTOFF,:);
    u = movingmeanxy(Ti.TSS1,Ti.Rd1,20,unique(Ti.TSS1));
    plot(unique(Ti.TSS1),u)
    title('WT Inside')
    
%% WT Behaviors outside csgA behaivors inside
AllDataTable = WT_T(WT_T.rho1 <= DENSITY_CUTOFF,:);
AllDataTable = [AllDataTable; CSGA_T(CSGA_T.rho1 > DENSITY_CUTOFF,:)];
    
save_folder = [PROJECT_BASENAME '/data/processed/Development/Chimeras/WTOutsideCsgAInside/'];
mkdir(save_folder)
save([save_folder '/AllDataTable'],'AllDataTable')

%Data Check
T = AllDataTable;
T = T(T.state1 < 3,:);

figure,
    subplot(2,1,1)
    Ti = T(T.rho1 <= DENSITY_CUTOFF,:);
    u = movingmeanxy(Ti.TSS1,Ti.Rd1,20,unique(Ti.TSS1));
    plot(unique(Ti.TSS1),u)
    title('WT Outisde')
    
    subplot(2,1,2)
    Ti = T(T.rho1 > DENSITY_CUTOFF,:);
    u = movingmeanxy(Ti.TSS1,Ti.Rd1,20,unique(Ti.TSS1));
    plot(unique(Ti.TSS1),u)
    title('CsgA Inside')
