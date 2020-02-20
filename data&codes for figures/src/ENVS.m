%%%
% These varaible allow to be overwritten by the script loading ENVS.m
if(~exist('OVERWRITE','var'))        OVERWRITE = false;        end
if(~exist('PROJECT_BASENAME','var')) PROJECT_BASENAME = '/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/';  end
if(~exist('HEADLESS','var'))         HEADLESS = 1;             end
if(~exist('RAW_DATA_LOC','var'))
    RAW_DATA_LOC = [PROJECT_BASENAME, '/data/raw/']; 
end

FOV_X = 1440; % X size of microscope field of view (in pixels)
FOV_Y = 1920; % Y size of microscope filed of view (in pixels)
FOV_MU_CONV = 0.5136; % micrometers per pixel

%Images are shrunk during normalization to save space and reduce
% computation time
Kum_X = 2^10; % X size of microscope images after normalization
Kum_Y = 2^10; % Y size of microscope images after normalization
Kum_MU_CONV = [FOV_X/Kum_X * FOV_MU_CONV, ...
               FOV_Y/Kum_Y * FOV_MU_CONV]; %witdh of each pixel in Kum in 
                                           %  micrometers

%%%
% Constants defining data locations and names

WT_RAW_DEV_MOVIES = containers.Map();
WT_RAW_DEV_MOVIES('10112015') = [RAW_DATA_LOC, '/Development/WT in WT/10112015/10112015_1_MMStack.ome.tif'];
WT_RAW_DEV_MOVIES('10152015') = [RAW_DATA_LOC, '/Development/WT in WT/10152015/10112015_2_MMStack.ome.tif'];
WT_RAW_DEV_MOVIES('10282015') = [RAW_DATA_LOC, '/Development/WT in WT/10282015/10282015_1_MMStack.ome.tif'];
                 
CSGA_RAW_DEV_MOVIES = containers.Map();
CSGA_RAW_DEV_MOVIES('A_30')    = [RAW_DATA_LOC, '/Development/csgA in WT/A_30/A_30_MMStack_Pos0.ome.tif'];
CSGA_RAW_DEV_MOVIES('DEV_1')   = [RAW_DATA_LOC, '/Development/csgA in WT/DEV_1/DEV_1_MMStack_Pos0.ome.tif'];
CSGA_RAW_DEV_MOVIES('Exp1_35') = [RAW_DATA_LOC, '/Development/csgA in WT/Exp1_35/Exp1_35_MMStack_Pos0.ome.tif'];

CSGA_IN_CSGA_RAW_DEV_MOVIES = containers.Map();
CSGA_IN_CSGA_RAW_DEV_MOVIES('Exp1_31') = [RAW_DATA_LOC, '/Development/csgA in csgA/Exp1_31/Exp1_31_MMStack_Pos0.ome.tif'];
CSGA_IN_CSGA_RAW_DEV_MOVIES('Exp1_32') = [RAW_DATA_LOC, '/Development/csgA in csgA/Exp1_32/Exp1_32_MMStack_Pos0.ome.tif'];
CSGA_IN_CSGA_RAW_DEV_MOVIES('Exp1_34') = [RAW_DATA_LOC, '/Development/csgA in csgA/Exp1_34/Exp1_34_MMStack_Pos0.ome.tif'];

WT_RAW_RIPPLE_MOVIES = containers.Map();
WT_RAW_RIPPLE_MOVIES('RIPPLE_3') =  [RAW_DATA_LOC, '/Rippleing/WT in WT/RIPPLE_3/RIPPLE_3_MMStack_Pos0.ome.tif'];
WT_RAW_RIPPLE_MOVIES('RIPPLE_24') = [RAW_DATA_LOC, '/Rippleing/WT in WT/RIPPLE_24/RIPPLE_24_MMStack_Pos0.ome.tif'];
WT_RAW_RIPPLE_MOVIES('RIPPLE_28') = [RAW_DATA_LOC, '/Rippleing/WT in WT/RIPPLE_28/RIPPLE_28_MMStack_Pos0.ome.tif'];

CSGA_RAW_RIPPLE_MOVIES = containers.Map();
CSGA_RAW_RIPPLE_MOVIES('RIPPLE_25') = [RAW_DATA_LOC, '/Rippleing/csgA in WT/RIPPLE_25/RIPPLE_25_MMStack_Pos0.ome.tif'];
CSGA_RAW_RIPPLE_MOVIES('RIPPLE_26') = [RAW_DATA_LOC, '/Rippleing/csgA in WT/RIPPLE_26/RIPPLE_26_MMStack_Pos0.ome.tif'];
CSGA_RAW_RIPPLE_MOVIES('RIPPLE_27') = [RAW_DATA_LOC, '/Rippleing/csgA in WT/RIPPLE_27/RIPPLE_27_MMStack_Pos0.ome.tif'];

DIFA_RAW_RIPPLE_MOVIES = containers.Map();
DIFA_RAW_RIPPLE_MOVIES('RIPPLE_6') = [RAW_DATA_LOC, '/Rippleing/difA in WT/RIPPLE_6/RIPPLE_6_MMStack_Pos0.ome.tif'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auto Generated Variables
%%%

%%%
% DEV_FOLDERS provide the path locations to put data inside
% project folders such as: data/processed, reports and data/interim
%%%

%%% DEV Variables 
    ALL_RAW_DEV_MOVIES = [WT_RAW_DEV_MOVIES; CSGA_RAW_DEV_MOVIES];

    WT_DEV_FOLDERS = containers.Map();
    for i = keys(WT_RAW_DEV_MOVIES)
        k = i{1};
        fname = split(WT_RAW_DEV_MOVIES(k),'/');
        s = join(fname(end-3:end-1),'/');
        WT_DEV_FOLDERS(k) = ['/', s{1}, '/'];
    end

    CSGA_DEV_FOLDERS = containers.Map();
    for i = keys(CSGA_RAW_DEV_MOVIES)
        k = i{1};
        fname = split(CSGA_RAW_DEV_MOVIES(k),'/');
        s = join(fname(end-3:end-1),'/');
        CSGA_DEV_FOLDERS(k) = ['/', s{1}, '/'];
    end
    
    CSGA_IN_CSGA_DEV_FOLDERS = containers.Map();
    for i = keys(CSGA_IN_CSGA_RAW_DEV_MOVIES)
        k = i{1};
        fname = split(CSGA_IN_CSGA_RAW_DEV_MOVIES(k),'/');
        s = join(fname(end-3:end-1),'/');
        CSGA_IN_CSGA_DEV_FOLDERS(k) = ['/', s{1}, '/'];
    end

    All_DEV_FOLDERS = [WT_DEV_FOLDERS; CSGA_DEV_FOLDERS];

%%% RIPPLE Variables
    ALL_RAW_RIPPLE_MOVIES = [WT_RAW_RIPPLE_MOVIES; CSGA_RAW_RIPPLE_MOVIES; DIFA_RAW_RIPPLE_MOVIES];

    WT_RIPPLE_FOLDERS = containers.Map();
    for i = keys(WT_RAW_RIPPLE_MOVIES)
        k = i{1};
        fname = split(WT_RAW_RIPPLE_MOVIES(k),'/');
        s = join(fname(end-3:end-1),'/');
        WT_RIPPLE_FOLDERS(k) = ['/', s{1}, '/'];
    end

    CSGA_RIPPLE_FOLDERS = containers.Map();
    for i = keys(CSGA_RAW_RIPPLE_MOVIES)
        k = i{1};
        fname = split(CSGA_RAW_RIPPLE_MOVIES(k),'/');
        s = join(fname(end-3:end-1),'/');
        CSGA_RIPPLE_FOLDERS(k) = ['/', s{1}, '/'];
    end
    
    DIFA_RIPPLE_FOLDERS = containers.Map();
    for i = keys(DIFA_RAW_RIPPLE_MOVIES)
        k = i{1};
        fname = split(DIFA_RAW_RIPPLE_MOVIES(k),'/');
        s = join(fname(end-3:end-1),'/');
        DIFA_RIPPLE_FOLDERS(k) = ['/', s{1}, '/'];
    end

    ALL_RIPPLE_FOLDERS = [WT_RIPPLE_FOLDERS; CSGA_RIPPLE_FOLDERS; DIFA_RIPPLE_FOLDERS];
    
%%%
% Cleanup
%%%exi
    clear i k ans fname s
