%% Calculate Simulation Fraction of cells inside aggreagte curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Libraries/Utils')
addpath('Models/ModelBase/Analysis')
addpath('Libraries/')
addpath('Libraries/LAPTracker/')
run('ENVS.m')

%Simulation properties
SIM_NAME = 'PUB_NEW_12';     %Used to name output files
Nruns = 3;            %Number of repliate runs performed
prerun_length = 180;  %Length of prerun
bandwidth = 15;       %bandwith to use for estimating deisity
AGG_DENSITY_CUTOFF = 2.3226; %Cutoff density for detecting aggreagtes
                             %  in cells/um^2
AVE_CELL_DENSITY = 1.12;    %Average cell density within the field of view
                            %  in cells/um^2
                        
%AGG_DENSITY_CUTOFF = 5.5; % OLD PNAS VALUE

%Where sim data is located
SIM_DATA = [PROJECT_BASENAME, '/data/raw/Simulations/ClosedLoop/' SIM_NAME '/Results/' SIM_NAME '/'];

% Where and if to save resluting figure
SAVE = true;
SAVE_FOLDER = [PROJECT_BASENAME '/data/processed/simulations/development/' SIM_NAME '/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

%Loacation of real data 
MOVIE_DATA = [PROJECT_BASENAME '/data/processed/Development/frac_curves.mat'];

colormap = [[0 0 0];[1 1 1] * 0.5];

%% Calculate Fraction of cells in aggreagtes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simFracCurves = zeros(600 + prerun_length - 1,Nruns);

p = Progress(Nruns)
for r = 1:Nruns  
    p.d(r)
    %%
    load([SIM_DATA '/sim_tracks-' num2str(r) '.mat']);
    sim_tracks = CreateSimTracks(sim_tracks);
    sim_tracks.in_agg = sim_tracks.density > AGG_DENSITY_CUTOFF; %in_out_agg(sim_tracks,agg_tracks);

    sim_inside  = subStruct(sim_tracks,sim_tracks.in_agg > 0);
    sim_inside_count  = histc(sim_inside.frame,min(sim_tracks.frame):max(sim_tracks.frame));
    sim_total_count   = histc(sim_tracks.frame,min(sim_tracks.frame):max(sim_tracks.frame));
    %simFracCurves(:,set,run) = sim_inside_count ./ (sim_total_count);
    simFracCurves(:,r) = sim_inside_count ./ (sim_total_count);
end
p.done()

save([ SAVE_FOLDER '/simFracCurves'],'simFracCurves')

%% Calculate aggreagte properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[AggProps] = extractAggregateProperties(SIM_DATA,Nruns,bandwidth,prerun_length,AGG_DENSITY_CUTOFF,AVE_CELL_DENSITY);

save([SAVE_FOLDER '/AggProps'],'AggProps')
