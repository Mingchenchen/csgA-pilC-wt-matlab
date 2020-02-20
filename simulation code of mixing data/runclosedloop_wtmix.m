clear all
close all
addpath('./MxanthusCellTracking-master/src/Models/ModelBase')
addpath('./MxanthusCellTracking-master/src/Models/ModelBase/Analysis')
addpath('./MxanthusCellTracking-master/src/Models/ClosedLoopModel/ClosedLoop')
addpath('./MxanthusCellTracking-master/src/Models/ModelBase/SimulationCreation')
addpath('./MxanthusCellTracking-master/src/Libraries/Utils')
addpath('./MxanthusCellTracking-master/src/Models/ClosedLoopModel/SimulationCreation')
prerun = 180;
rng('shuffle');
%load pilc data table
load('AllDataTable.mat')
AllDataTable2=AllDataTable;
[Mdlb2, Mdlr2, Mdlst2, Init2] = create_search_trees(AllDataTable2);
AllDataTable_goin=AllDataTable2(cos(AllDataTable2.beta1)>=0,:);
[Mdlbin, Mdlrin, Mdlstin, Initin] = create_search_trees(AllDataTable_goin);
AllDataTable_goout=AllDataTable2(cos(AllDataTable2.beta1)<0,:);
[Mdlbout, Mdlrout, Mdlstout, Initout] = create_search_trees(AllDataTable_goout);

%load WT data table
load('AllDataWStopsum.mat')
%%% Generate Probs
%%%%%%%%%%%%%%%%%%  
[Mdlb, Mdlr, Mdlst, Init] = create_search_trees(AllDataTable);
probs = ClosedLoopProbs(Mdlb,Mdlr,Mdlst,Init,prerun);

probs2 = ClosedLoopProbsmx2(Mdlb, Mdlr, Mdlst, Mdlb2,Mdlr2,Mdlst2,Init,prerun);

%probs3 = ClosedLoopProbs3(Mdlbout_goin,Mdlbout_goout,Mdlrout,Mdlstout,Initout,prerun);

%%% Run Simulations
%%%%%%%%%%%%%%%%%%%
E = ClosedLoopmx2(probs,10000,probs2,8000,'bandwidth',14);
%E.view('on')
E.loop(600+prerun );

%% Save Results
[data,data2 ]= E.exportData();
[sim_tracks,sim_tracks2] =  E.exportTracks();
sim_tracks=CreateSimTracks(sim_tracks);
sim_tracks2=CreateSimTracks(sim_tracks2);
save('sim_tracks_wt_csgawtbiaslongdur2','sim_tracks');
save('sim_tracks_pilc_csgawtbiaslongdur2','sim_tracks2');
save('datawt_csgawtbiaslongdur2','data');
save('datacsga_csgawtbiaslongdur2','data2')
den=E.field.density;
save('density_wt_csgawtbiaslongdur2','den')
in_agg=E.inagg;
save('in_agg_csgawtbiaslongdur2','in_agg')
