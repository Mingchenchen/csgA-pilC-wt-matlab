%% Create dataset with WT behaviors inside aggreagte, complemented csgA 
%    behaviors outside aggreates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/AllDataTable.mat';
csgAADT = AllDataTable;

load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/AllDataTable.mat';
wtADT = AllDataTable;

clear AllDataTable

%% Combined the two
DENSITY_CUTOFF = 2.3226;
    
combinedADT = csgAADT(csgAADT.rho1 < DENSITY_CUTOFF ,:);

combinedADT = [combinedADT; wtADT(wtADT.rho1 >= DENSITY_CUTOFF,:)];

%% Plot the behaviors to see if the combination worked
% combined;
%%%
addpath('util_functions')
NBOOT = 10;
T = combinedADT;
T = T(T.state1 < 3,:);

%%% Weight by the number of runs in each movie
weights = zeros(size(T,1),1);
for curMovie = 1:length(unique(T.set))
    f = T.set == curMovie;
    weights(f) = sum(f);
end
weights = 1 - (weights / sum(weights));
[CO] = bootstrap_runs_in_time(T,weights,NBOOT);

%%  WT in WT
% WT in WT
%%%
load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/AllDataTable.mat'
T = AllDataTable;
T = T(T.state1 < 3,:);
%%% Weight by the number of runs in each movie
weights = zeros(size(T,1),1);
for curMovie = 1:length(unique(T.set))
    f = T.set == curMovie;
    weights(f) = sum(f);
end
weights = 1 - (weights / sum(weights));

[WW] = bootstrap_runs_in_time(T,weights,NBOOT);

%% CsgA in WT
% CsgA in WT
%%%
load '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/AllDataTable.mat'
T = AllDataTable;
T = T(T.state1 < 3,:);

%%% Weight by the number of runs in each movie
weights = zeros(size(T,1),1);
for curMovie = 1:length(unique(T.set))
    f = T.set == curMovie;
    weights(f) = sum(f);
end
weights = 1 - (weights / sum(weights));
[CW] = bootstrap_runs_in_time(T,weights,NBOOT);

%% Pub Figure
figure,
    subplot(1,3,1)
        hold on
        plot_boot(WW.outside.xi / 2 / 60,WW.outside.Rs,'b')
        plot_boot(WW.inside.xi / 2 / 60,WW.inside.Rs,'k')
        plot_boot(CW.outside.xi / 2 / 60,CW.outside.Rs,'--b')
        plot_boot(CW.inside.xi / 2 / 60,CW.inside.Rs,'--k')
        plot_boot(CO.outside.xi / 2 / 60,CO.outside.Rs,':b')
        plot_boot(CO.inside.xi / 2 / 60,CO.inside.Rs,':k')
        hold off
        ax = gca;
        ylim([0 ax.YLim(2)])
        ylabel('Run Speed (\mum/min)')
        xlabel('Time (hr)')
    subplot(1,3,2)
        hold on
        plot_boot(WW.outside.xi / 2 / 60,WW.outside.Rt,'b')
        plot_boot(WW.inside.xi / 2 / 60,WW.inside.Rt,'k')
        plot_boot(CW.outside.xi / 2 / 60,CW.outside.Rt,'--b')
        plot_boot(CW.inside.xi / 2 / 60,CW.inside.Rt,'--k')
        plot_boot(CO.outside.xi / 2 / 60,CO.outside.Rt,':b')
        plot_boot(CO.inside.xi / 2 / 60,CO.inside.Rt,':k')
        hold off
        ax = gca;
        ylim([0 ax.YLim(2)])
        ylabel('Run Duration (min)')
        xlabel('Time (hr)')
    subplot(1,3,3)
        hold on
        plot_boot(WW.outside.xi / 2 / 60,WW.outside.Rd,'b')
        plot_boot(WW.inside.xi / 2 / 60,WW.inside.Rd,'k')
        plot_boot(CW.outside.xi / 2 / 60,CW.outside.Rd,'--b')
        plot_boot(CW.inside.xi / 2 / 60,CW.inside.Rd,'--k')
        plot_boot(CO.outside.xi / 2 / 60,CO.outside.Rd,':b')
        plot_boot(CO.inside.xi / 2 / 60,CO.inside.Rd,':k')
        hold off
        ax = gca;
        ylim([0 ax.YLim(2)])
        ylabel('Run Distance (\mum)')
        xlabel('Time (hr)')
        
%% Save the combined AllDataTable;

AllDataTable = combinedADT;
save '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in CsgA out/AllDataTable' AllDataTable