%% Estimate bandwidth to use for aggreagte detection using the mean bandwith form the last
% frame from the simulations
folder_name = 'S7';
Nruns = 2;
nCells = 5000;
p = Progress(length(fnames) * Nruns);
bandwidths = zeros(2,p.max_value);
for set = 1:length(fnames)
    for run = 1:Nruns
        p.d(run + (Nruns * (set - 1)))
        load(['/Users/cotter/Desktop/NewData/Good/Results/OpenLoop/' folder_name '/Set' num2str(set) '_Run' num2str(run) ,'_nCells' num2str(nCells), '_sim_tracks.mat'])
        sim_tracks =  TrajectoryAnalysis.CreateSimTracks(sim_tracks);

        a = subStruct(sim_tracks,sim_tracks.frame == max(sim_tracks.frame));
        [bdw, d, ~] = kde2d([a.x a.y],2^10,[0 0],[986 740]);
        bandwidths(:,run + (Nruns * (set - 1))) = bdw;
    end
end

%%
A = imresize(d,[739 986]);
B = A; %imfilter(E,fspecial('gaussian',20./0.5)); 
bw = B > 5.5;
bounds = bwboundaries(bw);

figure,
imagesc(B)

hold on;
for i = 1:length(bounds)
    bb = bounds{i};

    plot(bb(:,2),bb(:,1),'-','Color',color_chooser(i));
end

colormap gray;
