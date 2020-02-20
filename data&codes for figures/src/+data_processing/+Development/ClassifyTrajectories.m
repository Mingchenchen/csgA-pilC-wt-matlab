% %%%
% % CsgA in WT
% %%%
% fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/A_30/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/A_31/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/DEV_1/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/Exp1_35/' , ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/Exp1_36/' , ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/Exp1_38/'}
                  
%%% 
% WT in WT
%%%
fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10282015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10152015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10112015/'};
      
      
addpath('Libraries/TrajectoryAnalysis')

for set = 1:length(fnames)
    set
    load([fnames{set} 'Kum.mat']);
    load([fnames{set} 'tracks_new.mat']);
    
    %%
    m_tracks = createMTracks(tracks,Kum); 
    save([fnames{set} 'm_tracks_new'],'m_tracks');
    
    %%
    clear Kum
end

%%
for set = 1:length(fnames)
   load([fnames{set} 'm_tracks.mat'])
   tracks = m_tracks;
   
   x = tracks.x;
   y = tracks.y;
   o = tracks.o;
   id = tracks.id;
   frames = tracks.frame;
    
   num_frames = zeros(length(unique(id)),1);
   last_pos = zeros(length(unique(id)),2);
    
    for i = unique(id)'
        xi = x(id == i);
        yi = y(id == i);
        frame = frames(id == i);

        num_frames(i) = max(frame) - min(frame);
        last_pos(i,:) = [xi(end) yi(end)];
    end

    num_cells = histc(frames,min(frames):max(frames));
    
    figure(1),
        hold on
        cdfplot(num_frames ./ 2);
        title('Length of Trajectories');
        xlabel('Trajectory Length (min)');
        ylabel('CDF');

    figure(2),
        hold on
        plot(num_cells);
        title('Number of cells tracked');
        xlabel('Frame');
        ylabel('Count');
        
    drawnow
end

figure(1),
    legend({'A_30','A_31','Dev_1','Exp1_35','Exp1_36','Exp1_38'})
    
figure(2),
    legend({'A_30','A_31','Dev_1','Exp1_35','Exp1_36','Exp1_38'})