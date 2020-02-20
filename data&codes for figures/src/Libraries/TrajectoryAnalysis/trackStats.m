function [ num_frames, num_cells, last_pos ] = trackStats( tracks, varargin)
% [ num_frames, num_cells, last_pos ] = trackStats( tracks, varargin)
%
% Output statistics about the tracks
%
% INPUTS
% tracks    tracks from LAPtracking
% varargin  output figure for any combiniation of:
%               NumFrames
%               NumCells
%               LastPos
%
%Author: Chris Cotter (cotter@uga.edu)
    %%DO Work
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
    
    for opt = varargin
        figure
        switch opt{1}
            case 'NumFrames'
                cdfplot(num_frames ./ 2);
                title('Length of Trajectories');
                xlabel('Trajectory Length (min)');
                ylabel('CDF');

            case 'NumCells'
                plot(num_cells);
                title('Number of cells tracked');
                xlabel('Frame');
                ylabel('Count');

            case 'LastPos'
                scatter(last_pos(:,1),last_pos(:,2)); 
        end
    end
end

