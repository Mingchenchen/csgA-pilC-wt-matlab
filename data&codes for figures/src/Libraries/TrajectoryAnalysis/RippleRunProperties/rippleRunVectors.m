function [runVectors,state_begin,state_end] = rippleRunVectors(m_tracks,varargin)
%   [runVectors,in_ripple_list] = rippleRunVectors(m_tracks,varargin)
%
%   opt-in:
%        threshold
%        window
%        weight
%

    addpath('Libraries/Utils');
    %% Extract data
    p = inputParser;
    
    %addOptional(p,'DistanceCutoff', 40, @isnumeric);
    addOptional(p,'threshold', 2.4e-3, @isnumeric);
    addOptional(p,'mean_window', 5, @isnumeric);
    addOptional(p,'std_window', 2, @isnumeric);
    addOptional(p,'weight', 0.8, @isnumeric);
    addRequired(p,'tracks');
    
    parse(p,m_tracks, varargin{:} );
    
    threshold = p.Results.threshold;
    mean_window    = p.Results.mean_window;
    std_window  = p.Results.std_window
    weight      = p.Results.weight;
    tracks      = p.Results.tracks;
    
    %% Find when cell is in ripple
    state_begin = [];
    state_end   = [];
    
    p = Progress(max(tracks.id));
    tracks.index = [1:length(tracks.x)]';
    
    for id = 1:max(tracks.id)
        p.d(id);
        t = subStruct(tracks,tracks.id == id);
       
        [wmmean,xi,c,s] = movingmeanxy(t.frame,t.density,mean_window);
        [u,xi,c,wmstd] = movingmeanxy(t.frame,t.density,std_window);
        wmstd = max(wmstd,threshold);
        
        in_ripple = t.density > (wmmean + wmstd .* weight);
        
        props = regionprops(in_ripple,'Centroid');
        a = [props.Centroid];
        a = round(a(2:2:end))';
%         same_ripple = diff(t.frame(a)) <= 4;
%         same_r_start = find(diff([0; same_ripple]) > 0);
%         same_r_stop  = find(diff([same_ripple; 0]) < 0) + 1;
%         different_ripple = diff(t.frame(a)) > 4;
% 
%         a = a(sort([find(different_ripple) + 1; round((same_r_stop + same_r_start)/2)]));

        a = union(a,[1 length(t.density)]');
        state_begin = [state_begin; a(1:end-1) + t.index(1) - 1];
        state_end   = [state_end; a(2:end) + t.index(1) - 1];
    end
    p.done();
    
    %% Generate ripple run vectors
    runVectors.state = tracks.state(state_end);
    runVectors.id = tracks.id(state_begin);

    runVectors.count = tracks.count(state_end);
    runVectors.start.x = tracks.x(state_begin);
    runVectors.start.y = tracks.y(state_begin);
    runVectors.start.density = tracks.density(state_begin);
    runVectors.start.frame = tracks.frame(state_begin);


    runVectors.stop.x = tracks.x(state_end);
    runVectors.stop.y = tracks.y(state_end);
    runVectors.stop.density = tracks.density(state_end);
    runVectors.stop.frame = tracks.frame(state_end);
    
    runVectors.distance = sqrt((runVectors.start.x - runVectors.stop.x).^2 ...
                                + (runVectors.start.y - runVectors.stop.y).^2);
                            
    runVectors.length = (runVectors.stop.frame - runVectors.start.frame) / 2;
    runVectors.speed = runVectors.distance ./ runVectors.length;
end