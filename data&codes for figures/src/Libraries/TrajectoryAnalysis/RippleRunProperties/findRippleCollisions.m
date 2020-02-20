function [tracks] = findRippleCollisions(m_tracks,varargin)
%   [tracks] = findRippleCollisions(m_tracks,opt_in)
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
    std_window = p.Results.std_window
    weight    = p.Results.weight;
    tracks    = p.Results.tracks;
    
    tracks.in_crest = false(length(tracks.x),1);
    p = Progress(max(tracks.id));
    for id = 1:max(tracks.id)
        p.d(id);
        t = subStruct(tracks,tracks.id == id);
        
        [wmmean,xi,c,s] = movingmeanxy(t.frame,t.density,mean_window);
        [u,xi,c,wmstd] = movingmeanxy(t.frame,t.density,std_window);
        wmstd = max(wmstd,threshold);
        tracks.in_crest(tracks.id == id) = t.density > (wmmean + wmstd .* weight);    
    end
    p.done();
end