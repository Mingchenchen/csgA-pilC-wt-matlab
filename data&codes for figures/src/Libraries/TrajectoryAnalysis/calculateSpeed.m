function [spd, distance, state, dens] = calculateSpeed(tracks)
    %[spd id] = calculateSpeed(tracks,interactive)
    %
    % Calculates instentainous cell speeds. 
    %
    %
    %
    %PARAMS
    % tracks      m_tracks struct. That is, tracks must contain:
    %                tracks.x         x coordiante of cell in um from origin
    %                tracks.y         y coordiante of cell in um from origin
    %                tracks.frame
    %                tracks.id
    %                tracks.state     state of movement
    %
    %RETRUNS
    % spd(i)      speed of cell between each frame
    
    %
    % Private Variables
    % 
    x = tracks.x;
    y = tracks.y;
    ids = tracks.id;
    frame = tracks.frame;
    state_ = tracks.state;
    dens_    = tracks.density; 
    %
    % Code
    %
    spd = [];
    distance = [];
    state = [];
    dens = [];
    for i = unique(ids)'
        curX = x(ids == i);
        curY = y(ids == i);
        curF = frame(ids == i);
        curS = state_(ids == i);
        curD = dens_(ids == i);
        
        if(length(curX) < 2)
            continue
        end
        
        dist = sqrt( (curX(1:end - 1) - curX(2:end)).^2 + (curY(1:end - 1) - curY(2:end)).^2 );        
        distance = [distance; dist];
        
        %Adjust for missing points in the track
        time = diff(curF);
        
        spd = [spd; dist ./ (time / 2)]; 
        state = [state; curS(1:end-1)];
        dens  = [dens; curD(1:end-1)];
    end
end