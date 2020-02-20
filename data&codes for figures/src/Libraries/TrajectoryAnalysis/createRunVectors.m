function [runVectors, vectorSpeeds] = createRunVectors(tracks,varargin)
%
% runVectors
% {
%   state           movement state of the vector
%   id              cell id
%   distance        vector distance (\mum)
%   speed           average cell speed during vector
%   length          time vector took (min)
%   start {         vector start
%        x 
%        y 
%        density
%        frame
%   }
%   start {         vector stop
%        x 
%        y 
%        density
%        frame
%   }
% }

% NOTE:
% the "WrapStops" function wraps a stop into the previous run, but
% does not wrap concestive same-direction runs toegather.
% 
% for example 
%   cells state transtions "1 to 3 to 2" become "1 to 2"
%   but state transitions "1 to 3 to 1" become "1 to 1"
%
% changeing the behavior to " "1 to 3 to 1" becomes "1" is a TODO!

    p = inputParser;

    addOptional(p,'WrapStops', false);
    addRequired(p,'tracks',@(x) isstruct(x));
    parse(p,tracks,varargin{:} );
    
    if(p.Results.WrapStops)
        cur_id = tracks.id(1);
        moving_state = 0;
        if(tracks.state(1) == 3)
            moving_index = find(tracks.state < 3,1);
            end_run_index = find(tracks.id ~= cur_id,1);

            if(moving_index < end_run_index)
                moving_state = tracks.state(moving_index);
            else
                moving_state = 1;
            end
        end

        for i = 1:length(tracks.id)
            if(cur_id == tracks.id(i))
                if(tracks.state(i) == 3)
                    tracks.state(i) = moving_state;
                else
                    moving_state = tracks.state(i);
                end
            else
                cur_id = tracks.id(i);

                if(tracks.state(i) == 3)
                    moving_index = find(tracks.state(i:end) < 3,1);
                    end_run_index = find(tracks.id(i:end) ~= cur_id,1);

                    if(moving_index < end_run_index)
                      moving_state = tracks.state(i + moving_index - 1);
                    else
                      moving_state = 1;
                    end

                    tracks.state(i) = moving_state;
                end
            end
        end
    end
 
    state_begin = find([true; [diff(tracks.id(1:end-1)) ~= 0 | (diff(tracks.id(2:end)) == 0 & diff(tracks.state(2:end)) ~= 0)]; false]);
    state_end   = find([false;[diff(tracks.id(2:end)) ~= 0   | (diff(tracks.state(2:end)) ~= 0 & diff(tracks.id(1:end-1)) == 0)]; true]);

    runVectors.state = tracks.state(state_end);
    runVectors.id = tracks.id(state_begin);

    runVectors.count = tracks.count(state_end);
    runVectors.start.x = tracks.x(state_begin);
    runVectors.start.y = tracks.y(state_begin);
    runVectors.start.density = tracks.density(state_begin);
    runVectors.start.frame = tracks.frame(state_begin);
    if(isfield(tracks,'in_crest'))
        runVectors.start.in_crest = tracks.in_crest(state_begin);
    end
    if(isfield(tracks,'in_agg'))
        runVectors.start.in_agg = tracks.in_agg(state_begin);
    end

    runVectors.stop.x = tracks.x(state_end);
    runVectors.stop.y = tracks.y(state_end);
    runVectors.stop.density = tracks.density(state_end);
    runVectors.stop.frame = tracks.frame(state_end);
    if(isfield(tracks,'in_crest'))
        runVectors.stop.in_crest = tracks.in_crest(state_end);
    end
    if(isfield(tracks,'in_agg'))
        runVectors.stop.in_agg = tracks.in_agg(state_end);
    end
    
    runVectors.distance = sqrt((runVectors.start.x - runVectors.stop.x).^2 ...
                                + (runVectors.start.y - runVectors.stop.y).^2);
                            
    runVectors.length = (runVectors.stop.frame - runVectors.start.frame) / 2;
    runVectors.speed = runVectors.distance ./ runVectors.length;
    
    %vectorSpeeds.speed = [];
    %vectorSpeeds.density = [];
    vectorSpeeds.time = [];
    vectorSpeeds.Rt = [];
    vectorSpeeds.speed1 = [];
    vectorSpeeds.speed2 = [];
    vectorSpeeds.rho = [];
    vectorSpeeds.state = [];
%     for i = 1:length(state_begin);
%         begin = state_begin(i);
%         tend = state_end(i);
%         
%         a = subStruct(tracks,begin:tend);
%         
%         
%         run_dist = runVectors.distance(i);
% 
%         dist_per_step = sqrt( ...
%                         (a.x(1:end-1) - a.x(2:end)).^2 + ...
%                         (a.y(1:end-1) - a.y(2:end)).^2 ...
%                      );
%         dist_traveled = sum(dist_per_step);
%         dt = a.frame(2:end) - a.frame(1:end-1);
% 
%         frac_dist_traveled = (dist_per_step ./ dt) / dist_traveled;
%         spd = frac_dist_traveled .* runVectors.distance(i)  * 2;
%         %vectorSpeeds.speed   = [vectorSpeeds.speed; spd];
%         vectorSpeeds.speed1  = [vectorSpeeds.speed1; spd(1:end-2)];
%         vectorSpeeds.speed2  = [vectorSpeeds.speed2; spd(2:end-1)];
%         vectorSpeeds.rho     = [vectorSpeeds.rho; a.density(2:end-2)];
%         vectorSpeeds.time    = [vectorSpeeds.time; repmat(runVectors.start.frame(i),length(a.frame) - 3,1)];
%         vectorSpeeds.Rt = [vectorSpeeds.Rt; repmat(runVectors.length(i),length( a.frame) - 3,1)];
%         vectorSpeeds.state = [vectorSpeeds.state; repmat(runVectors.state(i),length( a.frame) - 3,1)];
%         %vectorSpeeds.density = [vectorSpeeds.density; a.density(1:end-1)];
% 
%         %vectorSpeeds.Rt = [vectorSpeeds.Rt; repmat(runVectors.length(i),length( a.frame) - 1,1)];
%     end
    
    %runVectors.meanJumpSpeed = runVectors.distanceTraveled ./ runVectors.length;
end