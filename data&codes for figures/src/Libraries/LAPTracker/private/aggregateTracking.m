function [tracks] = aggregateTracking(posList,END_THRESHOLD)
    
    %% Extract data
    %     p = inputParser;
    %     
    %     %addOptional(p,'DistanceCutoff', 40, @isnumeric);
    %     addRequired(p,'z',@(x) isnumeric(x));
    % 	addRequired(p,'Q',@(x) isnumeric(x) && (size(x,1) == size(x,2)));
    % 	addRequired(p,'G');
    % 	addRequired(p,'H');
    % 	addRequired(p,'R');
    % 	addRequired(p,'P')
    %     
    %     parse(p,z,Q,G,H,R,Pint,varargin{:} );
    
    %%
    %Initialize Variables
	steps = max(posList.frame);
    posList.index = (1:length(posList.x))'; 
    aggs_in_slice = histc(posList.frame,1:max(posList.frame));
    END_COST = 20;
    END_COST_BACKUP= 15;
    
    
    lastFrame = subStruct(posList,posList.frame == 1);
    n_tracks  = length(lastFrame.x);
	lastFrame.id = (1:n_tracks)';
    tracks.x = [];
    tracks.y = [];
    tracks.majorAxis = [];
    tracks.minorAxis = [];
    tracks.id = [];
    tracks.index = [];
    tracks.frame = [];
    tracks.propagated = [];
    tracks = add_tracks(tracks,lastFrame,lastFrame.id,1);
    p = Progress(steps)
	for k = 2:steps
        p.d(k);
        curFrame = subStruct(posList,posList.frame == k);
        
        %% Link Aggs
        if(~isempty(lastFrame.x) & ~isempty(curFrame.x))
            %Calculate costs
            dCost = pdist2([lastFrame.x lastFrame.y],[curFrame.x curFrame.y]);
            aCost = pdist2(lastFrame.majorAxis,curFrame.majorAxis);
            bCost = pdist2(lastFrame.minorAxis,curFrame.minorAxis);

            %Combine the costs
            sigmaXY = 25;
            sigmaA = 50;
            sigmaB = 50;
            costs = Px(dCost,sigmaXY) + Px(aCost,sigmaA) + Px(bCost,sigmaB);

            [p12,p21] = lap(costs,-1,1,1,END_COST);

            i = p21(p21(1:aggs_in_slice(k)) <= length(lastFrame.x));
            j = p12(i);
        else
            i = [];
            j = [];
            costs = [];
        end
        
        %% Add new positions to tracks
        tt.x = [];
        tt.y = [];
        tt.majorAxis = [];
        tt.minorAxis = [];
        tt.id = [];
        tt.index = [];
        tt.frame = [];
        tt.propagated = [];
        
        if(~isempty(i))
            tt = add_tracks(tt, subStruct(curFrame,j), lastFrame.id(i), k);
        end

        %Add new aggregates
        new = setxor(j,1:aggs_in_slice(k));
        if(~isempty(new))
            tt = add_tracks(tt, subStruct(curFrame,new), ((n_tracks + 1):(n_tracks+length(new)))', k);
            n_tracks = n_tracks + length(new);
        end

        %Propagate old aggreagtes
        gone = setxor(i,1:length(lastFrame.x));
        keep = gone(lastFrame.frame(gone) > k - END_THRESHOLD);
        if(~isempty(keep))
            tt = add_tracks(tt, ...
                            subStruct(lastFrame,keep), ...
                            lastFrame.id(keep), ...
                            lastFrame.frame(keep), ...
                            true(length(keep),1));
        end
 
        %Add frame tracks to track list
        tracks = add_tracks(tracks,tt,tt.id,k,tt.propagated);

        %Remove propagated aggreagtes that were not relinekd before END_THRESHOLD
        tracks = subStruct(tracks,~(tracks.frame > k - END_THRESHOLD & ...
                                          ismember(tracks.id,lastFrame.id(gone(lastFrame.frame(gone) == k - END_THRESHOLD)))));
        
        
        %Update last Frame
        lastFrame = tt;


        END_COST = max(costs(sub2ind(size(costs),i,j))) * 1.05;

        if(isnan(END_COST))
            END_COST = END_COST_BACKUP;
        else
            END_COST_BACKUP = END_COST;
        end
            
    end
    
    tracks = tracks; %concatenateTracks(tracks,20);
    p.done();
end

function tracks = add_tracks(tracks, add, id, frame, propagated)
    tracks.x = [tracks.x; add.x];
    tracks.y = [tracks.y; add.y];
    tracks.majorAxis = [tracks.majorAxis; add.majorAxis];
    tracks.minorAxis = [tracks.minorAxis; add.minorAxis];
    tracks.index = [tracks.index; add.index];
    
    if(nargin > 4)
        tracks.propagated = [tracks.propagated; propagated];
    else
        tracks.propagated = [tracks.propagated; false(length(add.x),1)];
    end
    
    if(length(frame) == 1)
        tracks.frame = [tracks.frame; repmat(frame,length(add.x),1)];
    else
        tracks.frame = [tracks.frame; frame];
    end
    
    assert(length(id) == length(add.x),'Id length no equal to struct length');
    tracks.id = [tracks.id; id];
end

function [p] = Px(x,sigma)
    p = (x).^2./(2*sigma.^2) - log(1/(sigma.*sqrt(2*pi)));
end

function [tt] = concatenateTracks(tt,min_track_length)
    if(min_track_length == 0)
        return
    end
    
    unv = unique(tt.id);
    frameCount = histc(tt.id,unv);
    short = find(frameCount <= min_track_length);
    tt = subStruct(tt,~ismember(tt.id,unv(short)));
    count = 1;
    for i = unique(tt.id)'
        tt.id(tt.id == i) = count;
        count = count + 1;
    end
end