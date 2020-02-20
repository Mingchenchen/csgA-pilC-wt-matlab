function [ tracks, NO_LINKING_COST, CUTOFFS, deltaD,deltaO,deltaX,deltaY,deltaV,model, displacement] = LAPtracking( posList, varargin )
% [ tracks, NO_LINKING_COST, CUTOFFS, deltaD,deltaO, displacement] = LAPtracking( posList, varargin )
%
% Uses the methods described in [1] for mutiple object tracking. Optimized
% for linking trajectories of fluorescently labeled M. xanthus cells by
% penalizing links that make sharp turns (i.e. assumes cells move along
% their long axis). 
%
% Default values are also chosen for M. xanthus tracking from 200x
% maginification and 30 sec/frame. 
%
% INPUTS:
% positionlist: an array listing the scrambled coordinates and data 
%     of the different particles at different times, such that:
%       positionlist([1:2,:): contains the d coordinates and
%           data for all the particles, at the different times. must be positve
%       positionlist(3,:): contains the time t that the position 
%           was determined, must be integers (e.g. frame number.  These values must 
%           be monotonically increasing and uniformly gridded in time.
% OPTIONAL INPUTS
%                    %%%% FRAME-TO-FRAME LINKING %%%%
%  NoLinkingMinCost: Lower bound for the cost to end a trajectory. (Default
%        = 0)
%  EndCostMutiplier: The mutiplier used to calculate the no linking cost in
%       the frame to frame linking matrix. See Jaqman et al supplementals
%       (Default = 1.2)
%  SearchRadiusMutiplier: Mutiplier used to calculate the search radius in
%       the frame to frame linking matrix. See Jaqman et al supplementals
%       (Default = 3)
%  EndCostSeed: The initial cost for no linking in the frame to frame cost
%       matrix. Only used for linking of frame 2 to frame 1. Value is
%       dynamically caculated after that, see EndCostMutiplier. 
%       (Default = 15)
%                       %%%% GAP CLOSING %%%%
%  NoGapClosing: Skips gap closing. WARNING: this ingores the
%       MinTrackLength flag and always includes the PosListIndex column in
%       tracks. Meant for debugging purposes only. (Default = false)
%  MaxSliceGap: Maximum number of frames the beginning of 1 trajectory and
%       the end of another can be to be considered for gap closing (Default
%       is 0)
%  NoiseCutoff: Minimum number of frames a trajectory must be to be
%       considered for gap closing. Useful if the input data was noisy to reduce
%       the complexity of solving the GAP closing matrix (Default is 0)
%  SliceGapScaler: Scaling factor for determining the cost of linking gaps
%       between trajectorys that are more than 1 frame apart. See Jaqman et
%       al supplementals for more info (Default = 2)
%
%                        %%%% OUTPUT %%%%
%  Debug: If set to true, prints out progress information. (Default is false)
%  MinTrackLength: Minimum number of frames a trajectory must be to be
%       included in the final trajectory list (Default is 0)
%  IncludePosListIndex: adds a column to tracks which linkes the row in
%       tracks to the row it came from in posList (Default = false)
%  
% OUTPUTS:
% result:  a list containing the original data rows sorted 
%     into a series of trajectories.  To the original input 
%     data structure there is appended an additional column 
%     containing a unique 'id number' for each identified 
%     particle trajectory.  The result array is sorted so 
%     rows with corresponding id numbers are in contiguous 
%     blocks, with the time variable a monotonically
%     increasing function inside each block.  For example:
%     
%     For the input data structure (posList is a struct with):
%         (posList.x)  (posList.y) (posList.o)  (posList.frame)
%           3.60000      5.00000      50.8        0.00000
%           15.1000      22.6000      45.0        0.00000
%           4.10000      5.50000     -36.2        1.00000 
%           15.9000      20.7000     -66          2.00000
%           6.20000      4.30000     -85          2.00000
%
%     will return the result as a struct tracks:
%         (tracks.x)   (tracks.y) (tracks.o)  (tracks.frame)     (tracks.id)     (tracks.index)
%           3.60000      5.00000   50.8       0.00000            1.0000              1
%           15.1000      22.6000   45.0       0.00000            1.0000              2
%           4.10000      5.50000  -36.2       1.00000            1.0000              3
%           15.9000      20.7000  -66         2.00000            2.0000              4
%           6.20000      4.30000  -85         2.00000            2.0000              5
%
% NO_LINKING_COST: a nx1 vector containing the cost associated with 
%     terminating a trajectory in the linking in frame n 
%
% CUTOFFS: cutoffs is a nx3 where row n contains the mean, max and min
%     search radius cutoffs for linking in frame n
%
% DEPENDS:
%   Statistics Toolbox
%
% REFERENCES
%   [1] 1. Jaqaman K, et al. (2008) Robust single-particle tracking in 
%           live-cell time-lapse sequences. Nat Methods 5: 695-702.
% 
% Author: Chris Cotter cotter@sciencesundries.com
%

%%TODO (in to order):
% 1) Make more efficient use of the sparse matrix (or replace it)
%       ie. remove all the uses of find.
% 2) Add Options:
%      GapClosingSearchRadiusCutoff: The maximum distance that can be closed
%         duing gap closing. (Default = 25, Set to Inf for no cutoff)

%Log at bottom

    %% Extract data
    p = inputParser;
    
    %addOptional(p,'DistanceCutoff', 40, @isnumeric);
    addOptional(p,'Debug', 0, @islogical);
    addOptional(p,'MaxSliceGap', 0, @isnumeric);
    addOptional(p,'MinTrackLength', 0, @isnumeric);
    addOptional(p,'NoiseCutoff', 0, @isnumeric);
    addOptional(p,'SliceGapScaler',2,@isnumeric);
    addOptional(p,'IncludePosListIndex',1,@islogical);
    addOptional(p,'NoLinkingMinCost',0,@isnumeric);
    addOptional(p,'EndCostMutiplier',1.2,@isnumeric);
    %addOptional(p,'GapCostMutiplier',0.9,@isnumeric);
    addOptional(p,'GapClosingSearchRadiusCutoff',25,@isnumeric);
    addOptional(p,'EndCostSeed',10,@isnumeric);
    addOptional(p,'NoGapClosing',false,@islogical);
    addRequired(p,'posList');
    
    
    parse(p,posList, varargin{:} );
     
    %Input Variables
                 x = posList.x;
                 y = posList.y;
                 o = posList.o;
               pos = (1:size(x,1))';
            slices = posList.frame;
    cells_in_slice = histc(slices,1:max(slices));

    %Output Variables
    CUTOFFS                          = NaN(size(slices,2),3);
    NO_LINKING_COST                  = NaN(size(slices,2),1);
    % deltaD{\tau} contains the \delta_(\tau,\tau+1) for all valid links between images \tau and tau+1
    deltaD                           = cell(size(slices,2),1); 
    % deltaO{\tau} contains the \theta_(\tau,\tau+1) for all valid links between images \tau and tau+1
    deltaO                           = cell(size(slices,2),1);
    deltaX                           = cell(size(slices,2),1); 
    deltaY                           = cell(size(slices,2),1); 
    deltaV                           = cell(size(slices,2),1);
    model                            = cell(size(slices,1),1);
    displacement                     = cell(size(slices,2),1); 
    
    %Constants
    DEBUG                            = p.Results.Debug;
    MIN_TRACK_LENGTH                 = p.Results.MinTrackLength;
    END_COST                         = p.Results.EndCostSeed;
    END_COST_MUTIPLIER               = p.Results.EndCostMutiplier;  
    NOISE_CUTOFF                     = p.Results.NoiseCutoff;
    %INCLUDE_INDEX                    = p.Results.IncludePosListIndex;
    %GAP_COST_MUTIPLIER               = p.Results.GapCostMutiplier;
    NO_GAP_CLOSING                   = p.Results.NoGapClosing;
    MAX_SLICE_GAP                    = p.Results.MaxSliceGap;
    MAX_GAP                          = p.Results.GapClosingSearchRadiusCutoff;
    
    %% Link Frames

    if(DEBUG)
        fprintf('Linking Frames: ');
        fprintf('%04d/%04d', 0, max(slices));
        tic;
    end
    
    nTracks = cells_in_slice(1);
    prev_pos = [x(slices == min(slices)) y(slices == min(slices)) o(slices == min(slices)) (1:nTracks)'];
    dx = NaN(size(prev_pos,1),6); %Last 6 displacements of cell i in the x plane
    dy = NaN(size(prev_pos,1),6); %Last 6 displacements of cell i in the y plane
    tracks = cell(max(slices)-1,1);
    tracks{1} = [prev_pos(:,1:3) ones(cells_in_slice(1),1) prev_pos(:,4) pos(slices == min(slices))];       
    
    for slice = 1:max(slices) - 1;
        if(DEBUG)
            fprintf(repmat('\b', 1, 9));
            fprintf('%04d/%04d', slice, max(slices));
        end
        
        %Calculate Distance Matrix for next slice
        next_pos = [x(slices == slice + 1) y(slices == slice + 1) o(slices == slice+1)];
        [D,disD,disO,disX,disY,disp,modl] = linearMotionVelocityPropagationCostMatrix(prev_pos(:,1:2),next_pos(:,1:2),dx,dy,o(slices == slice),o(slices == slice + 1));
        
        if(isempty(D))
            continue
        end
        
        %Link Frames
        [p12,p21] = lap(D,-1,1,1,END_COST);

        %Update the states
        positions =  zeros(cells_in_slice(slice + 1),5);
        end_cost = 0;
        
        prev_dx = dx;
        prev_dy = dy;
        dx = NaN(cells_in_slice(slice+1),6);
        dy = NaN(cells_in_slice(slice+1),6);
        pos_i = pos(slices == slice + 1);
		for r = 1:cells_in_slice(slice + 1)
			if(p12(r) <= cells_in_slice(slice))
                                %    x,y                i
                positions(r,:) = [next_pos(r,:), prev_pos(p12(r),4), pos_i(r)];
                
                %Reshuffel displacements accroding to the new tracks and
                % rotate the window      
                dx(r,1:end-1) = prev_dx(p12(r),2:end);
                dy(r,1:end-1) = prev_dy(p12(r),2:end);
                
                dx(r,end) = next_pos(r,1) - prev_pos(p12(r),1);
                dy(r,end) = next_pos(r,2) - prev_pos(p12(r),2);
                
                end_cost = max(end_cost,D(r,p12(r)));
			else
				%If not linked to a state, create a new one
                nTracks = nTracks + 1;
                positions(r,:) = [next_pos(r,:), nTracks, pos_i(r)];
			end
        end

        %                   [x y o frame id]
        tracks{slice + 1} = [positions(:,1:3) repmat(slice + 1,cells_in_slice(slice + 1),1) positions(:,4:5)];
       
        %%Extract the distributions of \delat_(\tau,\tau+1) and \theta_(\tau,\tau+1)
        i = p21(p21(1:cells_in_slice(slice)) <= cells_in_slice(slice+1));
        j = p12(i);
        deltaD{slice} = disD(sub2ind(size(disD),i,j));
        deltaO{slice} = disO(sub2ind(size(disO),i,j));
        displacement{slice} = disp(:,j)';
        deltaX{slice} = disX(sub2ind(size(disD),i,j));
        deltaY{slice} = disY(sub2ind(size(disO),i,j));
        model{slice} = modl(sub2ind(size(disO),i,j));
        V = sqrt( (next_pos(i,1) - prev_pos(j,1)).^2 +...
                  (next_pos(i,2) - prev_pos(j,2)).^2);
        deltaV{slice} = V - sqrt(sum(disp(:,j).^2))';
        
        %%Update no linking cost
        END_COST = end_cost * END_COST_MUTIPLIER;
        if(END_COST == 0)
            END_COST = NO_LINKING_COST(slice - 1);
        end
        NO_LINKING_COST(slice) = END_COST;
        
        if(range(NO_LINKING_COST(1:slice)) > 50 && NO_LINKING_COST(slice) >  38)
            warning(['Suspected no-linking-cost overflow. Range: ', num2str(range(NO_LINKING_COST(1:slice))), ', Current: ', num2str(END_COST)]);
        end
        
        %Prepare for next slice 
        prev_pos = positions;
    end
    
    if(DEBUG)
        fprintf('\n    Finished in %03.2f Seconds\n',toc);
    end

    %% Generate a list of trajectories
    if(DEBUG)
        fprintf('Generating Preliminary Tracks \n');
        tic
    end
    

    tracks = concatenateTracks(cell2mat(tracks),NOISE_CUTOFF);

    
    if(DEBUG)
        fprintf('   Finished in %03.2f Seconds\n',toc);
    end
    
    if(NO_GAP_CLOSING)
        t = concatenateTracks(tracks,MIN_TRACK_LENGTH);
        
        clear tracks;
        tracks.x = t(:,1);
        tracks.y = t(:,2);
        tracks.o = t(:,3);
        tracks.frame = t(:,4);
        tracks.id = t(:,5);
        tracks.index = t(:,6); 
        return
    end
        
    %% Populate Gap Closing Matrix
    %Extract Strating and ending times
    if(DEBUG)
        fprintf('Populating Gap Closing Matrix \n');
        tic
    end
    
    [track_starts,track_ends,dist_traveled,o_change,track_end_angle] = get_track_data(tracks);
    
    [pdf_dist, dist_bins, pdf_angle, angle_bins, pdf_o, o_bins] = calculate_pdfs(MAX_SLICE_GAP,dist_traveled,o_change,track_end_angle);
    
    [dist_pdfs,  dist_bin_width,  dist_max_edge,  dist_min_edge, dist_n_bins,  dist_time_offset]  = combine_cell_pdfs(pdf_dist,dist_bins);
    [angle_pdfs, angle_bin_width, angle_max_edge, angle_min_edge, angle_n_bins, angle_time_offset] = combine_cell_pdfs(pdf_angle,angle_bins);
    [o_pdfs,     o_bin_width,     o_max_edge,     o_min_edge, o_n_bins, o_time_offset]     = combine_cell_pdfs(pdf_o,o_bins);
    
    timeMatrix = pdist2(track_ends(:,4),track_starts(:,4),@(x,y) x - y);
    oMatrix    = pdist2(track_ends(:,3),track_starts(:,3),@(x,y) atan2d(sind(2*x - 2*y),cosd(2*x-2*y)) / 2);
    dMatrix    = pdist2(track_ends(:,1:2),track_starts(:,1:2));
    pMatrix    = pdist2(track_ends,track_starts,@find_angle); 
    
    D = ones(size(timeMatrix)) * -1;
    for j = 1:size(timeMatrix,2)
        Trow = timeMatrix(:,j);
        Drow = dMatrix(:,j);
        Prow = pMatrix(:,j);
        Orow = oMatrix(:,j);
        
        Tfilter = Trow > 0 & Trow < MAX_SLICE_GAP + 1;
        
        Dcosts = Inf(size(timeMatrix,1),1);
        Pcosts = Inf(size(timeMatrix,1),1);
        Tcosts = Inf(size(timeMatrix,1),1);
        Ocosts = Inf(size(timeMatrix,1),1);
        Dfilter = false(size(timeMatrix,1),1);
        Pfilter = false(size(timeMatrix,1),1);
        Ofilter = false(size(timeMatrix,1),1);
        
        %-- Filter out values that are outside the the acceptable linking values
        if(~any(Tfilter))
            continue
        end
        
        Dfilter(Tfilter) = Drow(Tfilter) < min(dist_max_edge(Trow(Tfilter))',MAX_GAP);
        Pfilter(Tfilter) = Prow(Tfilter) < angle_max_edge(Trow(Tfilter))' & ...
                           Prow(Tfilter) > angle_min_edge(Trow(Tfilter))'; 
        Ofilter(Tfilter) = Orow(Tfilter) < o_max_edge(Trow(Tfilter))' & ...
                           Orow(Tfilter) > o_min_edge(Trow(Tfilter))';
        
        if(~any(Ofilter) || ~any(Pfilter) | ~any(Dfilter))
            continue
        end
        
        %-- Calculate the remaning gap closing costs using the PDFs
        Dcosts(Dfilter) = dist_pdfs(dist_time_offset(Trow(Dfilter))'   + ...
            c_bin(Drow(Dfilter),dist_min_edge(Trow(Dfilter))',dist_n_bins(Trow(Dfilter))',dist_bin_width(Trow(Dfilter))'));
        
        Pcosts(Pfilter) = angle_pdfs(angle_time_offset(Trow(Pfilter))' + ...
            c_bin(Prow(Pfilter),angle_min_edge(Trow(Pfilter))',angle_n_bins(Trow(Pfilter))',angle_bin_width(Trow(Pfilter))'));
                       
        Ocosts(Ofilter) = o_pdfs(o_time_offset(Trow(Ofilter))' + ... %Which pdf to use
                                 c_bin(Orow(Ofilter),o_min_edge(Trow(Ofilter))',o_n_bins(Trow(Ofilter))',o_bin_width(Trow(Ofilter))')); %which bin in that pdf
        Tcosts(Tfilter) = poisspdf(Trow(Tfilter),0.5);
        
        D(:,j) = Dcosts + Pcosts + Ocosts + Tcosts;
    end
    D(isinf(D)) = -1;
    
    no_link_cost = prctile(D(D ~= -1),90);
    
    if(DEBUG)
        fprintf('   90th %%tile no link cost: %03.2f \n',no_link_cost);
    end
    
    if(DEBUG)
        AJ = D;
        AJ(D == -1) = 0;
        fprintf('   # of possible gap closures %03.2f \n',sum(sum(sum(AJ > 0,1) > 0)));
        fprintf('   max closure cost %03.2f \n',max(D(:)));
    end
    
    if(DEBUG)
        fprintf('   Finished in %03.2f Seconds\n',toc);
    end
    
    %% Solve the LAP for gap closing
    if(DEBUG)
        fprintf('Solving Gap Closing \n');
        fprintf('   Number of preliminary tracks: %u\n',max(tracks(:,5)));
        tic;
    end
    
    %[target] = lapjv(C,1);
    [p12,~] = lap(D,-1,1,1,no_link_cost);
   
    if(DEBUG)
        fprintf('   Finished in %03.2f Seconds\n',toc);
    end
        
    %% Apply the changes to the sparse matrix
    if(DEBUG)
        fprintf('Apply Gap Closing \n');
        tic
    end
    
    id = tracks(:,5);
    for r = 1:max(id)
        if(p12(r) <= max(id))
            tracks(id == p12(r),5) = r;
        else
            %No valid gap to fill found
        end
    end
    
     if(DEBUG)
        fprintf('   %u Trajectories After Gap Closing \n',length(unique(tracks(:,5))));
    end
    

    t = concatenateTracks(tracks,MIN_TRACK_LENGTH);
    
    clear tracks;
    tracks.x = t(:,1);
    tracks.y = t(:,2);
    tracks.o = t(:,3);
    tracks.frame = t(:,4);
    tracks.id = t(:,5);
    tracks.index = t(:,6); 
    
    if(DEBUG)
        fprintf('   Finished in %03.2f Seconds\n',toc);
        fprintf('DONE \n');
        fprintf('   %u Trajectories \n',length(unique(tracks.id)));
    end
end

%Calculates cost matrix by assuming the error in predicted cells state [x,y,theta] and measured
% states are well represented by a distribution with a mean of zero.
function [D,disD,disO,disX,disY,displacement, model] = linearMotionVelocityPropagationCostMatrix(cur_slice,next_slice,dx,dy,o,o_next) 
    sigma_o = 0.328035371823715; %radians
    sigma_d = 1.936155917204507; %pixels
    fallback_disp = 3.577341157940233; %pixels
    
    %%Movement Propagation
    %Possible positions
    if all(isnan(dx(:))) %initilization if nothing can be approximated
        displacement = repmat(fallback_disp,2,size(dx,1)) .* [cosd(o) sind(o)]';
    else
        displacement = [dx(:,end) dy(:,end)]';
        displacement(1,isnan(dx(:,end))) = nanmean(sqrt(dx(:,end).^2)) * cosd(o(isnan(dx(:,end))));
        displacement(2,isnan(dy(:,end))) = nanmean(sqrt(dy(:,end)).^2) * sind(o(isnan(dx(:,end))));
    end
   
    forward = cur_slice' + displacement;
    reverse = cur_slice' + displacement * -1;
    stopped = cur_slice'; 

    %%Displacements between cells
    dx = zeros(size(next_slice,1),size(reverse,2),3);
    dx(:,:,1) = pdist2(next_slice(:,1),forward(1,:)',@(x,y) (x-y));
    dx(:,:,2) = pdist2(next_slice(:,1),reverse(1,:)',@(x,y) (x-y));
    dx(:,:,3) = pdist2(next_slice(:,1),stopped(1,:)',@(x,y) (x-y));

    dy = zeros(size(next_slice,1),size(reverse,2),3);
    dy(:,:,1) = pdist2(next_slice(:,2),forward(2,:)',@(x,y) (x-y));
    dy(:,:,2) = pdist2(next_slice(:,2),reverse(2,:)',@(x,y) (x-y));
    dy(:,:,3) = pdist2(next_slice(:,2),stopped(2,:)',@(x,y) (x-y));
    
    %Total displacement between cells for each movement model
    costMatrix = sqrt(dx.^2 + dy.^2);
    
    %%Cost Calculations
    oMatrix = pdist2(o_next,o,@(x,y) atan2(sind(2*x - 2*y),cosd(2*x-2*y)) / 2);
    
    Po = (oMatrix).^2./(2*sigma_o^2) - log(1/(sigma_o*sqrt(2*pi)));
    Pd = (costMatrix).^2./(2*sigma_d^2) - log(1/(sigma_d*sqrt(2*pi))); %bsxfun(@(x,y) x.^2 ./ (2 * y.^2),costMatrix,standard_dev') - log(1/(50*sqrt(2*pi)));
    
    %%Choose the min cost
    [D,model] = min(bsxfun(@plus,Pd,Po),[],3);
    
    %%Calculate Stats
    [X,Y] = meshgrid(1:size(D,1),1:size(D,2));
    X = X'; Y = Y';
    disD = reshape(costMatrix(sub2ind(size(costMatrix),X(:),Y(:),model(:))),size(D));
    disO = oMatrix;
    disX = reshape(dx(sub2ind(size(dx),X(:),Y(:),model(:))),size(D));
    disY = reshape(dy(sub2ind(size(dy),X(:),Y(:),model(:))),size(D));
end

%% Helper Functions %%

%Calculates the bin x would fall in given that bins start at min_edge
% and are of width bin_width
function [bin] = c_bin(x,min_edge,n_bins,bin_width)
    %the min is required to catch floating point roundoff error that can
    % cause bin > n_bins. 
    bin = min(floor((x - min_edge) ./ bin_width),n_bins-1);
end

%Pulls the start, end, distance traveled, change in bearing, and final orientation
% of each trajecotry from the vector of trajectories
function [track_starts,track_ends,dist_traveled,o_change,track_end_angle] = get_track_data(tracks)
    id = tracks(:,5);
    
    dist_traveled = cell(max(id),1);
    o_change      = cell(size(dist_traveled));
    track_end_angle = cell(size(dist_traveled)); %Angle of the track at its end
    
    track_starts  = NaN(max(id),6); % [X Y o START_SLICE]
    track_ends    = NaN(max(id),6); % [X Y o END_SLICE sxi syi]
    for i = 1:max(id) %Populates the variables above
        t = tracks(id == i,1:4); %[X Y o START_SLICE]
        track_starts(i,1:4) = t(1,:);
        track_ends(i,1:4)   = t(end,:);  
        
        track_a = NaN(max(size(t,1) - 5,1),5);
        dist_t  = NaN(size(track_a));
        o_t     = NaN(size(track_a));

        for j = 1:max(size(t,1) - 5,1)
            
            range1 = j:min(j+5,size(t,1) - 1);
            range2 = j + 1:min(j+6,size(t,1));
            
            if(size(t,1) > 3)
                sxi = t(range1,1) - t(range2,1);
                syi = t(range1,2) - t(range2,2);
                track_starts(i,5) = sxi(end);
                track_starts(i,6) = syi(end);

                lvl = sqrt(sxi.^2 + syi.^2);
                xy = sxi(1) .* sxi(2:end) + syi(1) .* syi(2:end);

                track_a(j,1:length(range1) - 1) = real(acosd( xy ./ ( lvl(1) .* lvl(2:end) ) ) .* ...
                        sign( sxi(1) .* syi(2:end) + syi(1) .* sxi(2:end) ));
            end

            range2 = j + 1:min(j+5,size(t,1));   
            dist_t(j,1:min(size(t,1) - 1,5)) = sqrt((t(j,1) - t(range2,1)).^2 + (t(j,2) - t(range2,2)).^2);
            o_t(j,1:min(size(t,1) - 1,5)) = atan2d(sind(2*t(j,3) - 2*t(range2,3)),cosd(2*t(j,3) - 2*t(range2,3))) / 2;
            
        end %end window loop
        
        dist_traveled{i} = dist_t;
        o_change{i} = o_t;
        track_end_angle{i} = track_a;
    end
    
    dist_traveled   = cell2mat(dist_traveled);
    o_change        = cell2mat(o_change);
    track_end_angle = cell2mat(track_end_angle);
end

%Calculates the probability density functions for
%   pdf_dist{t} is the distance between the ending of track I and the beginning of track j at time t
%   pdf_angel{t} is angular bearing change between end of track I and begging of track j at time t
%   pdf_o{t} is the difference in the orientation between the cell at the end of track I and the beginning of track j at time t
function [pdf_dist, dist_bins, pdf_angle, angle_bins, pdf_o, o_bins] = calculate_pdfs(max_slice_gap,dist_traveled,o_change,track_end_angle)
    pdf_dist = cell(1,max_slice_gap);
    pdf_o = cell(1,max_slice_gap);
    pdf_angle = cell(1,max_slice_gap);
    dist_bins = cell(1,max_slice_gap);
    o_bins = cell(1,max_slice_gap);
    angle_bins = cell(1,max_slice_gap);
    for ti = 1:max_slice_gap
        [pdf_dist{ti}, dist_bins{ti}] = histcounts(dist_traveled(~isnan(dist_traveled(:,ti)),ti),calcnbins(dist_traveled(~isnan(dist_traveled(:,ti)),ti)));
        [pdf_o{ti}, o_bins{ti}] = histcounts(o_change(~isnan(o_change(:,ti)),ti),calcnbins(o_change(~isnan(o_change(:,ti)),ti)));
        [pdf_angle{ti}, angle_bins{ti}] = histcounts(track_end_angle(~isnan(track_end_angle(:,ti)),ti),calcnbins(track_end_angle(~isnan(track_end_angle(:,ti)),ti)));
        
        pdf_dist{ti}(pdf_dist{ti} == 0) = min(pdf_dist{ti}(pdf_dist{ti} > 0));
        pdf_o{ti}(pdf_o{ti} == 0) = min(pdf_o{ti}(pdf_o{ti} > 0));
        pdf_angle{ti}(pdf_angle{ti} == 0) = min(pdf_angle{ti}(pdf_angle{ti} > 0));
        
        pdf_dist{ti} = -log(pdf_dist{ti} ./ sum(pdf_dist{ti}));
        pdf_o{ti} = -log(pdf_o{ti} ./ sum(pdf_o{ti}));
        pdf_angle{ti} = -log(pdf_angle{ti} ./ sum(pdf_angle{ti}));       
    end
end

%Flattens the cell pdfs from calculate_pdfs to a vector by concatonating the cells togeather
function [pdfs, bin_width, max_edge, min_edge, n_bins, time_offset] = combine_cell_pdfs(pdf,bins)
    pdfs  = cell2mat(pdf);
    bin_width = cellfun(@(x) diff(x(1:2)),bins);
    max_edge  = cellfun(@(x) x(end),bins);
    min_edge  = cellfun(@(x) x(1),bins);
    n_bins = cellfun(@(x) length(x),bins) - 1;
    time_offset = cumsum(cellfun('length',{1, pdf{1:end-1}}));
end

function [D] = find_angle(ts,te)

        tx_start = ts(1,5);
        ty_start = ts(1,6);
        lvl_start = sqrt(tx_start^2 + ty_start^2);
        
        tx_end = te(:,1) - ts(:,1);
        ty_end = te(:,2) - ts(:,1);
        lvl_end = sqrt(tx_end.^2 + ty_end.^2);
        
        xy = tx_start .* tx_end + ty_start .* tx_end;

        D = real(acosd( xy ./ ( lvl_start .* lvl_end ) ) .* ...
                sign( tx_start .* ty_end + ty_start .* tx_end ));
end

function [tt] = concatenateTracks(tt,min_track_length)
    if(min_track_length == 0)
        return
    end
    
    unv = unique(tt(:,5));
    frameCount = histc(tt(:,5),unv);
    short = find(frameCount <= min_track_length);
    tt(ismember(tt(:,5),unv(short)),:) = [];
    count = 1;
    for i = unique(tt(:,5))'
        i;
        tt(tt(:,5) == i,5) = count;
        count = count + 1;
        count;
    end
end
%% Log
%$Log: LAPtracking.m,v $
%Revision 1.8  2014/05/20 17:03:31  cotter
%linearMotionVelocityPropagationCostMatrix now applies the linear motion cost prior to taking the min of the propagations. linear mostion cost calcuation cleaned up for all lin motion cost calculations.
%
%Revision 1.7  2014/05/09 12:34:08  cotter
%Removed scaler from no_linking_cost/adjusted angle penalty for gap closing
%
%Revision 1.6  2014/05/08 17:48:12  cotter
%Fixed MANY Bugs found in the gap closing code, Implemented a costMatrix that accounts for angle for frame to frame linking and gap closing. Add Function Heading comments. Added more optional set variables.
%
%Revision 1.5  2014/04/19 14:51:24  cotter
%Rewrote The Gap linking matrx to better fit what is decribed in Jaqman et al.
%
%Revision 1.4  2014/03/26 02:45:56  cotter
%Added data driven cutoff and no-linking costs. Fixed Bug in determining the matrix for the GAP linking.
%
%Revision 1.3  2014/03/25 18:01:24  cotter
%Added Kalman Next Position Estimations for linking stage. Fixed indexing of final tracks
%
%Revision 1.2  2014/03/24 20:24:37  cotter
%Added GAP Closing Added optional input params Moved to LAP solver from u-track (MUCH Faster)
%
%Revision 1.1  2014/03/23 22:06:35  cotter
%Frame to Frame Tracking finished and tested.
%


