function [Tab ] = CompileVariables( AllData )
    % [AllDataTable ] = CompileVariables( AllData )
    % Generates  A table combining all the run data from AllData into one
    % structure. This was used to run the simulations and generate all figures.
    % The resulting table of runs is temporally aligned according to the values
    % of AllData{i}.AlignedStart and AllData{i}.AlignedStop
    % 
    % Each row represents one cell run, with columns (variables) representing the
    % quantified data relating to that run. Most variables are followed by
    % a number (i.e. Rd1, Rd0). The number 1 indicates the value for that run
    % (i.e. AllDataTable.Rd1(1) is the run duration for the first run in the table).
    % Variables that end in a 0 (i.e. Rd0) is the value for the run the same cell
    % performed prior to the current run (Rd1). In the same way Variables that end
    % in a 2 indicate variables for the run the the same cell performed after the
    % current run. See figure S1 in Cotter et al 2017 PNAS for a visual explanation.
    % 
    % Unnumbered variables relate to the run represented by that row
    % (the current run).
    % 
    % The table consists of the variables:
    % 
    % Rd (float): Run distance (um)
    % 
    % Rt (float): Run duration (min)
    % 
    % Rs (float): Run speed (um/min)
    % 
    % rho (float): Local Cell density at start of the run (cells/um^2)
    % 
    % drho (float): Change in cell density from the start to the end of the run (cells/um^2)
    % 
    % TSS (float): Time since the beginning of the movie (Frames [1 frame = 30 seconds])
    % 
    % Dn (float): Distance to the nearest aggregate boundary (um)
    % 
    % phi (float): Angle (radians 0 <= beta < pi/2)
    % 
    % beta (float): Angle between the run and a vector pointing to the center of the neareast aggregate (radians 0 <= beta < pi/2)
    % 
    % startInside (bool): Run start inside of an aggregate 
    % 
    % theta (float) : Angle enclosed between the previous run of a cell and the current run
    % 
    % xstart (float): x coordinate of the start of the run (um)
    % 
    % ystart (float): y coordinate of the start of the run (um)
    % 
    % xstop (float): x coordinate of the end of the run (um)
    % 
    % ystop (float): y coordinate of the end of the run (um)
    % 
    % orientation (float): Angle enclosed between the run vector and the x axis (radians 0 <= beta < pi/2)
    % 
    % state (int): State of the cell during the run (1: Persistent forward, 2: Persistent forward [opposite direction], 3: Nonpersistent )
    % 
    % movie (string): The name of the fluorescent movie the run came from 
    % 
    % set (int): An id unique to the movie the run came from
    % 
    % ksi (float): Mean Nematic alignment orientation at cell location (radians 0 <= beta < pi/2)
    % 
    % ksiN (float): Number of runs that started without the alignment window of the current run, 
    %         0 indicates no other runs were close enough to calculate alignment. The size of the
    %         window is set by:
    %            TIME: the surroudning runs must occur within +/- TIME frames (Set to 10)
    %          RADIUS: the surrounding runs must occur within +/- RADIUS um   (Set to 14)
    % 
    % chi (float): ksi1 - orientation1
    % 
    % neighbor_alignment (float):  mean(cos(2 * (orientation1 - surrounding))) where surrounding
    %    is the cells within the alignment calculation window. 
    % 
    % mean_alignment (float): cos(2 * (orientation1 - T.ksi1))

    Tab = table();
    addpath('Libraries/Utils')
    for j = 1:length(AllData)
        all_runs = AllData{j}.all_runs;
        name = AllData{j}.movie;

        HAS_AGG_DATA = false;
        if(isfield(AllData{j},'agg_tracks'))
            HAS_AGG_DATA = true;

            agg_tracks = AllData{j}.agg_tracks;
            agg_props  = AllData{j}.AggProps;
            first_agg_seen = min(agg_tracks.frame);
            early_agg_cutoff = first_agg_seen + 30;
        end

        p = Progress(max(all_runs.id));
        T = table();
        for i = unique(all_runs.id)'
            p.d(i);
            runs = subStruct(all_runs,all_runs.id == i,'r');
            if(length(runs.start.x) < 3)
                continue
            end

            %%%
            %Calculate angles. Since we do not have a reversal angle for the first
            %run it is dropped from any further analysis
            %%%
            V1 = [runs.stop.x(1:end-1) - runs.start.x(1:end-1), runs.stop.y(1:end-1) - runs.start.y(1:end-1)]';
            V2 = [runs.stop.x(2:end)   - runs.start.x(2:end),   runs.stop.y(2:end)   - runs.start.y(2:end)]';
            cross_product12 = V1(1,:) .* V2(2,:) - V1(2,:).*V2(1,:);
            costheta=dot(V1,V2) ./ (sqrt(sum(V1.^2)) .* sqrt(sum(V2.^2)));
            
            theta = acos(costheta)' .* sign(cross_product12)';
            theta(imag(theta) > 0) =0; %Catches a rounding error
            theta(imag(theta) < 0) =acos(-1);
            if(any(imag(theta)))
                    error('Imaginary reversal angle detected')
            end
            direc = atan(V1(2,:) ./ V1(1,:));

            %%%
            % If aggreagte data exits, populate variables realted
            % to aggreagtes
            %%%
            if(HAS_AGG_DATA)
                %%%
                % Find the nearest aggregate to the START of the run
                %%%
                D = zeros(length(runs.start.x)-1,1);
                centroids = NaN(size(D,1),2);
                AggProps = array2table(nan(size(D,1),width(agg_props)));
                AggProps.Properties.VariableNames = agg_props.Properties.VariableNames;
                for k = 1:(length(runs.start.x)-1)
                    start_frame = runs.start.frame(k+1);

                    %if(start_frame < early_agg_cutoff)
                    %   at = subStruct(agg_tracks,agg_tracks.frame <= early_agg_cutoff,'r');
                    %else
                       at = subStruct(agg_tracks,agg_tracks.frame == start_frame,'r');
                    %end

                    if(isempty(at.x))
                        D(k) = NaN; 
                    else
                        [D(k),I] = pdist2([at.x at.y],[runs.start.x(k+1) runs.start.y(k+1)],'euclidean','Smallest',1);
                        AggProps(k,:) = agg_props(find(agg_props.index == at.index(I),1),:);
                        centroids(k,:) = [at.x(I) at.y(I)];
                    end
                end

                %%%
                %Calculate angles related to aggregates
                %%%
                V3 = [centroids(:,1) - runs.start.x(2:end),centroids(:,2) - runs.start.y(2:end)]';

                cross_product13 = V1(1,:) .* V3(2,:) - V1(2,:).*V3(1,:);
                phi = acos(dot(V1,V3) ./ (sqrt(sum(V1.^2)) .* sqrt(sum(V3.^2))))' .* sign(cross_product13)';
                cross_product23 = V2(1,:) .* V3(2,:) - V2(2,:).*V3(1,:);
                beta = acos(dot(V2,V3) ./ (sqrt(sum(V2.^2)) .* sqrt(sum(V3.^2))))' .* sign(cross_product23)';


        %             figure
        %             for k = 1:length(V1)
        %                 clf
        %                 line([0 V1(1,k)],[0 V1(2,k)],'Color','b');
        %                 line([0 cos(direc(k))],[0 sin(direc(k))],'Color','m','LineStyle','--')
        %                 title(['orientation = ' num2str(rad2deg(direc(k)))]);
        %                 w = waitforbuttonpress;
        %             end
        %             
                %%%
                %Find the radius of the aggreagte along the vector between
                %the beginning of the run and the aggreate centroid
                %%%
                start_dist_to_agg_norm = NaN(length(D),1);
                for k = find(~isnan(D))'
                    ph = AggProps.orientation(k);
                    a = AggProps.majorAxis(k)/2;
                    b = AggProps.minorAxis(k)/2;
                    x = centroids(k,1);
                    y = centroids(k,2);

                    R = [ cos(ph)   -sin(ph)
                          sin(ph)   cos(ph)];
                    xy = R*[runs.start.x(k+1)-x; runs.start.y(k+1)-y];
                    xrotated = xy(1);
                    yrotated = xy(2);
                    t = atan2(yrotated,xrotated);
                    start_dist_to_agg_norm(k) = D(k) - (a*b / sqrt((b*cos(t)).^2+(a*sin(t)).^2));

            %                 clf
            %                 orien0 = atan2(V1(2,k),V1(1,k));
            %                 orien1 = atan2(V2(2,k),V2(1,k));
            %                 line([0 V3(1,k)],[0 V3(2,k)],'Color','r');
            %                 line([0 V1(1,k)],[0 V1(2,k)],'Color','b');
            %                 line([0 V2(1,k)],[0 V2(2,k)],'Color','g');
            %                 
            %                 line([0 cos(orien0)],[0 sin(orien0)],'Color','w','LineStyle','--')
            %                 line([0 cos(orien1)],[0 sin(orien1)],'Color','y','LineStyle','--')  
            %                 
            %                 line([0 cos(orien0 + phi(k))] * 0.5,[0 sin(orien0 + phi(k))] * 0.5,'Color','g');
            %                 line([0 cos(orien1 + beta(k))] * 0.25,[0 sin(orien1 + beta(k))] * 0.25,'Color','b');
            %                 legend('Aggregate', ...
            %                        'Run 0', ...
            %                        'Run 1', ...
            %                        'Orinetation 0', ...
            %                        'Orientation 1', ...
            %                        'Phi 0', ...
            %                        'Beta 1')
            %                 axis equal
            %                 xlim([-1 1]);
            %                 ylim([-1 1]);
            %                 title(['\phi_{n-1} = ' num2str(rad2deg(phi(k))) ...
            %                        ', \beta_n = ' num2str(rad2deg(beta(k))) ...
            %                        ', \theta_n = ' num2str(rad2deg(theta(k)))]);
            %                 w = waitforbuttonpress;
                end

                runs = subStruct(runs,1:(length(runs.start.x)-1),'r');
                runs.theta = theta;
                
                runs.orientation = direc';
                
                runs.phi = phi;
                runs.beta = beta;
                runs.start_dist_to_agg_norm = start_dist_to_agg_norm;

                T = [T; table(... 
                      runs.distance(1:end-1), ... %1
                      runs.distance(2:end), ...
                      runs.length(1:end-1), ...
                      runs.length(2:end), ...
                      runs.speed(1:end-1), ...    %5 
                      runs.speed(2:end), ...     
                      runs.start.density(1:end-1), ...
                      runs.start.density(2:end), ...
                      runs.stop.density(2:end), ...
                      (runs.stop.density(1:end-1) - runs.start.density(1:end-1)) ./ runs.length(1:end-1), ... %10
                      (runs.stop.density(2:end) - runs.start.density(2:end)) ./ runs.length(2:end), ...
                      runs.start.frame(1:end-1), ...
                      runs.start.frame(2:end), ...
                      runs.start_dist_to_agg_norm(1:end-1), ...
                      runs.start_dist_to_agg_norm(2:end), ... %15
                      runs.phi(1:end-1), ...
                      runs.phi(2:end), ...
                      runs.beta(1:end-1), ...
                      runs.beta(2:end), ...
                      runs.theta(1:end-1), ...
                      runs.theta(2:end), ...
                      runs.start.x(1:end-1), ...
                      runs.start.x(2:end), ...
                      runs.start.y(1:end-1), ...
                      runs.start.y(2:end), ...
                      runs.stop.x(1:end-1), ...
                      runs.stop.x(2:end), ...
                      runs.stop.y(1:end-1), ...
                      runs.stop.y(2:end), ...
                      runs.orientation(1:end-1), ...
                      runs.orientation(2:end), ...
                      runs.state(1:end-1), ...
                      runs.state(2:end), ...
                      runs.id(1:end-1)) AggProps(1:end-1,:)];
            else
                runs = subStruct(runs,1:(length(runs.start.x)-1),'r');
                runs.theta = theta;
                if(any(imag(theta)))
                    error('Imaginary reversal angle detected')
                end
                runs.orientation = direc';
                
                T = [T; table(... 
                      runs.distance(1:end-1), ... %1
                      runs.distance(2:end), ...
                      runs.length(1:end-1), ...
                      runs.length(2:end), ...
                      runs.speed(1:end-1), ...    %5 
                      runs.speed(2:end), ...     
                      runs.start.density(1:end-1), ...
                      runs.start.density(2:end), ...
                      runs.stop.density(2:end), ...
                      (runs.stop.density(1:end-1) - runs.start.density(1:end-1)) ./ runs.length(1:end-1), ... %10
                      (runs.stop.density(2:end) - runs.start.density(2:end)) ./ runs.length(2:end), ...
                      runs.start.frame(1:end-1), ...
                      runs.start.frame(2:end), ...
                      runs.theta(1:end-1), ...
                      runs.theta(2:end), ...
                      runs.start.x(1:end-1), ...
                      runs.start.x(2:end), ...
                      runs.start.y(1:end-1), ...
                      runs.start.y(2:end), ...
                      runs.stop.x(1:end-1), ...
                      runs.stop.x(2:end), ...
                      runs.stop.y(1:end-1), ...
                      runs.stop.y(2:end), ...
                      runs.orientation(1:end-1), ...
                      runs.orientation(2:end), ...
                      runs.state(1:end-1), ...
                      runs.state(2:end), ...
                      runs.id(1:end-1))];
            end
        end
        
        T.movie = repmat({name},size(T,1),1);
        T.set   = repmat(j,size(T,1),1);
        
        %Align movie and normailze timepoints
        StartFrame = AllData{j}.AlignedStart;
        StopFrame = AllData{j}.AlignedStop;
        T(T{:,12} < StartFrame,:) = [];
        T(T{:,13} > StopFrame,:) = [];
        T{:,12} = T{:,12} - StartFrame;
        T{:,13} = T{:,13} - StartFrame;
        
        Tab = [Tab; T];
        p.done();
    end

    if(HAS_AGG_DATA)
        Tab.Properties.VariableNames = [{
                                    'Rd0','Rd1', ... %1
                                    'Rt0','Rt1', ... %3
                                    'Rs0','Rs1', ... %5
                                    'rho0','rho1','rho2', ... %7
                                    'drho0','drho1', ... %10
                                    'TSS0','TSS1', ... %12
                                    'Dn1' ,'Dn2' , ... %14
                                    'phi0','phi1', ... %16
                                    'beta1','beta2', ... %18
                                    'theta1','theta2', ...
                                    'xstart0','xstart1', ...
                                    'ystart0','ystart1', ...
                                    'xstop0','xstop1', ...
                                    'ystop0','ystop1', ...
                                    'orientation0', ...
                                    'orientation1' ,...
                                    'state0', ...
                                    'state1', ...
                                    'id', ...
                                    }, ...
                                    agg_props.Properties.VariableNames, ...
                                    {
                                    'movie', ...
                                    'set'}
                                    ];
    else
        Tab.Properties.VariableNames = [{
                            'Rd0','Rd1', ... %1
                            'Rt0','Rt1', ... %3
                            'Rs0','Rs1', ... %5
                            'rho0','rho1','rho2', ... %7
                            'drho0','drho1', ... %10
                            'TSS0','TSS1', ... %12
                            'theta1','theta2', ...
                            'xstart0','xstart1', ...
                            'ystart0','ystart1', ...
                            'xstop0','xstop1', ...
                            'ystop0','ystop1', ...
                            'orientation0', ...
                            'orientation1' ,...
                            'state0', ...
                            'state1', ...
                            'id', ...
                            'movie', ...
                            'set'}
                            ];
    end
    
    %% Add Neighbor Alignment Data
    RADIUS = 14;
    TIME = 10;
    
    Tab.ksi = zeros(height(Tab),1);
    Tab.ksiN = zeros(height(Tab),1);
    Tab.chi = zeros(height(Tab),1);
    
    Tab.neighbor_alignment = zeros(height(Tab),1);
    Tab.mean_alignment = zeros(height(Tab),1);
    
    for curMovie = 1:length(unique(Tab.set))
        filt = Tab.set == curMovie;
        T = Tab(filt,:);
        Mdl = createns([T.xstart1 T.ystart1]);
        [result] = rangesearch(Mdl,[T.xstart1 T.ystart1],RADIUS);

        T.mean_alignment = zeros(height(T),1);
        T.chi = zeros(height(T),1);
        T.mean_alignment = zeros(height(T),1);
        T.ksiN = zeros(height(T),1);
        
        p = Progress(length(result));
        for j = 1:length(result)
            idx = result{j};
            idx = idx(T.TSS1(idx) > T.TSS1(idx(1)) - TIME & T.TSS1(idx) <= T.TSS1(idx(1)) + TIME);
            p.d(j)
            if(length(idx) == 1)
                continue
            end

            cur_cell = T.orientation1(idx(1));
            surrounding = T.orientation1(idx(2:end));

            T.ksi(j) = atan2(sum(sin(2 * surrounding)),sum(cos(2 * surrounding))) / 2; %Mean Nematic alignment orientation at cell location
            T.chi(j) = atan(sin(cur_cell - T.ksi(j))./cos(cur_cell - T.ksi(j)));       % chi = ksi1 - ori1
            T.mean_alignment(j) = cos(2 * (cur_cell - T.ksi(j) ));
            T.neighbor_alignment(j) = mean(cos(2 * (cur_cell - surrounding)));
            T.ksiN(j) = length(surrounding);
        end
        p.done();

        Tab.ksi(filt) = T.ksi;
        Tab.chi(filt) = T.chi;
        Tab.mean_alignment(filt) = T.mean_alignment;
        Tab.neighbor_alignment(filt) = T.neighbor_alignment;
        Tab.ksiN(filt) = T.ksiN;
    end
    
end

