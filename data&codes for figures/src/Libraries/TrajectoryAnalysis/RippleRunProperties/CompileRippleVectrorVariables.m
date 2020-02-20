function [Tab ] = CompileRippleVectrorVariables( AllData )
    % [Tab ] = CompileVariables( AllData )
    %
    % Tab.
    %   

    Tab = table();

    for j = 1:length(AllData)
    all_runs = AllData{j}.all_ripple_vectors;
    name = AllData{j}.set;
    
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
        theta = acos(dot(V1,V2) ./ (sqrt(sum(V1.^2)) .* sqrt(sum(V2.^2))))' .* sign(cross_product12)';
        theta(imag(theta) > 0) = acos(1); %Catches a rounding error
        direc = atan(V1(2,:) ./ V1(1,:));

        runs = subStruct(runs,1:(length(runs.start.x)-1),'r');
        runs.theta = theta;
        runs.orientation = direc'; 
        
        T = [T; table(... 
                  runs.distance(1:end-1), ... %1
                  runs.distance(2:end), ...
                  runs.length(1:end-1), ...
                  runs.length(2:end), ...
                  runs.speed(1:end-1), ...    %5 
                  runs.speed(2:end), ...
                  runs.start.frame(1:end-1), ...
                  runs.start.frame(2:end), ...
                  runs.theta(1:end-1), ...
                  runs.theta(2:end), ...      %10
                  runs.start.x(1:end-1), ...
                  runs.start.x(2:end), ...
                  runs.start.y(1:end-1), ...
                  runs.start.y(2:end), ...
                  runs.stop.x(1:end-1), ...   %15
                  runs.stop.x(2:end), ...
                  runs.stop.y(1:end-1), ...
                  runs.stop.y(2:end), ...
                  runs.orientation(1:end-1), ...
                  runs.orientation(2:end))];  
    end
    T.movie = repmat({name},size(T,1),1);
    T.set = repmat(j,size(T,1),1);

    %Align movie and normailze timepoints
    StartFrame = AllData{j}.AlignedStart;
    StopFrame = AllData{j}.AlignedStop;
    T(T{:,7} < StartFrame,:) = [];
    T(T{:,8} > StopFrame,:) = [];
    T{:,7} = T{:,7} - StartFrame;
    T{:,8} = T{:,8} - StartFrame;
    %        
    Tab = [Tab; T];
    p.done();
    end

    Tab.Properties.VariableNames = {
                                'Rd0','Rd1', ... %1
                                'Rt0','Rt1', ... %3
                                'Rs0','Rs1', ... %5
                                'TSS0','TSS1', ... %7
                                'theta1','theta2', ... %9
                                'xstart0','xstart1', ... %11
                                'ystart0','ystart1', ... %13
                                'xstop0','xstop1', ... %15
                                'ystop0','ystop1', ... %17
                                'orientation0', ... %19
                                'orientation1' ,... %20
                                'movie', ... %21
                                'set', ...
                                };
    
end

