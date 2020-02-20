function h = visRuns(tracks,varargin)
% h = visRuns(tracks, ...)
%
% Plots all trajectories in tracks
%
% params
% tracks
%
% ...
%   'NoTracks' supresses plotting tracks
%   'NoRunVectors' supresses plotting run vectors
%
    PLOT_TRACKS = true;
    PLOT_RUN_VECTORS = true;
    
    for i = varargin
        if(strcmp(i{1},'NoTracks'))
            PLOT_TRACKS = false;
        elseif strcmp(i{1},'NoRunVectors')
            PLOT_RUN_VECTORS = false;
        end
    end
    
    if(isfield(tracks,'units') & strcmp(tracks.units,'um'))
        DISTCONV = 1;
    else
        DISTCONV = 0.5136;
    end
    
    if(nargout == 0)
        figure
    else
        h = 1;
    end
    
    hold on;
    axis([0 1920 0 1440] .* 0.5136)
    xlabel('x (\mum)');
    ylabel('y (\mum)');
    set(gca,'YDir','reverse');

    if PLOT_RUN_VECTORS && isfield(tracks,'runVectors')
        runVectors = tracks.runVectors;
        id = runVectors.id;
        for i = unique(id)'
            cellVectors = subStruct(runVectors,id == i & runVectors.length > 3.5,'r');
            
            stopped = cellVectors.state == 3;
            line([cellVectors.start.x(stopped) cellVectors.stop.x(stopped)]',[cellVectors.start.y(stopped) cellVectors.stop.y(stopped)]','Color',color_chooser(i),'LineStyle',':');
            line([cellVectors.start.x(~stopped) cellVectors.stop.x(~stopped)]',[cellVectors.start.y(~stopped) cellVectors.stop.y(~stopped)]','Color',color_chooser(i));

            %text(Vi(1,1),Vi(1,2),num2str(i),'FontSize',12,'Color',color_chooser(i));
        end
    end

    if PLOT_TRACKS
        for i = unique(tracks.id)'
            xi = tracks.x(tracks.id == i) .* DISTCONV;
            yi = tracks.y(tracks.id == i) .* DISTCONV;
            plot(xi,yi,'-','Color',color_chooser(i))
            %text(xi(1),yi(1),num2str(i),'FontSize',12,'Color',color_chooser(i));
        end
    end

%     if nargin > 1
%         %Aggreagte Edges
%         for i = 1:length(boundaries)
%           plot(boundaries{i}(:,2) .* DISTCONV,boundaries{i}(:,1) .* DISTCONV,'LineWidth',2,'Color','b');
%         end
%     end
end
