function [] = PlotCell( tracks , cellID )
    %[] = PlotCell( tracks , cellID )
    %
    %Plots the trajectory of cell cellID from tracks 
    %
    %PARAMS
    %   tracks   output from track
    %   cellID   cell id from tracks to plot
    %
    %Author: $Author: cotter $
    %Revision: $Revision: 1.1 $
    
    %
    % Constants
    %
    if(isfield(tracks,'units') & strcmp(tracks.units,'um'))
        DISTCONV = 1;
    else
        DISTCONV = 0.5136;
    end
   
    xcell = tracks.x(tracks.id == cellID);
    ycell = tracks.y(tracks.id == cellID);
    if(isfield(tracks,'state')); mcell = tracks.state(tracks.id == cellID); end
    ocell = tracks.o(tracks.id == cellID);

    figure; hold on; %axcells([0 2000 0 1500]); 
    
    %%plot(xcell .* DISTCONV,ycell .* DISTCONV,':','Color','k');
    
    cmap(1,:) = [0 0 1];
    cmap(2,:) = [1 0 0];
	cmap(3,:) = [1 .8 0];
    
    if(isfield(tracks,'state'))
        multiColorLine(xcell .* DISTCONV,ycell .* DISTCONV,[mcell(2:end); mcell(end)],cmap);
    else
        plot(xcell .* DISTCONV,ycell .* DISTCONV,'o');
    end
    
    if(isfield(tracks,'runVectors'))
        cellVectors = subStruct(tracks.runVectors,tracks.runVectors.id == cellID,'r')
        
        line([cellVectors.start.x(cellVectors.state == 3) cellVectors.stop.x(cellVectors.state == 3)]' .* DISTCONV, ...
             [cellVectors.start.y(cellVectors.state == 3) cellVectors.stop.y(cellVectors.state == 3)]' .* DISTCONV, ...
             'Color','y','LineStyle','-');
        line([cellVectors.start.x(cellVectors.state == 2) cellVectors.stop.x(cellVectors.state == 2)]' .* DISTCONV, ...
            [cellVectors.start.y(cellVectors.state == 2) cellVectors.stop.y(cellVectors.state == 2)]' .* DISTCONV, ...
            'Color','r','LineStyle','-');
        line([cellVectors.start.x(cellVectors.state == 1) cellVectors.stop.x(cellVectors.state == 1)]' .* DISTCONV, ...
            [cellVectors.start.y(cellVectors.state == 1) cellVectors.stop.y(cellVectors.state == 1)]' .* DISTCONV, ...
            'Color','b','LineStyle','-');
    end
    
    %\plot(xcell(I),ycell(I),'color','green');
    axis equal
    xlabel('y (\mum))');
    ylabel('x (\mum)');
    grid on;
end

% $Log: PlotCell.m,v $
% Revision 1.1  2013/08/16 01:10:40  cotter
% _
%