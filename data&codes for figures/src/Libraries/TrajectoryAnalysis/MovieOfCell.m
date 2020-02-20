function [F] = MovieOfCell( tracks , cellID, in_fname, varargin)
    %[F] = MovieOfCell( tracks , cellID, in_fname, varargin)
    
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
    FRAMERATE = 30; % 30 seconds / frame
    DISTCONV = 0.5136; % um / pixel
        
    f = tracks.frame(tracks.id == cellID);
    xcell = tracks.x(tracks.id == cellID);
    ycell = tracks.y(tracks.id == cellID);
    axis_limits = [min(xcell) - 50 max(xcell) + 50 min(ycell) - 50 max(ycell) + 50];
        
    hfig = figure('Visible','off','Color','w');
        
    start_frame = min(f);
    F(length(f)-1).cdata = [];
    F(length(f)-1).colormap= [];
    index = 1;
    for fi = f(2:end)'
        fi
        
        clf, hold on;
        
        if(~isempty(in_fname))
            A = imread(in_fname,fi);
            A = A(:,:,1);
            %imagesc(bpass(A,1,10));
            imagesc(A);
            colormap('gray');
        end
        
        plot(xcell(f < fi),ycell(f < fi),'w');
        
        axis(axis_limits);
        xlabel('x (\mum)');
        ylabel('y (\mum)');
        set(gca,'YDir','Reverse')
        %set(gca,'YTickLabel',round(get(gca,'YTick') .* DISTCONV));
        %set(gca,'XTickLabel',round(get(gca,'XTick') .* DISTCONV));
        grid on;
        
        drawnow;
        F(index) = getframe(hfig);
        
        index = index + 1;
    end   
end

function circle(x,y,r,c)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,c);
end

% $Log: PlotCell.m,v $
% Revision 1.1  2013/08/16 01:10:40  cotter
% _
%