function [F] = createTrackingMovie(tiff_stack,fname_out,tracks)
% [] = createTrackingMovie(tiff_stack,fname_out,tracks)
%
% creates a tiff stack with the cells in the floresent channel
% labeled with their ids
%
% INPUT
% tiff_stack an RGB stack containing the floresent channel images in the R channel
% fname_out  the filename of the output stack
% tracks     the tracks from LAPTracking function where
%             1 <= min(tracks(:,3)). If it is > 1, min(tracks(:,3) is
%             assigned to frame 1 of tiff_stack and so on.
%             It is also possible to pass posList instead of tracks.
%             In that case the cell ids are left off the output tiff stack
% 
% tiff_stack must have the same number of frames as tracks
% requires ExtractTracks()
%
% Author: Chris Cotter (cotter@uga.edu)

    %Setup a figure to build the image in
    hf = figure('visible','off','Renderer','opengl','Position',[0 0 1920 1440]);
    
    x = tracks.x;
    y = tracks.y;
    o = tracks.o;
    frame = tracks.frame;
    showIds = false;
   
    if(isfield(tracks,'units'))
        if(strcmp(tracks.units,'um'));
            x = x ./ 0.5136;
            y = y ./ 0.5136;
        end
    end
    
    if(isfield(tracks,'id'))
        showIds = true;
        id = tracks.id;
    end
    
    %frame_offset = min(frame);
    
    %Initilize frame struct
    F(max(frame)).cdata = [];
    F(max(frame)).colormap= [];
    info = imfinfo(tiff_stack);
    count = 1;
    p = Progress(max(frame));
    
    v = VideoWriter('fname_out');
    for j = 1:max(frame)
        p.d(j);
        x_frame = x(frame == j);
        y_frame = y(frame == j);
        o_frame = o(frame == j);

        A = imread(tiff_stack,j,'info',info);
        %A = imresize(A,[740 986]);
        %

        hold on;
        colormap('gray'), imagesc(bpass(A,1,10));
        %colormap('gray'), imagesc(A);
        plot(x_frame,y_frame,'ob');
        %for i = 1:length(x_frame)
        %    [x2, y2] = rotate(x_frame(i),y_frame(i),o_frame(i),5);
        %    [x1r, y1r] = rotate(x_frame(i),y_frame(i),o_frame(i),-2);
        %    line([x1r x2]',[y1r y2]','Color','r');
        %end

        if(showIds) 
            id_frame = id(frame == j);
            for k = 1:length(x_frame)
               y_cur = y_frame(k);
               x_cur = x_frame(k);
               id_cur = id_frame(k);

               text(x_cur+4,y_cur+4,num2str(id_cur),'FontSize',7,'Color',[1 1 1]);
            end;
        end
        
        drawnow;
        F(count) = getframe(hf);
        %cdata = print(hf,'-RGBImage','-r100');
        %imwrite(cdata,fname_out,'Compression','none','WriteMode','Append','Resolution',100);
        clf(hf)
        count = count + 1;
    end
    p.done();
end

function [x1, y1] = rotate(x,y,theta,length)
    theta = pi*theta/180;
    R = [ cos(theta)   sin(theta);
         -sin(theta)   cos(theta)];
    xy = R * [length; 0];
    x1 = xy(1) + x;
    y1 = xy(2) + y;
end