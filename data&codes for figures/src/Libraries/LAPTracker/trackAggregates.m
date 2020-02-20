function [agg_tracks,stable_aggs,unstable_aggs,posList] = trackAggregates(Knormalized,cutoff)

    %% Detect aggreagte positions
    bw = false([739 986 size(Knormalized,3)]);
    for i = 1:size(Knormalized,3)
        IM = imresize(Knormalized(:,:,i),[739 986]);
        A = false([739 986]);
        
        A(IM > cutoff) = true;
        bw(:,:,i) = A;
    end
    
    %% Get information about each aggreagte
    posList.x = [];
    posList.y = [];
    posList.majorAxis = [];
    posList.minorAxis = [];
    posList.frame = [];
    posList.orientation = [];
    posList.frameIndex = [];
    posList.area = [];
    posList.eccentricity = [];
    posList.meanIntensity = [];
    posList.maxIntensity = [];
    for i = 1:size(bw,3)
        A = Knormalized(:,:,i);
%         % Note: 9/2016 (Chris Cotter)
%         % A recalculation of the mean cell density was perfomed, requiring a change of the mean from 2.66 to 1.12.
%         % instead of going back and redoing all the density estimations and simulations, the change was hacked in here
%         IM = imresize(A .* 0.4227,[739 986]);
        IM = imresize(A,[739 986]);

        props = regionprops(bw(:,:,i),IM, ...
            'MajorAxis','MinorAxis','Centroid','BoundingBox', ...
            'Orientation','Area','Eccentricity','MeanIntensity','MaxIntensity');

       xy = reshape([props.Centroid],[2,length([props.MajorAxisLength])])';
       posList.x = [posList.x; xy(:,1)];
       posList.y = [posList.y; xy(:,2)];
       posList.majorAxis = [posList.majorAxis; [props.MajorAxisLength]'];
       posList.minorAxis = [posList.minorAxis; [props.MinorAxisLength]'];
       posList.frame = [posList.frame; repmat(i,length([props.MajorAxisLength]),1)];
       posList.orientation = [posList.orientation; deg2rad([props.Orientation]')];
       posList.frameIndex = [posList.frameIndex; [1:length([props.MajorAxisLength]')]'];
       posList.area = [posList.area; [props.Area]'];
       posList.eccentricity = [posList.eccentricity; [props.Eccentricity]'];
       posList.meanIntensity = [posList.meanIntensity; [props.MeanIntensity]'];
       posList.maxIntensity = [posList.maxIntensity; [props.MaxIntensity]'];
    end
    
    %% Do the tracking
    agg_tracks = aggregateTracking(posList,20);
    
    %% Tack on varaibles that were not used in the tracking
    agg_tracks.orientation = posList.orientation(agg_tracks.index);
    agg_tracks.frameIndex = posList.frameIndex(agg_tracks.index);
    agg_tracks.area = posList.area(agg_tracks.index);
    agg_tracks.eccentricity = posList.eccentricity(agg_tracks.index);
    agg_tracks.meanIntensity = posList.meanIntensity(agg_tracks.index);
    agg_tracks.maxIntensity = posList.maxIntensity(agg_tracks.index);
    
    %% Makes a list of "Stable" aggreagtes, these are aggreagtes
    % that are present in the last from of the video or merged into
    % an aggreagte that is in the last frame (or merged into one of these and so on)
    stable_aggs = agg_tracks.id(agg_tracks.frame == max(agg_tracks.frame));
    ids = unique(agg_tracks.id(agg_tracks.frame ~= max(agg_tracks.frame)))';
    ids = sort(ids,'descend'); %Work backwords from the last aggreagtes to merge
    for k = ids
        tt = subStruct(agg_tracks,agg_tracks.id == k);
        last_frame = tt.frame(end);
        
        aggs_in_frame = subStruct(agg_tracks,agg_tracks.frame == last_frame + 1 ...
                                             & ismember(agg_tracks.id,stable_aggs));   
        x = tt.x(end);
        y = tt.y(end);
        
        for j = 1:length(aggs_in_frame.id)
            if(inEllipse(aggs_in_frame.x(j), ...
                         aggs_in_frame.y(j), ...
                         x,y, ...
                         aggs_in_frame.majorAxis(j)/2, ...
                         aggs_in_frame.minorAxis(j)/2, ...
                         -aggs_in_frame.orientation(j)))
                      
                     stable_aggs = [stable_aggs; k];
                     break
            end
        end
    end
    
    %% Makes a list of unstable aggreagtes
    % To be considered an aggreagte at all, the blob must be at least
    % 1000px in area and exist for 10 frames
    all_unstable_aggs = setdiff(agg_tracks.id,stable_aggs);
    unstable_aggs = [];
    for j = all_unstable_aggs'
        a = subStruct(agg_tracks,agg_tracks.id == j);

        if(max(a.area) > 1000 & ...
                a.frame(end) - a.frame(1) > 10)
            unstable_aggs = [unstable_aggs; j];
        end
    end
    
    agg_tracks = subStruct(agg_tracks,ismember(agg_tracks.id,union(stable_aggs,unstable_aggs)),'r');
    agg_tracks.stable = ismember(agg_tracks.id,stable_aggs);
end

function [inside] = inEllipse(x,y,xp,yp,a,b,alpha)
    %xp,yp are point coordinates
    %x,y is teh center of the ellipse
    %alpha is positive in the COUNTERCLOCKWISE direction
    %   relative to the x axis. 
    
    inside = (((cos(alpha).*(xp - x) + sin(alpha).*(yp - y)).^2 / a.^2 + ...
               (sin(alpha).*(xp - x) - cos(alpha).*(yp - y)).^2 / b.^2) <= 1);
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
% hf = figure('Visible','on'), hold on;
% frames = agg_tracks.frame;
% 
% Kmin = min(Knormalized(:));
% Kmax = max(Knormalized(:));
% 
% xc = cos(0:0.1:2*pi);
% yc = sin(0:0.1:2*pi);
% %s = {'/Users/cotter/Desktop/LS3629 0_04 LS3909 07172014 7hr 1-250.tif', ...
% %    '/Users/cotter/Desktop/LS3629 0_04 LS3909 07172014 7hr 251-500.tif', ...
% %   };
% for i = 1:max(frames)
%     i
%     clf
%     xi = agg_tracks.x(frames == i);
%     yi = agg_tracks.y(frames == i);
%     id = agg_tracks.id(frames == i);
%     %fi = tracks(frames == i,8);
%     idxi = agg_tracks.index(frames == i);
%     maja = agg_tracks.majorAxis(frames == i);
%     mina = agg_tracks.minorAxis(frames == i);
%     ori = agg_tracks.orientation(frames == i);
%     id;
%     
%    % 
%     %A = imread('/Volumes/USB/06192015_1/06192015_1_MMStack.ome.tif',i);
%     %imagesc(imresize(A,[739 986]));
%     %imagesc(imresize(K(:,:,i),size(bw(:,:,1))));
%     imagesc(imresize(Knormalized(:,:,i),[740 986]));
%     %imagesc(bw(:,:,i));
%     hold on;
%     
% 
%     
%     for j = find(ismember(id,stable_aggs))'
%          text(xi(j) + 10,yi(j) - 30,num2str(id(j)),'Color','r');
%          plot(xi(j),yi(j),'+');
% %         rectangle('Position',bb(j,:),'EdgeColor',color_chooser(id(j)),'LineWidth',2);
% %         
%         %plot(B{idxi(j)}(:,2),B{idxi(j)}(:,1),'Color',color_chooser(id(j)),'LineWidth',2)
% 
%       
%         
%         ph = ori(j);
%         a = maja(j)/2;
%         b = mina(j)/2;
% 
%         R = [ cos(ph)   sin(ph)
%              -sin(ph)   cos(ph)];
%         xy = R*[xc .* a; yc .* b];
%         xrotated = xy(1,:) + xi(j);
%         yrotated = xy(2,:) + yi(j);
%         plot(xrotated,yrotated,'Color',color_chooser(id(j)),'LineWidth',2);
%         
%     end
% 
%     caxis([Kmin Kmax])
%     colormap jet
%     title(['Frame: ',num2str(i)]);
%     xlim([0 986]);
%     ylim([0 739]);
%     drawnow;
%     %cdata = print(hf,'-RGBImage','-r300');
%     %imwrite(cdata,'~/Desktop/aggtracking.tiff','Compression','none','WriteMode','Append','Resolution',300);
% end