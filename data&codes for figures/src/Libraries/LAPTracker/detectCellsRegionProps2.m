function [posList] = detectCellsRegionProps2(fname,varargin)
% [posList] = detectCells(fname)
%
% Creates of a list of cell positions in each frame as contigous sets of pixels
% with intensities above MinLevel. 
%
%
% INPUTS
% fname  filename of multipage tiff of images to detect cells in
%
% OPTIONAL INPUTS
% check     if set, only identifies in the frame from the tiff stack identifed
%               by check (max(frames) <= check <= min(frames)). The indentifed
%               cells are plotted and shown. Use this function to confirm
%               the values of pkfnd
%                 [posList] = detectCells(fname,'check',1);
% MinLevel  the minimum brightness of a pixel that might be local maxima.
%               see input th of pkfnd() for more info (Default: 1.3)
%                 [posList] = detectCells(fanem,'MinLevel',2);
%
%
% OUTPUTS
% posList  a struct with values:
%         (posList.x)  (posList.y)  (posList.t)
%           3.60000      5.00000      0.00000
%           15.1000      22.6000      0.00000
%           4.10000      5.50000      1.00000
%           15.9000      20.7000      2.00000
%           6.20000      4.30000      2.00000
%
% Author: Chris Cotter (cotter@uga.edu)

    %%
    % TODO:
    % 1) Work find a fix for the warning off hack
    %%
    info = imfinfo(fname); %Info about the tiff file

    p = inputParser;
    addOptional(p,'check',[],@(x) all(isnumeric(x)) && (0 < min(x)) && (max(x) <= length(info)));
    addOptional(p,'EccCutoff',0);
    addOptional(p,'MinLevel',1.3,@(x) isnumeric(x) && x > 0);
    addOptional(p,'large_cutoff',10,@(x) isnumeric(x) && x > 0);
    addOptional(p,'small_cutoff',1,@(x) isnumeric(x) && x > 0);
    addOptional(p,'startframe',0);
    parse(p,varargin{:});

    ECC_CUTOFF = p.Results.EccCutoff;

    import trackingUtils.*;

    %s = warning('off','MATLAB:colon:nonIntegerIndex');
    warning off;

    fprintf('Analyzing Frame:   ');


    %preallocating
    A = zeros(info(1).Height,info(1).Width);
    b = zeros(info(1).Height,info(1).Width);
    PADDING = 1200;

    posList.x = NaN(PADDING,1);
    posList.y = NaN(PADDING,1);
    posList.o = NaN(PADDING,1);
    posList.frame = NaN(PADDING,1);
    posList.ecc = NaN(PADDING,1);
    posList.inten = NaN(PADDING,1);
    
    if(isempty(p.Results.check))
    	frames = 1:length(info);
    else
        frames = p.Results.check;
    end

    posList.x = [];
    start = 1;
    stop = 0;

    for k = frames
        fprintf([repmat(8,1,4) '%4u'],k);

        %For multipage tiffs
        A = imread(fname,k,'info',info);

        %red = A(:,:,1);
        %green = A(:,:,2);
        %blue = A(:,:,3);

        %Using imadjust actucally adds
        %    noise to the results of bpass
        %A = imadjust(A(:,:,1));
        A = A(:,:,1); %Red image

        %For tracking in phase images
        %A = 255-red;
        %

        %Remove Noise
        %b = bpass(A,1,10);
        b = bpass(A,p.Results.small_cutoff,p.Results.large_cutoff); %Removes FB noise better

        %colormap('gray'),imagesc(b)

        %Find cells
        %pk = pkfnd(b,3,15);
        c = false(size(b));
        c(b > p.Results.MinLevel) = 1;
        %c = bwareaopen(c,8);
        imbw=im2bw(b);
        im_label=bwlabel(imbw);
        prp= regionprops(c,b,'MeanIntensity');
        props = regionprops(c,'Centroid','Orientation','Eccentricity');
        cnt = [reshape([props.Centroid],2,[])]';
        start = stop + 1;
        stop = start + size(cnt,1) - 1;
        spots = start:stop;

        if stop > size(posList.x,1)
            ncells = length(cnt);
            padding = max(ncells * (max(frames) - k),ncells) + ncells;
            posList.x = [posList.x; NaN(padding,1)];
            posList.y = [posList.y; NaN(padding,1)];
            posList.o = [posList.o; NaN(padding,1)];
            posList.frame = [posList.frame; NaN(padding,1)];
            posList.ecc = [posList.ecc; NaN(padding,1)];
            posList.inten = [posList.inten; NaN(padding,1)];
        end

        posList.x(spots) = cnt(:,1);
        posList.y(spots) = cnt(:,2);
        posList.o(spots) = [props.Orientation]';
        posList.ecc(spots) = [props.Eccentricity]';
        posList.inten(spots) = [prp.MeanIntensity]';
        posList.frame(spots) = ones(length(cnt(:,1)),1) * (k+p.Results.startframe); %Frame

    %    hf = figure('visible','off');
    %    set(hf,'Renderer','painters');	% Fast animation, many objects
    %       figure
    %       colormap('gray'),imagesc(b);
    %       hold on;
    %       plot(pk(1:end,1),pk(1:end,2),'o');
    %       drawnow;
    %    print('-dpng',sprintf('img%d.png',k));
    %    clf(hf);
    end



    posList.x(stop+1:end) = [];
    posList.y(stop+1:end) = [];
    posList.o(stop+1:end) = [];
    posList.frame(stop+1:end) = [];
    posList.ecc(stop+1:end) = [];
    posList.inten(stop+1:end) = [];
    posList = subStruct(posList,posList.ecc >= ECC_CUTOFF);

    if(~isempty(p.Results.check))
    	imagesc(b);
    	hold on;
    	plot(posList.x,posList.y,'bo');
    	colormap gray;
        x1 = posList.x;
        y1 = posList.y;
        theta = posList.o;
        %for i = 1:length(x1)
            %[x2, y2] = rotate(x1(i),y1(i),theta(i),5);
            %[x1r, y1r] = rotate(x1(i),y1(i),theta(i),1);
            %line([x1r x2]',[y1r y2]','Color','r');
        %end
    end

    warning on;

end

function [x1, y1] = rotate(x,y,theta,length)
    theta = pi*theta/180;
    R = [ cos(theta)   sin(theta);
         -sin(theta)   cos(theta)];
    xy = R * [length; 0];
    x1 = xy(1) + x;
    y1 = xy(2) + y;
end
