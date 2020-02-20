function [posList] = detectCells(fname,varargin)
% [posList] = detectCells(fname)
%
% creates of a list of cell positions in each frame by finding
% local maxima in the images, see private/pkfnd.m for more info on how cells
% are detected.
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
% posList  with structure:
%             (x)          (y)          (t)
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
    addOptional(p,'check',[],@(x) isnumeric(x) && (0 < x) && (x <= length(info)));
    addOptional(p,'MinLevel',1.3,@(x) isnumeric(x) && x > 0);
    parse(p,varargin{:});

    import trackingUtils.*;

    %s = warning('off','MATLAB:colon:nonIntegerIndex');
    warning off;

    fprintf('Analyzing Frame:   ');

    %preallocating
    A = zeros(info(1).Height,info(1).Width);
    b = zeros(info(1).Height,info(1).Width);
    PADDING = 1200;

    if(isempty(p.Results.check))
    	frames = 1:length(info);
    else
        frames = p.Results.check;
    end

    posListTemp = NaN(PADDING * max(frames),3);  %Position of cells [x y frame]
    start = 1;
    stop = 0;

    for k = frames
        fprintf([repmat(8,1,4) '%4u'],k);

        %For multipage tiffs
        A = imread(fname,k);

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
        b = bpass(A,1,5); %Removes FB noise better

        %colormap('gray'),imagesc(b)

        %Find cells
        %pk = pkfnd(b,3,15);
        pk = pkfnd(b,p.Results.MinLevel,15);

    %     imagesc(b);
    %     hold on;
    %     scatter(pk(:,1),pk(:,2));
    %     colormap gray;

        %pk = pkfnd(b,3,15);
        %cnt = pk;
        cnt = cntrd(b,pk,18);

        start = stop + 1;
        stop = start + size(cnt,1) - 1;
        spots = start:stop;

        if stop > size(posListTemp,1)
            posListTemp = [posListTemp; NaN((max(frames) - k) * PADDING,3)];
        end

        posListTemp(spots,1) = cnt(1:end,1); %X pos
        posListTemp(spots,2) = cnt(1:end,2); %Y pos
        posListTemp(spots,3) = ones(length(cnt),1) * (k); %Frame

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

    posList = posListTemp(1:stop,:);

    if(~isempty(p.Results.check))
    	imagesc(b);
    	hold on;
    	scatter(pk(:,1),pk(:,2));
    	colormap gray;
    end

    warning on;

end
