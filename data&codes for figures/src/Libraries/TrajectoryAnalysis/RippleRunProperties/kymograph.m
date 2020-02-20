function [kymo,xmin,xpoints,YTick,YTickLabel,frame_px_ratio,mutiplier] = kymograph(x,t,I)
    SLICE_HEIGHT = 15;

    p = Progress(length(t));
    xmin = max(ceil(min(x)) - 50,1);
    xmax = min(ceil(max(x)) + 50,987);
    kymo = zeros(SLICE_HEIGHT * 2 * length(t),xmax-xmin+1);
    for i = 1:size(I,3)
        p.d(i);
        A = I(:,:,i); %imread(fname1,t.frame(i),'info',info);
        A = imresize(A,[1440 * 0.5136 1920 * 0.5136]);

        ymin = max(ceil(t(i)) - SLICE_HEIGHT,1);
        ymax = min(ceil(t(i)) + SLICE_HEIGHT,size(A,1));

        %imagesc(A(ymin:ymax,xmin:xmax))
        kymo((((i-1) * (ymax-ymin+1)) + 1):(i * (ymax-ymin+1)),:) = A(ymin:ymax,xmin:xmax);
        %drawnow
    end
    p.done()
    
    mutiplier = (SLICE_HEIGHT + 1) * 2;
    xpoints = ([1:length(t)]'-1) .* (SLICE_HEIGHT + 1) * 2;
    frame_px_ratio = size(kymo,1) / (t(end) - t(1));
    YTick = (1:10:length(t)) .* frame_px_ratio;
    YTickLabel = (1:10:length(t));
    
    figure, hold on;
        imagesc(kymo)

        ax = gca;
        ax.YTick = YTick;
        ax.YTickLabel = YTickLabel;
        ax.YLim = [0,size(kymo,1)];
        colormap gray

        xlabel('Distance (\mum)');
        ylabel('Time (frames)');
end