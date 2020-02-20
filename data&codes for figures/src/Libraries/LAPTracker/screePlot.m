function [cfh] = screePlot(fname,range)

info = imfinfo(fname);
Ncells_first = countCells(fname,range,1,info);
Ncells_last  = countCells(fname,range,length(info),info);

cfh = figure;
hold on;
plot(range,Ncells_first,'-o');
plot(range,Ncells_last,'-o');
legend('First Frame','Last Frame')
hold off
end

function [Ncells] = countCells(fname,range,frame,info)
    Ncells = zeros(length(range),1); 
    
    %p = Progress(length(range));
    for i = 1:length(range)
        %p.d(i);
        A = imread(fname,frame,'info',info);

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
        b = bpass(A,3,10);

        %colormap('gray'),imagesc(b)

        %Find cells
        %pk = pkfnd(b,3,15);
        c = false(size(b));
        c(b > range(i)) = 1;
        %c = bwareaopen(c,8);

        props = regionprops(c,'Centroid');
        Ncells(i) = length(props);
    end
    %p.done();
end