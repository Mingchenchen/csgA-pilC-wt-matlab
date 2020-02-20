function [K] = normalizeImages_(fname1,calibrate,BINS)
    info = imfinfo(fname1);
    L = length(info);
    K = zeros(BINS,BINS,L);
    
    p = Progress(L)
    for i = 1:L
        p.d(i)
        %Normalizes for:
        % 1) differences in illumination in different areas. Necessary
        %     because microscope flourscence becomes less intensce
        %     the further away from the center of the image
        %
        %  3)  Loss of flourscence intensity over time most likely due to 
        %       photobleaching
        %  1 is  addressed by dividing by calibrate
        %  3 is addresed by subtracting the mean
        %  The order of these operations, as done below, is important!
        
        A = double(imread(fname1,i,'Info',info)) ./ calibrate;
        B = bpass(A,1,25);
        A(B > prctile(B(:),99)) = NaN;
        A = inpaint_nans(A);
        Afilt = imfilter(double(A),fspecial('gaussian',100,20),'replicate');
        
        Acalibrated = double(Afilt) - mean(Afilt(:));
        
        resized = imresize(Acalibrated,[BINS BINS]);
        K(:,:,i) = resized;
    end
    p.done();
end
