function calibrate = createCalibrationImage(fname1)
    info = imfinfo(fname1);
    L = length(info);
    p = Progress(15)
    for i = 1:15
        p.d(i)
        A = double(imread(fname1,i,'Info',info));
        calibrate_array(:,:,i) = A;
    end
    p.done()

    calibrate = imfilter(mean(double(calibrate_array(:,:,:)),3),fspecial('average', 500),'replicate');
    calibrate = calibrate .* mean(calibrate(:));
end