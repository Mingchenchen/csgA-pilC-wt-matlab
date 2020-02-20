function [K] = normalizeImages(img_fname,save_folder)
    % Loads the raw TIFF files from the microscope software, reduces the size of,
    % performs filters, tries to adjust for uneven illumination, and subtracts the global mean
    % from the images to produce normalized readings for calulcating cell density.
    % saves the resulting datafile to the save_folder
    %
    % K: Filtered and normalized images
    Nbins = 2^10;

    calibrate = createCalibrationImage(img_fname);
    %save('calibrate1','calibrate');
    %load('calibrate1');
    K = normalizeImages_(img_fname,calibrate,Nbins);
    Kmin = min(K(:));
    Kmax = max(K(:));
    save([save_folder 'Knormalized.mat'],'K','-v7.3');
end
