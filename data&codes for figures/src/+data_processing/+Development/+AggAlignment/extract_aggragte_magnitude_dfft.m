function [M,S] = extract_aggragte_magnitude_dfft(Kum,low_cutoff,high_cutoff)
% Extracts the sum of the magnitudes of the spatial wavelengths for each
% image in Kum between low_cutoff and high_cutoff using dfft. 
% High and low cutoffs are only  approximate as the dft bin nearest to the  
% cutoff value is chosen.
%
% :param Kum: image sequence to extract magnidues from
% :param low_cutoff: lowest wavelength (in micrometers) to include in magnitude
% :param high_cutoff: highest wavelength (in micrometers) to include in magnitude
%
% :returns M: Where M(i) is the sum of magnitudes of image Kum(:,:,i)

    %run('ENVS.m') %Get FOV sizes and distance converstion constants
    FOV_X = 1440; % X size of microscope field of view (in pixels)
FOV_Y = 1920; % Y size of microscope filed of view (in pixels)
FOV_MU_CONV = 0.5136; % micrometers per pixel

%Images are shrunk during normalization to save space and reduce
% computation time
Kum_X = 2^10; % X size of microscope images after normalization
Kum_Y = 2^10; % Y size of microscope images after normalization
Kum_MU_CONV = [FOV_X/Kum_X * FOV_MU_CONV, ...
               FOV_Y/Kum_Y * FOV_MU_CONV]; %witdh of each pixel in Kum in 
                                           %  micrometers
    Nframes = size(Kum,3);
    Nx = size(Kum,1);
    Ny = size(Kum,2);
    
    low_cutoff_bins  = ceil( Nframes ./ (low_cutoff .* Kum_MU_CONV) );
    high_cutoff_bins = floor( Nframes ./ (high_cutoff .* Kum_MU_CONV) );
     
    % window chosen based on http://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf
    [hx,hy] = meshgrid(hanning(size(Nx,1)),hanning(size(Ny,2)));
    w = hx.*hy;

    M = zeros(Nframes,1);
    S = zeros(Nframes,1);
    for i = 1:Nframes
        cz = Kum(:,:,i);
        cz = cz .* w;
        
        fcz = fft2(cz);
        A = abs(fcz);
        
        B = A(1:(Ny/2 + 1),1:(Ny/2 + 1));
        B(:,2:end) = B(:,2:end) + fliplr(A(1:(Ny/2 + 1),(Ny/2 + 1):end));
        
        values = B(high_cutoff_bins(1):low_cutoff_bins(1),high_cutoff_bins(2):low_cutoff_bins(1));
        M(i) = mean(values(:));
        S(i) = std(values(:));
    end
end