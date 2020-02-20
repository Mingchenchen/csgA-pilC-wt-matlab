function [weighted_dfft_bin,window_centers] = ripplePeriod(K,WINDOW_SIZE,STEP)
    %%
    assert(WINDOW_SIZE > 0 & mod(WINDOW_SIZE,2) == 0,'Window size must be > 0 and even')
    
    Nframes = size(K,3);
    half_window = WINDOW_SIZE / 2;
    window_centers = half_window:STEP:(size(K,3) - half_window);
    Nwindows = length( window_centers);

    weighted_dfft_bin = zeros(Nwindows,1);
    i = 1;
    p = Progress(Nwindows);
    for r = 1:Nwindows
        p.d(r)
        cz = K(:,:,(STEP * (r - 1) + 1):(WINDOW_SIZE + (STEP * (r-1))));
        
        % window chosen based on http://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf
        w = hanning(size(cz,3));
        cz = bsxfun(@times,cz,reshape(w,1,1,length(w)));

        cz = fft(cz,[],3);
        cz = abs(cz);

        Nt = size(cz,3);
 
        cz = cz(:,:,2:(Nt/2 + 1));
        
        cz = reshape(cz,[],Nt/2);
        
        magnitude = sum(cz,1);
        weighted_dfft_bin(r) = sum( (1:size(cz,2)) .* (magnitude / sum(magnitude(:))));
    end
end
