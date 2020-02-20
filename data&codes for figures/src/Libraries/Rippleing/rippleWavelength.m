function [C] = rippleWavelength(K)
    % [C] = rippleWavelength(K)
    % 
    % Performs the 2D DFT on each image t in K ( image = K(:,:,t) )
    % Returns the total magnitude of the rotation invarent frequency.
    %
    % params:
    %  K mxnxt a stack of t image frames of szie mxn
    % returns;
    %  C C(i,t) is the ith DFT mangitude for image t
    %%
    NPoints = size(K,1); %Number of points to generate
    Nframes = size(K,3);
    
    %
    [x,y] = meshgrid(0:NPoints);
    lambda = sqrt(x.^2 + y.^2);

    C = zeros(NPoints/2 + 2,Nframes);
    i = 1;
    %p = Progress(Nframes)
    for r = 1:Nframes
        %p.d(r)
        cz = K(:,:,r);

        % window chosen based on http://download.ni.com/evaluation/pxi/Understanding%20FFTs%20and%20Windowing.pdf
        [hx,hy] = meshgrid(hanning(size(cz,1)),hanning(size(cz,2)));
        w = hx.*hy;
        cz = cz .* w;

        fcz = fft2(cz);
        A = abs(fcz);

        Nx = size(A,1);
        Ny = size(A,2);

        %Combine the top left and bottom right quadrents of the 2d FFT
        % to capture all possible roatations of the pattern.
        %NOTE: the top right and bottom left quadrents are mirror images 
        %   of the top left and bottom right quadrents, respectively
        B = A(1:(Ny/2 + 1),1:(Ny/2 + 1));
        B(:,2:end) = B(:,2:end) + fliplr(A(1:(Ny/2 + 1),(Ny/2 + 1):end));

        %Extract the bin centers for B
        l = lambda(1:(Ny/2 + 1),1:(Nx/2 + 1));
        
        %
        [~, idx] = histc(l(:),l(:,1));
        
        % accumulates elements of the vector B using the subscripts in idx.
        binsums = accumarray(idx + 1,B(:));

        C(:,i) = binsums;
        i = i + 1;
    end
    %p.done()
end
