function [K,fK,mfK,mask] = fftRippleFilter(Kripple,cutoff,DEBUG)
%Filters 3D matrix Kripple along the 3rd dimension, including only
%frequencies in fft bins [2,cutoff].
%
%Params:
% Kripple: 3D matrix of values to filter
% cutoff: the maximum frequency bin to include
%Optional:
% DEBUG: Add's debugging info
%   DEBUG = 0: Skip all debugging
%   DEBUG = 1: Show figures to check work
%   DEBUG = 2: VERY memory intesive for large matricies!!
%               DEBUG=1 + Do not clear intermideate files 
%               if DEBUG < 2, fK, mfK, and mask are empty, otherwise
%               fK = fft(Kripple,[],3)
%               mfK = fK with filter mask applied
%               mask = mask applied to fK
%Returns:
% K: filtered version of Kripple
% fk,mfk,mask: see DEBUG = 3

    %% Param parsing
    if(nargin < 3)
        DEBUG = 0;
    end

    %% Apply FFT
    fK = fft(Kripple,[],3);

    %% Mask
    N = size(fK,3);
    maskStart = 2;
    maskStop = floor(N/cutoff);
    mask = zeros(1,1,N);
    mask(1,1,maskStart:maskStop) = 1;
    mask(1,1,(N - maskStop + 2):(N - maskStart + 2)) = 1;

    mfK = bsxfun(@times,fK,mask);

    if(DEBUG < 3)
        fK = [];
    end

    % Check Result %
    if(DEBUG > 0)
        %%
        figure
        Nplot = 10;
        f = mfK(1:Nplot,1:Nplot,:);
        plot(abs(reshape(f,[size(f,3),Nplot*Nplot]))')
    end

    %% Apply ifft
    K = ifft(mfK,[],3);

    if(DEBUG < 3)
        mfK = [];
    end

    assert(isreal(K),'K is not real, is the fft mask correct?')

    if(DEBUG > 0)
        %%
        figure
        %for t = 1:size(K,3)
        t = 1;
            subplot(2,3,1)
                imagesc(K(:,:,t))
                title('Filtered')
            subplot(2,3,2)
                imagesc(Kripple(:,:,t))
                title('Origional')
            subplot(2,3,3)
                cla
                hold on
                f1 = K(10,:,t);
                f2 = Kripple(10,:,t);
                plot(f1(:))
                plot(f2(:))
                legend('Filtered','Original')
            subplot(2,3,4)
                imagesc(K(:,:,t) > 0)
                title('Filtered')
            subplot(2,3,5)
                imagesc(Kripple(:,:,t) > 0)
                title('Origional')

            pause(0.2)
            drawnow
        %end
    end
end        