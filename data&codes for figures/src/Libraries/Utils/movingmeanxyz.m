function [u,xi,c,s] = movingmeanxyz(x,z,y,window,varargin)
    % [u,xi, zi] = movingmeanxyz(x,y,window)
    % [u,xi, zi] = movingmeanxyz(x,y,window,xi,zi)
    % u moving mean
    if(length(x) ~= length(y))
        error('X and Y must be same length');
    end
    if(nargin == 3)
        xi = sort(unique(x));
        zi = sort(unique(z));
    else
        xi = varargin{1};
        zi = varargin{2};
    end
    
    Xwindow = window(1);
    Zwindow = window(2);
    
    u = zeros(length(xi),length(zi),1);
    c = zeros(length(xi),length(zi),1);
    s = zeros(length(xi),length(zi),1);
    for i = 1:length(xi)
        Xfilter = x >  (xi(i) - Xwindow) & ...
                 x <= (xi(min(i+1,length(xi))) + Xwindow);
             
        for j = 1:length(zi)
            Zfilter = z >  (zi(j) - Zwindow) & ...
                    z <= (zi(min(j+1,length(zi))) + Zwindow);
                
            u(i,j) = mean(y(Xfilter & Zfilter));
            c(i,j) = sum(Xfilter & Zfilter);
            s(i,j) = std(y(Xfilter & Zfilter));
        end
    end
end