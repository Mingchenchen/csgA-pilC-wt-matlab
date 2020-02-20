function [ ] = multiColorLine(x,y,c,cmap)
%MULTICOLORLINE Summary of this function goes here
%   Detailed explanation goes here
    
    if(nargin < 4)
        cmap = colormap('jet');
    else
        cmap = cmap;
    end
   
    if(size(x,1) > 1)
        x = x';
    end
    
    if(size(y,1) > 1)
        y = y';
    end
    
    if(size(c,1) > 1)
        c = c';
    end
   
    z = zeros(1,length(x));
    
%     c = floor((c - min(c)) / (max(c - min(c))) * (size(cmap,1) - 1) + 1)
%     c = (c - 1) / 2 * size(cmap,1) + 1
    
    %surface([x;x],[y;y],[z;z],[c;c],'edgecol','flat','linew',2,'CDataMapping','direct')
    surface([x;x],[y;y],[z;z],[c;c],'edgecol','flat','linew',1)
    colormap(cmap);
end

