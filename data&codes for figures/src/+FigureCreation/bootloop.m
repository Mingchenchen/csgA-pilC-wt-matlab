function [R] = bootloop(data,points,bootfun)
    %Used to perform a loop in a bootfunction
    % 
    % Eg:
    % points = 1:100 %poitns to loop over
    % fun = @(x,i) x * i;
    % bootfun = @(x) bootloop(x,points,fun)
    % bootstat(rand(1000,1),bootfun)
    %
    R = [];
    c = 1;
    
    if(size(points,1) ~= 1)
        error('bootloop: points must be a 1 x N array')
    end
    
    for i = points
        R(c,:) = bootfun(data,i);
        c = c+1;
    end
end