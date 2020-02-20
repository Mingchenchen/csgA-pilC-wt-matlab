function [x1, y1] = rotatePoint(x,y,theta,length)
    theta = pi*theta/180;
    R = [ cos(theta)   sin(theta);
         -sin(theta)   cos(theta)];
    xy = R * [length; 0];
    x1 = xy(1) + x;
    y1 = xy(2) + y;
end