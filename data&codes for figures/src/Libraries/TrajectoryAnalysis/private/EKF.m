function [u,P,p, z_u,S] = EKF(u,P,z,G,H,R,Q,g)
% u variable ans at time t
% P covariance matrix at time t
% z measurment at time t
    
    %%%
    %EKF Prediction
    %%%
    u =   u + g;
          
    P = G * P * G' + Q;
    
    %%%
    %EKF Update
    %%%
    
    S = H * P * H' + R;
    K = P * H' * inv(S);
    z_u = H * u;
    u = u + K * k(z,z_u);
    %P = (eye(size(P)) - K * H) * P;

    P = (eye(size(P)) - K * H) * P;

    p = likelyhood(S,z_u,z);
end
       
function z = k(z,z_u)
      z(1) = z(1) - z_u(1);
      z(2) = z(2) - z_u(2);
      
      y = z_u(3);
      x = z(3);
      theta = [atan2d(sind(x-y),cosd(x-y)),atan2d(sind(x-y-180),cosd(x-y-180))]; [~,I] = min(abs(theta)); 
      z(3) = theta(I);
end

function p = likelyhood(S,u,z)
    p = exp(-0.5 * k(z,u)' * inv(S) * k(z,u)) / sqrt(det(2 * pi * S));
    %p = log((z - u)' * inv(S) * (z - u));
end

