function [I] = floro2dens(x)
%     a2 = 33.371401888256216;
%     a1 = 2.979142424499991;
%     c1 = 0.00293;
%     I = a2*(I./(c1+I))+a1;

%Fit params for 8.2e5 cells
%      Linear model Poly3:
%      fit3(x) = p1*x^3 + p2*x^2 + p3*x + p4
%      Coefficients (with 95% confidence bounds):
%        p1 =   5.855e+07  (5.802e+07, 5.907e+07)
%        p2 =  -7.617e+05  (-7.647e+05, -7.586e+05)
%        p3 =        4517  (4513, 4520)
%        p4 =       1.204  (1.204, 1.205)

%Fit params for 1.9400e+06 cells
    p1 = 1.482e+08;
    p2 = -2.01e+06;
    p3 = 1.069e+04;
    p4 = 2.851;
    I = p1.*x.^3 + p2.*x.^2 + p3.*x + p4;
    %I = 9.28527e3 * I + 0.002641e3;
    I(I < 0) = 0;
end

%%
%G =  
%            sse: 5.4249e+06
%        rsquare: 0.8581
%            dfe: 3145724
%     adjrsquare: 0.8581
%           rmse: 1.3132
%
% fit3 = 
% 
%      Linear model Poly3:
%      fit3(x) = p1*x^3 + p2*x^2 + p3*x + p4
%      Coefficients (with 95% confidence bounds):
%        p1 =   1.482e+08  (1.47e+08, 1.495e+08)
%        p2 =   -2.01e+06  (-2.018e+06, -2.003e+06)
%        p3 =   1.069e+04  (1.068e+04, 1.07e+04)
%        p4 =       2.851  (2.85, 2.853)