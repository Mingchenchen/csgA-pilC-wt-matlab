function [K] = floro2dens(K,background)
    if(nargin < 2)
        background = 'WT'
    end
%     I=x;
%     I(I<0)=0;
   % I=ones(size(x));
    if(strcmp(background,'WT'))
        %Fit params for WT in WT, 1.9400e+06 cells
        p1 = 1.482e+08;
        p2 = -2.01e+06;
        p3 = 1.069e+04;
        p4 = 2.851;
       
        for i=1:length(K(1,1,:))
             A=K(:,:,i);
        K (:,:,i)= p1.*A.^3 + p2.*A.^2 + p3.*A + p4;
        
        end
        %I(I < 0) = 0;
%     elseif(strcmp(background,'CSGA'))
%         %Fit params for CsgA in WT, 1.9400e+06 cells
%         p1 = -4.411e+05;
%         p2 = 5307;
%         p3 = 1.115;
%         I = p1.*x.^2 + p2.*x + p3;
% 
%         I(I < 0) = 0;
%     else
%         error('Not a valid background, choose either WT or CSGA')
    end
end

%%
%Fit params for WTinWT 8.2e5 cells
%      Linear model Poly3:
%      fit3(x) = p1*x^3 + p2*x^2 + p3*x + p4
%      Coefficients (with 95% confidence bounds):
%        p1 =   5.855e+07  (5.802e+07, 5.907e+07)
%        p2 =  -7.617e+05  (-7.647e+05, -7.586e+05)
%        p3 =        4517  (4513, 4520)
%        p4 =       1.204  (1.204, 1.205)