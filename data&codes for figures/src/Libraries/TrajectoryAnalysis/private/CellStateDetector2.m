function [P,S] = CellStateDetector2(x,y,o)
%[P] = CellStateDetector(x,y,o)
%
%INPUTS
% x(i) x position in micrometers at step i
% y(i) y position in micrometers at step i
% o(i) orientation (in degrees) at step i
%       Note, the orientation output of ExtractTracks has the oppisate
%       sign of the input needed here. ie. muiply o by -1 before inputting
%
%OUTPUTS
% P(i,m) probability of model m and step i + 1
% S(i,:) cell state [x y v o m] at step i + 1

    delta_t = 30; %(Seconds) This variable is used in nested function!
    S = zeros(length(x)-1,5);
    P = zeros(length(x)-1,3);

    %% Matrix Setup
    %%%
    %Process Noise
%     Q = [ 0.7230   -0.0191    0.0001   -0.0404
%          -0.0191    0.7175    0.0004   -0.1447
%           0.0001    0.0004    0.0020   -0.0039
%          -0.0404   -0.1447   -0.0039  277.6630];
    Qfr = [ 0.8707    0.0392    0.0002   -0.0070
            0.0392    0.9231   -0.0005    0.0195
            0.0002   -0.0005    0.0014   -0.0019
           -0.0070    0.0195   -0.0019  357.3595];


    Qs =[ 0.4327    0.0214   -0.0003   -0.0094
          0.0214    0.4690    0.0003   -0.0131
         -0.0003    0.0003    0.0027    0.0024
         -0.0094   -0.0131    0.0024  338.0039];

%     Qfr   = [ 0.8707    0.0392    0.0060   -0.0070
%             0.0392    0.9231   -0.0141    0.0195
%             0.0060   -0.0141    1.2476   -0.0560
%            -0.0070    0.0195   -0.0560  357.3595];
%        
%     Qs = Qfr; 

    %p(x) covairance matrix
    P_f = Qfr;
    P_r = P_f;
    P_s = Qs;
    
    %State conversion matrix
    H =[1 0 0 0;
        0 1 0 0
        0 0 0 1];

    %Measurment Noise
    R = 0;

    %Model tranisition probabilities
    %T = [0.773 0.067 0.16; 0.067 0.773 0.16; 0.3 0.3 0.4];
%      T = [ 0.9029    0.0467    0.0504
%            0.0467    0.9029    0.0504
%            0.1017    0.1017    0.7966];
%     T =[ 0.8847    0.0806    0.0348
%          0.0789    0.8881    0.0330
%          0.0420    0.0235    0.9345];
    T =[ 0.9194    0.0458    0.0174 0.0174
         0.0458    0.9194    0.0174 0.0174
         0.0655    0.0655    0.9345      0
         0.0655    0.0655    0      0.9345];
      
    assert(1 - sum(T(:))/size(T,1) < 1e10)

    %% State initiation
    %%%
    v_init = mean(sqrt((x(1) - x(2)).^2 + (y(1) - y(2)).^2) / delta_t); %um/sec
    u = [x(1), y(1), v_init, mod(o(1),360)]';

    %Model Probabilities
    model_prob = [1 1 1]'./3;
    active_model = 1;
    p = zeros(3,3);
    p(:,1) = [1,1,1];

    %% Hypothesis Testing
    %%%
    for i = 2:length(x)-1

        %% Update measurments
        %%%
        u(:,1) = [x(i-1)   y(i-1)   abs(u(3)) u(4)]';
        u(:,2) = [x(i) y(i) abs(u(3)) u(4)]';
        z(:,1) = [x(i) y(i) o(i)]';
        z(:,2) = [x(i+1) y(i+1) o(i+1)]';
        
        %% RUN EKFs
        %%%
        [p(:,2),u,P_f,P_r,P_s,~,~] = run_ekfs(u(:,1),z(:,1),P_f,P_r,P_s,Qfr,Qs,H,R); %T estimate
        [p(:,3),~,~,~,~,~,~]       = run_ekfs(u(:,2),z(:,2),P_f,P_r,P_s,Qfr,Qs,H,R); %T + 1 estimate

        %% Cacluate Model probabilidies
        %%%
        prior = T * model_prob;

        p_f = zeros(3,1);
        for h = 1:3
            for j = 1:3
                for k = 1:3
                    p_f(h) = p_f(h) + prior(h) * p(j,1) * T(j,h) * p(h,2) * T(h,k) * p(k,3);
                end
            end
        end
        p(:,1) = p(:,2);
        c = sum(p_f);

        model_prob(1) = p_f(1) ./ c;
        model_prob(2) = p_f(2) ./ c;
        model_prob(3) = p_f(3) ./ c;

        active_model_prob(i,:) = model_prob';

        %Determine Active Model based on model probabilities
        [~, active_model] = max(model_prob);
        switch(active_model)
            case 1
                u = u(:,1);
            case 2
                u = u(:,2);
            case 3
                u = u(:,3);
        end

        %% Save State
        %%%
        assert(1 - abs(cosd(u(4)-o(i))) < .00001,num2str(1 - cosd(u(4)-o(i))));
        P(i-1,:) = model_prob';
        S(i-1,:) = [u' active_model];
    end
    %% Make estimate for last point
    % (Which does not have a t+1 point to help in the estimate)
    %%%
    u = [x(end-1)   y(end-1)   abs(u(3)) u(4)]';
    z(:,1) = [x(end) y(end) o(end)]';
    [p(:,2),u,P_f,P_r,P_s,~,~] = run_ekfs(u(:,1),z(:,1),P_f,P_r,P_s,Qfr,Qs,H,R); %T estimate
    prior = T * model_prob;
    p_f = zeros(3,1);
    for h = 1:3
        for j = 1:3
            for k = 1:3
                p_f(h) = p_f(h) + prior(h) * p(j,1) * T(j,h) * p(h,2);
            end
        end
    end
    c = sum(p_f);
    model_prob(1) = p_f(1) ./ c;
    model_prob(2) = p_f(2) ./ c;
    model_prob(3) = p_f(3) ./ c;
    active_model_prob(end,:) = model_prob';
    [~, CURRENT_MODEL] = max(model_prob);
    switch(CURRENT_MODEL)
        case 1
            u = u(:,1);
        case 2
            u = u(:,2);
        case 3
            u = u(:,3);
    end
    P(end,:) = model_prob';
    S(end,:) = [u' active_model];
    
    %%
    %Does the actucal ekf stuff
    %%%
    function [p,u_r,P_rf,P_rr,P_rs,z_r,S] = run_ekfs(u,z,P_f,P_r,P_s,Qfr,Qs,H,R)
        
        %% EKF Matricies
        %%%
        v = u(3);
        theta = u(4);
        G_f = [1 0  cosd(theta) * delta_t   -v * sind(theta) * delta_t; 
               0 1  sind(theta) * delta_t    v * cosd(theta) * delta_t;
               0 0              1                         0          ;
               0 0              0                         1          ;];

        G_r = [1 0 -cosd(theta) * delta_t    v * sind(theta) * delta_t; 
               0 1 -sind(theta) * delta_t   -v * cosd(theta) * delta_t;
               0 0              1                         0          ;
               0 0              0                         1          ;];

        G_s = [1 0              0                         0
               0 1              0                         0
               0 0             -1                         0          ;
               0 0              0                         1          ;];

        g_f = [ v * cosd(theta) * delta_t  v * sind(theta) * delta_t 0  0]';
        g_r = [-v * cosd(theta) * delta_t -v * sind(theta) * delta_t 0  0]';
        g_s = [0 0 -v 0]';

        %% RUN EKFs
        %%%
        [u_r(:,1),P_rf,p(1),z_r(:,1),S(:,:,1)] = EKF(u,P_f,z,G_f,H,R,Qfr,g_f);
        [u_r(:,2),P_rr,p(2),z_r(:,2),S(:,:,2)] = EKF(u,P_r,z,G_r,H,R,Qfr,g_r);
        [u_r(:,3),P_rs,p(3),z_r(:,3),S(:,:,3)] = EKF(u,P_s,z,G_s,H,R,Qs,g_s);

        %Negative velocities are impossible in forward and reverse models. 
        if(u_r(3,1) < 0)
            p(1) = 0;
        end
        if(u_r(3,2) < 0)
            p(2) = 0;
        end
    end
end
