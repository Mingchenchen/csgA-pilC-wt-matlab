classdef ClosedLoopProbsmx2<Probabilities
	%
	% Extends the behvaor logic in the Probabilites class to 
	% include the cell-cell orientation alignemnt included in 
	% the closed loop models
	%
    
    properties
        Init;
        InitFrames; 
        Mdlb2;
        Mdlr2;
        Mdlst2;
	rt_tow_m;
	rt_tow_s;
	rt_away_m;
	rt_away_s;
	mu_rt_tow_m;
        mu_rt_tow_s;
        mu_rt_away_m;
        mu_rt_away_s;
	rs_tow_m;
          rs_tow_s;
          rs_away_m;
          rs_away_s;
          mu_rs_tow_m;
          mu_rs_tow_s;
          mu_rs_away_m;
          mu_rs_away_s;
    end
    
    methods
        function obj = ClosedLoopProbsmx2( Mdlb2,Mdlr2,Mdlst2,Mdlb, Mdlr, Mdlst,Init2,InitFrames)
	        %
	        % Constructor
	        %
	        %params:
	        %   Mdlb  KNNSampler for sampling next run paramaters
	        %   Mdlr  KNNSampler for sampling reversal angle
	        %   Mdlst KNNSampler for sampling next state
			%   Init  Struct containing:
			%          Init.Mdlr   KNNSampler for sampling next reversal angle during the prerun
			%          Init.Mdlrst KNNSampler for sampling next run state during the prerun
			%          Init.Mdlrb KNNSampler for sampling next run behaviors during the prerun
			%   InitFrames number of frames in the prerun
	        %
            
            obj@Probabilities(Mdlb,Mdlr,Mdlst);
            obj.Init = Init2;
            obj.InitFrames = InitFrames;
            obj.Mdlb2=Mdlb2;
            obj.Mdlr2=Mdlr2;
            obj.Mdlst2=Mdlst2;
            load('mu_r.mat');
            load('wt_r.mat');
            obj.rt_tow_m=rt_tow_m;
            obj.rt_tow_s=rt_tow_s;
            obj.rt_away_m=rt_away_m;
            obj.rt_away_s=rt_away_s;
            obj.mu_rt_tow_m=mu_rt_tow_m;      
            obj.mu_rt_tow_s=mu_rt_tow_s;       
            obj.mu_rt_away_m=mu_rt_away_m;
            obj.mu_rt_away_s=mu_rt_away_s;
            obj.rs_tow_m=rs_tow_m;
              obj.rs_tow_s=rs_tow_s;
              obj.rs_away_m=rs_away_m;
              obj.rs_away_s=rs_away_s;
              obj.mu_rs_tow_m=mu_rs_tow_m;
              obj.mu_rs_tow_s=mu_rs_tow_s;
              obj.mu_rs_away_m=mu_rs_away_m;
              obj.mu_rs_away_s=mu_rs_away_s;
        end
        
        function [time,spd,theta1,beta1,dist,state1] = getTransitions(obj,TSS1,phi0,chi,Dn1,rho1,state0)   
	        %
	        % Get the transition probilities for a cell
	        %
	        % Variables:
	        %   TSS1    time since beginning of simulation
	        %   phi0    angle to the nearest aggreagte
	        %   Dn1     distnace to the nearest aggreagte
	        %   rho1    local cell density
	        %   state0  current cell state
			%   chi
	        %
	        % returns
	        %   time    duration of next run
	        %   spd     speed of next run
	        %   theta1  reversal angle prior to next run
	        %   beta1   angle between next run and nearest aggreagte
	        %   dist    distance between chosen run from run database and the
	        %            current cell state for (1) Mdlst, (2) Mdlr, and (3) Mdlb
	        %   state1  state of next run
	        %
			
            if(TSS1 < obj.InitFrames)
				% Use different behavior logic during the prerun
                [time,spd,theta1,beta1,dist,state1] = obj.getTransitions__(phi0,Dn1,chi,rho1,state0);
                return
            else
                TSS1 = TSS1 - obj.InitFrames - 1;
            end
            
            dist = zeros(3,1);
            %Get new state
            if(state0 < 3)
                s0 = 1;
            else
                s0 = 2;
            end

			[n, dist(1)] = obj.Mdlst2(s0).KNN.sample([TSS1, ...
			                                   NaN, ...rho1
			                                   cos(phi0), ...
			                                   Dn1]);
											   
            state1 = obj.Mdlst2(s0).state1(n);

            if (state0 < 3 && state1 < 3)
                s = 1;
            elseif(state0 < 3 && state1 == 3)
                s = 2;
            elseif(state0 == 3 && state1 < 3)
                s = 3;
            else
                error('Bad Transition');
            end
            
            %Get new reversal angle
            [neighbor, dist(2)] = obj.Mdlr(s).KNN.sample([chi,...
                TSS1,...
                Dn1,...
                phi0]);
			%neighbor = randi(length(obj.Mdlr(s).theta1),1,1); % Uncomment to make reversal angle random
            theta1 = obj.Mdlr(s).theta1(neighbor);
            beta1 = atan2(sin(theta1 - phi0),cos(theta1 - phi0));
            
            %Get new beahvior
            if state1==3
            [neighbors, dist(3)] = obj.Mdlb(s).KNN.sample([TSS1, ...
                                               NaN, ...rho1
                                               cos(beta1), ...
                                          Dn1]);
            time = obj.Mdlb(s).T.Rt1(neighbors) ;
            spd =  obj.Mdlb(s).T.Rs1(neighbors) ; 
            dis=ceil(Dn1+50);
            if dis<1 
		dis=1; end;
            if dis>201||isnan(dis) 
		dis=201; end;
                if cos(beta1)>0%toward
                    time=ceil(time*obj.rt_tow_m(dis)/obj.mu_rt_tow_m(dis));
                    %time=(time-obj.mu_rt_tow_m(dis))/obj.mu_rt_tow_s(dis)*obj.rt_tow_s(dis)+obj.rt_tow_m(dis);
                    spd=(spd*obj.rs_tow_m(dis))/obj.mu_rs_tow_m(dis);
                    %spd=(spd-obj.mu_rs_tow_m(dis))/obj.mu_rs_tow_s(dis)*obj.rs_tow_s(dis)+obj.rs_tow_m(dis);
                else %away
                    time=ceil(time/obj.rt_away_m(dis)*obj.mu_rt_away_m(dis));
                    %time=(time-obj.mu_rt_away_m(dis))/obj.mu_rt_away_s(dis)*obj.rt_away_s(dis)+obj.rt_away_m(dis);
                    spd=(spd*obj.rs_away_m(dis))/obj.mu_rs_away_m(dis);
                    %spd=(spd-obj.mu_rs_away_m(dis))/obj.mu_rs_away_s(dis)*obj.rs_away_s(dis)+obj.rs_away_m(dis);
                end
                time=time*2;
                spd=spd/2;
            else 
                [neighbors, dist(3)] = obj.Mdlb(s).KNN.sample([TSS1, ...
                                               NaN, ...rho1
                                               cos(beta1), ...
                                              Dn1]);
				[neighbors2, dist(3)] = obj.Mdlb2(s).KNN.sample([TSS1, ...
                                               NaN, ...rho1
                                               cos(beta1), ...
                                              Dn1]);							  
                time = obj.Mdlb(s).T.Rt1(neighbors) * 2;
                spd =  obj.Mdlb(s).T.Rs1(neighbors) / 2;                         
             end

            
        end
        
        function [time,spd,theta1,beta1,dist,state1] = getTransitions__(obj,phi0,Dn1,chi,rho1,state0)   
			%
			% Seperate reversal logic used during the prerun. Uses reveral and beahvoir
			% databases contained within obj.Init
			%
			
            dist = zeros(3,1);

            %Get new state
            if(state0 < 3)
                s0 = 1;
            else
                s0 = 2;
            end
            
            [n, dist(1)] = obj.Init.Mdlst(s0).KNN.sample([rho1]);
            state1 = obj.Init.Mdlst(s0).state1(n);

            if (state0 < 3 && state1 < 3)
                s = 1;
            elseif(state0 < 3 && state1 == 3)
                s = 2;
            elseif(state0 == 3 && state1 < 3)
                s = 3;
            else
                error('Bad Transition');
            end
            
            %Get new reversal angle
            [neighbor, dist(2)] = obj.Init.Mdlr(s).KNN.sample([chi]);
            %neighbor = randi(length(obj.Init.Mdlr(s).theta1),1,1); % Uncomment to make reversal angle random
            theta1 = obj.Init.Mdlr(s).theta1(neighbor);
            beta1 = atan2(sin(theta1 - phi0),cos(theta1 - phi0));
            
            %Get new beahvior
            %[neighbors, dist(3)] = obj.Init.Mdlb(s).KNN.sample([NaN]);
            neighbors = randi(length(obj.Init.Mdlb(s).T.Rt1),1,1); % Uncomment to make behavior random
            time = obj.Init.Mdlb(s).T.Rt1(neighbors) * 2;
            spd =  obj.Init.Mdlb(s).T.Rs1(neighbors) / 2;
        end
    end
    
end

