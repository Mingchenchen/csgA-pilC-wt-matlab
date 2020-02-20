classdef ClosedLoop3<Engine3
    %The Closed-Loop model esitamtes the biofilm cell desntiy from the
	%  agent positions. Aggregate locations are then extracted from
	%  the estimated cell density. Cell run behaviors are then chosen
	%  based on the estimated cell density and dected aggregate positions
    
    properties
        start_density;
    end
    
    methods        
        function obj = ClosedLoop3(probs,nCells,probspilin,probspilout,nCellspil,varargin)
	        %
	        %Construct a model using DensityField.m
	        %
	        %params
			% probs         Probabilites() class that drives cell beahvior
			% nCells        number of cells in simulation
			%optional:
	        % start_density A 2D probability matrix of size (Field.xSize,Field.ySize)
			%               used to sample initial cell positions
	        % bandw         KDE bandwith used to calculate cell density
			%
			
            p = inputParser;
            addRequired(p,'probs',@(x) isa(x,'ClosedLoopProbs'));
            addRequired(p,'nCells',@isnumeric);
            addRequired(p,'probspilin',@(x) isa(x,'ClosedLoopProbs2'));
            addRequired(p,'probspilout',@(x) isa(x,'ClosedLoopProbs2'));
            addRequired(p,'nCellspil',@isnumeric);
            addOptional(p,'StartDensity',[],@isnumeric);
            addOptional(p,'bandwidth',15,@isnumeric)
            parse(p,probs,nCells,probspilin,probspilout,nCellspil,varargin{:});  
            
            obj@Engine3( ...
                        DensityField3(p.Results.bandwidth, ...
                                     p.Results.StartDensity, ...
                                     p.Results.nCells,...
                                     p.Results.nCellspil), ...
                        p.Results.probs, ...
                        p.Results.probspilin, ...
                        p.Results.probspilout ...
                       );
        end
        
        function reset(obj)
            obj.stepsElapsed = 0;
            obj.modelStats = ClosedLoopModelStats(); 
            obj.modelStatspil = ClosedLoopModelStats(); 
            obj.field.reset(obj.modelStats,obj.modelStatspil,obj.probs,obj.probspilin,obj.probspilout);

            obj.fieldView.update();
            obj.sanityCheck();
        end
    end
end

