classdef Engine2<handle
    %Performs the actucal stepping of the model though time and provides
    % helper functions for viewing the model and extracting model data

    properties
        field;             %The field in which the cell exist
        fieldView;         %Draws a GUI view of the model
        modelStats;        %Model stats
        modelStatspil;     %Model stats of pilc
        stepsElapsed;      %number of steps that have occured
        probs;             %Cell behavior model
        probspil;          %Cell behavior model of pilc
        inagg=logical([]);
        dencutoff=3;
    end

    methods
        %
        %Construct a simulation
        %
        function obj = Engine2(field,probs,probspil)
            obj.field = field;
            obj.probs = probs;
            obj.probspil = probspil;
            obj.fieldView = FieldView(obj.field);
            obj.reset(); %Populate the simulation
            fprintf('\n');
        end

        %
        %Reset the simulation
        %
        function reset(obj)
            obj.stepsElapsed = 0;
            obj.modelStats = ModelStats();
            obj.modelStatspil = ModelStats(); 
            obj.field.reset(obj.modelStats,obj.modelStatspil,obj.probs,obj.probspil);

            obj.fieldView.update();
            obj.sanityCheck();
        end

        %
        % Step the model forwrad N steps
        %
        % N the number of steps forward
        function loop(obj,N)
            fprintf('\n\n');

            p = Progress(N);
            for i = 1:N
                p.d(i);
                obj.step();
            end
            p.done();
        end

        %
        % Step the model one step forward in time
        %
        function step(obj)
            obj.stepsElapsed = obj.stepsElapsed + 1;
            obj.field.step(obj.stepsElapsed);

            for i = 1:obj.field.nCells+obj.field.nCellspil
                obj.field.cellList{i}.step(obj.stepsElapsed);
            end
            obj.inagg(:,:,end+1)=obj.field.density>obj.dencutoff;
            obj.fieldView.update();
        end

        %
        % Toggle the view
        %
        % value:
        %     'on'    Turn view on
        %     'off'   Turn view off
        function view(obj,value)
            switch value
                case 'off'
                    obj.fieldView.viewOFF();
                case 'on'
                    obj.fieldView.viewON();
            end
        end

        %
        % Return model data
        %
        % Without "Flusing" simTracks will return the tracks for each agent
        % up to its last state change. "Flushing" will "flush" the data after
        % the last state change into the simTracks, but the trajectories will
        % be corrupted if the model is steped forward in time after a flush
        %
        % simTracks =
        %       exportTracks()
        %           Exports a struct simular to m_tracks
        %       exportTracks(noflush)
        %           Passing in any value skips flushing data from cells
        %
        function [simTracks,simTracks2] = exportTracks(obj)
            %Force remaning track data into modelStats
            if(nargin == 1)
                for i = 1:obj.field.nCells+obj.field.nCellspil;
                    obj.field.cellList{i}.flush();
                end
            end
            simTracks2=obj.modelStatspil.tracks;
            simTracks = obj.modelStats.tracks;
        end

        %
        % Return model data
        %
        % See exportTracks for details about noflush. In the case of runs
        % a flushed agent adds a partial run to the run database
        % that ends at the agents current position
        %
        % modelData =
        %       exportData()
        %           Exports a table of run variables simular to AllDataTable
        %       exportData(noflush)
        %           Passing in anything skips flushing data from cells
        %
        function [modelData,modelData2] = exportData(obj,varargin)
            if(nargin == 1)
                for i = 1:obj.field.nCells+obj.field.nCellspil 
                    obj.field.cellList{i}.flush();
                end
            end
            modelData2=obj.modelStatspil.exportData();
            modelData = obj.modelStats.exportData();
        end

        %
        %Destructor
        %
        function delete(obj)
            obj.fieldView.viewOFF();
        end

        %
        % Checks run after model reset or initlization
        %
        function sanityCheck(obj)
            % Check on a few things and throw warnings
            if(obj.stepsElapsed ~= 0)
                warning('Steps Elapsed not starting at 0')
            end
            if(obj.field.currentStep < 1 || obj.field.currentStep > size(obj.field.density,3))
                Error('Field currentStep is starting outside the range 1 , size(density,3)')
            end
        end
    end

end
