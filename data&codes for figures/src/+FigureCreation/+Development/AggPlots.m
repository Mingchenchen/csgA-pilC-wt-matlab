classdef AggPlots
   methods (Static)
       function[AggProps] = combine_agg_props(AllData)
           % Combines the AllData{set}.AggProps tables from each set
           % into one AggProps table and aligns to AligedStart and AlignedStop
           
           AggProps = table();
           for set = 1:length(AllData)
               props = AllData{set}.AggProps;
               t1 = AllData{set}.AlignedStart;
               t2 = AllData{set}.AlignedStop;
               
               props = props(props.frame > t1 & props.frame <= t2,:);
               props.frame = props.frame - t1;
               props.movie = repmat({AllData{set}.movie},height(props),1);
               
               AggProps = [AggProps; props];
           end
       end
       
        function [] = plot_property(G1,G2,G3,points,groupvar,varargin)
            G1 = G1(ismember(G1.frame,points),:);
            grps = findgroups(G1.frame);
            G1grps = splitapply(@(x) mat2cell(x,length(x),1),G1(:,groupvar),grps);

            G2 = G2(ismember(G2.frame,points),:);
            grps = findgroups(G2.frame);
            G2grps = splitapply(@(x) mat2cell(x,length(x),1),G2(:,groupvar),grps);

            G3 = G3(ismember(G3.frame,points),:);
            grps = findgroups(G3.frame);
            G3grps = splitapply(@(x) mat2cell(x,length(x),1),G3(:,groupvar),grps);    
            
            combined = {};
            for i = 1:length(points)
                combined{i} = {G1grps{i}, G2grps{i},G3grps{i}};
            end
           
            gboxplot(combined,[],varargin{:})
            h=findobj(gca,'tag','Outliers');
            delete(h)
            ax = gca;
            ax.FontSize = 10;
            ylim([0 ax.YLim(2)]);
            ax.XTickLabel = points / 2 / 60;
            ylabel(groupvar)
            xlabel('Time (hr)')
        end
        
        function [] = agg_count(G1,G2,G3,points,varargin)
            grps = findgroups(G1.movie_id);
            G1Count = cell2mat(splitapply(@(x) mat2cell(histc(x,1:max(G1.frame)),max(G1.frame),1),G1.frame,grps)');
            grps = findgroups(G2.movie_id);
            G2Count = cell2mat(splitapply(@(x) mat2cell(histc(x,1:max(G2.frame)),max(G2.frame),1),G2.frame,grps)');
            grps = findgroups(G3.movie_id);
            G3Count = cell2mat(splitapply(@(x) mat2cell(histc(x,1:max(G3.frame)),max(G3.frame),1),G3.frame,grps)');
            
            combined = {};
            for i = 1:length(points)
                combined{i} = {G1Count(points(i),:), G2Count(points(i),:),G3Count(points(i),:)};
            end
            
            gboxplot(combined,{},varargin{:})
            h=findobj(gca,'tag','Outliers');
            delete(h)
            ax = gca;
            ax.FontSize = 10;
            ylim([0 ax.YLim(2)]);
            ax.XTickLabel = points / 2 / 60;
            ylabel('Aggregate Count')
            xlabel('Time (hr)')
        end
   end
end