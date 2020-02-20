function [ Data] = CalculateData(Data, m_tracks, agg_tracks, stable_aggs, unstable_aggs, varargin)
    %Takes in Data = AllData{i} and adds
    %   data.all_runs
    %   data.agg_tracks
    %   data.stable_aggs -- list of stable agg ids
    %   data.tracks

    p = inputParser;
    addOptional(p,'WrapStops', false);
    addOptional(p,'IncludeUnstable', false);
    addOptional(p,'AggDensityCutoff',2.32);
    parse(p,varargin{:});
    
    m_tracks.in_agg = m_tracks.density > p.Results.AggDensityCutoff;
    all_runs = createRunVectors(m_tracks,p.Results.WrapStops);

    agg_tracks.stable = ismember(agg_tracks.id,stable_aggs);
    if(p.Results.IncludeUnstable)
         %Use both stable and unstable aggreagtes
        agg_tracks = subStruct(agg_tracks,ismember(agg_tracks.id,union(stable_aggs,unstable_aggs)),'r');
    else
        % Sort out only the filtered stable aggreagtes
        % This is the version used by Cotter et al. PNAS 2017 
        agg_tracks = subStruct(agg_tracks,ismember(agg_tracks.id,stable_aggs),'r');
    end


    Data.agg_tracks = agg_tracks;
    Data.all_runs = all_runs;
    Data.stable_aggs = stable_aggs;
    Data.tracks = m_tracks;

end

