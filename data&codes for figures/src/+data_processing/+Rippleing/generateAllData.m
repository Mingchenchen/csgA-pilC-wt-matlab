function [AllData] = generateAllData(MAP,PROJECT_BASENAME,WRAP_STOPS,MAX_FRAME,tracks_fname)
    AllData = {};

    assert(length(unique(keys(MAP))) == length(keys(MAP)),'All keys in MAP must be unique')
    
    k = keys(MAP);
    p = Progress(length(MAP));
    for i = 1:length(MAP)
        p.d(i)

        folder = [PROJECT_BASENAME, '/data/processed/', MAP(k{i})];

        load([folder '/' tracks_fname]);

        %%
        AllData{i}.movie = k{i};
        AllData{i}.all_runs = createRunVectors(ripple_tracks,WRAP_STOPS);
        AllData{i}.m_tracks = ripple_tracks;
        AllData{i}.AlignedStart = 1;
        AllData{i}.AlignedStop = MAX_FRAME;

        AllData{i}.STOPS_WRAPPED = WRAP_STOPS;
    end
end