function [R] = revs_in_ripple(AllDataTable,AllData,SETS,RIPPLE_START,NBOOT,DENSITY_CUTOFF)
    T = AllDataTable;
    T = T(ismember(T.set,SETS),:);
    T = T(T.TSS1 >= RIPPLE_START,:);

    tracks.density = [];
    for i = SETS
        tracks.density = [tracks.density; AllData{i}.m_tracks.density];
    end
    tracks.density = tracks.density(~isnan(tracks.density));

    R = FigureCreation.Rippleing.bootstrap_revs_in_ripple(T,tracks,NBOOT,DENSITY_CUTOFF); 
end