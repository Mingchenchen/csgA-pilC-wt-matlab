function [AggProps] = extractAggregateProperties(folder_name,Nruns,bandwidth,prerun_length,AGG_DENSITY_CUTOFF,AVE_CELL_DENSITY)
    bwd = ([bandwidth bandwidth] ./ [986 740]).^2;
    p = Progress(Nruns);
    AggProps = table();
    for run = 1:Nruns
        p.d(run);
        %% Load simulation results
        load([folder_name '/sim_tracks-' num2str(run) '.mat'])
        sim_tracks = CreateSimTracks(sim_tracks);

        %% Calculate density from agent locations
        K = zeros(2^10,2^10,max(sim_tracks.frame) - prerun_length);
        p = Progress(max(sim_tracks.frame)-prerun_length);
        for ti = prerun_length:max(sim_tracks.frame)
           t = ti - prerun_length + 1;
           p.d(t);
           a = subStruct(sim_tracks,sim_tracks.frame == ti);

           [~, d, ~] = kde2d([a.x a.y],2^10,[0 0],[986 740],bwd(1),bwd(2));
           K(:,:,t) = imresize(d * AVE_CELL_DENSITY * 2^10 * 2^10,[2^10 2^10]);   
        end
        p.done()

        %% Track aggreagtes
        [agg_tracks,stable_aggs,unstable_aggs] = trackAggregates(K,AGG_DENSITY_CUTOFF);
        save([folder_name '/agg_tracks_with_area-' num2str(run)],'agg_tracks','stable_aggs','unstable_aggs');

        %% Calculate agg properties
        real_aggs = subStruct(agg_tracks,ismember(agg_tracks.id,stable_aggs));

        T = table();

        T.majorAxis    = real_aggs.majorAxis;
        T.minorAxis    = real_aggs.minorAxis;
        T.orientation  = real_aggs.orientation;
        T.id           = real_aggs.id;
        T.area         = real_aggs.area;
        T.eccentricity = real_aggs.eccentricity;
        T.mean_density = real_aggs.meanIntensity;
        T.frame        = real_aggs.frame;
        T.movie        = repmat(run,length(real_aggs.frame),1);

        AggProps = [AggProps; T];

    end %run loop
end