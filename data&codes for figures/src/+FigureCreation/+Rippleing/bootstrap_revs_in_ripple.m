function [R] = bootstrap_revs_in_ripple(T,tracks,NBOOT,DENSITY_CUTOFF)    
    %%% Weight by the number of runs in each movie
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.set)'
        f = T.set == curMovie;
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
        
    tic
    [R.revs.ci,R.revs.samples] = bootci(NBOOT,{@bootfun_RunVectors,T.rho1,DENSITY_CUTOFF},'type','cper','Weights',weights);
    toc   
    
    tic
    [R.trajectories.ci,R.trajectories.samples] = bootci(NBOOT,{@bootfun_Trajectories,tracks.density,DENSITY_CUTOFF},'type','cper');
    toc   
end

function [prct] = bootfun_RunVectors(rhos,DENSITY_CUTOFF) 
    N = length(rhos);
    prct = sum(rhos > DENSITY_CUTOFF) / N;
end

function [prct] = bootfun_Trajectories(rhos,DENSITY_CUTOFF)
    N = length(rhos);
    prct = sum(rhos > DENSITY_CUTOFF) / N;
end
    