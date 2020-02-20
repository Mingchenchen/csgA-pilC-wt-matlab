function [R] = bootstrap_angle_to_agg(xi,T,WINDOW,NBOOT)
    %%% Weight by the number of runs in each movie
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.set)'
        f = T.set == curMovie;
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
    
    bootfun = @(x,y) movingmeanxy(x,y,WINDOW,xi);

    %Real
    tic 
        [R.ci,R.samples] = bootci(NBOOT,{bootfun,T.Dn1,cos(2 * T.beta1)},'type','per','Weight',weights);
    toc

    %Random
    x = T.beta1(randperm(height(T))); %This is dummy data. 

    tic 
        [R.rand.ci,R.rand.samples] = bootci(NBOOT,{bootfun,T.Dn1,cos(2 * x)},'type','per','Weight',weights);
    toc
end