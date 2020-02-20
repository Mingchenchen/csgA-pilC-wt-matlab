function [R] = bootstrap_runs_in_time(T,NBOOT,STABLE_AGG_CUTOFF,WINDOW,DENSITY_CUTOFF)    
    %%% Weight by the number of runs in each movie
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.set)'
        f = T.set == curMovie;
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
    
    %%%
    % Outside
    outsideFilt = T.rho1 < DENSITY_CUTOFF;
    outsideTSS  = T.TSS1(outsideFilt);
    xi = sort(unique(outsideTSS));
    bootfun = @(x,y) movingmeanxy(x,y,WINDOW,xi);
    filt_weights = weights(outsideFilt) / sum(weights(outsideFilt));
    
    tic
    [R.outside.Rd.ci,R.outside.Rd.samples] = bootci(NBOOT,{bootfun,outsideTSS,T.Rd1(outsideFilt)},'type','cper','Weights',filt_weights);
    toc
    [R.outside.Rs.ci,R.outside.Rs.samples] = bootci(NBOOT,{bootfun,outsideTSS,T.Rs1(outsideFilt)},'type','cper','Weights',filt_weights);
    [R.outside.Rt.ci,R.outside.Rt.samples] = bootci(NBOOT,{bootfun,outsideTSS,T.Rt1(outsideFilt)},'type','cper','Weights',filt_weights);
    R.outside.xi = xi;
    
    %%%
    % Inside
    insideFilt = T.rho1 >= DENSITY_CUTOFF & T.TSS1 > STABLE_AGG_CUTOFF;
    insideTSS  = T.TSS1(insideFilt);
    xi = sort(unique(insideTSS));
    bootfun = @(x,y) movingmeanxy(x,y,WINDOW,xi);
    filt_weights = weights(insideFilt) / sum(weights(insideFilt));
    
    [R.inside.Rd.ci,R.inside.Rd.samples] = bootci(NBOOT,{bootfun,insideTSS,T.Rd1(insideFilt)},'type','cper','Weights',filt_weights);
    [R.inside.Rs.ci,R.inside.Rs.samples] = bootci(NBOOT,{bootfun,insideTSS,T.Rs1(insideFilt)},'type','cper','Weights',filt_weights);
    [R.inside.Rt.ci,R.inside.Rt.samples] = bootci(NBOOT,{bootfun,insideTSS,T.Rt1(insideFilt)},'type','cper','Weights',filt_weights);
    R.inside.xi = xi;
    
end