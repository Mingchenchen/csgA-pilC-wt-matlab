function [R] = bootstrap_runs_in_time(T,NBOOT,WINDOW)    
    %%% Weight by the number of runs in each movie
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.set)'
        f = T.set == curMovie;
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
    
    xi = sort(unique(T.TSS1));
    bootfun = @(x,y) movingmeanxy(x,y,WINDOW,xi);
    filt_weights = weights / sum(weights);
    
    tic
    [R.Rd.ci,R.Rd.samples] = bootci(NBOOT,{bootfun,T.TSS1,T.Rd1},'type','cper','Weights',filt_weights);
    toc
    [R.Rs.ci,R.Rs.samples] = bootci(NBOOT,{bootfun,T.TSS1,T.Rs1},'type','cper','Weights',filt_weights);
    [R.Rt.ci,R.Rt.samples] = bootci(NBOOT,{bootfun,T.TSS1,T.Rt1},'type','cper','Weights',filt_weights);
    R.xi = xi;    
end