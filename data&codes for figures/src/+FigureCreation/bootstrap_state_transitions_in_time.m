function [P3] = bootstrap_state_transitions_in_time(T,NBOOT,WINDOW)  
addpath('./');
    xi = sort(unique(T.TSS1))';
    
    bootfun1 = @(Ti) sum(Ti.state0 < 3 & Ti.state1 == 3) / sum(Ti.state0 < 3);
    bootfun2 = @(Tb,x) bootfun1(Tb(Tb.TSS1 > x - WINDOW & Tb.TSS1 < x + WINDOW,:));
    bootfun3 = @(Tb) bootloop(Tb,xi,bootfun2);

    %%% Weight by the number of runs in each movie
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.movie)'
        f = strcmp(curMovie,T.movie);
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
    filt_weights = weights / sum(weights);
    
    tic
    [P3.ci,P3.samples] = bootci(NBOOT,{bootfun3,T},'type','cper','Weights',filt_weights);
    toc
    
    P3.xi = xi;
end