function [R] = bootstrap_neighbor_alignment(T,WINDOW,NBOOT)
    %%% Weight by the number of runs in each movie
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.set)'
        f = T.set == curMovie;
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
    
    xi = 1:max(T.TSS1);
    
    bootfun = @(x,y) movingmeanxy(x,y,WINDOW,xi);

    %Real
    tic 
        [R.ci,R.samples] = bootci(NBOOT,{bootfun,T.TSS1,T.neighbor_alignment},'type','per','Weight',weights);
    toc

    %Random
    Trand = zeros(height(T),1);

    for i = 1:height(T)
        Trand(i) = mean(cos(2 * (T.orientation1(i) - T.orientation1(randi(height(T),3,1)))));
    end

    tic 
        [R.rand.ci,R.rand.samples] = bootci(NBOOT,{bootfun,T.TSS1,Trand},'type','per','Weight',weights);
    toc
    R.xi = xi;
end