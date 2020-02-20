function [R,xi] = bootstrap_runs_in_beta(xi,T,WINDOW,NBOOT)
    %%% Weight by the number of runs in each movie
    CUTOFF = cos(pi/4);
    
    weights = zeros(size(T,1),1);
    for curMovie = unique(T.set)'
        f = T.set == curMovie;
        weights(f) = sum(f);
    end
    weights = 1 - (weights / sum(weights));
    
    bootfun = @(x,y) movingmeanxy(x,y,WINDOW,xi);
    
    tic
        [R.towards.Rd.ci,R.towards.Rd.samples] = bootci(NBOOT,{bootfun,T.Dn1(cos(T.beta1) > CUTOFF),T.Rd1(cos(T.beta1) > CUTOFF)},'type','cper','Weights',weights(cos(T.beta1) > CUTOFF)/sum(weights(cos(T.beta1) > CUTOFF)));
    toc
        [R.away.Rd.ci,R.away.Rd.samples]       = bootci(NBOOT,{bootfun,T.Dn1(cos(T.beta1) < -CUTOFF),T.Rd1(cos(T.beta1) < -CUTOFF)},'type','cper','Weights',weights(cos(T.beta1) < -CUTOFF)/sum(weights(cos(T.beta1) < -CUTOFF)));
        
        [R.towards.Rs.ci,R.towards.Rs.samples] = bootci(NBOOT,{bootfun,T.Dn1(cos(T.beta1) > CUTOFF),T.Rs1(cos(T.beta1) > CUTOFF)},'type','cper','Weights',weights(cos(T.beta1) > CUTOFF)/sum(weights(cos(T.beta1) > CUTOFF)));
        [R.away.Rs.ci,R.away.Rs.samples]       = bootci(NBOOT,{bootfun,T.Dn1(cos(T.beta1) < -CUTOFF),T.Rs1(cos(T.beta1) < -CUTOFF)},'type','cper','Weights',weights(cos(T.beta1) < -CUTOFF)/sum(weights(cos(T.beta1) < -CUTOFF)));

        [R.towards.Rt.ci,R.towards.Rt.samples] = bootci(NBOOT,{bootfun,T.Dn1(cos(T.beta1) > CUTOFF),T.Rt1(cos(T.beta1) > CUTOFF)},'type','cper','Weights',weights(cos(T.beta1) > CUTOFF)/sum(weights(cos(T.beta1) > CUTOFF)));
        [R.away.Rt.ci,R.away.Rt.samples]       = bootci(NBOOT,{bootfun,T.Dn1(cos(T.beta1) < -CUTOFF),T.Rt1(cos(T.beta1) < -CUTOFF)},'type','cper','Weights',weights(cos(T.beta1) < -CUTOFF)/sum(weights(cos(T.beta1) < -CUTOFF)));
end