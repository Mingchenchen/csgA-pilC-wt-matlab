function [outside,inside] = bootstrap_state_transitions_in_time(T,NBOOT,STABLE_AGG_CUTOFF,WINDOW,DENSITY_CUTOFF)

    outsideFilt = T.rho1 < DENSITY_CUTOFF;
    tic
    outside = FigureCreation.bootstrap_state_transitions_in_time(T(outsideFilt,:),NBOOT,WINDOW);
    toc
    
    insideFilt = T.rho1 >= DENSITY_CUTOFF & T.TSS1 > STABLE_AGG_CUTOFF;
    inside =  FigureCreation.bootstrap_state_transitions_in_time(T(insideFilt,:),NBOOT,WINDOW);
end