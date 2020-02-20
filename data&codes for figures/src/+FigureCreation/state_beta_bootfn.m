function [res] = state_beta_bootfn(d,s0,s1,xi,WINDOW)
% boot function used in bootstrap_runs_in_beta.m to bootstrap the P(state1 = non-persistent| state0 = persistent)
% d = distance run started from aggreagte (T.Dn1)
% s0 = state0 (T.state0)
% s1 = state1 (T.state1)
% xi d points at which to calculate the probability
% WINDOW window around xi to calcalute the probaiblty in.
    res = NaN(length(xi),1);
    
    i = 1;
    for x = xi
        filt = d > x - WINDOW & d < x + WINDOW; 
        res(i) = sum(s0(filt) < 3 & s1(filt) == 3) ./ sum(s0(filt) < 3);
        i = i + 1;
    end
end