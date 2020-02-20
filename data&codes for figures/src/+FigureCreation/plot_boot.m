function [h] = plot_boot(xi,R,varargin)
    hold on
    lower = mean(R.samples)' - R.ci(1,:)';
    upper = R.ci(2,:)' - mean(R.samples)';
    h = boundedline(xi,mean(R.samples)',[lower,upper],varargin{:})
    hold off
end