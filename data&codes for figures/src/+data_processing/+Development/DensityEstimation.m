%%%
% Fitts an equation between the KDE of the last frame with the flourscence of the last frame
% to estimate cell desnity from flourscence. Adds:pwd
%    AllData{set}.coeffs:
%           the coefficients for converting Knormilzed values to 
%           cells/um^2; cells/um^2 = x * AllData{set}.coeffs(1) + AllData{set}.coeffs(2)
%    AllData{set}.dens_est(x)
%            converst x = Knormzlied into cells/um^2 using equation: 
%               x * AllData{set}.coeffs(1) + AllData{set}
%
% This was required because the amount of background noise in the flourscence images
% made direct calculation impossible and the dark image required to subtrack background
% was not created. 
%%%
run('ENVS.m')

%% Create Friast and Last Frame arrays
for set = 1:length(fnames)
    set
    load([fnames{set} 'Kum.mat']);
     
    first_frame = Kum(:,:,AllData{set}.AlignedStart);
    last_frame = Kum(:,:,AllData{set}.AlignedStop);
    
    save([fnames{set} 'first_frame.mat'],'first_frame')
    save([fnames{set} 'last_frame.mat'],'last_frame')
end

%%
SIZE = 1024 * 1024;

floro = zeros(SIZE * 3,1);
kdee  = zeros(size(floro)); 

xSIZE = 1920 * 0.5136; %um
ySIZE = 1440 * 0.5136; %um
CELLS_IN_FRAME = 8.2e5; %Cells
%CELLS_IN_FRAME = 2.66 * xSIZE * ySIZE;

bw = zeros(3,2);
%p = Progress(length(fnames));
coeff = [];
%for bandw = 20:1:30
    %bandw
    %bandw = 22 %Linear Correaltion Coefficient;
    %bandw = 29.8 %Spearman Correlation Coefficient;
    bandw = 23.7674 %kde2d Estimated Correaltion Coefficient
    %bandw = 11;
    for set = 1:length(fnames)
        %p.d(set);
        bwd = ([bandw bandw] ./ [986 740]).^2;
        load([fnames{set} 'last_frame.mat']);
        load([fnames{set} 'm_tracks.mat']);
        a = subStruct(m_tracks,m_tracks.frame == AllData{set}.AlignedStop);
        [bandwidth, d, ~, ~, ~, t_x, t_y] = kde2d([a.x a.y],2^10,[0 0],[986 740],bwd(1),bwd(2));
        %[bandwidth, d, ~, ~, ~, t_x, t_y] = kde2d([a.x a.y],2^10,[0 0],[986 740]);
        bw(set,:) = bandwidth;

        d = d * CELLS_IN_FRAME * 2^10 / xSIZE * 2^10 / ySIZE;

        floro((SIZE * (set - 1) + 1):(SIZE * set)) = last_frame(:);
        kdee((SIZE * (set - 1) + 1):(SIZE * set)) = d(:);
    end
    
    %coeff = [coeff; corr(floro(:),kdee(:),'Type',)];
%end

%%
figure, hold on;
    plot(floro(1:200:end),kdee(1:200:end),'.');
    xlabel('Fluorescence intensity (A.U.)');
    ylabel('KDE (cell/\mum^2)');
    kmin = min(floro(:));
    kmax = max(floro(:));
    x = linspace(-2.2e-04,kmax,100);
    plot(x,floro2dens(x),'.-r','LineWidth',2)


%%
%c1 = 0.00293;
%ft = floro(:) ./ (c1 + floro(:));
%fit = polyfit(floro(:),kdee(:),3);
f = fittype('poly3');

%Scale floro
%fmin = min(floro(:));
%fmax = max(floro(:));
%frange = (fmax - fmin);
%floro1 = (floro - fmin) / (fmax - fmin);
[fit3, G] = fit(floro(:),kdee(:),f,'Robust','bisquare');
% X = floro(1:100:end);
% Y = kdee(1:100:end);
% [X,I] = sort(X);
% Y = Y(I);
% tic
% yy = smooth(X,Y,0.3,'loess');
% toc
%floro2dens = @(x) fit(1) * (x ./ (c1 + x)) + fit(2);
 %[brob, STATS] = robustfit(ft(:),kdee(:));
 %fun =  @(x) brob(2) * x + brob(1);
 %floro2dens = @(x) brob(2) * (x ./ (c1 + x)) + brob(1);
 %STATS.robust_s

%%
%x = linspace(min(floro(:)),max(floro(:)),100);
%y = polyval(fit,x);

figure, hold on
    %plot(X,Y,'.');
    %;
    %plot(x,y,'r');
    %y = fit3(((floro(1:100:end) .* 0.4227) - fmin) / (fmax - fmin)) .* frange + fmin;
    %scatterDensity(floro(1:100:end) .* 0.4227,kdee(1:100:end));
    %plot(floro(1:100:end),y,'-r');
    %plot(fit3,'r-',);
    x = sort(floro(1:100:end));
    yy = fit3(x) .* 0.4227;
    plot(fit3,floro(1:100:end),kdee(1:100:end))
    xlabel('')
    ylabel('')
    %plot(floro(1:100:end),kdee(1:100:end),'.')
    plot(x,yy,'g-')
    
    
%%
figure, hold on
    xx = floro2dens(floro(1:100:end)) .* 0.4039;
    plot(xx,kdee(1:100:end),'.')
    %plot(x,yy,'r-')
 %%
lf_kernel_density = zeros(71,3);
ff_kernel_density = zeros(21,3);

figure, hold on;
    for set = 1:3
        r = (SIZE * (set - 1) + 1):(SIZE * set);
        plot3(floro(r),kdee(r),repmat(4 - set,numel(kdee(r)),1),'.','color',color_chooser(1))
    end
    title(num2str(corr(floro2dens(floro(:)),kdee(:))))
    kmin = min(floro(:));
    kmax = max(floro(:));
    x = linspace(kmin,kmax,25);
    plot3(x,floro2dens(x),repmat(set+1,25,1),'.-r','LineWidth',5)
    
    xlabel('Floro (A.U.)');
    ylabel('KDE (cell/\mum^2)');
   
    
    
    
%%
for set = 1:length(fnames)
    load([fnames{set} 'last_frame.mat']);
    load([fnames{set} 'first_frame.mat']);
    
    figure,
        imagesc(floro2dens(last_frame)) 
        title(num2str(set))
        colorbar
        %caxis([0 50])

   %[lf_kernel_density(:,set)]=ksdensity(fun(last_frame(:)),0:0.5:35,'Bandwidth', 0.0298);
   %[ff_kernel_density(:,set)]=ksdensity(fun(first_frame(:)),0:0.5:10,'Bandwidth', 0.0477);
end

% figure
%     plot(0:0.5:35,lf_kernel_density);
% figure
%     plot(0:0.5:10,ff_kernel_density);

%%
frame = 'first_frame.mat';
load([fnames{3} frame]);
[r3,xi3]=ksdensity(floro2dens(first_frame(:)));
mean(floro2dens(first_frame(:)))

load([fnames{2} frame]);
[r2,xi2]=ksdensity(floro2dens(first_frame(:)));
mean(floro2dens(first_frame(:)))

load([fnames{1} frame]);
[r1,xi1]=ksdensity(floro2dens(first_frame(:)));
mean(floro2dens(first_frame(:)))

figure, hold on
    plot(xi1,r1)
    plot(xi2,r2)
    plot(xi3,r3)
    legend('1','2','3')

%%
frame = 'last_frame.mat';
load([fnames{3} frame]);
[r3,xi3]=ksdensity(floro2dens(last_frame(:)));
mean(floro2dens(last_frame(:)))

load([fnames{2} frame]);
[r2,xi2]=ksdensity(floro2dens(last_frame(:)));
mean(floro2dens(last_frame(:)))

load([fnames{1} frame]);
[r1,xi1]=ksdensity(floro2dens(last_frame(:)));
mean(floro2dens(last_frame(:)))

figure, hold on
    plot(xi1,r1)
    plot(xi2,r2)
    plot(xi3,r3)
    legend('1','2','3')