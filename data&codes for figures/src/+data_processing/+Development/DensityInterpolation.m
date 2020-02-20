%%%
% Fitts an equation between the KDE of the last frame with the flourscence of the last frame
% to estimate cell desnity from flourscence. Adds:
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

%% 
fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/Exp1_35/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/DEV_1/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/A_30/'};

% fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT New/10282015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT New/10152015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT New/10112015/'};

%% Create First and Last Frame arrays
for set = 1:length(fnames)
    set
    load([fnames{set} 'Knormalized.mat']);
     
    first_frame = K(:,:,end);
    last_frame = K(:,:,end);
    
    save([fnames{set} 'first_frame.mat'],'first_frame')
    save([fnames{set} 'last_frame.mat'],'last_frame')
end

%%
addpath('Libraries/Utils')
SIZE = 1024 * 1024;

floro = zeros(SIZE * 3,1);
kdee  = zeros(size(floro)); 

xSIZE = 1920 * 0.5136; %um
ySIZE = 1440 * 0.5136; %um
AVE_CELL_DENSITY = 1.12; %Average cell density within the field of view
                        %in cells/um^2

bw = zeros(3,2);
%p = Progress(length(fnames));
coeff = [];
% for bandw = 20:1:40
    %bandw
    %bandw = 15 %Linear Correaltion Coefficient;
    %bandw = 29.8 %Spearman Correlation Coefficient;
    %bandw = 23.7674 %kde2d Estimated Correaltion Coefficient for Cotter et al PNAS 2017
    bandw = 32;
    %bandw = 11;
    for set = 1:length(fnames)
        %p.d(set);
        bwd = ([bandw bandw] ./ [986 740]).^2;
        load([fnames{set} 'last_frame.mat']);
        load([fnames{set} 'tracks.mat']);
        a = subStruct(tracks,tracks.frame == max(tracks.frame - 25));
        [bandwidth, d, ~, ~, ~, t_x, t_y] = kde2d([a.x a.y].* 0.5136,2^10,[0 0],[986 740],bwd(1),bwd(2));
        %[bandwidth, d, ~, ~, ~, t_x, t_y] = kde2d([a.x a.y],2^10,[0 0],[986 740]);
        bw(set,:) = bandwidth;

        d = d * AVE_CELL_DENSITY * 2^10 * 2^10;

        floro((SIZE * (set - 1) + 1):(SIZE * set)) = last_frame(:);
        kdee((SIZE * (set - 1) + 1):(SIZE * set)) = d(:);
    end
    
%     coeff = [coeff; corr(floro(:),kdee(:),'Type','Spearman')];
% end

%%
addpath('Libraries/DensityEstimation')
figure, hold on;
    plot(floro(1:200:end),kdee(1:200:end),'.');
    xlabel('Fluorescence intensity (A.U.)');
    ylabel('KDE (cell/\mum^2)');
    kmin = min(floro(:));
    kmax = max(floro(:));
    x = linspace(-2.2e-04,kmax,100);
    plot(x,floro2dens(x,'CSGA'),'.-r','LineWidth',2)


%%
%c1 = 0.00293;
%ft = floro(:) ./ (c1 + floro(:));
%fit = polyfit(floro(:),kdee(:),3);
f = fittype('poly3');

%Scale floro
filt = floro < 3e-3;
f1 = floro(filt);
k = kdee(filt);
[fit3, G] = fit(f1(:),k(:),f,'Robust','bisquare');


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
    x1 = floro(1:100:end)
    k = kdee(1:100:end);
    filt = x1 < 3e-3;
    x1 = x1(filt);
    k = k(filt);
    yy = fit3(x);
    plot(fit3,x1,k)
    xlabel('')
    ylabel('')
    %plot(floro(1:100:end),kdee(1:100:end),'.')
    plot(x,yy,'g-')
    plot(x,floro2dens(x,'CSGA'),'.-r','LineWidth',2)
    plot(x,floro2dens(x,'WT').* 0.4227,'.-k','LineWidth',2)
    
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