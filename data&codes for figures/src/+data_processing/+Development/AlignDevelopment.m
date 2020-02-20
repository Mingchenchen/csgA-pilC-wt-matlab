% %%%
% % CsgA in WT
% %%%
% fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/A_30/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/A_31/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/DEV_1/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/Exp1_35/' , ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/Exp1_36/' , ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/csgA in WT/Exp1_38/'}
                  
% %%% 
% % WT in WT
% %%%
% fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10282015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10152015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10112015/'};
      
%%%
% Combined
%%%
fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT New/10112015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT New/10152015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT New/10282015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/A_30/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/DEV_1/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/Exp1_35/'};
      
save_folder = '/Users/cotter/GitHub/CotterPNAS2017/Data/';
      
filenames = {};
for set = 1:length(fnames)
    s = strsplit(fnames{set},'/');
    filenames{set} = s{9};
end

addpath('Libraries/Utils')

%%
w1 = hamming(1024,'periodic');
w2 = w1(:) * w1(:).';

agg_curves = {};

for set = 1:length(fnames)
    set
    load([fnames{set} 'Kum.mat']);

    %%
    a = zeros(size(Kum,3),1);
    p = Progress(length(a))
    for i = 1:size(Kum,3)
        p.d(i)
        A = Kum(:,:,i);
        fA = abs(fft2(A .* w2));
        a(i) = sum(sum(fA(3:6,3:6)));
    end
    p.done()
    
    %%
    agg_curves{set} = a;
    
    clear Kum
end

save([save_folder, 'agg_curves'],'agg_curves')

%%
figure, hold on
for i = 1:length(agg_curves)
    plot(agg_curves{i})
end

ylabel('Total Amplitude');
xlabel('Aligned Frame')
%line([maxAgg maxAgg],get(gca,'YLim'))
grid on;
legend(filenames,'Location','SouthEast')

%%
nCurves = {}
for i = 1:length(agg_curves);
    c = (agg_curves{i} - min(agg_curves{i})) / range(agg_curves{i});
    nCurves{i} = c(:);
end

%% Plot Results
figure, hold on;

for i = 1:length(nCurves)
    c = nCurves{i}
    plot(c);
end
ylabel('Normailzed Amplitude');
xlabel('Frame')
grid on;
legend(filenames,'Location','SouthEast')

agg_start = zeros(length(nCurves),1);

%% Find centerpoint in aggregate growth by finding the min MSE between
% the aggregate growth plot and a simple step function
%

agg_start = zeros(length(nCurves),1);

% Detect where aggreagtion begins
for j = 1:length(nCurves)
    c1 = nCurves{j};
    
    agg_start(j) = find(c1 > 0.1,1,'first');
end

% % agg_start = zeros(length(nCurves),1);
% % 
% % % Detect where aggreagtion begins
% % p = Progress(length(nCurves))
% % for j = 1:length(nCurves)
% %     p.d(j)
% %     c1 = nCurves{j};
% %     mse = zeros(size(c1));
% %     for i = 0:length(c1)-1
% %         nc = (c1 - min(c1)) / range(c1);
% %         %fx = [zeros(i,1); [0.5]; linspace(1,nc(end),length(nc) - 1 - i)'];
% %         fx = [zeros(i,1); [0.5]; ones(length(nc) - 1 - i,1)];
% %         mse(i+1) = nanmean((nc - fx).^2);
% %     end
% %     [~,I] = min(mse);
% %     agg_start(j) = I;
% % end
% % p.done()

%% Plot the aligned graphs by padding the data
padding = max(agg_start) - agg_start;
plot_types = {'-','--',':','.-'};

figure, hold on
for i =  1:length(nCurves)
    c1 = agg_curves{i};
    c1 = c1-c1(1)
    plot([zeros(padding(i),1); c1]);
end
ylabel('Total Amplitude');
xlabel('Aligned Frame')
grid on;
legend(filenames,'Location','SouthEast')

%% Crop the data so that they all start at the same time
leading_crop = agg_start - min(agg_start) + 1;
len = min(cellfun(@length,agg_curves)' - (agg_start - min(agg_start)));

plot_types = {'-','--',':','.-'};

figure, hold on
for i =  1:length(agg_curves)
    c1 = nCurves{i}
    plot(c1(leading_crop(i):end));
end
ylabel('Total Amplitude');
xlabel('Aligned Frame')
grid on;
legend(filenames,'Location','SouthEast')

%% Find the maximum poin in the curve and crop the data to that point
% This is assumed to be just before aggreagtes begin to dissapear

start = agg_start - min(agg_start) + 1;
len = min(cellfun(@length,agg_curves)' - (agg_start - min(agg_start)));

conCurves = zeros(len,length(filenames));

for i =  1:length(agg_curves)
    c1 = agg_curves{i}; % - min(curves{i});
    c1 = c1 - min(c1)
    conCurves(:,i) = c1(start(i):(start(i) + len-1));
end


figure, hold on
plot(conCurves)
plot(mean(conCurves,2),'--k')
ylabel('Total Amplitude');
xlabel('Aligned Frame')
[~,maxAgg] = max(mean(conCurves,2));
%line([maxAgg maxAgg],get(gca,'YLim'))
grid on;
legend(filenames,'Location','SouthEast')

%% Create the structure that will hold the useful data during the analysis
AllData = {};

for set = 1:length(fnames)
    AllData{set}.set = fnames{set};
    AllData{set}.AlignedStart = start(set);
    AllData{set}.AlignedStop = start(set) + len-1;
end

save([save_folder, 'AllData_Combined_Init.mat'],'AllData')


