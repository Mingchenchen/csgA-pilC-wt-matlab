clear all
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src')
load('csgainoutcell')
load('pilcinoutcell')
load('wtinoutcell')
csgameannumcellin=csgameannumcellin(300:end);
csgameannumcellout=csgameannumcellout(300:end);
pilcmeannumcellin=pilcmeannumcellin(300:779);
pilcmeannumcellout=pilcmeannumcellout(300:779);
wtmeannumcellin=wtmeannumcellin(300:end);
wtmeannumcellout=wtmeannumcellout(300:end);

figure

xi = (1:length(csgameannumcellin)) / 2 / 60+1;
m=movingmean(csgameannumcellin,20);
s=movingstd(csgameannumcellin,20)/sqrt(20);
CI95 = tinv([0.025 0.975], 20-1);                    
yCI95 = bsxfun(@times, s', CI95(:));             
boundedline(xi,m,[yCI95(2,:)',yCI95(2,:)'],'-','cmap',[1 0 0],'alpha');
m=movingmean(csgameannumcellout,20);
s=movingstd(csgameannumcellout,20)/sqrt(20);
CI95 = tinv([0.025 0.975], 20-1);                    
yCI95 = bsxfun(@times, s', CI95(:));             
boundedline(xi,m,[yCI95(2,:)',yCI95(2,:)'],'-','cmap',[0 0 0],'alpha');
xlim([1,5]);
ylim([0,6e-3]);
hold on
ratio=csgameannumcellout./csgameannumcellin;
plotyy(xi,m,xi,movingmean(ratio,20));
hold off
xlim([1,5]);
figure
xi = (1:length(wtmeannumcellin)) / 2 / 60+1;
m=movingmean(wtmeannumcellin,20);
s=movingstd(wtmeannumcellin,20)/sqrt(20);
CI95 = tinv([0.025 0.975], 20-1);                    
yCI95 = bsxfun(@times, s', CI95(:));             
boundedline(xi,m,[yCI95(2,:)',yCI95(2,:)'],'-','cmap',[1 0 0],'alpha');
m=movingmean(wtmeannumcellout,20);
s=movingstd(wtmeannumcellout,20)/sqrt(20);
CI95 = tinv([0.025 0.975], 20-1);                    
yCI95 = bsxfun(@times, s', CI95(:));             
boundedline(xi,m,[yCI95(2,:)',yCI95(2,:)'],'-','cmap',[0 0 0],'alpha');
xlim([1,5]);
ylim([0,6e-3]);
hold on
ratio=wtmeannumcellout./wtmeannumcellin;
plotyy(xi,m,xi,movingmean(ratio,20));
hold off
xlim([1,5]);
figure
xi = (1:length(pilcmeannumcellin)) / 2 / 60+1;
m=movingmean(pilcmeannumcellin,20);
s=movingstd(pilcmeannumcellin,20)/sqrt(20);
CI95 = tinv([0.025 0.975], 20-1);                    
yCI95 = bsxfun(@times, s', CI95(:));             
boundedline(xi,m,[yCI95(2,:)',yCI95(2,:)'],'-','cmap',[1 0 0],'alpha');
m=movingmean(pilcmeannumcellout,20);
s=movingstd(pilcmeannumcellout,20)/sqrt(20);
CI95 = tinv([0.025 0.975], 20-1);                    
yCI95 = bsxfun(@times, s', CI95(:));             
boundedline(xi,m,[yCI95(2,:)',yCI95(2,:)'],'-','cmap',[0 0 0],'alpha');

ylim([0,6e-3]);

hold on
ratio=pilcmeannumcellout./pilcmeannumcellin;
plotyy(xi,m,xi,movingmean(ratio,20));
xlim([1,5]);
hold off
