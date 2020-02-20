clear all
close all
addpath('..\src\Libraries\LAPTracker');
addpath('..\src\Libraries\Utils');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\DensityEstimation');
addpath('..\src\Libraries\TrajectoryAnalysis');
addpath('..\src');

load('../pilC/pilCDataTable')
load('../csgA/csgADataTable')
load('../WT/WTDataTable')
Nboot=10;
csgAin=csgADataTable(csgADataTable.Dn1<0,:);
WTin=WTDataTable(WTDataTable.Dn1<0,:);
pilCin=pilCDataTable(pilCDataTable.Dn1<0,:);
csgAout=csgADataTable(csgADataTable.Dn1>0,:);
WTout=WTDataTable(WTDataTable.Dn1>0,:);
pilCout=pilCDataTable(pilCDataTable.Dn1>0,:);
WTinrun=WTin(WTin.state0~=3,:);
WTinstop=WTin(WTin.state0==3,:);
pilCinrun=pilCin(pilCin.state0~=3,:);
pilCinstop=pilCin(pilCin.state0==3,:);
csgAinrun=csgAin(csgAin.state0~=3,:);
csgAinstop=csgAin(csgAin.state0==3,:);
WToutrun=WTout(WTout.state0~=3,:);
WToutstop=WTout(WTout.state0==3,:);
pilCoutrun=pilCout(pilCout.state0~=3,:);
pilCoutstop=pilCout(pilCout.state0==3,:);
csgAoutrun=csgAout(csgAout.state0~=3,:);
csgAoutstop=csgAout(csgAout.state0==3,:);
%% Panel A
y=[mean(WTinrun.Rs0),mean(pilCinrun.Rs0),mean(csgAinrun.Rs0);mean(WToutrun.Rs0),mean(pilCoutrun.Rs0),mean(csgAoutrun.Rs0)];

ngroups = size(y, 1);
nbars = size(y, 2);
bootfun = @(x)mean(x);
ci=[bootci(Nboot,bootfun,WTinrun.Rs0),bootci(Nboot,bootfun,pilCinrun.Rs0),bootci(Nboot,bootfun,csgAinrun.Rs0) ];
co=[bootci(Nboot,bootfun,WToutrun.Rs0),bootci(Nboot,bootfun,pilCoutrun.Rs0),bootci(Nboot,bootfun,csgAoutrun.Rs0) ];
ypos=[ci(2,1), ci(2,2),ci(2,3);co(2,1),co(2,2),co(2,3)]-y;
yneg=y-[ci(1,1), ci(1,2),ci(1,3);co(1,1),co(1,2),co(1,3)];
% Calculating the width for each bar group+
figure
bar(y)
hold on
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), yneg(:,i),ypos(:,i), '.');
end
hold off
%% Panel B
y=[mean(WTinrun.Rt0),mean(pilCinrun.Rt0),mean(csgAinrun.Rt0);mean(WToutrun.Rt0),mean(pilCoutrun.Rt0),mean(csgAoutrun.Rt0)];

ngroups = size(y, 1);
nbars = size(y, 2);
bootfun = @(x)mean(x);
ci=[bootci(Nboot,bootfun,WTinrun.Rt0),bootci(Nboot,bootfun,pilCinrun.Rt0),bootci(Nboot,bootfun,csgAinrun.Rt0) ];
co=[bootci(Nboot,bootfun,WToutrun.Rt0),bootci(Nboot,bootfun,pilCoutrun.Rt0),bootci(Nboot,bootfun,csgAoutrun.Rt0) ];
ypos=[ci(2,1), ci(2,2),ci(2,3);co(2,1),co(2,2),co(2,3)]-y;
yneg=y-[ci(1,1), ci(1,2),ci(1,3);co(1,1),co(1,2),co(1,3)];
% Calculating the width for each bar group
figure
bar(y)
hold on
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), yneg(:,i),ypos(:,i), '.');
end
hold off

%% Panel C
WTruns=WTDataTable(WTDataTable.state1~=3,:);
pilCruns=pilCDataTable(pilCDataTable.state1~=3,:);
csgAruns=csgADataTable(csgADataTable.state1~=3,:);
WTtow=WTruns(cos(WTruns.beta1)>0,:);
WTawy=WTruns(cos(WTruns.beta1)<0,:);
pilCtow=pilCruns(cos(pilCruns.beta1)>0,:);
pilCawy=pilCruns(cos(pilCruns.beta1)<0,:);
csgAtow=csgAruns(cos(csgAruns.beta1)>0,:);
csgAawy=csgAruns(cos(csgAruns.beta1)<0,:);
figure
bias=[(mean(WTtow.Rt1)-mean(WTawy.Rt1))/mean(WTDataTable.Rt1),(mean(pilCtow.Rt1)-mean(pilCawy.Rt1))/mean(pilCDataTable.Rt1),(mean(csgAtow.Rt1)-mean(csgAawy.Rt1))/mean(csgADataTable.Rt1)];
hold on
h=bar(1,bias(1));
set(h,'FaceColor','b');
h=bar(2,bias(2));
set(h,'FaceColor','g');
h=bar(3,bias(3));
set(h,'FaceColor','y');
hold off

%% Panel D
y=[mean(WTinstop.Rt0),mean(pilCinstop.Rt0),mean(csgAinstop.Rt0);mean(WToutstop.Rt0),mean(pilCoutstop.Rt0),mean(csgAoutstop.Rt0)];

ngroups = size(y, 1);
nbars = size(y, 2);
bootfun = @(x)mean(x);
ci=[bootci(Nboot,bootfun,WTinstop.Rt0),bootci(Nboot,bootfun,pilCinstop.Rt0),bootci(Nboot,bootfun,csgAinstop.Rt0) ];
co=[bootci(Nboot,bootfun,WToutstop.Rt0),bootci(Nboot,bootfun,pilCoutstop.Rt0),bootci(Nboot,bootfun,csgAoutstop.Rt0) ];
ypos=[ci(2,1), ci(2,2),ci(2,3);co(2,1),co(2,2),co(2,3)]-y;
yneg=y-[ci(1,1), ci(1,2),ci(1,3);co(1,1),co(1,2),co(1,3)];
% Calculating the width for each bar group
figure
bar(y)
hold on
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:,i), yneg(:,i),ypos(:,i), '.');
end
hold off

%% Panel E
WTinrun2stop=WTinrun(WTinrun.state1==3,:);
WToutrun2stop=WToutrun(WToutrun.state1==3,:);
pilCinrun2stop=pilCinrun(pilCinrun.state1==3,:);
pilCoutrun2stop=pilCoutrun(pilCoutrun.state1==3,:);
csgAinrun2stop=csgAinrun(csgAinrun.state1==3,:);
csgAoutrun2stop=csgAoutrun(csgAoutrun.state1==3,:);
prob=[length(WTinrun2stop.Rd0)/length(WTinrun.Rd0),length(pilCinrun2stop.Rd0)/length(pilCinrun.Rd0),length(csgAinrun2stop.Rd0)/length(csgAinrun.Rd0);length(WToutrun2stop.Rd0)/length(WToutrun.Rd0),length(pilCoutrun2stop.Rd0)/length(pilCoutrun.Rd0),length(csgAoutrun2stop.Rd0)/length(csgAoutrun.Rd0)];
figure
bar(prob)