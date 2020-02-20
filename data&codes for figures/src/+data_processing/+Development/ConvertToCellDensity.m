% %%%
% % CsgA in WT
% %%%
fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/Exp1_35/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/DEV_1/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/csgA in WT/A_30/'};
background = 'CSGA';

%%% 
% WT in WT
%%%
% fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10282015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10152015/', ...
%           '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Csga Development/WT in WT/10112015/'};
% background = 'WT';
      
addpath('Libraries/DensityEstimation')

for set = 2:length(fnames)
    set
    load([fnames{set} 'Knormalized.mat']);
    
    if(strcmp(background,'WT'))
        %
        % 0.4227 is a workaround. Floro2dens assumes a 2.66 cell/um^2
        % average cell density. The density used was 1.12 cell/um^2
        Kum = floro2dens(K,background) .* 0.4227;
    elseif(strcmp(background,'CSGA'))
        Kum = floro2dens(K,background);
    else
        error('Not a valid background, choose either WT or CSGA')
    end
    
    save([fnames{set} 'Kum.mat'],'Kum','-v7.3');
    clear K Kum
end