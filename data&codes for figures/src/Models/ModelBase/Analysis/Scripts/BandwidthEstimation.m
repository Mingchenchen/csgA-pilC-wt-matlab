%%
fnames = {'/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/10282015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/10152015/', ...
          '/Volumes/Scratch/Chris/Paper2_Rippling/Results/Development/WT in WT/10112015/'};

%%
for set = 1:length(fnames)
    set
    load([fnames{set} 'Kum.mat']);
     
    first_frame = Kum(:,:,end);
    last_frame = Kum(:,:,end);
    
    save([fnames{set} 'first_frame_um.mat'],'first_frame')
    save([fnames{set} 'last_frame_um.mat'],'last_frame')
end
% %%
% N = 10000;
% 
% for set = 1:length(fnames)
%         %p.d(set);
%         bwd = ([bandw bandw] ./ [986 740]).^2;
%         load([fnames{set} 'last_frame_um.mat']);
% end
%%

N = 10000;
vals=zeros(2,N);
for i=1:N
    [vals(1,i),vals(2,i)]=pinky(linspace(0,986,2^10),linspace(0,740,2^10),first_frame);
end

%%
bandw = 15
bwd = ([bandw bandw] ./ [986 740]).^2;
[~, d, ~] = kde2d(vals',2^10,[0 0],[986 740],bwd(1),bwd(2));
d = d * 1.12 * 2^10 * 2^10;

sum((first_frame(:) - d(:)).^2) / numel(first_frame)