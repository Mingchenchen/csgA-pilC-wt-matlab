%%%
% Creates Figures on cell bias towards aggreagtes
%%%
addpath('Libraries/Utils')
PROJECT_BASENAME = ['/Volumes/Scratch/Chris/Paper2_Rippling/Paper2/'];
run('ENVS.m')

% Publication: setting to true supresses some labeling during plotting
% to make plot publication ready, but less visually readable. 
PUB_READY = false;

SAVE = false;
SAVE_FOLDER = [PROJECT_BASENAME '/reports/figures/rippleing/Wavelengths/'];
if(SAVE)
    mkdir(SAVE_FOLDER)
end

% Figure Properties
FONT_SIZE = 12;
FIG_SIZE = [20, 20];

%% Check ripple wavelength estimation results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure, hold on
NPoints = 2^10;

MAP = ALL_RIPPLE_FOLDERS;
% MAP = WT_RIPPLE_FOLDERS;
% MAP = CSGA_RIPPLE_FOLDERS;
% MAP = DIFA_RIPPLE_FOLDERS;

c = 0;
for i = keys(MAP)
    i
    k = i{1};
    
    folder = [PROJECT_BASENAME, '/data/processed/', MAP(k)];

    load([folder 'wavelengths.mat'])
    Nframes = size(wavelengths,2);
    
    wl = NPoints ./ (0:(NPoints / 2 + 1));
    
    %
    % wavelengths is of dimensions: (time,DFT bin)
    total_magnitude = sum(wavelengths,1);
    D = wavelengths ./ total_magnitude;
    D = wl' .* D;
    major_wavelength = sum(D(2:end,:));
    
    total_magnitude  = (total_magnitude - min(total_magnitude)) ./ range(total_magnitude);

    %Crude conversion to um
    major_wavelength = major_wavelength .* (986 / 2^10 + 739 / 2^10) / 2;
    
    multiColorLine((1:Nframes) / 2 / 60,major_wavelength,total_magnitude,1 - colormap('gray'))
    c = c+1;
end
ax = gca;

ylim([0,ax.YLim(2)])
xlim([0,3.75])
grid on

ylabel('Major Wavelength (\mum)')
xlabel('Time (hours)')
% 
% if (SAVE)
%     saveFigures('SaveAs',[SAVE_FOLDER 'Wavelength'],'Style','none','Formats',{'svg','png','fig'},'Size',FIG_SIZE);
% end

    