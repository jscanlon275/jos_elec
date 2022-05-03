clear all
close all


EXPERIMENT = 'jos_sync';
MAINPATH = 'O:\projects\jos_sync1\sync_walking\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATH = [MAINPATH, 'data\R1_ana4_processed\'];
PATHOUT = [PATH, 'graphs\Noise Analysis\'];

SUBJ = { '002', '004', '005', '006', '007', '009', '012', '013', '014',...
    '015', '018', '020', '021', '022', '023', '024', '025', '026' };
allCONDS = {'standard_standing','deviant_standing',...
    'standard_walkingA','deviant_walkingA',...
    'standard_walkingT','deviant_walkingT'};
WALKCONDS = {'standing', 'walkingA', 'walkingT'};
STIMULI = { 'deviant', 'standard'};
ETYPE = {'Active', 'Passive'};
CONDNAMES = {'Standing'; 'Walking'; 'Walking Together' };

cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; ...
    1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4; .7 0 0 ; 0 0 .7];

cd(MAINPATH);

trialrej = {'rej';'norej'};
for whichtrials = 1:length(trialrej);;  % 1 or 2 for rej or norej
    
    %whichtrials = 2
    % Data for graphs
    ERPrej(whichtrials) = load([PATH,trialrej{whichtrials}, 'ERP.mat'])
    STD_st_allstim(whichtrials) = load([PATH,trialrej{whichtrials} ,'STD_singletrial1_allstim.mat'])
    ERP_st1(whichtrials) = load([PATH,trialrej{whichtrials} ,'ERP_singletrial1.mat'])
    prestim_std_erp_analyze(whichtrials) = load([PATH,trialrej{whichtrials} ,'prestim_std_erp_analyze.mat'])
    STDnoise_analyze(whichtrials) = load([PATH,trialrej{whichtrials} ,'STDnoise_analyze.mat'])
    SNR_analyze2(whichtrials) = load([PATH,trialrej{whichtrials} ,'SNR_analyze2.mat'])
    
end
% additional data
load([PATH, trialrej{1},'P3TWINDOW_Gave.mat'])
load([PATH,trialrej{1}, 'EEG_1sub.mat']) % to get eeg infos like EEG.times and EEG.chanlocs


%%
%% SNR Calculation
% electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));
SNR_chron= nan(18,2, 3,100,2);
for i_wt = 1:length(trialrej)
    for i_sub = 1:length(SUBJ)
        
        i_stim = 1; % just targets
        panel = 1;
        %                figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
        
        for i_type = 1:2;
            for i_cond = 1:3;
                
                curr_trials = 1;
                for i_trials = 1:100
                    subplot(2,3,panel)
                    ERP_odd =  squeeze(mean(ERP_st1(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, 1:2:curr_trials), 7));
                    ERP_even = squeeze(mean(ERP_st1(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, 2:2:curr_trials), 7));
                    signal = nanmean([ERP_odd ERP_even], 2);
                    noise =    abs(ERP_odd-ERP_even);
                    SNR_chron(i_sub,i_type, i_cond, i_trials, i_wt) = mean(signal)/mean(noise);
                    
                    curr_trials = curr_trials+1;
                end
                
            end
        end
      
    end
end
%%   SNR comparison graph
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};

f8  = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);

for i_type = 1:2;
    i_cond = 1;
    subplot(2,2,i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
    b=   boundedline([1:100] , SNR_means1, SNR_std1, 'b',[1:100] , SNR_means2, SNR_std2, 'r' );
    legend({'Trial Rejected','Non Rejected'})
    
    ylim([0 4])
    title([WALKCONDS{i_cond}, ' ' ETYPE{i_type}, ' '])
    
    i_cond = 2;
    subplot(2,2,2+i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :, 1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2 ),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
    b=   boundedline([1:100] , SNR_means1, SNR_std1, 'b',[1:100] , SNR_means2, SNR_std2, 'r' );
    legend({'Trial Rejected','Non Rejected'})
    
    ylim([0 4])
    title([WALKCONDS{i_cond}, ' ' ETYPE{i_type}, ' '])
    
end

thetitle = [ 'SNR with increasing chronological trial numbers'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')

%%  SNR raincloud plots
colours = [8 10 7 11];
%cl = [cb(13,:); cb(8,:); cb(7,:); cb(14,:)]
cl = [[1 0 0]; [0 0 1];[1 0 1] ;[0 1 1]]

%     W1E1 W1E2 W2E1 W2E2

% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = SNR_analyze2(1).SNR_analyze2(:,((2*i_cond-2)+i_type))
    end
end

% make figure
f8  = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
subplot(2,1,1);
%xlabel('Standard deviation (\muV)')
xlim([-4 12])
%ylabel('Density')
h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['Data Noise Distribution (Rejected)']);
% set(gca, 'YTickLabel', {' ', '0', '0.2', '0.4', '0', '0.2', '0.4' })% opposite than usual I guess because the axis is transposed??? anyway it works.
set(gca, 'YTickLabel', {' ', 'Walking', ' ', ' ', 'Standing', ' ', ' ' })% opposite than usual I guess because the axis is transposed??? anyway it works.

% set(gca, 'XTick', [])
% set(gca, 'YTick', [])

legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)
% save

% Figure 2: SNR raincloud plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = SNR_analyze2(2).SNR_analyze2(:,((2*i_cond-2)+i_type))
    end
end

% make figure
subplot(2,1,2);
%xlabel('Standard deviation (\muV)')
xlim([-4 12])
%ylabel('Density')
h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['Data Noise Distribution (Non-Rejected)']);
% set(gca, 'YTickLabel', {' ', '0', '0.2', '0.4', '0', '0.2', '0.4' })% opposite than usual I guess because the axis is transposed??? anyway it works.
set(gca, 'YTickLabel', {' ', 'Walking', ' ', ' ', 'Standing', ' ', ' ' })% opposite than usual I guess because the axis is transposed??? anyway it works.


legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)
% save


thetitle = [ 'SNR distribution before and after trial rejection'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')


%% ERP Data noise raincloud plots
colours = [8 10 7 11];
%cl = [cb(13,:); cb(8,:); cb(7,:); cb(14,:)]
cl = [[1 0 0]; [0 0 1];[1 0 1] ;[0 1 1]]

%     W1E1 W1E2 W2E1 W2E2

% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = prestim_std_erp_analyze(1).prestim_std_erp_analyze(:,((2*i_cond-2)+i_type))
    end
end

% make figure
f8  = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
subplot(2,1,1);
%xlabel('Standard deviation (\muV)')
xlim([0 .8])
%ylabel('Density')
h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['Data Noise Distribution (Rejected)']);
% set(gca, 'YTickLabel', {' ', '0', '0.2', '0.4', '0', '0.2', '0.4' })% opposite than usual I guess because the axis is transposed??? anyway it works.
set(gca, 'YTickLabel', {' ', 'Walking', ' ', ' ', 'Standing', ' ', ' ' })% opposite than usual I guess because the axis is transposed??? anyway it works.

% set(gca, 'XTick', [])
% set(gca, 'YTick', [])

legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)
% save

% Figure 2: SNR raincloud plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = prestim_std_erp_analyze(2).prestim_std_erp_analyze(:,((2*i_cond-2)+i_type))
    end
end

% make figure
subplot(2,1,2);
%xlabel('Standard deviation (\muV)')
xlim([0 .8])
%ylabel('Density')
h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['Data Noise Distribution (Non-Rejected)']);
% set(gca, 'YTickLabel', {' ', '0', '0.2', '0.4', '0', '0.2', '0.4' })% opposite than usual I guess because the axis is transposed??? anyway it works.
set(gca, 'YTickLabel', {' ', 'Walking', ' ', ' ', 'Standing', ' ', ' ' })% opposite than usual I guess because the axis is transposed??? anyway it works.

legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)
% save


thetitle = [ 'Prestimulus ERP data noise before and after trial rejection'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')

%% SINGLE TRIAL Data noise raincloud plots
colours = [8 10 7 11];
%cl = [cb(13,:); cb(8,:); cb(7,:); cb(14,:)]
cl = [[1 0 0]; [0 0 1];[1 0 1] ;[0 1 1]]

%     W1E1 W1E2 W2E1 W2E2

% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = STDnoise_analyze(1).STDnoise_analyze(:,((2*i_cond-2)+i_type))
    end
end

% make figure
f8  = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
subplot(2,1,1);
%xlabel('Standard deviation (\muV)')
xlim([2 16])
%ylabel('Density')
h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['Data Noise Distribution (Rejected)']);
% set(gca, 'YTickLabel', {' ', '0', '0.2', '0.4', '0', '0.2', '0.4' })% opposite than usual I guess because the axis is transposed??? anyway it works.
set(gca, 'YTickLabel', {' ', 'Walking', ' ', ' ', 'Standing', ' ', ' ' })% opposite than usual I guess because the axis is transposed??? anyway it works.

% set(gca, 'XTick', [])
% set(gca, 'YTick', [])

legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)
% save

% Figure 2: SNR raincloud plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = STDnoise_analyze(2).STDnoise_analyze(:,((2*i_cond-2)+i_type))
    end
end

% make figure
subplot(2,1,2);
%xlabel('Standard deviation (\muV)')
xlim([2 16])
%ylabel('Density')
h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['Data Noise Distribution (Non-Rejected)']);
% set(gca, 'YTickLabel', {' ', '0', '0.2', '0.4', '0', '0.2', '0.4' })% opposite than usual I guess because the axis is transposed??? anyway it works.
set(gca, 'YTickLabel', {' ', 'Walking', ' ', ' ', 'Standing', ' ', ' ' })% opposite than usual I guess because the axis is transposed??? anyway it works.

% set(gca, 'XTick', [])
% set(gca, 'YTick', [])

legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)
% save


thetitle = [ 'Prestimulus single-trial data noise before and after trial rejection'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')
