clear all
close all

trialrej = {'rej';'norej'};
whichtrials = 1; % 1 or 2 for rej or norej

EXPERIMENT = 'jos_sync';
MAINPATH = 'O:\projects\jos_sync1\sync_walking\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATH = [MAINPATH, 'data\data_R1\R1_ana4_processed\',  trialrej{whichtrials}];
PATHOUT = [MAINPATH, 'data\data_R1\R1_ana4_processed\graphs\',  trialrej{whichtrials}];

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

% Data for graphs
load([PATH, 'ERP.mat'])
load([PATH, 'STD_singletrial1_allstim.mat'])

% data for 5 main research questions
% structured in a table with the following headers: W1E1 W1E2 W2E1 W2E2
% W3E1 W3E2 (may include another 6 columns for S2)
load([PATH, 'Trialnums_analyze.mat'])
load([PATH(1:66), 'Accdata_analyze.mat'])
load([PATH, 'P3_analyze.mat'])
load([PATH, 'STDnoise_analyze.mat'])
load([PATH, 'SNR_analyze2.mat'])
load([PATH, 'prestim_std_erp_analyze.mat'])

% additional data
load([PATH, 'P3TWINDOW_Gave.mat'])
load([PATH, 'EEG_1sub.mat']) % to get eeg infos like EEG.times and EEG.chanlocs



for whichtrials = 1:length(trialrej);;  % 1 or 2 for rej or norej
    
    %whichtrials = 2
    % Data for graphs
    %     STD_st_allstim(whichtrials) = load([PATH,trialrej{whichtrials} ,'STD_singletrial1_allstim.mat'])
    ERP_st1(whichtrials) = load([PATH(1:66),trialrej{whichtrials} ,'ERP_singletrial1.mat'])
    prestim_std_erp_analyze2(whichtrials) = load([PATH(1:66),trialrej{whichtrials} ,'prestim_std_erp_analyze.mat'])
    %     STDnoise_analyze(whichtrials) = load([PATH,trialrej{whichtrials} ,'STDnoise_analyze.mat'])
    SNR_analyze22(whichtrials) = load([PATH(1:66),trialrej{whichtrials} ,'SNR_analyze2.mat'])
    
end

% %Save as csv files for JASP [only needed to do this once]
% W1E1, W1E2, W2E1, W2E2, W3E1, W3E2
% csvwrite([PATH, 'Trialnums_analyze.csv'], Trialnums_analyze)
% csvwrite([PATH, 'Accdata_analyze.csv'], Accdata_analyze)
% csvwrite([PATH, 'P3_analyze.csv'], P3_analyze)
% csvwrite([PATH, 'STDnoise_analyze.csv'], STDnoise_analyze)
% csvwrite([PATH, 'SNR_analyze2.csv'], SNR_analyze2)
% csvwrite([PATH, 'SNR_end_analyze.csv'], SNR_end_analyze)

%% Figure 2a: ERP grand average plots
electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));

figure('Renderer', 'painters', 'Position', [10 10 800 600]);
for i_cond = 1:2
    amp_lim = [-6 10];
    i_type = 1;
    switch i_cond
        case 1
            colour = 'r';
        case 2
            colour = 'm';
    end
    
    %standard and target waveforms
    % ERP = walkconds, frames, electrodes,stim, type, subj
    
    subplot(2,2,i_cond);
    
    f = fill( [   EEG.times(P3TWINDOW_Gave(1));EEG.times(P3TWINDOW_Gave(end));
        EEG.times(P3TWINDOW_Gave(end)); EEG.times(P3TWINDOW_Gave(1))],...
        [ amp_lim(1); amp_lim(1);amp_lim(end); amp_lim(end)], 'y')
    
    set(f,'edgecolor','none');
    uistack(f, 'bottom')
    
    Targetmeans = double(squeeze(mean(ERP(i_cond, :,electrode,1, i_type, :),6)));
    Targetstds = double(squeeze(           std(ERP(i_cond, :,electrode,1, i_type, :),[],6)            )./sqrt(length(SUBJ)));
    Standardmeans = double(squeeze(mean    (ERP  (i_cond, :,electrode,2, i_type, :)   ,6)));
    Standardstds = double(squeeze(std(ERP(i_cond, :,electrode,2, i_type, :),[],6))./sqrt(length(SUBJ)));
    
    b=   boundedline( EEG.times , Targetmeans, Targetstds, colour, EEG.times , Standardmeans , Standardstds , 'k'  );
    
    
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    
    axis tight; ylim(amp_lim);
    
    L1 =     line([-200 800],[0 0],'color','k');
    L2 =     line([0 0],amp_lim,'color','k');
    %   P3TWINDOW_subj(i_sub, :,  i_cond, i_type)
    
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    set(L1,'LineWidth',1.5);
    set(L2,'LineWidth',1.5);
    set(gca, 'Box', 'off')
    
    legend(b,{'Target', 'Standard'},'Location','NorthEast', 'LineWidth', 1);
    % Second electrode type
    %standard and target waveforms
    i_type = 2;
    switch i_cond
        case 1
            colour = 'b';
        case 2
            colour = 'c';
    end
    subplot(2,2,2+i_cond);
    
    f = fill( [   EEG.times(P3TWINDOW_Gave(1));EEG.times(P3TWINDOW_Gave(end));
        EEG.times(P3TWINDOW_Gave(end)); EEG.times(P3TWINDOW_Gave(1))],...
        [ amp_lim(1); amp_lim(1);amp_lim(end); amp_lim(end)], 'y')
    
    set(f,'edgecolor','none');
    uistack(f, 'bottom')
    
    Targetmeans = double(squeeze(mean(ERP(i_cond, :,electrode,1, i_type, :),6)));
    Targetstds = double(squeeze(           std(ERP(i_cond, :,electrode,1, i_type, :),[],6)            )./sqrt(length(SUBJ)));
    Standardmeans = double(squeeze(mean    (ERP  (i_cond, :,electrode,2, i_type, :)   ,6)));
    Standardstds = double(squeeze(std(ERP(i_cond, :,electrode,2, i_type, :),[],6))./sqrt(length(SUBJ)));
    
    b=   boundedline( EEG.times , Targetmeans, Targetstds, colour,  EEG.times , Standardmeans , Standardstds , 'k'  );
    
    
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    set(gca,'Color',[1 1 1]);
    set(gca,'YDir','reverse');
    
    
    axis tight; ylim(amp_lim);
    L1 =     line([-200 800],[0 0],'color','k');
    L2 =     line([0 0],amp_lim,'color','k');
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    set(L1,'LineWidth',1.5);
    set(L2,'LineWidth',1.5);
    set(gca, 'Box', 'off')
    legend(b,{'Target', 'Standard'},'Location','NorthEast', 'LineWidth', 1);
    
end

thetitle = ['All subs mean ERPs at ',EEG.chanlocs(electrode).labels ]
suptitle({thetitle ; ' '});
print([PATHOUT, thetitle], '-dpng')

%% Figure 2b: average topographies

% ERP = walkconds, frames, electrodes,stim, type, subj
%TWINDOW = find(EEG.times>300,1)-1:find(EEG.times>430,1)-2;
for i_stim = 1:2;
    
    figure('Renderer', 'painters', 'Position', [10 10 800 600])
    i_type = 1;
    for i_cond = 1:3;
        subplot(2,3,i_cond)
        for i_sub = 1:length(SUBJ)
            temp(i_sub, i_cond,:, i_type)=  mean(mean(ERP(i_cond, P3TWINDOW_Gave ,:, i_stim, i_type, :), 6),2);
        end
        
        topoplot(mean(temp(i_sub, i_cond,:, i_type),1), EEG.chanlocs, 'maplimits', [-5 5]);
        colorbar;
        title([CONDNAMES{i_cond},' (', ETYPE{i_type}, ')']);
        %title('Targets')
    end
    
    i_type = 2;
    for i_cond = 1:3;
        subplot(2,3,3+i_cond)
        for i_sub = 1:length(SUBJ)
            temp(i_sub, i_cond,:, i_type)=  mean(mean(ERP(i_cond, P3TWINDOW_Gave ,:, i_stim, i_type, :), 6),2);
        end
        topoplot(mean(temp(i_sub, i_cond,:, i_type),1), EEG.chanlocs, 'maplimits', [-5 5]);
        colorbar;
        title([CONDNAMES{i_cond},' (', ETYPE{i_type}, ')']);
        %title('Targets')
        
    end
    %
    %
    thetitle = [ 'All subjects average topo ', STIMULI{i_stim}];
    suptitle(thetitle);
    print([PATHOUT, thetitle], '-dpng')
    
end


%% Figure 3a: Data noise rainclouds: ERP noise & SNR
colours = [8 10 7 11];
cl = [[1 0 0]; [0 0 1];[1 0 1] ;[0 1 1]]

%     W1E1 W1E2 W2E1 W2E2
% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = prestim_std_erp_analyze(:,((2*i_cond-2)+i_type))
    end
end

% make figure
f8  = figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
subplot(2,1,1);
h   = rm_raincloud_2(tempdat, cl, 1);
title(['ERP Data Noise Distribution']);
xlim([-.1 .8])
% ylim([-2 6])
set(gca, 'FontSize', 11)

legendcolours = nan(4,10);
legendcoloursx = 1:10;

%just an invisible plot so I can control my darn legend
p = plot(legendcoloursx,legendcolours(1,:), 'r', legendcoloursx, legendcolours(2,:),  'b', legendcoloursx, legendcolours(3,:), 'm',...
    legendcoloursx, legendcolours(4,:), 'c', 'LineWidth', 4);
lgd = legend(p, {'Standing - Active', 'Standing - Passive', 'Walking - Active', 'Walking - Passive'}, 'Location', 'Northeast', 'LineWidth', 1);
set(gca, 'LineWidth', 2)

% Figure 4b: SNR raincloud plots
colours = [8 10 7 11];
cl = [[1 0 0]; [0 0 1];[1 0 1] ;[0 1 1]]

%     W1E1 W1E2 W2E1 W2E2
% read into cell array of the appropriate dimensions
for i_cond = 1:2
    for i_type = 1:2
        tempdat{i_cond, i_type} = SNR_analyze2(:,((2*i_cond-2)+i_type))
    end
end

% make figure
subplot(2,1,2);

h   = rm_raincloud_2(tempdat, cl, 1);
%set(gca, 'YLim', [2 16]);
title(['SNR distribution']);
xlim([-4 11])
ylim([-.2 .9])
set(gca, 'FontSize', 11)
legendcolours = nan(4,10);
legendcoloursx = 1:10;

set(gca, 'LineWidth', 2)
% save
thetitle = ['ERP and single trial data noise' ]
suptitle({thetitle });
print([PATHOUT, thetitle], '-dpng')


%% Figure 5: VERSION 3 correlation plots

figure('Renderer', 'painters', 'Position', [10 10 950 500]);


subplot(1, 3, 1);
%Shepherd's pi
[stats] = ScatterOutliers_2(SNR_analyze2(:,3), SNR_analyze2(:,4), 1000, false, '-winter',[-2 8],[-2 8]);
title({'Walking SNR'}, 'FontSize', 14)
set(gca, 'LineWidth', 2)
xlabel('Active')
ylabel('Passive')


subplot(1, 3, 2);
%Shepherd's pi
[stats] = ScatterOutliers_2(P3_analyze(:, 1),P3_analyze(:,2) , 1000, false, '-winter',[-5 20],[-1 20]);
title({'Standing P3 Amplitude (\muV)'}, 'FontSize', 14)
set(gca, 'LineWidth', 2)
xlabel('Active')
ylabel('Passive')

subplot(1, 3, 3);
[stats] = ScatterOutliers_2(P3_analyze(:, 3),P3_analyze(:,4), 1000, false, '-winter', [-5 20],[-1 15]);
title({'Walking P3 Amplitude (\muV)'}, 'FontSize', 14)
set(gca, 'LineWidth', 2)
xlabel('Active')
ylabel('Passive')

% %Print to folder
% thetitle = [ 'Correlations ERP data noise Accdata'];
% suptitle(thetitle);
% print([PATHOUT, thetitle], '-dpng')

%% SNR Calculation (chronological)
% electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));
SNR_chron= nan(18,2, 3,100,2);
for i_wt = 1:length(trialrej)
    for i_sub = 1:length(SUBJ)
        
        i_stim = 1; % just targets
        panel = 1;
        
        for i_type = 1:2;
            for i_cond = 1:3;
                
                curr_trials = 1;
                
                
                for i_trials = 1:100
                    
                    ERP_odd =  squeeze(mean(ERP_st1(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, 1:2:curr_trials), 7));
                    ERP_even = squeeze(mean(ERP_st1(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, 2:2:curr_trials), 7));
                    signal = nanmean([ERP_odd ERP_even], 2);
                    noise =    abs(ERP_odd-ERP_even);
                    SNR_chron(i_sub,i_type, i_cond, i_trials, i_wt) = nanmean(signal)/nanmean(noise);
                    
                    if i_trials > 1;
                        if isnan(SNR_chron(i_sub,i_type, i_cond, i_trials, i_wt));
                            if ~isnan(SNR_chron(i_sub,i_type, i_cond, i_trials-1, i_wt));
                                SNR_end(i_sub, i_type, i_cond, i_wt) = SNR_chron(i_sub,i_type, i_cond, i_trials-1, i_wt);
                            end
                        end
                    end
                    curr_trials = curr_trials+1;
                end
                
            end
        end
        
    end
end
%%  SNR comparison graph (chron)
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};
mymap = [ 0.1328    0.5430    0.1328;  0.5977    0.1953    0.7969];
f8  = figure('Renderer', 'painters', 'Position', [10 10 800 600]);

for i_type = 1:2;
    i_cond = 1;
    switch i_type
        case 1
            colour = 'r';
        case 2
            colour = 'b';
    end
    subplot(2,2,i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
    b=   boundedline( [1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour,  'nan', 'gap', 'alpha','transparency', 0.4 );
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNRc_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                    SNRc_end_all(:, i_type, i_cond, 1) = SNR_chron(:,i_type, i_cond, i_snr-1,1);

                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNRc_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                    SNRc_end_all(:, i_type, i_cond, 2) = SNR_chron(:,i_type, i_cond, i_snr-1,2);

                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNRc_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNRc_end( i_type, i_cond, 2), 1))])
    
    %note: non-trial-rejected is first here so that trial rejected is on top
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','SouthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Trials over time');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 12)
    
    i_cond = 2;
    switch i_type
        case 1
            colour = 'm';
        case 2
            colour = 'c';
    end
    subplot(2,2,2+i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :, 1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2 ),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
    b=   boundedline([1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour, 'nan', 'gap', 'alpha','transparency', 0.4  );
    %   b=   boundedline([1:100] , SNR_means1, SNR_std1,[1:100] , SNR_means2, SNR_std2, 'cmap', mymap(1), 'nan', 'gap', 'alpha','transparency', 0.4 );
    
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNRc_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                    SNRc_end_all(:, i_type, i_cond, 1) = SNR_chron(:,i_type, i_cond, i_snr-1,1);

                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNRc_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                    SNRc_end_all(:, i_type, i_cond, 2) = SNR_chron(:,i_type, i_cond, i_snr-1,2);

                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNRc_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNRc_end( i_type, i_cond, 2), 1))])
    
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','NorthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Trials over time');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 12)
    
end

thetitle = [ 'SNR with increasing chronological trial numbers'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')

%SNRncc_end_all(:, i_type, i_cond, wt)
SNR_end_chron_ana = [SNRc_end_all(:, 1, 1, 1) SNRc_end_all(:, 2, 1, 1) SNRc_end_all(:, 1, 2, 1) SNRc_end_all(:, 2, 2, 1)...
                         SNRc_end_all(:, 1, 1, 2) SNRc_end_all(:, 2, 1, 2) SNRc_end_all(:, 1, 2, 2) SNRc_end_all(:, 2, 2, 2)];


save('O:\projects\jos_sync1\sync_walking\data\R1_data\R1_ana4_processed\SNR_end_chron_ana', 'SNR_end_chron_ana')


%% Polyfit calculation for supplemental SNR plot
% i_cond = 1;
i_wt = 2;
i_type = 1;

SNR_means11 = squeeze(mean(SNR_chron(:,i_type, 1, :,i_wt ),1));
SNR_means12 = squeeze(mean(SNR_chron(:,i_type, 2, :,i_wt),1));

i_type = 2;

SNR_means21 = squeeze(mean(SNR_chron(:,i_type, 1, :,i_wt ),1));
SNR_means22 = squeeze(mean(SNR_chron(:,i_type, 2, :,i_wt ),1));

minlength = min([length(SNR_means11(~isnan(SNR_means11))) length(SNR_means21(~isnan(SNR_means21))) ...
    length(SNR_means22(~isnan(SNR_means22))) length(SNR_means12(~isnan(SNR_means12)))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standing
x = SNR_means11(2:minlength+1);
y = SNR_means21(2:minlength+1);

p=polyfit(x, y, 1);
yfit = polyval(p,x);
% makes equation: yfit = p(1) * x[SNR_means...] + p(2)
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1- SSresid/SStotal;

disp(['Active vs. Passive (Standing) R^2 = ', num2str(rsq)])
disp([' yfit = ', num2str(p(1)),'x + ', num2str(p(2))])
% This demonstrates how much of the equation yfit = p(1) * x[SNR_means...]
% + p(2) predicts y (percentage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Walking
x = SNR_means12(2:minlength+1);
y = SNR_means22(2:minlength+1);

p=polyfit(x, y, 1);
yfit = polyval(p,x);
% makes equation: yfit = p(1) * x[SNR_means...] + p(2)
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1- SSresid/SStotal;

disp(['Active vs. Passive (Walking) R^2 = ', num2str(rsq)])
disp([' yfit = ', num2str(p(1)),'x + ', num2str(p(2))])

% This demonstrates how much of the equation yfit = p(1) * x[SNR_means...]
% + p(2) predicts y (percentage)
%%  SUPPLEMENTAL/REVIEWER ONLY: SNR rej comparison graph
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};
f8  = figure('Renderer', 'painters', 'Position', [10 10 800 500]);

subplot(1,2,1)
%Standing
i_cond = 1;
i_type =1;
i_wt = 2;

switch i_wt
    case 1;
        mytitle = ' (trial-rejected)';
    case 2;
        mytitle = ' (non trial-rejected)';
end
%     subplot(2,2,2+i_type)

SNR_means20 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,i_wt),1));
SNR_std20 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, i_wt)))./sqrt(length(SUBJ));

i_type = 2;

SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,i_wt),1));
SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, i_wt)))./sqrt(length(SUBJ));
b=   boundedline([1:100] , SNR_means20, SNR_std20, 'r',[1:100]  , SNR_means2, SNR_std2, 'b','nan', 'gap', 'alpha','transparency', 0.4  );
legend({'Active','Passive'}, 'Location','NorthEast')

set(b,'LineWidth',2);
set(gca,'LineWidth',2);
ylim([0 4])
title([ CONDNAMES{i_cond}, mytitle] );
xlabel('Trials over time');
ylabel('Signal-to-noise ratio');
set(gca, 'FontSize', 12)

subplot(1,2,2)
% Walking
i_cond = 2;
i_type =1;
%     subplot(2,2,2+i_t
SNR_means20 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,i_wt),1));
SNR_std20 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, i_wt)))./sqrt(length(SUBJ));
i_type = 2;

SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,i_wt ),1));
SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, i_wt)))./sqrt(length(SUBJ));
b=   boundedline([1:100] , SNR_means20, SNR_std20, 'm',[1:100]  , SNR_means2, SNR_std2, 'c','nan', 'gap', 'alpha','transparency', 0.4  );
legend({'Active','Passive'}, 'Location','NorthEast')

set(b,'LineWidth',2);
set(gca,'LineWidth',2);
ylim([0 4])
title([ CONDNAMES{i_cond}, mytitle ] );
xlabel('Trials over time');
ylabel('Signal-to-noise ratio');
set(gca, 'FontSize', 12)

thetitle = [ 'SNR with increasing chronological trial numbers'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')

%% SNR Calculation Permutation (non chron)
% electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));

SNR_nonchron= nan(18,2, 3,100,2);
ERP_stl_perms = ERP_st1;

for i_wt = 1:length(trialrej)
    for i_sub = 1:length(SUBJ)
        
        i_stim = 1; % just targets
        
        for i_type = 1:2;
            for i_cond = 1:3;
                for i_perm = 1:100;
                    
                    % shuffle: just shuffle the ones that are not nans
                    ERP_calc = ERP_stl_perms(i_wt).ERP_singletrial1(i_sub,1, 1,i_cond,i_type, 1, 1:100);
                    nonnan =   1:length(find(~isnan(ERP_calc)));
                    shuf_nonnan = shuffle(nonnan);
                    ynan = [length(shuf_nonnan)+1:600];
                    neworder = [shuf_nonnan ynan];
                    
                    %TODO fix this. the indices isn't working and the left side
                    %is mostly zeros. might need a specific function for this.
                    
                    ERP_stl_perms(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, :)= ERP_stl_perms(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, neworder);
                    curr_trials = 1;
                    for i_trials = 1:100
                        ERP_odd =  squeeze(mean(ERP_stl_perms(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, 1:2:curr_trials), 7));
                        ERP_even = squeeze(mean(ERP_stl_perms(i_wt).ERP_singletrial1(i_sub,1, P3TWINDOW_Gave,i_cond,i_type, i_stim, 2:2:curr_trials), 7));
                        signal = nanmean([ERP_odd ERP_even], 2);
                        noise =    abs(ERP_odd-ERP_even);
                        SNR_nonchron(i_sub,i_type, i_cond,i_perm, i_trials, i_wt) = nanmean(signal)/nanmean(noise);
                        
                        if i_trials > 1;
                            if isnan(SNR_nonchron(i_sub,i_type, i_cond,i_perm, i_trials, i_wt));
                                if ~isnan(SNR_nonchron(i_sub,i_type, i_cond,i_perm, i_trials-1, i_wt));
                                    SNRnc_end(i_sub, i_type, i_cond, i_perm, i_wt) = SNR_nonchron(i_sub, i_type, i_cond, i_perm, i_trials-1, i_wt);
                                end
                            end
                        end
                        curr_trials = curr_trials+1;
                    end
                end
            end
        end
        
    end
end



%% (version 2): SNR comparison graph (NONCHRON)
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};
mymap = [ 0.1328    0.5430    0.1328;  0.5977    0.1953    0.7969];
f8  = figure('Renderer', 'painters', 'Position', [10 10 800 600]);

for i_type = 1:2;
    i_cond = 1;
    switch i_type
        case 1
            colour = 'r';
        case 2
            colour = 'b';
    end
    subplot(2,2,i_type)
    SNR_means1 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,1),4),1));
    SNR_std1 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond,:, :, 1),4)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,2),4),1));
    SNR_std2 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, 2),4)))./sqrt(length(SUBJ));
    b=   boundedline( [1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour,  'nan', 'gap', 'alpha','transparency', 0.4 );
    
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNRncc_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                    SNRncc_end_all(:, i_type, i_cond, 1) = nanmean(SNR_nonchron(:,i_type, i_cond,:, i_snr-1,1),4);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNRncc_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                    SNRncc_end_all(:, i_type, i_cond, 2) = nanmean(SNR_nonchron(:,i_type, i_cond,:, i_snr-1,2),4);

                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNRncc_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNRncc_end( i_type, i_cond, 2), 1))])
    
    %note: non-trial-rejected is first here so that trial rejected is on top
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','NorthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Trials over time');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 12)
    
    i_cond = 2;
    switch i_type
        case 1
            colour = 'm';
        case 2
            colour = 'c';
    end
    subplot(2,2,2+i_type)
    SNR_means1 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,1),4),1));
    SNR_std1 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond,:, :, 1),4)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,2),4),1));
    SNR_std2 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, 2),4)))./sqrt(length(SUBJ));
    b=   boundedline([1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour, 'nan', 'gap', 'alpha','transparency', 0.4  );
    %   b=   boundedline([1:100] , SNR_means1, SNR_std1,[1:100] , SNR_means2, SNR_std2, 'cmap', mymap(1), 'nan', 'gap', 'alpha','transparency', 0.4 );
    
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNRncc_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                    SNRncc_end_all(:, i_type, i_cond, 1) = nanmean(SNR_nonchron(:,i_type, i_cond,:, i_snr-1,1),4);

                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNRncc_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                    SNRncc_end_all(:, i_type, i_cond, 2) = nanmean(SNR_nonchron(:,i_type, i_cond,:, i_snr-1,2),4);

                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNRncc_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNRncc_end( i_type, i_cond, 2), 1))])
    
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','NorthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Trials over time');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 12)
    
end

thetitle = [ 'SNR with increasing trial numbers (100 perms)'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')

%SNRncc_end_all(:, i_type, i_cond, wt)
SNR_end_nonchron_ana = [SNRncc_end_all(:, 1, 1, 1) SNRncc_end_all(:, 2, 1, 1) SNRncc_end_all(:, 1, 2, 1) SNRncc_end_all(:, 2, 2, 1)...
                         SNRncc_end_all(:, 1, 1, 2) SNRncc_end_all(:, 2, 1, 2) SNRncc_end_all(:, 1, 2, 2) SNRncc_end_all(:, 2, 2, 2)];


save('O:\projects\jos_sync1\sync_walking\data\R1_data\R1_ana4_processed\SNR_end_nonchron_ana.mat', 'SNR_end_nonchron_ana')

%%  SUPPLEMENTAL/REVIEWER ONLY: SNR rej comparison graph
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};
f8  = figure('Renderer', 'painters', 'Position', [10 10 800 500]);

subplot(1,2,1)
%Standing
i_cond = 1;
i_type =1;
i_wt = 2;

switch i_wt
    case 1;
        mytitle = ' (trial-rejected)';
    case 2;
        mytitle = ' (non trial-rejected)';
end
%     subplot(2,2,2+i_type)

SNR_means20 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,i_wt),4),1));
SNR_std20 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, i_wt),4)))./sqrt(length(SUBJ));

i_type = 2;

SNR_means2 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,i_wt),4),1));
SNR_std2 =     squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, i_wt),4)))./sqrt(length(SUBJ));
b=   boundedline([1:100] , SNR_means20, SNR_std20, 'r',[1:100]  , SNR_means2, SNR_std2, 'b','nan', 'gap', 'alpha','transparency', 0.4  );
legend({'Active','Passive'}, 'Location','NorthEast')

set(b,'LineWidth',2);
set(gca,'LineWidth',2);
ylim([0 4])
title([ CONDNAMES{i_cond}, mytitle] );
xlabel('Number of trials');
ylabel('Signal-to-noise ratio');
set(gca, 'FontSize', 12)

subplot(1,2,2)
% Walking
i_cond = 2;
i_type =1;
%     subplot(2,2,2+i_t
SNR_means20 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,i_wt),4),1));
SNR_std20 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, i_wt),4)))./sqrt(length(SUBJ));

i_type = 2;

SNR_means2 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,i_wt),4),1));
SNR_std2 =     squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, i_wt),4)))./sqrt(length(SUBJ));
b=   boundedline([1:100] , SNR_means20, SNR_std20, 'm',[1:100]  , SNR_means2, SNR_std2, 'c','nan', 'gap', 'alpha','transparency', 0.4  );
legend({'Active','Passive'}, 'Location','NorthEast')

set(b,'LineWidth',2);
set(gca,'LineWidth',2);
ylim([0 4])
title([ CONDNAMES{i_cond}, mytitle ] );
xlabel('Number of trials');
ylabel('Signal-to-noise ratio');
set(gca, 'FontSize', 12)

thetitle = [ 'SNR with increasing trial numbers (100 perms)'];
suptitle(thetitle);
print([PATHOUT, thetitle], '-dpng')

%% Polyfit calculation for NONCHRON supplemental SNR plot
% i_cond = 1;
i_wt = 1;
i_type = 1;

SNR_means11 = squeeze(mean(mean(SNR_nonchron(:,i_type, 1,:, :,i_wt ),4),1));
SNR_means12 = squeeze(mean(mean(SNR_nonchron(:,i_type, 2,:, :,i_wt ),4),1));

i_type = 2;

SNR_means21 = squeeze(mean(mean(SNR_nonchron(:,i_type, 1,:, :,i_wt ),4),1));
SNR_means22 = squeeze(mean(mean(SNR_nonchron(:,i_type, 2,:, :,i_wt ),4),1));

minlength = min([length(SNR_means11(~isnan(SNR_means11))) length(SNR_means21(~isnan(SNR_means21))) ...
    length(SNR_means22(~isnan(SNR_means22))) length(SNR_means12(~isnan(SNR_means12)))]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standing
x = SNR_means11(2:minlength+1);
y = SNR_means21(2:minlength+1);

p=polyfit(x, y, 1);
yfit = polyval(p,x);
% makes equation: yfit = p(1) * x[SNR_means...] + p(2)
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1- SSresid/SStotal;

disp(['Active vs. Passive (Standing) R^2 = ', num2str(rsq)])
disp([' yfit = ', num2str(p(1)),'x + ', num2str(p(2))])
% This demonstrates how much of the equation yfit = p(1) * x[SNR_means...]
% + p(2) predicts y (percentage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Walking
x = SNR_means12(2:minlength+1);
y = SNR_means22(2:minlength+1);

p=polyfit(x, y, 1);
yfit = polyval(p,x);
% makes equation: yfit = p(1) * x[SNR_means...] + p(2)
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1- SSresid/SStotal;

disp(['Active vs. Passive (Walking) R^2 = ', num2str(rsq)])
disp([' yfit = ', num2str(p(1)),'x + ', num2str(p(2))])

% This demonstrates how much of the equation yfit = p(1) * x[SNR_means...]
% + p(2) predicts y (percentage)

%% Figure 4: Combined chron + nonchron graphs
% A:  Figure 5: SNR comparison graph (chron)
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};
mymap = [ 0.1328    0.5430    0.1328;  0.5977    0.1953    0.7969];
f8  = figure('Fi', 'painters', 'Position', [10 10 700 500]);

for i_type = 1:2;
    i_cond = 1;
    switch i_type
        case 1
            colour = 'r';
        case 2
            colour = 'b';
    end
    subplot(2,2,i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
    b=   boundedline( [1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour,  'nan', 'gap', 'alpha','transparency', 0.4 );
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNR_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNR_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNR_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNR_end( i_type, i_cond, 2), 1))])
    
    %note: non-trial-rejected is first here so that trial rejected is on top
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','SouthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Trials over time');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 10)
    
    i_cond = 2;
    switch i_type
        case 1
            colour = 'm';
        case 2
            colour = 'c';
    end
    subplot(2,2,2+i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :, 1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2 ),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
    b=   boundedline([1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour, 'nan', 'gap', 'alpha','transparency', 0.4  );
    %   b=   boundedline([1:100] , SNR_means1, SNR_std1,[1:100] , SNR_means2, SNR_std2, 'cmap', mymap(1), 'nan', 'gap', 'alpha','transparency', 0.4 );
    
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNR_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNR_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNR_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNR_end( i_type, i_cond, 2), 1))])
    
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','NorthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Trials over time');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 10)
    
end
thetitle = [ 'SNR with increasing chronological trial numbers'];
suptitle(thetitle);
% B: Figure 5 (version 2): SNR comparison graph
% WALKCONDS = {'standing', 'walkingA', 'walkingT'};
% ETYPE = {'Active', 'Passive'};
mymap = [ 0.1328    0.5430    0.1328;  0.5977    0.1953    0.7969];
f8  = figure('Renderer', 'painters', 'Position', [10 10 700 500]);

for i_type = 1:2;
    i_cond = 1;
    switch i_type
        case 1
            colour = 'r';
        case 2
            colour = 'b';
    end
    subplot(2,2,i_type)
    SNR_means1 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,1),4),1));
    SNR_std1 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond,:, :, 1),4)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,2),4),1));
    SNR_std2 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, 2),4)))./sqrt(length(SUBJ));
    b=   boundedline( [1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour,  'nan', 'gap', 'alpha','transparency', 0.4 );
    
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNR_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNR_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNR_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNR_end( i_type, i_cond, 2), 1))])
    
    %note: non-trial-rejected is first here so that trial rejected is on top
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','SouthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Number of trials');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 10)
    
    i_cond = 2;
    switch i_type
        case 1
            colour = 'm';
        case 2
            colour = 'c';
    end
    subplot(2,2,2+i_type)
    SNR_means1 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,1),4),1));
    SNR_std1 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond,:, :, 1),4)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(mean(SNR_nonchron(:,i_type, i_cond,:, :,2),4),1));
    SNR_std2 =      squeeze(std(mean(SNR_nonchron(:,i_type, i_cond, :,:, 2),4)))./sqrt(length(SUBJ));
    b=   boundedline([1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour, 'nan', 'gap', 'alpha','transparency', 0.4  );
    %   b=   boundedline([1:100] , SNR_means1, SNR_std1,[1:100] , SNR_means2, SNR_std2, 'cmap', mymap(1), 'nan', 'gap', 'alpha','transparency', 0.4 );
    
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNR_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNR_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNR_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNR_end( i_type, i_cond, 2), 1))])
    
    legend({'Non Trial-rejected','Trial-rejected'}, 'Location','NorthEast')
    
    set(b,'LineWidth',2);
    set(gca,'LineWidth',2);
    ylim([0 4])
    title([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}] );
    xlabel('Number of trials');
    ylabel('Signal-to-noise ratio');
    set(gca, 'FontSize', 10)
    
end


thetitle = [ 'SNR with increasing randomized trial numbers'];
suptitle(thetitle);
%print([PATHOUT, thetitle], '-dpng')

