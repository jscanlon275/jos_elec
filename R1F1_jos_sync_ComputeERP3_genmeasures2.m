clear all
close all

trialrej = {'rej';'norej'};
whichtrials = 1; % 1 or 2 for rej or norej

EXPERIMENT = 'jos_sync';
MAINPATH = 'D:\4. home office\otto backup\jos_sync1\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATH = [MAINPATH, 'data\R1_ana3_epochs_', trialrej{whichtrials}, '\'];
PATHOUT = [MAINPATH, 'data\R1_ana4_processed\',  trialrej{whichtrials}];

SUBJ = { '002'...
    , '004', '005', '006', '007', '009', '012', '013', '014',...
    '015', '018', '020', '021', '022', '023', '024', '025', '026'...
    };
allCONDS = {'standard_standing','deviant_standing',...
    'standard_walkingA','deviant_walkingA',...
    'standard_walkingT','deviant_walkingT'};
WALKCONDS = {'standing', 'walkingA', 'walkingT'};
STIMULI = { 'deviant', 'standard'};
ETYPE = {'act', 'pas'};
CONDNAMES = {'Standing'; 'Walking Alone'; 'Walking Together' };

cb = [0.5 0.8 0.9; 1 1 0.7; 0.7 0.8 0.9; 0.8 0.5 0.4; 0.5 0.7 0.8; 1 0.8 0.5; 0.7 1 0.4; 1 0.7 1; 0.6 0.6 0.6; 0.7 0.5 0.7; 0.8 0.9 0.8; 1 1 0.4];

cd(MAINPATH);

% ERP = walkconds, frames, electrodes,stim, type, subj
ERP = nan(length(WALKCONDS),250, 62,2,2, length(SUBJ));
NOISE = nan(length(WALKCONDS),62, 1,2,2, length(SUBJ));
tcounter = nan(18, 2, 3, 2);


%ERP_singletrial1 = nan(18,62,250,3,2,2,505); % too big: maybe not needed?
%maybe only one electrode, etc
ERP_singletrial1 = nan(18,1,250,3,2,2,600);
STD_singletrial1 = nan(18,3,2,2,600);
ERP_prestim1 = nan(18, 50,3,2, 2, 600);

STD_singletrial1_allstim = nan(18,3,2, 1200);
ERP_prestim2 = nan(18,50,3,2,1200);
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%% Load data
for i_sub = 1:length(SUBJ)
    for i_type = 1:length(ETYPE)
        for i_cond = 1:length(WALKCONDS)
            for i_stim = 1:length(STIMULI)
                clear data2
                
                EEG = pop_loadset('filename', [PATH,'jos_sync', SUBJ{i_sub}, '_', ETYPE{i_type}, '_oddball_epoched_', trialrej{whichtrials},'_',STIMULI{i_stim},'_' WALKCONDS{i_cond}, '.set']);
                
                %  ERP = walkconds, frames, electrodes,stim, type, subj
                ERP(i_cond,:,:, i_stim, i_type,i_sub) = mean(EEG.data ,3)'; %targets
                
                n_trials = EEG.trials;
                
                % now let's count trials
                tcounter(i_sub, i_type, i_cond, i_stim) = n_trials;
                
                for i_trial = 1:n_trials
                    
                    % single trial noise: for each trial, RMS value for each electrode, then averaged over electrodes, then average over  trials
                    % (for each participant), then over subjects
                    % here I want 1 value for each electrode, then to average over all electrodes
                    ERP_singletrial1(i_sub,:, :,i_cond,i_type, i_stim, i_trial) = ( EEG.data(12, :,i_trial) ); %select Pz
                    
                    % NEW revision 1: take only electrode Pz, then average
                    % over ALL trials, then STD of the timepoints for alltrials prestim ERP.
                    %ERP_prestim1(i_sub, frames,i_cond,i_type, i_stim, i_trial)
                    ERP_prestim1(i_sub, :,i_cond,i_type, i_stim, i_trial) = ( EEG.data(12, EEG.times<0,i_trial)); %select Pz
                    
                    data2= permute(EEG.data, [2, 1, 3]); % EEG.data(frames, elec, trials)
                    STD_singletrial1(i_sub, i_cond, i_type, i_stim, i_trial) = mean(std(data2(EEG.times<0, :, i_trial)));
                    
                end
                
            end
            
            % remember the size here should be preallocated to be big enough for both standards and targets
            STD_singletrial1_allstim(i_sub,i_cond,i_type, :) = [squeeze(STD_singletrial1(i_sub,i_cond,i_type, 1, :)); squeeze(STD_singletrial1(i_sub,i_cond,i_type, 2, :))];
            % ERP noise:  (i_sub, frames,i_cond,i_type, i_stim, i_trial)
            % concatonate all trials (targets and standards)
            ERP_prestim2(i_sub, :,i_cond,i_type, :) = [squeeze(ERP_prestim1(i_sub, :,i_cond,i_type, 1, :)) squeeze(ERP_prestim1(i_sub, :,i_cond,i_type, 2, :))];
        end
        chremove(i_sub, i_type) = EEG.rmChan;
        chansremov{i_sub,i_type} = {EEG.rmChanloc};
    end
end
% Make the ERP difference wave
erp_diff_out = squeeze(ERP(:,:,:,1,:,:)-ERP(:,:,:,2,:,:));


% ERP_prestim1(i_sub,elec, frames,i_cond,i_type, i_trial)

%STD_ERP_prestim(i_sub,frames, cond,type)
STD_ERP_prestim = nanmean(ERP_prestim2(:,:,:,:, :),5); % average over trials
STD_ERP_prestim2 = permute(STD_ERP_prestim, [2, 1,3,4]);% (frames,subs, cond,type)
STD_ERP_prestim3 = squeeze(std(STD_ERP_prestim2)); %averaged all trials, then std of frames (only Pz)
% STD_ERP_prestim3(subs, conds, type)


%%

for whichtrials = 1:length(trialrej);;  % 1 or 2 for rej or norej
    
    %whichtrials = 2
    % Data for graphs
    %     STD_st_allstim(whichtrials) = load([PATH,trialrej{whichtrials} ,'STD_singletrial1_allstim.mat'])
    ERP_st1(whichtrials) = load([PATH(1:63),trialrej{whichtrials} ,'ERP_singletrial1.mat'])
    prestim_std_erp_analyze2(whichtrials) = load([PATH(1:63),trialrej{whichtrials} ,'prestim_std_erp_analyze.mat'])
    %     STDnoise_analyze(whichtrials) = load([PATH,trialrej{whichtrials} ,'STDnoise_analyze.mat'])
    SNR_analyze22(whichtrials) = load([PATH(1:63),trialrej{whichtrials} ,'SNR_analyze2.mat'])
    
end
%% Save (don't need to be done every time)
save([PATHOUT, 'ERP'], 'ERP')
save([PATHOUT,'STD_singletrial1_allstim'], 'STD_singletrial1_allstim')
save([PATHOUT, 'EEG_1sub'], 'EEG')
save([PATHOUT, 'ERP_singletrial1'], 'ERP_singletrial1')
save([PATHOUT, 'STD_ERP_prestim3'], 'STD_ERP_prestim3')



eegtimes = EEG.times;
eegchanlocs = EEG.chanlocs;
save([PATHOUT(1:63), 'times'], 'eegtimes');
save([PATHOUT(1:63), 'chanlocs'], 'eegchanlocs');
save([PATHOUT(1:63), 'EEG_1sub'], 'EEG')


eeglab redraw;
%% Q0 removed channels
% chremove;
% sum(ch_remov);

%% Q1 electrode trial num comparisons
electype = 1;
stimu = 1;
walkcond = 2;

%
%act vs. pas; Condition 1, deviants
[h, p, ci, stats] = ttest(tcounter(:, 1, walkcond, stimu),tcounter(:, 2, walkcond, stimu) );
MeanAct = mean(tcounter(:, 1, walkcond, stimu));
MeanPas = mean(tcounter(:, 2, walkcond, stimu));
disp([WALKCONDS{walkcond}, '_', STIMULI{stimu}])
disp(['Compare means: Act: ', num2str(MeanAct), ' Pas: ', num2str(MeanPas)])
disp([WALKCONDS{walkcond}, '_', STIMULI{stimu}, ' p = ', num2str(p)])
% Trial numbers anova
% tcounter(subs, type, walkcond, stim)

% set up trial numbers for SPSS ( for reference, this is what I saved)

Trialnums_analyze = [tcounter(:, :, 1, 1) tcounter(:, :, 2, 1) tcounter(:, :, 3, 1)...
    tcounter(:, :, 1, 2) tcounter(:, :, 2, 2) tcounter(:, :, 3, 2)]

save([PATHOUT, 'Trialnums_analyze'], 'Trialnums_analyze');

%% Q3 P3 peaks
electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));
TWINDOW = find(EEG.times>300,1)-1:find(EEG.times>500,1)-2;
clear peaktime peakframe P3TWINDOW_subj

P3TWINDOW_subj = nan(18, 51, 3,2);

for i_sub = 1:length(SUBJ);
    for  i_cond = 1:length(WALKCONDS);
        for i_type = 1:length(ETYPE);
            Targetmeans = double(squeeze(mean(ERP(i_cond, TWINDOW,electrode,1, i_type, i_sub),6)));
            %     Standardmeans = double(squeeze(mean    (ERP  (i_cond, :,electrode,2, i_type, i_sub)   ,6)));
            
            peakframe(i_sub, i_cond, i_type) = 127+find(Targetmeans == max(Targetmeans));
            peaktime(i_sub, i_cond, i_type) = EEG.times(peakframe(i_sub, i_cond, i_type));
            peakamp(i_sub, i_cond, i_type) = max(Targetmeans); % amplitude without the time window
            
            %  TWINDOW_subj(i_sub, :) = find(EEG.times>(peaktime(i_sub)-100),1)-1:find(EEG.times>(peaktime(i_sub)+200),1)-2;
            P3TWINDOW_subj(i_sub, :,  i_cond, i_type) = find(EEG.times>(EEG.times(peakframe(i_sub)-25)),1)-1:find(EEG.times>(EEG.times(peakframe(i_sub)+25)),1)-1;
        end
    end
end
% ERP = walkconds, frames, electrodes,stim, electype, subj
for i_sub = 1:length(SUBJ)
    P3_analyze1(i_sub,:) = [squeeze(mean(ERP(1, P3TWINDOW_subj(i_sub, :,  1, 1), electrode, 1, 1, i_sub),2)), squeeze(mean(ERP(1, P3TWINDOW_subj(i_sub, :,  1, 2), electrode, 1, 2, i_sub),2)) ...
        squeeze(mean(ERP(2, P3TWINDOW_subj(i_sub, :,  2, 1), electrode, 1, 1, i_sub),2)), squeeze(mean(ERP(2, P3TWINDOW_subj(i_sub, :,  2, 2), electrode, 1, 2, i_sub),2)) ...
        squeeze(mean(ERP(3, P3TWINDOW_subj(i_sub, :,  3, 1), electrode, 1, 1, i_sub),2)), squeeze(mean(ERP(3, P3TWINDOW_subj(i_sub, :,  3, 2), electrode, 1, 2, i_sub),2)) ...
        squeeze(mean(ERP(1, P3TWINDOW_subj(i_sub, :,  1, 1), electrode, 2, 1, i_sub),2)), squeeze(mean(ERP(1, P3TWINDOW_subj(i_sub, :,  1, 2), electrode, 2, 2, i_sub),2)) ...
        squeeze(mean(ERP(2, P3TWINDOW_subj(i_sub, :,  2, 1), electrode, 2, 1, i_sub),2)), squeeze(mean(ERP(2, P3TWINDOW_subj(i_sub, :,  2, 2), electrode, 2, 2, i_sub),2)) ...
        squeeze(mean(ERP(3, P3TWINDOW_subj(i_sub, :,  3, 1), electrode, 2, 1, i_sub),2)), squeeze(mean(ERP(3, P3TWINDOW_subj(i_sub, :,  3, 2), electrode, 2, 2, i_sub),2))]
end

Peakanalyze = [peakframe(:, 1, 1) peakframe(:, 1, 2) peakframe(:, 2, 1) peakframe(:, 2, 2)  peakframe(:, 3, 1) peakframe(:, 3, 2)]
% save([PATHOUT, 'P3_analyze1'], 'P3_analyze1')
% save([PATHOUT, 'P3TWINDOW_subj'], 'P3TWINDOW_subj')

Peakgrandaverage = round(mean([peakframe(:, 1, 1); peakframe(:, 1, 2); peakframe(:, 2, 1); peakframe(:, 2, 2) ; peakframe(:, 3, 1) ;peakframe(:, 3, 2)]));
P3TWINDOW_subj = find(EEG.times>(EEG.times(Peakgrandaverage-25)),1)-1:find(EEG.times>(EEG.times(Peakgrandaverage+25)),1)-1;
P3TWINDOW_Gave = P3TWINDOW_subj;

% % ERP = walkconds, frames, electrodes,stim, electype, subj
P3_analyze =          [squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 1, 1, :),2)), squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 1, 2, :),2)) ...
    squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 1, 1, :),2)), squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 1, 2, :),2)) ...
    squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 1, 1, :),2)), squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 1, 2, :),2)) ...
    squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 2, 1, :),2)), squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 2, 2, :),2)) ...
    squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 2, 1, :),2)), squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 2, 2, :),2)) ...
    squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 2, 1, :),2)), squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 2, 2, :),2))];

save([PATHOUT, 'P3_analyze'], 'P3_analyze')
save([PATHOUT, 'P3TWINDOW_Gave'], 'P3TWINDOW_Gave')


%% Q4 Alternate prestim noise
%STD_ERP_prestim3 = squeeze(std(STD_ERP_prestim2)); %averaged over elecs, then all trials, then std of frames
% STD_ERP_prestim3(subs, conds, type)
% structured in a table with the following headers: W1E1 W1E2 W2E1 W2E2
prestim_std_erp_analyze = [STD_ERP_prestim3(:, 1, 1), STD_ERP_prestim3(:, 1, 2),...
    STD_ERP_prestim3(:, 2, 1), STD_ERP_prestim3(:, 2, 2), ...
    STD_ERP_prestim3(:, 3, 1), STD_ERP_prestim3(:, 3, 2)];
save([PATHOUT, 'prestim_std_erp_analyze'], 'prestim_std_erp_analyze')




%% Q5 SNR Calculation
% ERP_singletrial1_even(i_sub,elec, frames,i_cond,i_type, i_stim, i_trial/2)
% electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));
for i_sub = 1:length(SUBJ)
    for i_stim = 1:2
        panel = 1;
        %                figure('Renderer', 'painters', 'Position', [10 10 1000 800]);
        
        for i_type = 1:2;
            for i_cond = 1:3;
                
                
                if i_stim == 1;
                    amp_lim = [-4 9];
                else
                    amp_lim = [-4 2];
                end
                subplot(2,3,panel)
                ERP_odd =  squeeze(mean(nanmean(ERP_singletrial1(i_sub,1, :,i_cond,i_type, i_stim, 1:2:end), 7),1));
                ERP_even = squeeze(mean(nanmean(ERP_singletrial1(i_sub,1, :,i_cond,i_type, i_stim, 2:2:end), 7),1));
                signal = mean([ERP_odd ERP_even], 2);
                noise =    abs(ERP_odd-ERP_even);
                
                
                signal_mean(i_sub,i_type, i_cond, i_stim)= mean(signal);
                noise_mean(i_sub,i_type, i_cond, i_stim)= mean(noise);
                SNR_ERP(i_sub,i_type, i_cond, i_stim) = signal_mean(i_sub,i_type, i_cond, i_stim)/noise_mean(i_sub,i_type, i_cond, i_stim);
                
                signalP3_mean(i_sub,i_type, i_cond, i_stim)= mean(signal(P3TWINDOW_Gave));
                noiseP3_mean(i_sub,i_type, i_cond, i_stim)= mean(noise(P3TWINDOW_Gave));
                SNR_P3(i_sub,i_type, i_cond, i_stim) = signalP3_mean(i_sub,i_type, i_cond, i_stim)/noiseP3_mean(i_sub,i_type, i_cond, i_stim);
                
                %                        panel = panel +1 ;
                
            end
        end
    end
    %                thetitle = ['SNR ', STIMULI{i_stim}, ' at ', EEG.chanlocs(electrode).labels] ;
    %                supertitle({thetitle ; ' '});
    %                print([GRAPHPATH, 'Data_noise\', thetitle], '-dpng');
end

% SNR_P3(i_sub,i_type, i_cond, i_stim)
SNR_analyze = [SNR_P3(:,1, 1, 1), SNR_P3(:,1, 2, 1), SNR_P3(:,1, 3, 1), SNR_P3(:,2, 1, 1), SNR_P3(:,2, 2, 1), SNR_P3(:,2, 3, 1)...
    SNR_P3(:,1, 1, 2), SNR_P3(:,1, 2, 2), SNR_P3(:,1, 3, 2), SNR_P3(:,2, 1, 2), SNR_P3(:,2, 2, 2), SNR_P3(:,2, 3, 2)]

SNR_analyze2 = [SNR_P3(:,1, 1, 1), SNR_P3(:,2, 1, 1), SNR_P3(:,1, 2, 1), SNR_P3(:,2, 2, 1), SNR_P3(:,1, 3, 1), SNR_P3(:,2, 3, 1)...
    SNR_P3(:,1, 1, 2), SNR_P3(:,2, 1, 2), SNR_P3(:,1, 2, 2), SNR_P3(:,2, 2, 2), SNR_P3(:,1, 3, 2), SNR_P3(:,2, 3, 2)]

save([PATHOUT, 'SNR_analyze2'], 'SNR_analyze2')


%% SNR Calculation 2 (CHRON)
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
%%  Figure 5: SNR comparison graph (CHRON)
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
   % subplot(2,2,i_type)
    SNR_means1 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,1),1));
    SNR_std1 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 1)))./sqrt(length(SUBJ));
    SNR_means2 =    squeeze(mean(SNR_chron(:,i_type, i_cond, :,2),1));
    SNR_std2 =      squeeze(std(SNR_chron(:,i_type, i_cond, :, 2)))./sqrt(length(SUBJ));
   % b=   boundedline( [1:100] , SNR_means2, SNR_std2, 'k',[1:100] , SNR_means1, SNR_std1, colour,  'nan', 'gap', 'alpha','transparency', 0.4 );
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNR_end(:, i_type, i_cond, 1) = SNR_means1(i_snr-1);
                 SNR_end_all(:, i_type, i_cond, 1) = SNR_chron(:,i_type, i_cond, i_snr-1,1);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNR_end(:, i_type, i_cond, 2) = SNR_means2(i_snr-1);
                   SNR_end_all(:, i_type, i_cond, 2) = SNR_chron(:,i_type, i_cond, i_snr-1,2);

                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNR_end(:, i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNR_end(:, i_type, i_cond, 2), 1))])
    
      
    
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
 
    for i_snr = 1:length(SNR_means1)
        if i_snr > 1;
            if isnan(SNR_means1(i_snr));
                if ~isnan(SNR_means1(i_snr-1));
                    SNR_end( i_type, i_cond, 1) = SNR_means1(i_snr-1);
                 SNR_end_all(:, i_type, i_cond, 1) = SNR_chron(:,i_type, i_cond, i_snr-1,1);
                end
            end
            if isnan(SNR_means2(i_snr));
                if ~isnan(SNR_means2(i_snr-1));
                    SNR_end( i_type, i_cond, 2) = SNR_means2(i_snr-1);
                   SNR_end_all(:, i_type, i_cond, 2) = SNR_chron(:,i_type, i_cond, i_snr-1,2);

                end
            end
        end
    end
    disp([ CONDNAMES{i_cond}, ' - ', ETYPE{i_type}])
    disp(['final SNR (rej) = ', num2str(mean(SNR_end( i_type, i_cond, 1), 1))])
    disp(['final SNR (nonrej) = ', num2str(mean(SNR_end( i_type, i_cond, 2), 1))])
    
end
 %W1E1TR W1E2TR W2E1TR W2E2TR W3E1TR W3E2TR W1E1TN W1E2TN W2E1TN W2E2TN
    %W3E1TN W3E2TN
    SNR_end_analyze = [SNR_end_all(:, 1, 1, 1) SNR_end_all(:, 2, 1, 1) SNR_end_all(:, 1, 2, 1) SNR_end_all(:, 2, 2, 1) ...
        SNR_end_all(: ,1, 1, 2) SNR_end_all(:, 2, 1, 2) SNR_end_all(:, 1, 2, 2) SNR_end_all(:, 2, 2, 2)  ];
 save([PATHOUT, 'SNR_end_analyze'], 'SNR_end_analyze');
    close
    
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

%% SUPPLEMENTAL? Average Single trial STD
% MAY NEED TO CHANGE THIS (AND ALSO MARK CHANGE IN MANUSCRIPT)
% STD_singletrial1_allstim(i_sub,i_cond,i_type, trials)


% Data noise analysis
STDnoise_analyze = [squeeze(nanmean(STD_singletrial1_allstim(:,1,1, :),4)),squeeze(nanmean(STD_singletrial1_allstim(:,1,2, :),4)), ...
    squeeze(nanmean(STD_singletrial1_allstim(:,2,1, :),4)) ,squeeze(nanmean(STD_singletrial1_allstim(:,2,2, :),4)),...
    squeeze(nanmean(STD_singletrial1_allstim(:,3,1, :),4)) , squeeze(nanmean(STD_singletrial1_allstim(:,3,2, :),4))];
%
% thetitle = 'Mean Single trial STD Data noise' ;
% supertitle({thetitle ; ' '});
% print([GRAPHPATH, 'Data_noise\', thetitle], '-dpng');

save([PATHOUT, 'STDnoise_analyze'], 'STDnoise_analyze');