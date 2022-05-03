clear all
close all

trialrej = {'rej';'norej'};
whichtrials = 2; % 1 or 2 for rej or norej

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
RMS_single = nan(18, 3, 2, 2, 505);
tcounter = nan(18, 2, 3, 2);

ERPRMSNOISE = nan(3, 2, 2, 18);
ERPSTDNOISE = nan(3, 2, 2, 18);

%ERP_singletrial1 = nan(18,62,250,3,2,2,505); % too big: maybe not needed?
%maybe only one electrode, etc
ERP_singletrial1 = nan(18,1,250,3,2,2,600);
RMS_singletrial1 = nan(18,3,2,2,600);
STD_singletrial1 = nan(18,3,2,2,600);
% ERP_singletrial1_even = nan(18,62,250,3,2,2,253);
% ERP_singletrial1_odd  = nan(18,62,250,3,2,2,253);
% too big: maybe only use one electrode, shorter time period, only 2 conds, etc
%
 RMS_singletrial2 = nan(18,3,2,2); %(i_sub,i_cond,i_type, i_stim);
 STD_singletrial2 = nan(18,3,2,2); %(i_sub,i_cond,i_type, i_stim);
 
 RMS_singletrial1_allstim = nan(18,3,2, 1200);
 STD_singletrial1_allstim = nan(18,3,2, 1200);

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
%% Load data
for i_sub = 1:length(SUBJ)
    for i_type = 1:length(ETYPE)
        for i_cond = 1:length(WALKCONDS)
            for i_stim = 1:length(STIMULI)
                clear data2
                       
                EEG = pop_loadset('filename', [PATH,'jos_sync', SUBJ{i_sub}, '_', ETYPE{i_type}, '_oddball_epoched_', trialrej{whichtrials},'_',STIMULI{i_stim},'_' WALKCONDS{i_cond}, '.set']);

                %  ERP = walkconds, frames, electrodes,stim, type, subj
              %ERP = nan(length(WALKCONDS),250, 62,2,2, length(SUBJ));

                ERP(i_cond,:,:, i_stim, i_type,i_sub) = mean(EEG.data ,3)'; %targets
                NOISE(i_cond,:,:, i_stim, i_type,i_sub) = squeeze(rms(ERP(i_cond, find(EEG.times<0),:, i_stim, i_type, i_sub)));
                
        % ERP noise:average over trials, then RMS of <0 times, then average over electrodes, then over participants 
        ERPRMSNOISE(i_cond, i_stim, i_type, i_sub) =  nanmean(squeeze( rms(ERP(i_cond, find(EEG.times<0), :, i_stim, i_type, i_sub))));
        ERPSTDNOISE(i_cond, i_stim, i_type, i_sub) =  nanmean(squeeze( std(ERP(i_cond, find(EEG.times<0), :, i_stim, i_type, i_sub))));
        
              %  EEG.setname = ['jos_sync', SUBJ{i_sub}, '_', ETYPE{i_type}, '_oddball_epoched_rejected_',STIMULI{i_stim},'_' WALKCONDS{i_cond}]
              %  [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
                
                n_trials = EEG.trials;
                
                RMS(i_sub,i_cond,i_type, i_stim) = mean(mean(squeeze(rms( EEG.data(:, EEG.times<0,n_trials) ,2)),2));                
                rms_erp_out(i_sub,i_cond,i_type, i_stim) = mean(rms(mean(EEG.data(:,EEG.times<0,n_trials),3),2));                   
                
                % now let's count trials
                tcounter(i_sub, i_type, i_cond, i_stim) = n_trials;
                
                for i_trial = 1:n_trials
                    
% single trial noise: for each trial, RMS value for each electrode, then averaged over electrodes, then average over  trials (for each participant), then over subjects
% here I want 1 value for each electrode, then to average over all electrodes
                 
                 ERP_singletrial1(i_sub,:, :,i_cond,i_type, i_stim, i_trial) = ( EEG.data(12, :,i_trial) ); %select Pz
               %  RMS_singletrial1(i_sub,i_cond,i_type, i_stim, i_trial) = mean(rms( EEG.data(:, EEG.times<0,i_trial) ,2));
                     % EEG.data(elec, frames, trials)
                 data2= permute(EEG.data, [2, 1, 3]); % data(frames, elec, trials)
                 STD_singletrial1(i_sub, i_cond, i_type, i_stim, i_trial) = mean(std(data2(EEG.times<0, :, i_trial)));
             
                end
                % then average over  trials (for each participant), ...then
                % over subjects (later)
                 RMS_singletrial2(i_sub,i_cond,i_type, i_stim) = squeeze(nanmean(RMS_singletrial1(i_sub,i_cond,i_type, i_stim, :), 5));
                 STD_singletrial2(i_sub,i_cond,i_type, i_stim) = squeeze(nanmean(STD_singletrial1(i_sub,i_cond,i_type, i_stim, :), 5));

            end 
            %remember the size here should be preallocated to be big enough
%             for both standards and targets
                 RMS_singletrial1_allstim(i_sub,i_cond,i_type, :) = [squeeze(RMS_singletrial1(i_sub,i_cond,i_type, 1, :)); squeeze(RMS_singletrial1(i_sub,i_cond,i_type, 2, :))];
                 STD_singletrial1_allstim(i_sub,i_cond,i_type, :) = [squeeze(STD_singletrial1(i_sub,i_cond,i_type, 1, :)); squeeze(STD_singletrial1(i_sub,i_cond,i_type, 2, :))];

        end
    end
end
% Make the ERP difference wave
erp_diff_out = squeeze(ERP(:,:,:,1,:,:)-ERP(:,:,:,2,:,:));

% Average over standard and targets for the ERP noise (not sure if this is
% valid because the standards and targets have different numbers)
% ERPRMSNOISE(i_cond, i_stim, i_type, i_sub)
%ERPRMSNOISE_allstim(i_cond, i_type, i_sub)
  ERPRMSNOISE_alstim = squeeze(mean(ERPRMSNOISE(:, :, :, :), 2));

 %% Save (don't need to be done every time)
% save([PATHOUT, 'ERP_r'], 'ERP')
% save([PATHOUT,'STD_singletrial1_allstim_r'], 'STD_singletrial1_allstim')
% save([PATHOUT, 'EEG_1sub_r'], 'EEG')
% save([PATHOUT, 'ERP_singletrial1'], 'ERP_singletrial1')



% eegtimes = EEG.times;
% eegchanlocs = EEG.chanlocs;
% save([PATHOUT, 'times'], 'eegtimes');
% save([PATHOUT, 'chanlocs'], 'eegchanlocs');


eeglab redraw;

%% electrode trial num comparisons
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
Trialnums_r_analyze = [tcounter(:, :, 1, 1) tcounter(:, :, 2, 1) tcounter(:, :, 3, 1)...
                     tcounter(:, :, 1, 2) tcounter(:, :, 2, 2) tcounter(:, :, 3, 2)];
save([PATHOUT, 'Trialnums_norej'], 'Trialnums_r_analyze')

%% P3 peaks
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


Peakanalyze = [peakframe(:, 1, 1) peakframe(:, 1, 2) peakframe(:, 2, 1) peakframe(:, 2, 2)  peakframe(:, 3, 1) peakframe(:, 3, 2)];
% save('O:\projects\jos_sync1\sync_walking\data\ana4_processed\P3_analyze1', 'P3_analyze1')
% save('O:\projects\jos_sync1\sync_walking\data\ana4_processed\P3TWINDOW_subj', 'P3TWINDOW_subj')

Peakgrandaverage = round(mean([peakframe(:, 1, 1); peakframe(:, 1, 2); peakframe(:, 2, 1); peakframe(:, 2, 2) ; peakframe(:, 3, 1) ;peakframe(:, 3, 2)]));
P3TWINDOW_Gave = find(EEG.times>(EEG.times(Peakgrandaverage-25)),1)-1:find(EEG.times>(EEG.times(Peakgrandaverage+25)),1)-1;


% % ERP = walkconds, frames, electrodes,stim, electype, subj
P3_analyze =          [squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 1, 1, :),2)), squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 1, 2, :),2)) ...
                       squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 1, 1, :),2)), squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 1, 2, :),2)) ...
                       squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 1, 1, :),2)), squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 1, 2, :),2)) ...
                       squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 2, 1, :),2)), squeeze(mean(ERP(1, P3TWINDOW_Gave, electrode, 2, 2, :),2)) ...
                       squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 2, 1, :),2)), squeeze(mean(ERP(2, P3TWINDOW_Gave, electrode, 2, 2, :),2)) ...
                       squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 2, 1, :),2)), squeeze(mean(ERP(3, P3TWINDOW_Gave, electrode, 2, 2, :),2))];

 save([PATHOUT, 'P3_analyze'], 'P3_analyze');
 save([PATHOUT, 'P3TWINDOW_Gave'], 'P3TWINDOW_Gave');

%% grand average plots
electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));

figure('Renderer', 'painters', 'Position', [10 10 800 600]);
for i_cond = 1:2
    amp_lim = [-6 10];
    i_type = 1;
    switch i_cond
        case 1
            colour = 'm';
        case 2
            colour = 'g';
            %         case 3
            %             colour = 'b';
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
    
    if i_cond ==1;
        legend(b,{'Targets', 'Standards'},'Location','NorthEast');
    end
    
    axis tight; ylim(amp_lim);
    
    L1 =     line([-200 1000],[0 0],'color','k');
    L2 =     line([0 0],amp_lim,'color','k');
    %   P3TWINDOW_subj(i_sub, :,  i_cond, i_type)
    
    title([ CONDNAMES{i_cond}, ' (', ETYPE{i_type}, ')'] );
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    set(L1,'LineWidth',1.5);
    set(L2,'LineWidth',1.5);
    set(gca, 'Box', 'off')
    
    if i_type ==1;
        legend(b,{'Targets', 'Standards'},'Location','NorthEast', 'LineWidth', 1);
%         legend('boxoff')
    end
    % Second electrode type
    %standard and target waveforms
    i_type = 2;
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
    L1 =     line([-200 1000],[0 0],'color','k');
    L2 =     line([0 0],amp_lim,'color','k');
    title([CONDNAMES{i_cond}, ' (', ETYPE{i_type}, ')'] );
    xlabel('Time (ms)');
    ylabel('Voltage (\muV)');
    set(L1,'LineWidth',1.5);
    set(L2,'LineWidth',1.5);
    % legend({'Targets', 'Standards'},'Location','NorthEast');
    set(gca, 'Box', 'off')
    
end
%electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));

thetitle = ['meanERPs',EEG.chanlocs(electrode).labels ]
suptitle([{thetitle} ]);
print([PATHOUT, '\graphs\', thetitle], '-dpng') 

%% SNR Calculation
% ERP_singletrial1_even(i_sub,elec, frames,i_cond,i_type, i_stim, i_trial/2)
electrode = find(strcmp({EEG.chanlocs.labels}, 'Pz'));
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
                       ERP_odd =  squeeze(mean(nanmean(ERP_singletrial1(i_sub,electrode, :,i_cond,i_type, i_stim, 1:2:end), 7),1));
                       ERP_even = squeeze(mean(nanmean(ERP_singletrial1(i_sub,electrode, :,i_cond,i_type, i_stim, 2:2:end), 7),1));
                       signal = mean([ERP_odd ERP_even], 2);
                       noise =    abs(ERP_odd-ERP_even);
                       
%                        plot(EEG.times, ERP_odd, 'r')
%                        hold on
%                        plot(EEG.times, ERP_even, 'y')
%                        plot(EEG.times, signal, 'k', 'LineWidth', 2)
%                        plot(EEG.times, noise,'b', 'LineWidth', 2)
                       
                       
%                        
%                        axis tight; ylim(amp_lim);
%                        L1 =     line([-200 1000],[0 0],'color','k');
%                        L2 =     line([0 0],amp_lim,'color','k');
%                        ylabel('Amplitude(uV)')
%                        xlabel('Time (ms)')
%                        legend({'odd trials', 'even trials', 'Signal', 'Noise'})
%                        title([CONDNAMES{i_cond}, ' (', ETYPE{i_type}, ')'])
%                        set(L1,'LineWidth',2);
%                        set(L2,'LineWidth',2);
%                        
%                        set(gca,'LineWidth',2);
%                        set(gca,'Color',[1 1 1]);
%                        set(gca,'YDir','reverse');
%                        set(gca,'Box','off');
                       
%                        signal_mean(i_sub,i_type, i_cond, i_stim)= mean(signal);
%                        noise_mean(i_sub,i_type, i_cond, i_stim)= mean(noise);
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
%                print(['O:\projects\jos_sync1\sync_walking\graphs\Data_noise\', thetitle], '-dpng');
      end 
      
      % SNR_P3(i_sub,i_type, i_cond, i_stim)
      SNR_analyze = [SNR_P3(:,1, 1, 1), SNR_P3(:,1, 2, 1), SNR_P3(:,1, 3, 1), SNR_P3(:,2, 1, 1), SNR_P3(:,2, 2, 1), SNR_P3(:,2, 3, 1)...
                     SNR_P3(:,1, 1, 2), SNR_P3(:,1, 2, 2), SNR_P3(:,1, 3, 2), SNR_P3(:,2, 1, 2), SNR_P3(:,2, 2, 2), SNR_P3(:,2, 3, 2)]
                 
     SNR_analyze2 = [SNR_P3(:,1, 1, 1), SNR_P3(:,2, 1, 1), SNR_P3(:,1, 2, 1), SNR_P3(:,2, 2, 1), SNR_P3(:,1, 3, 1), SNR_P3(:,2, 3, 1)...
                     SNR_P3(:,1, 1, 2), SNR_P3(:,2, 1, 2), SNR_P3(:,1, 2, 2), SNR_P3(:,2, 2, 2), SNR_P3(:,1, 3, 2), SNR_P3(:,2, 3, 2)]
                 
    save([PATHOUT, 'SNR_analyze2'], 'SNR_analyze2')