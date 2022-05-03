clear all
close all

trialrej = {'rej';'norej'};
whichtrials = 1; % 1 or 2 for rej or norej

EXPERIMENT = 'jos_sync';
MAINPATH = 'D:\4. home office\otto backup\jos_sync1\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATH = [MAINPATH, 'data\R1_ana4_processed\',  trialrej{whichtrials}];

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
% structured in a 'table' with the following (imaginary) headers: W1E1 W1E2 W2E1 W2E2 
% W3E1 W3E2 (may include another 6 columns for S1 vs. S2 analyses)
load([PATH, 'Trialnums_analyze.mat'])
load([PATH(1:63), 'Accdata_analyze.mat'])
load([PATH, 'P3_analyze.mat'])
load([PATH, 'STDnoise_analyze.mat'])
load([PATH, 'SNR_analyze2.mat'])
load([PATH, 'prestim_std_erp_analyze.mat'])


% additional data
load([PATH, 'P3TWINDOW_Gave.mat'])
load([PATH, 'EEG_1sub.mat']) % to get eeg infos like EEG.times and EEG.chanlocs

%pool stimulus for trialnums
Trialnums2x2 = Trialnums_analyze(:, 1:6) + Trialnums_analyze(:, 7:12);
save([PATH, 'Trialnums_stimpool_analyze.mat'], 'Trialnums2x2')



% % Save as csv files for JASP [only needed to do this once]
% % W1E1, W1E2, W2E1, W2E2, W3E1, W3E2
csvwrite([PATH, 'Trialnums_analyze.csv'], Trialnums_analyze);
csvwrite([PATH, 'Accdata_analyze.csv'], Accdata_analyze);
csvwrite([PATH, 'P3_analyze.csv'], P3_analyze);
csvwrite([PATH, 'STDnoise_analyze.csv'], STDnoise_analyze);
csvwrite([PATH, 'SNR_analyze2.csv'], SNR_analyze2);
csvwrite([PATH, 'Trialnums2x2.csv'], Trialnums2x2 );
csvwrite([PATH, 'prestim_std_erp_analyze.csv'], prestim_std_erp_analyze );


%% Stats Q1a Trialnums
disp('Stats Q1a Trialnums -----------------------------------------');
% 2x2x2 anova is done in spss, and here are the paired t-tests to go with
% it

%ANOVA 
%F1: walkcond
%F2: electrode
measure = 'Trialnums';

sub_means_W1E1S1 = Trialnums_analyze(:,1);
sub_means_W1E2S1 = Trialnums_analyze(:,2);
sub_means_W2E1S1 = Trialnums_analyze(:,3);
sub_means_W2E2S1 = Trialnums_analyze(:,4);

sub_means_W1E1S2 = Trialnums_analyze(:,7);
sub_means_W1E2S2 = Trialnums_analyze(:,8);
sub_means_W2E1S2 = Trialnums_analyze(:,9);
sub_means_W2E2S2 = Trialnums_analyze(:,10);

% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS

%ttests!


%How many of each trial were kept after artifact rejection? (this pools
%walking conditions over electrodes. It's just for the methods
disp( ['Active (Mtarg = ',num2str(mean([sub_means_W1E1S1;sub_means_W2E1S1])),', rangetarg = ', num2str(min([sub_means_W1E1S1;sub_means_W2E1S1])), '-', num2str(max([sub_means_W1E1S1;sub_means_W2E1S1]))  ...
    '; Mstand = ', num2str(mean([sub_means_W1E1S2;sub_means_W2E1S2])), ', Rangestand = ' num2str(min([sub_means_W1E1S2;sub_means_W2E1S2])), '-', num2str(max([sub_means_W1E1S2;sub_means_W2E1S2]))  ')'])

disp( ['Passive (Mtarg = ',num2str(mean([sub_means_W1E2S1;sub_means_W2E2S1])),', rangetarg = ', num2str(min([sub_means_W1E2S1;sub_means_W2E2S1])), '-', num2str(max([sub_means_W1E2S1;sub_means_W2E2S1]))  ...
    '; Mstand = ', num2str(mean([sub_means_W1E2S2;sub_means_W2E2S2])), ', Rangestand = ' num2str(min([sub_means_W1E2S2;sub_means_W2E2S2])), '-', num2str(max([sub_means_W1E2S2;sub_means_W2E2S2]))  ')'])



%% Stats Q1b POOLED Trialnums
disp('Stats Q1b POOLED Trialnums------------------------------------------------');
% 2x2x2 anova is done in spss, but here we have a 2x2 anova. 
% Here we will ignore stimulus, because the prestimulus noise analysis also
% ignores stimulus, and having 3 factors is kind of ridiculously
% complicated. So we will start by adding these together (above)

%ANOVA 
%F1: walkcond
%F2: electrode
measure = 'Trialnums (pooled) anova';

sub_means_W1E1=  Trialnums2x2(:,1);
sub_means_W1E2=  Trialnums2x2(:,2);
sub_means_W2E1 = Trialnums2x2(:,3);
sub_means_W2E2 = Trialnums2x2(:,4);


disp( ['Standing (Mean = ',num2str(mean([sub_means_W1E1; sub_means_W1E2])),'; range = ', num2str(min([sub_means_W1E1; sub_means_W1E2])), '-', num2str(max([sub_means_W1E1; sub_means_W1E2]))  ')'])
disp( ['Walking (Mean = ',num2str(mean([sub_means_W2E1; sub_means_W2E2])),'; range = ', num2str(min([sub_means_W2E1; sub_means_W2E2])), '-', num2str(max([sub_means_W2E1; sub_means_W2E2]))  ')'])


% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS
alldat = [sub_means_W1E1 sub_means_W1E2...
            sub_means_W2E1 sub_means_W2E2];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}];

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = '])
disp(['Electrode: F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ' num2str(stats{3,6}), '; etap2 = '])
disp(['Mobility * Electrode: F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = '])


%ttests!

%Standing: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1,sub_means_W1E2);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W1E2);
cD = computeCohen_d(sub_means_W1E1, sub_means_W1E2, 'paired');
disp('Standing: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1,sub_means_W2E2);
Mdiff = mean(sub_means_W2E1)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W2E1, sub_means_W2E2, 'paired');
disp('Walking: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1, sub_means_W2E1);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W2E1);
disp('Active: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1, sub_means_W2E1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2, sub_means_W2E2);
Mdiff = mean(sub_means_W1E2)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W1E2, sub_means_W2E2, 'paired');
disp('Passive: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1) mean(sub_means_W2E1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2) mean(sub_means_W2E2)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title(measure)

disp('-------------------------------------------------------------------------------------');
%% Stats Q1c TARGET Trialnums
disp('Stats Q1c TARGET Trialnums --------------------------------------------')
% 2x2x2 anova is done in spss, but here we have a 2x2 anova. 
% Here we will ignore stimulus, because the prestimulus noise analysis also
% ignores stimulus, and having 3 factors is kind of ridiculously
% complicated. So we will start by adding these together (above)

%ANOVA 
%F1: walkcond
%F2: electrode
measure = 'Trialnums (targets) anova';

sub_means_W1E1=  Trialnums_analyze(:,1);
sub_means_W1E2=  Trialnums_analyze(:,2);
sub_means_W2E1 = Trialnums_analyze(:,3);
sub_means_W2E2 = Trialnums_analyze(:,4);
% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS


% % %Means & ranges:
% disp( ['Standing (Mtarg = ',num2str(mean([sub_means_W1E1;sub_means_W2E1])),', rangetarg = ', num2str(min([sub_means_W1E1S1;sub_means_W2E1S1])), '-', num2str(max([sub_means_W1E1S1;sub_means_W2E1S1]))  ...
%     '; Mstand = ', num2str(mean([sub_means_W1E1S2;sub_means_W2E1S2])), ', Rangestand = ' num2str(min([sub_means_W1E1S2;sub_means_W2E1S2])), '-', num2str(max([sub_means_W1E1S2;sub_means_W2E1S2]))  ')'])
% 
% disp( ['Walking (Mtarg = ',num2str(mean([sub_means_W1E2;sub_means_W2E2])),', rangetarg = ', num2str(min([sub_means_W1E2S1;sub_means_W2E2S1])), '-', num2str(max([sub_means_W1E2S1;sub_means_W2E2S1]))  ...
%     '; Mstand = ', num2str(mean([sub_means_W1E2S2;sub_means_W2E2S2])), ', Rangestand = ' num2str(min([sub_means_W1E2S2;sub_means_W2E2S2])), '-', num2str(max([sub_means_W1E2S2;sub_means_W2E2S2]))  ')'])


alldat = [sub_means_W1E1 sub_means_W1E2...
            sub_means_W2E1 sub_means_W2E2];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}];

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = '])
disp(['Electrode: F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ', num2str(stats{3,6}), '; etap2 = '])
disp(['Mobility * Electrode: F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = '])


%ttests!

%Standing: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1,sub_means_W1E2);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W1E2);
cD = computeCohen_d(sub_means_W1E1, sub_means_W1E2, 'paired');
disp('Standing: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1,sub_means_W2E2);
Mdiff = mean(sub_means_W2E1)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W2E1, sub_means_W2E2, 'paired');
disp('Walking: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1, sub_means_W2E1);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W2E1);
disp('Active: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1, sub_means_W2E1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2, sub_means_W2E2);
Mdiff = mean(sub_means_W1E2)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W1E2, sub_means_W2E2, 'paired');
disp('Passive: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1) mean(sub_means_W2E1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2) mean(sub_means_W2E2)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title(measure)
disp('-------------------------------------------------------------------------------------');
 %% Stats Q2 Accelerometer variance/ walking activity
disp('Stats Q2 Accelerometer variance/ walking activity------------------------------------')
% 2x2x2 anova is done in spss, but here we have a 2x2 anova. 
% Here we will ignore stimulus, because the prestimulus noise analysis also
% ignores stimulus, and having 3 factors is kind of ridiculously
% complicated. So we will start by adding these together (above)

%ANOVA 
%F1: walkcond
%F2: electrode

measure = 'Accelerometer variance';

sub_means_W1E1=  Accdata_analyze(:,1);
sub_means_W1E2=  Accdata_analyze(:,2);
sub_means_W2E1 = Accdata_analyze(:,3);
sub_means_W2E2 = Accdata_analyze(:,4);
% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS
alldat = [sub_means_W1E1 sub_means_W1E2...
            sub_means_W2E1 sub_means_W2E2];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}]

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = '])
disp(['Electrode: F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ' num2str(stats{3,6}), '; etap2 = '])
disp(['Mobility * Electrode: F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = '])


%ttests!

%Standing: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1,sub_means_W1E2);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W1E2);
cD = computeCohen_d(sub_means_W1E1, sub_means_W1E2, 'paired');
disp('Standing: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1,sub_means_W2E2);
Mdiff = mean(sub_means_W2E1)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W2E1, sub_means_W2E2, 'paired');
disp('Walking: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1, sub_means_W2E1);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W2E1);
disp('Active: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1, sub_means_W2E1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2, sub_means_W2E2);
Mdiff = mean(sub_means_W1E2)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W1E2, sub_means_W2E2, 'paired');
disp('Passive: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Check direction
mean(sub_means_W1E1);
std(sub_means_W1E1);
mean(sub_means_W1E2);
std(sub_means_W1E2);

mean(sub_means_W2E1);
std(sub_means_W2E1);
mean(sub_means_W2E2);
std(sub_means_W2E2);

%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1) mean(sub_means_W2E1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2) mean(sub_means_W2E2)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title(measure)

disp('-------------------------------------------------------------------------------------');
%% Stats Q3A P3 t-tests 2x2 anova
disp('Stats Q3A P3 t-tests 2x2 anova ------------------------------------------------------');

measure = 'P3 amplitude';

sub_means_W1E1S1=  P3_analyze(:,1);
sub_means_W1E2S1=  P3_analyze(:,2);
sub_means_W2E1S1 = P3_analyze(:,3);
sub_means_W2E2S1 = P3_analyze(:,4);


sub_means_W1E1S2=  P3_analyze(:,7);
sub_means_W1E2S2=  P3_analyze(:,8);
sub_means_W2E1S2 = P3_analyze(:,9);
sub_means_W2E2S2 = P3_analyze(:,10);
%ttests!

%Standing targs: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1S1,sub_means_W1E2S1);
Mdiff = mean(sub_means_W1E1S1)-mean(sub_means_W1E2S1);
cD = computeCohen_d(sub_means_W1E1S1, sub_means_W1E2S1, 'paired');
disp('Standing targs: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking targs: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1S1,sub_means_W2E2S1);
Mdiff = mean(sub_means_W2E1S1)-mean(sub_means_W2E2S1);
cD = computeCohen_d(sub_means_W2E1S1, sub_means_W2E2S1, 'paired');
disp('Walking targs: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active targs: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1S1, sub_means_W2E1S1);
Mdiff = mean(sub_means_W1E1S1)-mean(sub_means_W2E1S1);
disp('Active targs: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1S1, sub_means_W2E1S1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive targs: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2S1, sub_means_W2E2S1);
Mdiff = mean(sub_means_W1E2S1)-mean(sub_means_W2E2S1);
cD = computeCohen_d(sub_means_W1E2S1, sub_means_W2E2S1, 'paired');
disp('Passive targs: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

  
 
%Active standing: Targets vs. standards
[h p ci test] = ttest(sub_means_W1E1S1, sub_means_W1E1S2);
Mdiff = mean(sub_means_W1E1S1)-mean(sub_means_W1E1S2);
cD = computeCohen_d(sub_means_W1E1S1, sub_means_W1E1S2, 'paired');
disp('Active standing: Targets vs. standards')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Passive standing: Targets vs. standards
[h p ci test] = ttest(sub_means_W1E2S1, sub_means_W1E2S2);
Mdiff = mean(sub_means_W1E2S1)-mean(sub_means_W1E2S2);
cD = computeCohen_d(sub_means_W1E2S1, sub_means_W1E2S2, 'paired');
disp('Passive standing: Targets vs. standards')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

  
 
%Active walking: Targets vs. standards
[h p ci test] = ttest(sub_means_W2E1S1, sub_means_W2E1S2);
Mdiff = mean(sub_means_W2E1S1)-mean(sub_means_W2E1S2);
cD = computeCohen_d(sub_means_W2E1S1, sub_means_W2E1S2, 'paired');
disp('Active walking: Targets vs. standards')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Passive walking: Targets vs. standards
[h p ci test] = ttest(sub_means_W2E2S1, sub_means_W2E2S2);
Mdiff = mean(sub_means_W2E2S1)-mean(sub_means_W2E2S2);
cD = computeCohen_d(sub_means_W2E2S1, sub_means_W2E2S2, 'paired');
disp('Passive walking: Targets vs. standards')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])
  mean(sub_means_W1E1S1);
std(sub_means_W1E1S1);
mean(sub_means_W1E2S1);
std(sub_means_W1E2S1);

mean(sub_means_W2E1S1);
std(sub_means_W2E1S1);
mean(sub_means_W2E2S1);
std(sub_means_W2E2S1);

%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1S1) mean(sub_means_W2E1S1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2S1) mean(sub_means_W2E2S1)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title([measure, 'Target'])

disp('-------------------------------------------------------------------------------------');
%% Stats Q3B P3 TARGET t-tests
disp('Stats Q3B P3 TARGET t-tests ---------------------------------------------------------');
% no anova because it's 2x2x2 and was carried out on SPSS

measure = 'P3 amplitude';

sub_means_W1E1S1=  P3_analyze(:,1);
sub_means_W1E2S1=  P3_analyze(:,2);
sub_means_W2E1S1 = P3_analyze(:,3);
sub_means_W2E2S1 = P3_analyze(:,4);

% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS
alldat = [sub_means_W1E1S1 sub_means_W1E2S1...
            sub_means_W2E1S1 sub_means_W2E2S1];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}]

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = '])
disp(['Electrode: F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ' num2str(stats{3,6}), '; etap2 = '])
disp(['Mobility * Electrode: F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = '])

% t-tests

%Standing targs: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1S1,sub_means_W1E2S1);
Mdiff = mean(sub_means_W1E1S1)-mean(sub_means_W1E2S1);
cD = computeCohen_d(sub_means_W1E1S1, sub_means_W1E2S1, 'paired');
disp('Standing targs: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking targs: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1S1,sub_means_W2E2S1);
Mdiff = mean(sub_means_W2E1S1)-mean(sub_means_W2E2S1);
cD = computeCohen_d(sub_means_W2E1S1, sub_means_W2E2S1, 'paired');
disp('Walking targs: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active targs: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1S1, sub_means_W2E1S1);
Mdiff = mean(sub_means_W1E1S1)-mean(sub_means_W2E1S1);
disp('Active targs: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1S1, sub_means_W2E1S1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive targs: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2S1, sub_means_W2E2S1);
Mdiff = mean(sub_means_W1E2S1)-mean(sub_means_W2E2S1);
cD = computeCohen_d(sub_means_W1E2S1, sub_means_W2E2S1, 'paired');
disp('Passive targs: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

disp('-------------------------------------------------------------------------------------');
%% Stats Q4: ERP Baseline noise (NEW VERSION)
disp('Stats Q4 VERSION 2: ERP Baseline noise ----------------------------------------------');



%ANOVA 
%F1: walkcond
%F2: electrode

measure = 'ERP Baseline data noise';

sub_means_W1E1=  prestim_std_erp_analyze(:,1);
sub_means_W1E2=  prestim_std_erp_analyze(:,2);
sub_means_W2E1 = prestim_std_erp_analyze(:,3);
sub_means_W2E2 = prestim_std_erp_analyze(:,4);
% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS
alldat = [sub_means_W1E1 sub_means_W1E2...
            sub_means_W2E1 sub_means_W2E2];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}];

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: (F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = )'])
disp(['Electrode: (F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ' num2str(stats{3,6}), '; etap2 = )'])
disp(['Mobility * Electrode: (F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = )'])


%ttests!

%Standing: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1,sub_means_W1E2);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W1E2);
cD = computeCohen_d(sub_means_W1E1, sub_means_W1E2, 'paired');
disp('Standing: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1,sub_means_W2E2);
Mdiff = mean(sub_means_W2E1)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W2E1, sub_means_W2E2, 'paired');
disp('Walking: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1, sub_means_W2E1);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W2E1);
disp('Active: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1, sub_means_W2E1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2, sub_means_W2E2);
Mdiff = mean(sub_means_W1E2)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W1E2, sub_means_W2E2, 'paired');
disp('Passive: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Check direction
mean(sub_means_W1E1);
std(sub_means_W1E1);
mean(sub_means_W1E2);
std(sub_means_W1E2);

mean(sub_means_W2E1);
std(sub_means_W2E1);
mean(sub_means_W2E2);
std(sub_means_W2E2);

%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1) mean(sub_means_W2E1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2) mean(sub_means_W2E2)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title(measure)

disp('-------------------------------------------------------------------------------------');
%% Stats Q5 SNR 1
disp('Stats Q5 SNR ------------------------------------------------------------------------');
%2x2Anova on only targets, because the standards aren't realy the signal we
%are looking at


%ANOVA 
%F1: walkcond
%F2: electrode

measure = 'SNR';

sub_means_W1E1=  SNR_analyze2(:,1);
sub_means_W1E2=  SNR_analyze2(:,2);
sub_means_W2E1 = SNR_analyze2(:,3);
sub_means_W2E2 = SNR_analyze2(:,4);
% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS
alldat = [sub_means_W1E1 sub_means_W1E2...
            sub_means_W2E1 sub_means_W2E2];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}]

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: (F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = )'])
disp(['Electrode: (F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ' num2str(stats{3,6}), '; etap2 = )'])
disp(['Mobility * Electrode: (F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = )'])


%ttests!

%Standing: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1,sub_means_W1E2);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W1E2);
cD = computeCohen_d(sub_means_W1E1, sub_means_W1E2, 'paired');
disp('Standing: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1,sub_means_W2E2);
Mdiff = mean(sub_means_W2E1)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W2E1, sub_means_W2E2, 'paired');
disp('Walking: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1, sub_means_W2E1);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W2E1);
disp('Active: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1, sub_means_W2E1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2, sub_means_W2E2);
Mdiff = mean(sub_means_W1E2)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W1E2, sub_means_W2E2, 'paired');
disp('Passive: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Check direction
mean(sub_means_W1E1);
std(sub_means_W1E1);
mean(sub_means_W1E2);
std(sub_means_W1E2);

mean(sub_means_W2E1);
std(sub_means_W2E1);
mean(sub_means_W2E2);
std(sub_means_W2E2);

%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1) mean(sub_means_W2E1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2) mean(sub_means_W2E2)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title(measure)

disp('-------------------------------------------------------------------------------------');

%% Stats #6a: ERP Data noise Correlations
disp('Stats #6a: ERP Data noise Correlations ----------------------------------------------');
% STD_singletrial1_allstim(i_sub,i_cond,i_type, trials)
figure;
% Between electrodes correlations
% Standing cond
[rho, pval] = corr(prestim_std_erp_analyze(:,1),prestim_std_erp_analyze(:,2) );
disp('Prestim noise correlation between electrode types,  standing');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
%walkA
[rho, pval] = corr(prestim_std_erp_analyze(:,3),prestim_std_erp_analyze(:,4)  );
disp('Prestim noise correlation between electrode types, walking alone');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shepherd's pi
[stats] = ScatterOutliers(prestim_std_erp_analyze(:,1),prestim_std_erp_analyze(:,2) );
disp('Prestim noise (Shepherds pi) correlation between electrode types, standing');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(prestim_std_erp_analyze(:,3),prestim_std_erp_analyze(:,4) );
disp('Prestim noise (Shepherds pi) correlation between electrode types, walking ');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);

% %Print to folder
% thetitle = [ 'Correlation data noise'];
% supertitle(thetitle);
% print(['O:\projects\jos_sync1\sync_walking\graphs\Correlations\', thetitle], '-dpng')
disp('-------------------------------------------------------------------------------------');
%% Stats #6a2: SNR Correlations
disp('Stats #6a2: SNR Correlations ----------------------------------------------');
% STD_singletrial1_allstim(i_sub,i_cond,i_type, trials)
figure;
% Between electrodes correlations
% Standing cond
[rho, pval] = corr(SNR_analyze2(:,1),SNR_analyze2(:,2) );
disp('SNR correlation between electrode types,  standing');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
%walkA
[rho, pval] = corr(SNR_analyze2(:,3),SNR_analyze2(:,4)  );
disp('SNR correlation between electrode types, walking alone');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Shepherd's pi
[stats] = ScatterOutliers(SNR_analyze2(:,1),SNR_analyze2(:,2) );
disp('SNR (Shepherds pi) correlation between electrode types, standing');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(SNR_analyze2(:,3),SNR_analyze2(:,4) );
disp('SNR (Shepherds pi) correlation between electrode types, walking ');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);

% %Print to folder
% thetitle = [ 'Correlation data noise'];
% supertitle(thetitle);
% print(['O:\projects\jos_sync1\sync_walking\graphs\Correlations\', thetitle], '-dpng')
disp('-------------------------------------------------------------------------------------');
%% Stats #6b: Accelerometer data correlations
disp('Stats #6b: Accelerometer data correlations ------------------------------------------');
% Organized: W1E1 W1E2 W2E1 W2E2 W3E1 W3E2 

[rho, pval] = corr(Accdata_analyze(:,1), Accdata_analyze(:,2));
disp('Acc data correlation between electrode types, walking condition 1');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
[rho, pval] = corr(Accdata_analyze(:,3), Accdata_analyze(:,4));
disp('Acc data correlation between electrode types, walking condition 2');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shepherd's pi
[stats] = ScatterOutliers(Accdata_analyze(:,1), Accdata_analyze(:,2) );
disp('Acceleration noise (Shepherds pi) correlation between electrode types, standing');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(Accdata_analyze(:,3), Accdata_analyze(:,4));
disp('Acceleration noise (Shepherds pi) correlation between electrode types, walking ');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);

% %Print to folder
% thetitle = [ 'Correlation Acceleration walking patterns'];
% supertitle(thetitle);
% print(['O:\projects\jos_sync1\sync_walking\graphs\Correlations\', thetitle], '-dpng')

disp('-------------------------------------------------------------------------------------');
%% Stats #6c: Correlation between Acc and prestim data noise
disp('Stats #6c: Correlation between Acc and prestim data noise----------------------------');
% Organized: W1E1 W1E2 W2E1 W2E2 W3E1 W3E2 
% STD_singletrial1_allstim(i_sub,i_cond,i_type, trials)


[rho, pval] = corr(Accdata_analyze(:,1),prestim_std_erp_analyze(:,1));
disp('Correlation between Acc data and prestim data noise: standing, active elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
[rho, pval] = corr(Accdata_analyze(:,2),prestim_std_erp_analyze(:,2));
disp('Correlation between Acc data and prestim data noise: cond 1, passive elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);


[rho, pval] = corr(Accdata_analyze(:,3),prestim_std_erp_analyze(:,3));
disp('Correlation between Acc data and prestim data noise: cond 2, active elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
[rho, pval] = corr(Accdata_analyze(:,4),prestim_std_erp_analyze(:,4));
disp('Correlation between Acc data and prestim data noise: cond 2, passive elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shepherd's pi
[stats] = ScatterOutliers(Accdata_analyze(:,1),prestim_std_erp_analyze(:,1));
disp('Correlation (Shep) between Acc data and prestim data noise: standing, active elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(Accdata_analyze(:,2), prestim_std_erp_analyze(:,2));
disp('Correlation (Shep) between Acc data and prestim data noise: standing, passive elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);


%Shepherd's pi
[stats] = ScatterOutliers(Accdata_analyze(:,3),prestim_std_erp_analyze(:,3));
disp('Correlation (Shep) between Acc data and prestim data noise: walking, active elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(Accdata_analyze(:,4), prestim_std_erp_analyze(:,4));
disp('Correlation (Shep) between Acc data and prestim data noise: walking, passive elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);

disp('-------------------------------------------------------------------------------------');
%% Stats #6c2: Correlation between Acc and SNR
disp('Stats #6c2: Correlation between Acc and SNR----------------------------');
% Organized: W1E1 W1E2 W2E1 W2E2 W3E1 W3E2 
% STD_singletrial1_allstim(i_sub,i_cond,i_type, trials)


[rho, pval] = corr(Accdata_analyze(:,1),SNR_analyze2(:,1));
disp('Correlation between Acc data and SNR: standing, active elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
[rho, pval] = corr(Accdata_analyze(:,2),SNR_analyze2(:,2));
disp('Correlation between Acc data and SNR: cond 1, passive elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);


[rho, pval] = corr(Accdata_analyze(:,3),SNR_analyze2(:,3));
disp('Correlation between Acc data and SNR: cond 2, active elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);
[rho, pval] = corr(Accdata_analyze(:,4),SNR_analyze2(:,4));
disp('Correlation between Acc data and SNR: cond 2, passive elec');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shepherd's pi
[stats] = ScatterOutliers(Accdata_analyze(:,1),SNR_analyze2(:,1));
disp('Correlation (Shep) between Acc data anD SNR: standing, active elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(Accdata_analyze(:,2), SNR_analyze2(:,2));
disp('Correlation (Shep) between Acc data and SNR: standing, passive elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);


%Shepherd's pi
[stats] = ScatterOutliers(Accdata_analyze(:,3),SNR_analyze2(:,3));
disp('Correlation (Shep) between Acc data and SNR: walking, active elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(Accdata_analyze(:,4), SNR_analyze2(:,4));
disp('Correlation (Shep) between Acc data and SNR: walking, passive elec');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);

disp('-------------------------------------------------------------------------------------');
%% Stats #6d: P3 correlations
disp('Stats #6d: P3 correlations ----------------------------------------------------------');

[rho, pval] = corr(P3_analyze(:, 1),P3_analyze(:,2) );
disp('correlation in P3 deviant amplitude between active and passive electrodes: standing');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);

[rho, pval] = corr(P3_analyze(:, 3),P3_analyze(:,4) );
disp('correlation in P3 deviant amplitude between active and passive electrodes: walking');
disp([' r = ',  num2str(rho), '; p = ', num2str(pval)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Shepherd's pi
[stats] = ScatterOutliers(P3_analyze(:, 1),P3_analyze(:,2) );
disp('correlation (shep) in P3 deviant amplitude between active and passive electrodes: standing');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);
[stats] = ScatterOutliers(P3_analyze(:, 3),P3_analyze(:,4) );
disp('correlation (shep) in P3 deviant amplitude between active and passive electrodes: walking');
disp([' r = ',  num2str(stats(1)), '; p = ', num2str(stats(2))]);


%% SUPPLEMENTAL Single trial Baseline noise
disp('SUPPLEMENTAL STATS ----------------------------------------------------------------------')
disp('SUPPLEMENTAL Single trial Baseline noise ------------------------------------------------')
% 2x2x2 anova is done in spss, but here we have a 2x2 anova. 
% Here we will ignore stimulus, because the prestimulus noise analysis also
% ignores stimulus, and having 3 factors is kind of ridiculously
% complicated. So we will start by adding these together (above)

%ANOVA 
%F1: walkcond
%F2: electrode

measure = 'Baseline data noise';

sub_means_W1E1=  STDnoise_analyze(:,1);
sub_means_W1E2=  STDnoise_analyze(:,2);
sub_means_W2E1 = STDnoise_analyze(:,3);
sub_means_W2E2 = STDnoise_analyze(:,4);
% for SPSS anova + partia eta squared and GG correction, put the above
% matrix into SPSS
alldat = [sub_means_W1E1 sub_means_W1E2...
            sub_means_W2E1 sub_means_W2E2];
subs =     [1:length(SUBJ), 1:length(SUBJ), 1:length(SUBJ),1:length(SUBJ)]' ;
F1stim =   [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*1; 
           ones(length(SUBJ),1)*2; ones(length(SUBJ),1)*2]   ;
F2conds =  [ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2;
           ones(length(SUBJ),1)*1; ones(length(SUBJ),1)*2]   ;
Factnames = [{'Mobility'},{'Electrode'}]

% run the anova
stats = rm_anova2(alldat,subs,F1stim,F2conds,Factnames)
disp(measure)
disp(['Mobility: (F(1,17) = ', num2str(round(stats{2,5}*100)/100), '; p = ' num2str(stats{2,6}), '; etap2 = )'])
disp(['Electrode: (F(1,17) = ', num2str(round(stats{3,5}*100)/100), '; p = ' num2str(stats{3,6}), '; etap2 = )'])
disp(['Mobility * Electrode: (F(1,17) = ', num2str(round(stats{4,5}*100)/100), '; p = ' num2str(stats{4,6}), '; etap2 = )'])


%ttests!

%Standing: Active vs. Passive
[h p ci test] = ttest(sub_means_W1E1,sub_means_W1E2);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W1E2);
cD = computeCohen_d(sub_means_W1E1, sub_means_W1E2, 'paired');
disp('Standing: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Walking: Active vs. Passive
[h p ci test] = ttest(sub_means_W2E1,sub_means_W2E2);
Mdiff = mean(sub_means_W2E1)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W2E1, sub_means_W2E2, 'paired');
disp('Walking: Active vs. Passive')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])


%Active: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E1, sub_means_W2E1);
Mdiff = mean(sub_means_W1E1)-mean(sub_means_W2E1);
disp('Active: Standing vs. Walking')
cD = computeCohen_d(sub_means_W1E1, sub_means_W2E1, 'paired');
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')'])


%Passive: Standing vs. Walking
[h p ci test] = ttest(sub_means_W1E2, sub_means_W2E2);
Mdiff = mean(sub_means_W1E2)-mean(sub_means_W2E2);
cD = computeCohen_d(sub_means_W1E2, sub_means_W2E2, 'paired');
disp('Passive: Standing vs. Walking')
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), ' ; p = ', num2str(p), ...
      '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))  ')' ])

%Check direction
mean(sub_means_W1E1);
std(sub_means_W1E1);
mean(sub_means_W1E2);
std(sub_means_W1E2);

mean(sub_means_W2E1);
std(sub_means_W2E1);
mean(sub_means_W2E2);
std(sub_means_W2E2);

%QUICKPLOT
figure;
plot([1:2],[ mean(sub_means_W1E1) mean(sub_means_W2E1)], 'r'); hold on
plot([1:2],[ mean(sub_means_W1E2) mean(sub_means_W2E2)], 'b');
legend('Active', 'Passive')
xlabel(['Standing                                                Walking'])
title(measure)

disp('-------------------------------------------------------------------------------------');