% based on: std_exercise_8.m
% Oct 20, 2019
%
%
% INFO:
% This script is based on implementation of the "Preprocessing steps - advanced"
%    pipeline, as described on EEGLAB_03_slides.pdf slide 36:
%
% 1 	Import raw data
% 2 	Add channel location (and event) information
% 3 	Save dataset
% 4		Apply 1 Hz highpass filter
% 5		Epoch data into consecutive 1 sec segments
% 6		Reject epochs with atypical artifacts
% 7		Decompose data with independent component analysis (ICA)
% Save image and data then loop to next subject

clear all; close all; clc;
MAINPATH = 'D:\4. home office\otto backup\jos_sync1\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATHIN = [MAINPATH, 'rawdata\synchronized\'];
PATHOUT = [MAINPATH, 'data\R1_ana4_processed\'];
mkdir(PATHOUT);
CONDS = {'standard_standing','deviant_standing',...
    'standard_walkingA','deviant_walkingA',...
    'standard_walkingT','deviant_walkingT'};


cd(MAINPATH);

allsubs = {
          %  '001', 'act';...%%%crazy error
          %  '001', 'pas';...
            '002', 'act';...%manual fix
             '002', 'pas';...
         %   '003', 'pas'; % never did the second recording
            '004', 'act';
            '004', 'pas';
            '005', 'act'; %% Manual fix: Missing first EEG TTL. sync'd the accelerometers and clipped
          %  the EEG to size using TTL at the end. no time warping in EEG, only
          %  in acc. **interpolate elec 21 (?)
             '005', 'pas';  % check if this needs to be interpolated
           '006', 'act';
           '006', 'pas'; % manual fix 'assumed that an incomplete TTL was the right signal' - appears to work.
           '007', 'act';
           '007', 'pas';
          %  '008', 'act';  %%%didn't work at all...crazy error
        %   '008', 'pas';
           '009', 'act';
            '009', 'pas';   %note: interpolate PO10 (?)
          %  '010', 'act'; %%% crazy error again...
        %   '010', 'pas';
        %   '011', 'act';
          %  '011', 'pas' %%% there is no EEG data in this file...
           '012', 'act';
           '012', 'pas';
           '013', 'pas';
           '013', 'act';
            '014', 'act';
           '014', 'pas';
           '015', 'act'; % Missing second EEG TTL: sync'd the accelerometers and clipped
           %the EEG to size using TTL at the end. no time warping in EEG, only
           '015', 'pas'
          %  '016', 'act'; % crazy error
          % '016', 'pas';
          %  '017', 'act'; % new amp % crazy error again
         %  '017', 'pas'; %new amp
           '018', 'pas'
           '018', 'act'
          % '019' 'act' - % new amp % didn't work so we cancelled second session
           '020', 'act';
           '020', 'pas';
           '021', 'pas' % Manual fix: pruned out the extra TTL pulses
           '021', 'act';
           '022', 'act';
           '022', 'pas';
           '023', 'act';
           '023', 'pas';
           '024', 'act';
           '024', 'pas';
           '025', 'pas';
           '025', 'act';
           '026', 'act';
           '026', 'pas'
    };

% pairedsubs = {'002', '004', '005', '006', '007', '009', '012', '013',
% '014', '015', '018', '020', '021', '022', '023', '024', '025', '026'};

 % Interpolations: 
  % 005 act - elec 21
  % 005 pas - elec 51
  % 008 pas - elec 02
  % 009 pas - elec 51 & 52 
  % 014 act - elec 62
  % 022 pas - elec 53
  % 025 pas - elec 02

sub  = allsubs(:,1);
type = allsubs(:,2);
%type = {'act', 'pas'};

EVENTS = {'standard1', 'deviant1', 'standard2', 'deviant2', 'standard3', 'deviant3'};  % time-locking events

EP_from = -.2;                      % epoch beginning
EP_to = .8;                         % epoch end
BASE_from = -200;                   % baseline begin
AMP = 100;                          % rejection threshold
MYPLOT = 1;                         % 1 = own plot, 0 = pop_comperp
REJ_ICA = 2;                        % pruning for ICA
% SPEED_UP = 1;                       % 1 = run ICA with PCA 20, 0 = no ICA
testing = 0;

%%
for i_sub=1:length(sub)
    
    %i_sub = 1
   %% Step 1: Import raw data 
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = load([PATHIN,'jos_sync', sub{i_sub}, '_' type{i_sub} '.mat'])%
    EEG = EEG.EEG;
    ACC = EEG.data(65:end, :);
    EEG.Acclocs = EEG.chanlocs(65:end)
    EEG = pop_select(EEG,'channel', [1:64]);
    EEG.setname = ['jos_sync' sub{i_sub}  '_' type{i_sub}];
    
%     %channel removing (due to excessive noise likely caused by bad connection)
%     % see file 'B_channel-reject maps2.ppt' (in graphs folder) for details.
%     if      strcmp(sub{i_sub}, '002') && strcmp(type{i_sub},'pas')
%         EEG = pop_select(EEG, 'nochannel', 63)
%     elseif  strcmp(sub{i_sub}, '005') && strcmp(type{i_sub},'act')
%         EEG = pop_select(EEG, 'nochannel', 21)
%     elseif  strcmp(sub{i_sub}, '005') && strcmp(type{i_sub},'pas')
%         EEG = pop_select(EEG, 'nochannel', 51)
%     elseif  strcmp(sub{i_sub}, '008') && strcmp(type{i_sub},'pas')
%         EEG = pop_select(EEG, 'nochannel', 02)
%     elseif  strcmp(sub{i_sub}, '009') && strcmp(type{i_sub},'pas')
%         EEG = pop_select(EEG, 'nochannel', 52)
%         EEG = pop_select(EEG, 'nochannel', 51)
%     elseif  strcmp(sub{i_sub}, '014') && strcmp(type{i_sub},'act')
%         EEG = pop_select(EEG, 'nochannel', 62)
%     elseif  strcmp(sub{i_sub}, '022') && strcmp(type{i_sub},'pas')
%         EEG = pop_select(EEG, 'nochannel', 53)
%     elseif  strcmp(sub{i_sub}, '025') && strcmp(type{i_sub},'pas')
%         EEG = pop_select(EEG, 'nochannel', 02)
%     end
%     
    
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 1
   
%     % Find the first startpoint and endpoint
%     midpoint = [find(strcmp({EEG.event.type}, 'ASReyesopen'))  find(strcmp({EEG.event.type}, 'ASReyesclose'))];
%     startpoint1 = min(midpoint);
%     midpoint = max(midpoint);
%     if strcmp({EEG.event(midpoint+1).type}, 'end')
%         endpoint1 = midpoint+1;
%     elseif strcmp({EEG.event(midpoint+2).type}, 'end')
%         endpoint1 = midpoint+2;
%     end
%     
%     % find the second startpoint - the first countdown
%     startpoint =  [find(strcmp({EEG.event.type}, 'countdown_highdev'))];
%     if isempty(startpoint)
%         startpoint =  [find(strcmp({EEG.event.type}, 'countdown_lowdev'))];
%     end
%     startpoint2 = min(startpoint);
%     
%     % find the last endpoint- the first WalkingInstruct
%     if testing == 1; % shorten the file so ICA doesn't take forever during testing
%         endpoint =  [find(strcmp({EEG.event.type}, 'countdown_highdev'))];
%         if isempty(endpoint)
%             endpoint =  [find(strcmp({EEG.event.type}, 'countdown_lowdev'))];
%         end
%         endpoint2 = endpoint(4);
%     else
%         endpoint =  [find(strcmp({EEG.event.type}, 'WalkingInstruct1')) , find(strcmp({EEG.event.type}, 'WalkingInstruct2'))];
%         endpoint2 = min(endpoint);
%     end
%     % [startpoint1 endpoint1 startpoint2 endpoint2]
%     EEG.ACC = [ACC(:, EEG.event(startpoint1).latency: EEG.event(endpoint1).latency) ACC(:, EEG.event(startpoint2).latency:EEG.event(endpoint2).latency)];
%     EEG = pop_select(EEG, 'point', [EEG.event(startpoint1).latency EEG.event(endpoint1).latency; EEG.event(startpoint2).latency EEG.event(endpoint2).latency]);
%     %    Channel labels already included:
%     %     TMP = pop_chanedit(TMP, 'lookup', [MAINPATH, ...
%     %         'eeglab14_1_2b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp']);
%    
    %% Count original trial numbers
    %% Step 16:
    % EVENTS = {'standard1', 'deviant1', 'standard2', 'deviant2', 'standard3', 'deviant3'};  % time-locking events
for e = 1:length(EVENTS)

    if type{i_sub} == 'act'; 
    
       OG_Trialnums_act(i_sub, e) =  length(find(strcmp({EEG.event.type}, EVENTS(e))));
    else
       OG_Trialnums_pas(i_sub, e) =  length(find(strcmp({EEG.event.type}, EVENTS(e))));

    end
end
% these are organized like this: {'standard1', 'deviant1', 'standard2', 'deviant2', 'standard3', 'deviant3'};  % time-locking events
    if type{i_sub} == 'act'; 
OG_Trialnums_act(OG_Trialnums_act==0) = nan;
    save([PATHOUT, 'OGtrials_act.mat'], 'OG_Trialnums_act')

    else
OG_Trialnums_pas(OG_Trialnums_pas==0) = nan;
    save([PATHOUT, 'OGtrials_pas.mat'], 'OG_Trialnums_pas')

    end
    

 
   %%  
%     % start using TMP here
%     % filter the data
%     TMP = pop_eegfiltnew(EEG,[],1,1650,1,[],0);
%     TMP = pop_eegfiltnew(TMP, [],40,166,0,[],0);
%     %%
%     TMP = eeg_regepochs(TMP);
%     TMP = pop_jointprob(TMP,1,[1:length(TMP.chanlocs)],REJ_ICA,REJ_ICA,0,1,0,[],0);
%     
%     
%     TMP = pop_runica(TMP, 'extended',1,'interupt','on');
%     
%     %% Step 8 & 9:
%     EEG.icawinv = TMP.icawinv;
%     EEG.icasphere = TMP.icasphere;
%     EEG.icaweights = TMP.icaweights;
%     EEG.icachansind = TMP.icachansind;
%     clear TMP;
%     EEG = eeg_checkset(EEG);
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 2
%     
%     %% Step 10:
%     pop_topoplot(EEG,0, [1:32],EEG.setname,[6 6] ,0,'electrodes','on');
%     print([PATHOUT, 'jos_sync' sub{i_sub}  '_' type{i_sub}, '_icaweights'], '-dpng')
%     close
%     %size(EEG.icawinv,2)
%     pop_topoplot(EEG,0, [33:length(EEG.chanlocs)],EEG.setname,[6 6] ,0,'electrodes','on');   %size(EEG.icawinv,2)
%     print([PATHOUT, 'jos_sync' sub{i_sub}  '_' type{i_sub}, '_icaweights2'], '-dpng')
%     close
%     pop_eegplot(EEG, 0, 1, 1);
%     close
%     %save
%     EEG = pop_saveset(EEG, 'filename',[EEG.setname, '.set'],...
%         'filepath',PATHOUT);
%     %eeg_eventtypes(EEG)
end

eeglab redraw

