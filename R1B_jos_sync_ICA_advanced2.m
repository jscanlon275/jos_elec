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
PATHOUT = [MAINPATH, 'data\R1_ana1_3cons_noln\'];
mkdir(PATHOUT);
CONDS = {'standard_standing','deviant_standing',...
    'standard_walkingA','deviant_walkingA',...
    'standard_walkingT','deviant_walkingT'};

cd(MAINPATH);


 SUBJ = { '002'...
    , '004', '005', '006', '007', '009', '012', '013', '014',...
         '015', '018', '020', '021', '022', '023', '024', '025', '026'...
         };

%full pairs: {'002', '004', '005', '006', '007', '009', '012', '013',
%'014', '015', '018', '020', '021', '022', '023', '024', '025', '026'}


TYPE = {'act', 'pas'};

EVENTS = {'standard1', 'deviant1', 'standard2', 'deviant2', 'standard3', 'deviant3'};  % time-locking events
CONDS = {'odd Standing1','odd Standing2', 'odd WalkingA1','odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2'};

EP_from = -.2;                      % epoch beginning
EP_to = .8;                         % epoch end
BASE_from = -200;                   % baseline begin
AMP = 100;                          % rejection threshold
MYPLOT = 1;                         % 1 = own plot, 0 = pop_comperp
REJ_ICA = 2;                        % pruning for ICA
% SPEED_UP = 1;                       % 1 = run ICA with PCA 20, 0 = no ICA
testing = 0;

[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% Extra dataset for interpolation
ALLCHANS = pop_loadset(['D:\4. home office\otto backup\jos_sync1\rawdata\synchronized\jos_sync002_act.set']);

%%
for i_sub=1:length(SUBJ)
  for  i_type = 1:length(TYPE)
    %i_sub = 1
    %% Step 1: Import raw data
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    EEG = load([PATHIN,'jos_sync', SUBJ{i_sub}, '_' TYPE{i_type} '.mat'])%
    EEG = EEG.EEG;
    ACC = EEG.data(65:end, :);
    EEG.Acclocs = EEG.chanlocs(65:end)
    EEG = pop_select(EEG,'channel', [1:64]);
    EEG.setname = ['jos_sync' SUBJ{i_sub}  '_' TYPE{i_type}];
    
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
    
    % filter the data
    EEG = pop_eegfiltnew(EEG,[],1,1650,1,[],0);
    EEG = pop_eegfiltnew(EEG, [],40,166,0,[],0);
    EEG = pop_resample(EEG, 250) %downsample
    
    %% bad channel rejection: non default parameters
    nbChan = EEG.nbchan;
    %  EEG = clean_rawdata(EEG, 5, -1, 0.8, 4, -1, -1); % 1. default
    %  EEG = clean_rawdata(EEG, 10, -1, .5, 7, -1, -1); % 2. conservative
    %  with line noise
       EEG = clean_rawdata(EEG, 5, -1, .5, -1, -1, -1); % 3. cons. no line noise (doesn't need)

    % store number of removed channels
    numSubj.rmChan(i_sub) = nbChan-EEG.nbchan;
    numSubj.sub(i_sub) = str2num(SUBJ{i_sub});
    rm =1; removed = [];
    for i_chan = 1:length(ALLCHANS.chanlocs(1:64))
        rmelec = find(strcmp({EEG.chanlocs.labels}, {ALLCHANS.chanlocs(i_chan).labels}));
        if isempty(rmelec)
            removed{rm}= {ALLCHANS.chanlocs(i_chan).labels};
            rm = rm+1;
        end
        clear rmelec;
    end
    EEG.rmChan = numSubj.rmChan(i_sub);
    EEG.rmChanloc = removed;
    %save this to a file
    channels_removed{i_sub, i_type} =  [{EEG.rmChan}, EEG.rmChanloc ]
      chremove(i_sub, i_type) = EEG.rmChan;
        chansremov{i_sub, i_type} = EEG.rmChanloc;
        
            save([PATHOUT, 'chremove'],'chremove')
            save([PATHOUT, 'chansremov'],'chansremov')


    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 1
    %
    
    
    %% Sorting data: Make blocks
    % note: this segment is super long and I hate it. next time just include
    % these triggers in the data like a smart person.
    clear blocks
    stepcountexp = 1;
    stepcountpar = 1;
    blockcount =1; use_end = 0;
    latency_add = 3270; % approximate time of the countdown & button press
    for i_event = 1:length( EEG.event)
        
        eventadd = 0;
        
        if strcmp(EEG.event(i_event).type, 'Standing')
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'end_odd';
            blockcount = blockcount+1;
            if strcmp(EEG.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEG.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEG.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd Standing1')))
                blocks(blockcount).block = 'odd Standing1'
            else
                blocks(blockcount).block = 'odd Standing2'
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingA');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'end_odd';
            blockcount = blockcount+1;
            if strcmp(EEG.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEG.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEG.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd WalkingA1')))
                blocks(blockcount).block = 'odd WalkingA1';
            else
                blocks(blockcount).block = 'odd WalkingA2';
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingT');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'end_odd';
            blockcount = blockcount+1;
            if strcmp(EEG.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEG.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEG.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd WalkingT1')));
                blocks(blockcount).block = 'odd WalkingT1';
            else
                blocks(blockcount).block = 'odd WalkingT2';
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingInstruct1');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        
        if strcmp(EEG.event(i_event).type, 'WalkingInstruct2');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEG.event(i_event).type, 'WalkingInstruct3');
            blocks(blockcount).lat = EEG.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
    end
    %%
    for i_cond = 1:length(CONDS)
        START = blocks(strcmp({blocks.block}, [ CONDS{i_cond}])).lat;
        TO = blocks(find(strcmp({blocks.block}, [ CONDS{i_cond}]))+1).lat;
        %         TMP = pop_select(Acc, 'point', [START TO]);
        %         TMP2 = pop_select(Acc2, 'point', [START TO]);
        %         TMP.data = double(TMP.data);
        %         TMP2.data = double(TMP2.data);
        
        %now let's add these points as events
        EEG.event(end+1).type = ['start_',CONDS{i_cond}];
        EEG.event(end).latency = START;
        EEG.event(end).duration = .5;
        startpoint(i_cond) = START;
        EEG.event(end+1).type = ['end_',CONDS{i_cond}];
        EEG.event(end).latency = TO;
        EEG.event(end).duration = .5;
        endpoint(i_cond) = TO;
        
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');
    
    %CONDS = {'odd Standing1','odd Standing2', 'odd WalkingA1','odd WalkingA2', 'odd WalkingT1', 'odd WalkingT2'};
    %%% START HERE. SORT start and endpoints so they are in chronological order.
    %also check if all sections have identifzing triggers
    startpoint_sort = sort(startpoint);
    endpoint_sort = sort(endpoint);
    
    
    EEG.ACC = [ACC(:, startpoint_sort(1) : endpoint_sort(1)), ...
        ACC(:, startpoint_sort(2) : endpoint_sort(2)),...
        ACC(:, startpoint_sort(3):endpoint_sort(3)),...
        ACC(:, startpoint_sort(4): endpoint_sort(4)), ...
        ACC(:, startpoint_sort(5):endpoint_sort(5)),...
        ACC(:, startpoint_sort(6):endpoint_sort(6)) ];
    
    EEG = pop_select(EEG, 'point', [ startpoint_sort(1) endpoint_sort(1);startpoint_sort(2) endpoint_sort(2); startpoint_sort(3) endpoint_sort(3); ...
        startpoint_sort(4) endpoint_sort(4);startpoint_sort(5) endpoint_sort(5);  startpoint_sort(6) endpoint_sort(6)]);
    
    %%  start using TMP here
    TMP = eeg_regepochs(EEG);
    TMP = pop_jointprob(TMP,1,[1:length(TMP.chanlocs)],REJ_ICA,REJ_ICA,0,1,0,[],0);
    
    
    TMP = pop_runica(TMP, 'extended',1,'interupt','on');
    
    %% Step 8 & 9:
    EEG.icawinv = TMP.icawinv;
    EEG.icasphere = TMP.icasphere;
    EEG.icaweights = TMP.icaweights;
    EEG.icachansind = TMP.icachansind;
    clear TMP;
    EEG = eeg_checkset(EEG);
    [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 2
    
    %% Step 10:
    pop_topoplot(EEG,0, [1:32],EEG.setname,[6 6] ,0,'electrodes','on');
    print([PATHOUT, 'jos_sync' SUBJ{i_sub}  '_' TYPE{i_type}, '_icaweights'], '-dpng')
    close
    %size(EEG.icawinv,2)
    pop_topoplot(EEG,0, [33:length(EEG.chanlocs)],EEG.setname,[6 6] ,0,'electrodes','on');   %size(EEG.icawinv,2)
    print([PATHOUT, 'jos_sync' SUBJ{i_sub}  '_' TYPE{i_type}, '_icaweights2'], '-dpng')
    close
    pop_eegplot(EEG, 0, 1, 1);
    close
    %save
    EEG = pop_saveset(EEG, 'filename',[EEG.setname, 'ICA_alltrials.set'],...
        'filepath',PATHOUT);
    %eeg_eventtypes(EEG)
    save([PATHOUT, 'channels_removed'],'channels_removed')
end
end
save([PATHOUT, 'chremove'],'chremove')
save([PATHOUT, 'chansremov'],'chansremov')
eeglab redraw

