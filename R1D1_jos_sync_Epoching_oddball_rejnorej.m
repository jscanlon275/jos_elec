
%% in this file:
% 12	Apply (highpass) filter as appropriate for research question (e.g. 0.3 Hz)
% 13    Epoch the data for all events of interest (e.g., -200 to 800 ms)
% 14    Remove baseline
% 15	Reject epochs with residual artifacts not accounted for by ICA
% 16	Create and save datasets for each condition of interest





close all; clear all;
trialrej = {'rej';'norej'};
whichtrials = 1; % 1 or 2 for rej or norej

MAINPATH = 'D:\4. home office\otto backup\jos_sync1\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATHIN = [MAINPATH, 'data\R1_ana2_ica_3cons\'];
PATHOUT = [MAINPATH, 'data\R1_ana3_epochs_', trialrej{whichtrials},'\'];
mkdir(PATHOUT);

% %clear the folder
% which_dir = PATHOUT
% dinfo = dir(which_dir);
% dinfo([dinfo.isdir]) = [];   %skip directories
% filenames = fullfile(which_dir, {dinfo.name});
% if ~isempty(filenames)
%     delete( filenames{:} )
% end


cd(MAINPATH);

 allsubs = {
    %  '001', 'act';...%%%crazy error
    %     '001', 'pas';... %missing segment
    '002', 'act';...%manual fix
    '002', 'pas';...
    %    '003', 'pas'; % never did the second recording
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
    %       '008', 'pas';
    '009', 'act';
    '009', 'pas';   %note: interpolate PO10 (?)
    %  '010', 'act'; %%% crazy error again...
    %      '010', 'pas';
    %        '011', 'act';
    %  '011', 'pas' %%% there is no EEG data in this file...
    '012', 'act';
    '012', 'pas';
    '013', 'act';
    '013', 'pas';
    '014', 'act';
    '014', 'pas';
    '015', 'act'; % Missing second EEG TTL: sync'd the accelerometers and clipped
    %the EEG to size using TTL at the end. no time warping in EEG, only
    '015', 'pas'
    %  '016', 'act'; % crazy error
    %       '016', 'pas';
    %  '017', 'act'; % new amp % crazy error again
    %       '017', 'pas'; %new amp
    '018', 'act'
    '018', 'pas'
    % '019' 'act' - % new amp % didn't work so we cancelled second session
    '020', 'act';
    '020', 'pas';
    '021', 'act';
    '021', 'pas'; % Manual fix: pruned out the extra TTL pulses
    '022', 'act';
    '022', 'pas';
    '023', 'act';
    '023', 'pas';
    '024', 'act';
    '024', 'pas';
    '025', 'act';
    '025', 'pas';
    '026', 'act';
    '026', 'pas'...
    };
 % Interpolations: 
  % 005 act - elec 21
  % 005 pas - elec 51
  % 008 pas - elec 02
  % 009 pas - elec 51 & 52 
  % 022 pas - elec 53
  
sub  = allsubs(:,1);
type = allsubs(:,2);

CONDS = {'standard_standing','deviant_standing',...
    'standard_walkingA','deviant_walkingA',...
    'standard_walkingT','deviant_walkingT'};
EVENTS = {'standard1', 'deviant1', 'standard2', 'deviant2', 'standard3', 'deviant3'};  % time-locking events
% EVENTS2 = {'ASReyesclose', 'ASReyesopen'};  % time-locking events %0 30000

EP_from = -.2;                      % epoch beginning
EP_to = .8;                         % epoch end
BASE_from = -200;                   % baseline begin
AMP = 100;                          % rejection threshold
MYPLOT = 1;                         % 1 = own plot, 0 = pop_comperp
REJ_ICA = 2;                        % pruning for ICA
%testing = 0;

% Extra dataset for interpolation
ALLCHANS = pop_loadset(['D:\4. home office\otto backup\jos_sync1\rawdata\synchronized\jos_sync002_act.set']);
ALLCHANS = pop_select(ALLCHANS,'channel', [1:64]);


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% Load the data


for i_sub = 1:length(sub)
%i_sub = 1

EEG = pop_loadset([PATHIN, 'jos_sync', sub{i_sub}, '_' type{i_sub}, '_ica_processed.set']);
% [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 1




%% Step 12: Filter according to research question
EEG = pop_eegfiltnew(EEG, [],0.3,5500,1,[],1);
EEG = pop_eegfiltnew(EEG, [],40,166,0,[],1);

%% Interpolating bad electrodes
EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical');
% if      strcmp(sub{i_sub}, '002') && strcmp(type{i_sub},'pas')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '005') && strcmp(type{i_sub},'act')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '005') && strcmp(type{i_sub},'pas')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '008') && strcmp(type{i_sub},'pas')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '009') && strcmp(type{i_sub},'pas')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '014') && strcmp(type{i_sub},'act')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '022') && strcmp(type{i_sub},'pas')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% elseif  strcmp(sub{i_sub}, '025') && strcmp(type{i_sub},'pas')
%     EEG = pop_interp(EEG, ALLCHANS.chanlocs, 'spherical')
% end

%% Re-reference
EEG = pop_reref(EEG, [10 21]);

%% Step 13:
EEG = pop_epoch(EEG, EVENTS, [EP_from EP_to]);

%% Step 14 & 15:
EEG = pop_rmbase(EEG, [BASE_from 0]);
EEG.setname = ['jos_sync', sub{i_sub},'_', type{i_sub} '_oddball_epoched'];
EEG.OGtrialnums = EEG.trials;
% rejecting trials?
if strcmp(PATHOUT, [MAINPATH, 'data\R1_ana3_epochs_rej\'])
EEG = pop_jointprob(EEG,1,[1:length(EEG.chanlocs)],3,3,0,1,0,[],0);
EEG.setname = [EEG.setname , '_rej'];
else
    EEG.setname = [EEG.setname , '_norej'];
end

EEG = pop_saveset(EEG, 'filename',[EEG.setname, '.set'], 'filepath',PATHOUT);
EEG = eeg_checkset(EEG);
 [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % ALLEEG 3
idx = length(ALLEEG);

%% Step 16:
for e = 1:length(EVENTS)
    EEG = pop_selectevent(ALLEEG(idx), 'latency','-2<=2','type',{EVENTS{e}},...
        'deleteevents','off','deleteepochs','on','invertepochs','off');
    EEG.setname = [EEG.setname, '_', CONDS{e}];
    EEG = pop_saveset(EEG, 'filename',[EEG.setname, '.set'],...
        'filepath',PATHOUT);
    EEG = eeg_checkset(EEG);
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);
end

end