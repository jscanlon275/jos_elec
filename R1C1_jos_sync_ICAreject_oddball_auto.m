

%in the next installment:
% 8 	Load previously saved dataset
% 9 	Apply ICA weights
% 10 	Identify components representing artifact (blinks, lateral eye movement)
% 11	Back-project all other (non artifactual) components

%% TODO
% epoch the eyes open/ eyes closed data




close all; clear all


MAINPATH = 'D:\4. home office\otto backup\jos_sync1\';
addpath([MAINPATH, 'script\functions\']);
addpath([MAINPATH, 'eeglab14_1_2b\']);

PATHIN = [MAINPATH, 'data\R1_ana1_3cons_noln\'];
PATHOUT = [MAINPATH, 'data\R1_ana2_ica_3cons\'];
PATHOUT_image = [MAINPATH, 'data\R1_ana2_ica_3cons\ICLabel_images\'];

mkdir(PATHOUT);
CONDS = {'standard_standing','deviant_standing',...
    'standard_walkingA','deviant_walkingA',...
    'standard_walkingT','deviant_walkingT'};

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
    '026', 'pas'
    };


sub  = allsubs(:,1);
type = allsubs(:,2);

EVENTS = {'standard1', 'deviant1', 'standard2', 'deviant2', 'standard3', 'deviant3'};  % time-locking events

EP_from = -.2;                      % epoch beginning
EP_to = .8;                         % epoch end
BASE_from = -200;                   % baseline begin
AMP = 100;                          % rejection threshold
MYPLOT = 1;                         % 1 = own plot, 0 = pop_comperp
REJ_ICA = 2;                        % pruning for ICA
% SPEED_UP = 1;                       % 1 = run ICA with PCA 20, 0 = no ICA
testing = 0;


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


for i_sub = 1:length(sub)
    %% Step 8: Load previously saved dataset
    EEG = pop_loadset([PATHIN, 'jos_sync', sub{i_sub}, '_' type{i_sub}, ' resampledICA_alltrials.set']);
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 1
    
    %% Step 9 and 10: Apply ICA weights and Remove ICA components
    
    % INFO:
    % This script draws a threshold at 80% for eyes and 40% for heart, and
    % allows a visual inspection of these components, before allowing a yes/no decision on each removable component
    % It also marks muscle components, but we don't worry about removing them
    % here.
    
    %  Component removed if:
    % 1. It's highest label is Eye >85%, or Heart >40%
    % and
    % 2. It has a clear and dipolar map
    % 3. If these things conflict (i.e. good map with lower percentage chance of being eye or heart, or high percentage eye or heart with a close-but-not-so-clear map), check components in the component scroll.
    
    close all;
    clear inspect;
    
    pop_topoplot(EEG,0, [1:size(EEG.icawinv,2)],EEG.setname,[7 10] ,0,'electrodes','on');
    pop_eegplot(EEG, 0, 1, 1);
    EEG = iclabel(EEG);
    n = 1;
    figure(1)
    for i_comp = 1:size(EEG.icawinv, 2)
        subplot( 7, 10, i_comp);
        [perc, pos] = max(EEG.etc.ic_classification.ICLabel.classifications(i_comp, :));
        
        if strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Heart') || strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Eye')|| strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Muscle')
            ticolour = 'r';
            if strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Eye') && perc>0.85
                inspect(n) = i_comp;
                n = n+1;
            elseif strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Heart') && perc>0.4
                inspect(n) = i_comp;
                n = n+1;
                %              elseif strcmp(num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), 'Muscle') && perc>0.5
                %              inspect(n) = i_comp;
                %              n = n+1;
            end
        else
            ticolour = 'k';
        end
        
        title([num2str(i_comp), ': ', num2str(EEG.etc.ic_classification.ICLabel.classes{pos}), ' ', num2str(perc) ], 'Color', ticolour)%highest comp+ perc.;
    end
    
    figure
    %
    b = 1;
    for i_insp = 1:length(inspect);
        pop_prop_extended(EEG, 0, inspect(i_insp)); % for component properties
        
        %      remove = input(['Remove Component ', num2str(inspect(i_insp)), '? (1/0)']);
        %      if remove == 1;
        EEG.badcomps(b) = inspect(i_insp);
        b = b+1;
        %      end
    end
    
    
    %  EEG.notes = input('Questionable comps? [] : ');
    figure(1)
    subplot(7,10, i_comp+1);
    axis off
    title([ {' ' ; 'Removed:';'Questionable:'} ], 'Interpreter', 'none');
    subplot(7,10, i_comp+2);
    axis off
    % title([{EEG.setname ;  num2str(EEG.badcomps);  num2str(EEG.notes)} ], 'Interpreter', 'none');
    title([{EEG.setname ;  num2str(EEG.badcomps)} ], 'Interpreter', 'none');
    
    %% Save the png file to record which components were removed and what they looked like
    print([PATHOUT_image, 'jos_sync' sub{i_sub}  '_' type{i_sub}, '_icalabels'], '-dpng');
    close all
    
    %% Step 11: Back-project all other (non artifactual) components
    EEG = pop_subcomp(EEG, EEG.badcomps, 0);
    EEG.setname = [EEG.setname, '_icacorrected'];
    EEG = eeg_checkset(EEG);
%     [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG);%ALLEEG 2
    
    
    %% Save the data
    EEG.setname = ['jos_sync', sub{i_sub}, '_', type{i_sub}, ]
    EEG = pop_saveset(EEG, 'filename',[EEG.setname, '_ica_processed.set'], 'filepath',PATHOUT);
    
end