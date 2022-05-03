% Behaviour analysis



%% TODO
% 002_pas - had to start again, so triggers are mix. fix this.
% 007_act - 2 natural cnds
% normalization - between 0 and 333 (length of a standard half step)
% make a matrix of the stepping times so we dont have to re-load this
% every time!!
%%

clear all; close all; clc;
MAINPATH = 'D:\4. home office\otto backup\jos_sync2\';
addpath([MAINPATH, 'eeglab14_1_2b\']);
PATHIN = [MAINPATH, 'rawdata\ana1_A3_HS\'];
PATHOUT = [MAINPATH, 'Graphs\'];
mkdir(PATHOUT);


cd(MAINPATH);

%subs with non-int latencies: 002 act, 005 act, 015act

% ALLSUB = {
%      '002', 'act';...%manual fix.
%    '004', 'act';
%    '005', 'act'; %% Manual fix: Missing first EEG TTL. sync'd the accelerometers and clipped
% %      the EEG to size using TTL at the end. no time warping in EEG, only
% %      in acc. **interpolate elec 21 (?)
%     '006', 'act';
%   '007', 'act'; %fix: 2 walk nat segments
%     '009', 'act';
%     '011', 'act';
%     '012', 'act';
%     '013', 'act';
%     '014', 'act';
%     '015', 'act'; % Missing second EEG TTL: sync'd the accelerometers and clipped
% %     the EEG to size using TTL at the end. no time warping in EEG, only
%      '018', 'act' ;
%     '020', 'act';
%     '021', 'act';
%     '022', 'act';
%     '023', 'act';
%     '024', 'act';
%     '025', 'act';
%     '026', 'act'
%     
%         '001', 'pas';...   % Says 'EEG resampled'... fix
%    '002', 'pas';... %CHECK- starting again in walking cond may change something
%         '003', 'pas'; % never did the second recording
%         '004', 'pas';
%         '005', 'pas';  % check if this needs to be interpolated
%         '006', 'pas'; % manual fix 'assumed that an incomplete TTL was the right signal' - appears to work.
%         '007', 'pas';
%         '008', 'pas';
%         '009', 'pas';   %note: interpolate PO10 (?)
%         '010', 'pas';
%         '012', 'pas';
%         '013', 'pas';
%         '014', 'pas';
%         '015', 'pas'
%         '016', 'pas';
%         '017', 'pas'; %new amp
%         '018', 'pas'
%         '020', 'pas';
%         '021', 'pas' % Manual fix: pruned out the extra TTL pulses
%         '022', 'pas';
%         '023', 'pas';
%         '024', 'pas';
%         '025', 'pas';
%         '026', 'pas';
%     };
% 
% sub  = ALLSUB(:,1);
% type = ALLSUB(:,2);

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

%full pairs: {'002', '004', '005', '006', '007', '009', '012', '013',
%'014', '015', '018', '020', '021', '022', '023', '024', '025', '026'}



%%

% for
for i_sub = 1:length(SUBJ)
     for i_type = 1:length(ETYPE)
    % Step 1: Import raw data jos_sync002_actallAcc
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    FILENAME = ['jos_sync', SUBJ{i_sub}, '_' ETYPE{i_type}];
    load([PATHIN,FILENAME, '_steps.mat']);
    
    EEGstep = EEG;
    
    
    %% pull out and organize the important events
    clear blocks expstep parstep
    
    stepcountexp = 1;
    stepcountpar = 1;
    blockcount =1;
    latency_add = 3270; % approximate time of the countdown & button press
    for i_event = 1:length( EEGstep.event)
        
        eventadd = 0;
        if strcmp(EEGstep.event(i_event).type, 'end');%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency;
            blocks(blockcount).block = 'break';
            blockcount = blockcount +1;
        end
        if length(EEGstep.event(i_event).type) <5;
        % do nothing
         elseif strcmp(EEGstep.event(i_event).type(1:5), 'ExpHS')
            expstep(stepcountexp) = EEGstep.event(i_event).latency;
            stepcountexp = stepcountexp +1;
        elseif strcmp(EEGstep.event(i_event).type(1:5), 'ParHS')
            parstep(stepcountpar) = EEGstep.event(i_event).latency;
            stepcountpar = stepcountpar +1;
        end
        
        if strcmp(EEGstep.event(i_event).type, 'Standing')
            if strcmp(EEGstep.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEGstep.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEGstep.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd Standing1')))
                blocks(blockcount).block = 'odd Standing1'
            else
                blocks(blockcount).block = 'odd Standing2'
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEGstep.event(i_event).type, 'WalkingA');
            if strcmp(EEGstep.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEGstep.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEGstep.event(i_event+eventadd).latency+latency_add;
            
            
            if isempty(find(strcmp({blocks.block}, 'odd WalkingA1')))
                blocks(blockcount).block = 'odd WalkingA1';
            else
                blocks(blockcount).block = 'odd WalkingA2';
            end
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'WalkingT');
            if strcmp(EEGstep.event(i_event+1).type, '12');
                eventadd = 1;
            elseif strcmp(EEGstep.event(i_event+2).type, '12');
                eventadd = 2;
            end
            blocks(blockcount).lat = EEGstep.event(i_event+eventadd).latency+latency_add;
            if isempty(find(strcmp({blocks.block}, 'odd WalkingT1')));
                blocks(blockcount).block = 'odd WalkingT1';
            else
                blocks(blockcount).block = 'odd WalkingT2';
            end
            blockcount = blockcount +1;
        end
        
        if strcmp(EEGstep.event(i_event).type, 'WalkingInstruct1');
            blocks(blockcount).lat = EEGstep.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'WalkingInstruct2');
            blocks(blockcount).lat = EEGstep.event(i_event).latency;
            blocks(blockcount).block = 'Walkinginstruct';
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'walk1_start');%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk natural';
            blockcount = blockcount +1;
        end
        
        if strcmp(EEGstep.event(i_event).type, 'walk2_start') ;%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk control';
            blockcount = blockcount +1;
        end
        if strcmp(EEGstep.event(i_event).type, 'walk3_start');%make absolute sure that this is right.
            blocks(blockcount).lat = EEGstep.event(i_event).latency+latency_add;
            blocks(blockcount).block = 'walk sync';
            blockcount = blockcount +1;
        end
        
    end
    
    %% now organize into blocks
    clear walkingt1 walkingt2 walkingnat walkingcon walkingsyn pwalkingt1 pwalkingt2 pwalkingnat pwalkingcon pwalkingsyn;
  
    
    %% par steps
    wt1=1;
    wt2=1;
    wa1=1;
    wa2=1;
    wn =1;
    wc =1;
    ws =1;
    for i_step = 1:length(parstep) % remember for the participant we need the walking alone trial
        %   order =  find(strcmp({blocks.block}, 'odd WalkingT1'))
%         if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'odd WalkingT1'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'odd WalkingT1'))+1).lat;
%             pwalkingt1(wt1) = parstep(i_step);
%             wt1 = wt1+1;
%         end
%         if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'odd WalkingT2'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'odd WalkingT2'))+1).lat;
%             pwalkingt2(wt2) = parstep(i_step);
%             wt2 = wt2+1;
%         end
        if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'odd WalkingA1'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'odd WalkingA1'))+1).lat;
            pwalkinga1(wa1) = parstep(i_step);
            wa1 = wa1+1;
        end
        if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'odd WalkingA2'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'odd WalkingA2'))+1).lat;
            pwalkinga2(wa2) = parstep(i_step);
            wa2 = wa2+1;
        end
%         if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk natural'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk natural'))+1).lat;
%             pwalkingnat(wn) = parstep(i_step);
%             wn = wn +1;
%         end
%         if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk control'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk control'))+1).lat;
%             pwalkingcon(wc) = parstep(i_step);
%             wc = wc+1;
%         end
%         if parstep(i_step) > blocks(find(strcmp({blocks.block}, 'walk sync'))).lat && parstep(i_step) < blocks(find(strcmp({blocks.block}, 'walk sync'))+1).lat;
%             pwalkingsyn(ws) = parstep(i_step);
%             ws = ws+1;
%         end
    end
     
%         if strcmp(ALLSUB(i_sub, 2),'act');
%             typ = 1;
%         elseif strcmp(ALLSUB(i_sub, 2),'pas');
%             typ = 2;
%         end

% get the number of paces for each block
    numpaces(i_sub, i_type, 1) = wa1;
    numpaces(i_sub, i_type, 2) = wa2;
    
% get the time length for each block
   walktime(i_sub, i_type, 1) = blocks(find(strcmp({blocks.block}, 'odd WalkingA1'))).lat - blocks(find(strcmp({blocks.block}, 'odd WalkingA1'))+1).lat;
   walktime(i_sub, i_type, 2) = blocks(find(strcmp({blocks.block}, 'odd WalkingA2'))).lat - blocks(find(strcmp({blocks.block}, 'odd WalkingA2'))+1).lat;

% add the blocks together
   allpaces(i_sub, i_type) = numpaces(i_sub, i_type, 1)+numpaces(i_sub, i_type, 2);
   allwalktime(i_sub, i_type) = abs(walktime(i_sub, i_type, 1))./EEG.srate + abs(walktime(i_sub, i_type, 2))./EEG.srate;
   
% get paces per minute
   pacepersec(i_sub, i_type) = allpaces(i_sub, i_type)./(allwalktime(i_sub, i_type)./60);

    %% Question 1: how fast did the participant walk? 
     end
end

[h p ci test] = ttest(pacepersec(:, 1), pacepersec(:, 2));
Mdiff = mean(pacepersec(:, 1))-mean( pacepersec(:, 2));
disp('Paces per minute: Active vs. Passive')
cD = computeCohen_d( pacepersec(:, 1),  pacepersec(:, 2), 'paired');
disp(['Active: ', num2str(mean(pacepersec(:, 1))), '; Passive: ', num2str(mean( pacepersec(:, 2)))]);
disp(['(Mdiff = ', num2str(Mdiff), '; SDdiff = ',num2str(test.sd), '; t(' num2str(test.df) ') = ',num2str(test.tstat), '; p = ', num2str(p),  ...
    '; Cohen`s d = ', num2str(cD), '; CI95%: ', num2str(ci(1)), '-', num2str(ci(2))   ')']);
