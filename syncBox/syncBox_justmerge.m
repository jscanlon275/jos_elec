function [allData, allAcc] = syncBox(allAcc, varargin)
% pp_sync() - synchronzize data recorded by LiveAmp64, a Smarting amps and Faros accelerometers
%           by finding pulses sent by the syncBox with a template matching procedure
%
%--------------------------
% Attention: beta version!
%--------------------------
%
%Usage: >> [allData, allAcc] = pp_sync(allAcc);
%
% Inputs:
%   allAcc       - a cell with EEGLAB structures of EEG & accelerometer data with fields
%       data        - continous, chan x frames,
%       srate       - sampling rate,
%       TTL.chan    - index of channel on which the TTL pulse was recorded,
%       TTL.signal['Smarting'|'LiveAmp'|'accelerometer'] - type of stream
%
% Optional Inputs (with default parameters):
%   numTTL      - number of TTL pulses, min. 2, default is 2. Please adapt this
%               parameter, if indicated!
%   minDist     - minimal distance between two pulses in s. Default is 1s.
%   plotit[0|1] - Flag, default is 1. Plots alignment of data in the
%               beginning and end of the recording
%   chopit[0|1] - delete times at which any device was connected to the syncBox, default is 1,               
%               (only possible if connection to the syncBox was sucessfully detected)
%   timeChop    - time in s that are chopped away after the device was (dis-)connected from the syncBox,
%               default is 5 -> should eliminate the DC offset created by
%               this event
%   minDistPlug - minimal distance between two connectin devices to the syncBox and disconnecting them again
%               in s,  default is 15s, only relevant for fast detection of sync-pulses
%   warpstreams[0|1]   - Flag, warp accelerometer stream to same length as EEG
%               stream, default is 1.
%
% Outputs:
%   allAcc  - Cell with structures of EEG & accelerometer data with added
%            markers for the sync-pulses. For Smarting and Faros data the
%            correlation with the template is  stored at event.correlation and if
%            detected, there are also events for when trhe deviced were
%            (dis-)connected to the syncBox
%   allData - EEGLAB structure of synchronized data. TTL chan of accelerometers will be deleted,
%            resampled to EEG frequency
%
% naj, Jan 2019
% 
% 
% Caveats:
%----------------------------------------------------------------------------------
% Failed Fast Detection: (dis-) Connectiong to syncBox was not properly recognized
% - Recording was not started before connectiong to the syncBox
% - Parameters for detecting the DC offset produced by connection to the
% - syncBox need to be adjusted (s. timeChop & minDistPlug as well as l. 150)
%
% Failed template matching: 
% - Signal is to noisy (check data before submitting to function)
%
% Failed detection of LiveAmp64 sync-trigger
% - One LiveAmp did not receive a sync trigger (seems to happen from time
%   to time)
% - sync-triggers arrived more than one sample apart from each other and
%   were thus not recognized as these


%% parameters
% set dafault parameters
numTTL  = 2;
minDist = 1;
plotit  = 1;
chopit  = 0;
timeChop = 5;
minDistPlug = 15;
warpstreams = 1;

%% change them if additional input was given
if nargin < 1
    error('Provide EEG and accelerometer data'); end
for i=1:2:length(varargin)
    if strcmpi(varargin{i},'numTTL')
        numTTL= varargin{i+1};
        if numTTL < 2
            error('Only data with min. 2 TTL pulses can be processed!');
        end
    elseif strcmpi(varargin{i},'minDist')
        minDist=varargin{i+1};
    elseif strcmpi(varargin{i},'plotit')
        plotit=varargin{i+1};
    elseif strcmpi(varargin{i},'chopit')
        chopit=varargin{i+1};
    elseif strcmpi(varargin{i},'timeChop')
        timeChop=varargin{i+1};
    elseif strcmpi(varargin{i},'minDistPlug')
        minDistPlug=varargin{i+1};
        elseif strcmpi(varargin{i},'warpstreams')
        warpstreams=varargin{i+1};
    end
end

% get folder where function is stored
PATHIN = mfilename('fullpath');
PATHIN = PATHIN(1:end-7);
% PATHIN = 'O:\projects\all_gait\Scripts\syncBox\';

%% check input and set sampling frequency for synchronized streams
numAcc = length(allAcc); sratetmp = nan(length(allAcc),1);
idxEEG = []; idxLA = [];

for s = 1:numAcc
    Acc = allAcc{s};
    % check input structures
    if ~isfield(Acc, 'TTL')
        error(['Structure ', num2str(2), ' does not contain the field "TTL"!']);end
    if ~isfield(Acc.TTL, 'chan') && ~strcmpi(Acc.TTL.signal, 'LiveAmp')
        error(['Structure ', num2str(2), ' does not contain the field "TTL.chan"!']);end
    if ~isfield(Acc.TTL, 'signal')
        error(['Structure ', num2str(2), ' does not contain the field "TTL.signal"!']);end
    
    % get all sampling rates
    sratetmp(s) = Acc.srate;
    
    % set sampling rate to EEG sampling rate if one of the inputs is EEG
    if strcmp(Acc.TTL. signal, 'Smarting') %if there are multiple EEG datasets the last one is chosen
        idxEEG(end+1) = s;
        srate = Acc.srate;
        idxSmart = s;
    elseif strcmp(Acc.TTL. signal, 'LiveAmp') %if there are multiple EEG datasets the last one is chosen
        idxEEG(end+1) = s;
        srate = Acc.srate;
        idxLA(end+1) = s; 
    end
end

% otherwise use lowest sampling rate
if ~exist('srate', 'var')
    srate = min(sratetmp);
end

if warpstreams 
    if length(idxEEG)>1
        warning('Data can only be warped if only one EEG data set is present. Skipping warping.')
        warpstreams = 0;
    end
end

%% find sync pulse of each stream
numChopStreams = 0; durACC = nan(1,numAcc);

for s = 1:numAcc % loop through all structures
    
    Acc = allAcc{s};
    skip = 0; % flag whether template matching will be skipped (for LiveAmp)
        
    % load template
    if strcmp(Acc.TTL.signal, 'Smarting')
        load([PATHIN, 'tmpEEG'],'tmpAcc');
    elseif strcmp(Acc.TTL.signal, 'LiveAmp')
        skip = 1;
    elseif strcmp(Acc.TTL.signal, 'accelerometer')
        load([PATHIN, 'tmpAcc'],'tmpAcc');
    else
        error('Your TTL.signal is not of type "Smarting" , "LiveAmp" or " accelerometer".');
    end
    
    % check whether sampling rate of template and signal match the determined srate
       prob = 0;   % 0 or 1, 1= problems with TTL detection, run frame-by-frame comparison od signal and template
    if ~skip
        prob = 0;   % 0 or 1, 1= problems with TTL detection, run frame-by-frame comparison od signal and template
        if Acc.srate ~= srate
            disp(['Resample ', Acc.TTL.signal,' signal...']);
            evalc('Acc = pop_resample(Acc, srate);');% resample template
        end
        
        if tmpAcc.srate ~= srate
            disp('Resample sync-pulse template...');
            evalc('tmpAcc = pop_resample(tmpAcc, srate);');% resample template
        end
        sigtmpAcc = tmpAcc.data;
        
        % get data from channel with TTL pulse
        sigAcc = Acc.data(Acc.TTL.chan,:);
        
%     %% fast detection
%         
%         % Adapt find peaks setting if this does not work well
%         [~,locs]= findpeaks(diff(diff(diff(sigAcc))),'Threshold', 600,'MinPeakDistance' ,Acc.srate*minDistPlug, 'NPeaks', numTTL*2);
%         if length(locs) ~= numTTL*2
%             prob = 1; % skip fast detection if number of plugged in times does not match the number of pulses that should be detected
%             warning(['Fast TTL detection of ',Acc.setname,' failed. Initiating slow detection.']);
%         else
%             numChopStreams = numChopStreams+1;
%         end
%         
%         if ~prob
%             % find when amp was connected to the box with 3rd derivative
%             disp('Compare signal and template when the sensor was connected to the syncBox')
%             disp('Detected matches: ');
%             
%             for e = 1:2:length(locs)
%                 Acc.event(end+1).type = 'plug in';
%                 Acc.event(end).latency = locs(e);
%                 
%                 Acc.event(end+1).type = 'plug out';
%                 Acc.event(end).latency = locs(e+1);
%             end
%             
%             % only look for TTL when amp was connected to box
%             c = 1;       % count TTL pulses
%             
%             plugIn = [Acc.event(strcmp({Acc.event.type},'plug in')).latency];
%             plugOut = [Acc.event(strcmp({Acc.event.type},'plug out')).latency];
%             
%             for e = 1:length(plugIn) % loop through all times box was plugged in
%                 from = plugIn(e);
%                 to = plugOut(e);
%                 
%                 % calculate correlation of template and signal at each frame
%                 iterations = length(sigAcc(from:to))-length(sigtmpAcc);
%                 rACC = nan(1, iterations);
%                 for i = 1:iterations
%                     idx = from+i;
%                     rACCtmp = corrcoef(sigAcc(idx:idx+length(sigtmpAcc)-1),sigtmpAcc);
%                     rACC(i) = rACCtmp(1,2);
%                 end
%                 
%                 [~, i] = max(rACC);% get highest correlation
%                 
%                 % look for steepest incline @+- 10 samples
%                 idx = from+i-10;
%                 [~, iI] = max(diff(sigAcc(idx:idx+20)));
%                 
%                 % save its latency
%                 Acc.event(end+1).latency = idx+iI;
%                 Acc.event(end).type = 'sync1';
%                 
%                 % compute its correlation with the template
%                 rACCtmp = corrcoef(sigAcc(Acc.event(end).latency:Acc.event(end).latency+length(sigtmpAcc)-1),sigtmpAcc);
%                 
%                 % display warning and continue with frame-by-frame comparisons if r<.8
%                 if rACCtmp(1,2) < .8
%                     if c == 1
%                         warning('The first time you plugged the amp into the box, no TTL pulse was found. Probably the recording was only started after the amp was plugged in.');
%                     else
%                         warning(['The ', num2str(c),'. time you plugged the amp into the box, no TTL pulse was found'])
%                     end
%                     prob = 1; % set flag to enter slow detection
%                     warning(['Fast TTL detection  ',Acc.setname,' failed. Initiating slow detection.']);
%                     break
%                     
%                 else
%                     % otherwise store correlation in Acc.event
%                     Acc.event(end).correlation  = rACCtmp(1,2);
%                 end
%                 
%                 % counter
%                 disp(num2str(c)); % display
%                 if c == numTTL
%                     disp('done!');
%                 end
%                 c = c+1;% update
%             end
%         end
%         clear c e i idx iI iterations locs plugIn plugOut from to
%        
%        %% slow detection 
%         % (if previous identification failed at least one pulse)
%         if prob
%             % delete all previously addes events:
%             Acc.event = [];
%             
%             % calculate correlation of template and signal at each frame
%             iterations = length(sigAcc)-length(sigtmpAcc);
%             rACC = nan(1, iterations);
%             disp('Correlating whole signal and template...')
%             for i = 1:iterations
%                 rACCtmp = corrcoef(sigAcc(1,i:i+length(sigtmpAcc)-1),sigtmpAcc);
%                 rACC(i) = rACCtmp(1,2);
%             end
%             
%             TTLlat = nan(1, numTTL);
%             [~, lat] = sort(rACC, 'descend');% get highest correlation
%             % look for steepest incline @+- 10 samples
%             idx = lat(1)-10;
%             [~, iI] = max(diff(sigAcc(idx:idx+20)));% get steepest incline
%             TTLlattmp = idx+iI;
%             
%             TTLlat(1) = TTLlattmp ;                     % save latency of one pulse with highest correlation
%             
%             lat(abs(lat-TTLlat(1))<minDist*Acc.srate)=0;% set all latencies of correlations under the minDist to 0
%             for i = 2:numTTL                            % get the requested number of pulses by looking at the highest frequencies that are further apart than the minDist
%                 lat = lat(lat~=0);                      % delete latencies that were set to 0 (line 83)
%                 idx = lat(1)-10;                        % get latency of the highest remaining correlation and look from 10 samples before...
%                 [~, iI] = max(diff(sigAcc(idx:idx+20)));% .. to 10 samples after for the steepest incline
%                 TTLlat(i) = idx+iI;
%                 lat(abs(lat-TTLlat(i))<minDist*Acc.srate)=0;
%             end
%             
%             % add events 
%             
%             for e = 1:length(TTLlat)
%                 Acc.event(end+1).latency = TTLlat(e); % save sorted latencies...
%                 Acc.event(end).type = 'sync1';
%                 Acc.event(end).correlation = rACC(TTLlat(e)); % ...and their correlation with the template
%             end
%             
%             disp('done')
%         end
%      
%     %% LiveAmp
%     elseif skip
%         % determine wehther thera are alerady triggers called sync, if
%         % present skip the next lines
%         if ~sum(strcmp({Acc.event(:).type}, 'sync'))>2
%             % find markers that are not more than 1 sample apart
%             idxSync= find(abs(diff([Acc.event.latency]))<=1);
%             
%             % rename them
%             for e = idxSync
%                 Acc.event(e).type = 'sync';
%             end
%           
%             % get latency between 'sync' triggers
%             latSync = diff([Acc.event(idxSync).latency]);
%             
%             % rename them to 'sync1' if the following 3 'sync' markers are approx. 2s apart
%             diffSync = latSync>2*Acc.srate-4 & latSync<2*Acc.srate+2;
%             for e = 1:length(diffSync)-3
%                 if diffSync(e:e+3)
%                     Acc.event(idxSync(e)).type = 'sync1';
%                 end
%             end
%         end
%         
%         % resample to set srate
%         evalc('Acc = pop_resample(Acc,srate);');
%     end
%     
%     % duration of recording (samples)
%     evLat = geteventlat(Acc, 'sync1');
%     durACC(s) = evLat(end)-evLat(1);
%     Acc.TTL.dur = durACC(s);
%     
%     evalc('Acc = eeg_checkset(Acc,"eventconsistency");');
%     allAcc{s} = Acc;
%     
%     if prob && chopit
%         chopit = 0;
%         warning('Cannot chop data since detction of connection to the syncBox failed!')
%     end
% 
% end
% 
% 
% %% find differences in recording duration of EEG and accelerometers
% % if exist('idxEEG', 'var')
% %     offset = durACC-durACC(idxEEG(1));
% %     % if no EEG data is available, offset to first signal is calculated
% % else
% %     offset = durACC-durACC(1);
% % end
% 
% 
% %% determine legth of recording
% % choose longest one do all data ist there (for warping)
%     [dur, ~] = max(durACC);
%     dur = dur+length(sigtmpAcc);
% 
% 
% %% prune
% for a = 1:numAcc
%     % cut signal from first to last (including the pulse itself) TTL pulse and delete TTL pulse channel
%     idx = find(strcmp({allAcc{a}.event.type}, 'sync1'),1); % index of fist sync event
%     tmp = allAcc{a}.event(idx).latency; % its latency
%     evalc('allAcc{a} = pop_select( allAcc{a},"point",[tmp-1 tmp+dur])'); % cut signal
% end
% 
% %% timewarp
% % use interp1 since resample and timewarp cannot deal with so small dirrerences in long recordings
% if warpstreams
%     disp('Warping timeseries');
%     if isempty(idxEEG)
%             dur = durACC(1);
%         else
%             dur = durACC(idxEEG);
%              evalc('allAcc{idxEEG} = pop_select( allAcc{idxEEG},"point",[0 dur+length(sigtmpAcc)-1;])'); % cut sEEG signal
%     end
%         
%     for s = setdiff(1:numAcc,idxEEG) % only warp non-EEG streams
%         % extract data
%         DATAtmp = double(allAcc{s}.data(:,1:durACC(s)+length(sigtmpAcc)))'; 
%         
%         % create new time vector
%         newLat = linspace(1,size(DATAtmp,1),dur+length(sigtmpAcc));
%         
%         % interpolate
%         % TO-DO: evaluate whether spline or linear is better!
%          allAcc{s}.data= interp1([1:durACC(s)+length(sigtmpAcc)], DATAtmp, newLat, 'spline')'; 
%         
%          % adjust latency of events one by one
%         for e = 1:length(allAcc{s}.event)
%             % find new sample to which old latency was mapped  (smallest abs. difference)
%             [~, sample] = min(abs(newLat-allAcc{s}.event(e).latency)); 
%             allAcc{s}.event(e).latency = sample;
%         end
%         
%         evalc('allAcc{s} = pop_select( allAcc{s},"point",[0  dur+length(sigtmpAcc);])'); % cut signal to EEG length
%     end
%     disp('done!');
% end
%     
% 
% %% plot
% if plotit
%     figure()
%     set(gcf, 'Position', [1 1 1500 800]);
%      leg = cell(1, numAcc);
%     for a = 1:numAcc
%         Acc = allAcc{a};
%         if ismember(a,idxLA)
%             lat = geteventlat(Acc, 'sync1');
%             
%             subplot(211)
%             plot([lat(1), lat(1)], [-1000 1000]); hold on
%             
%             subplot(212)
%             plot([lat(end), lat(end)],[-1000 1000]); hold on
%         else
%             subplot(211)
%             plot(Acc.data(Acc.TTL.chan,:)); hold on
%             
%             subplot(212)
%             plot(Acc.data(Acc.TTL.chan,:)); hold on
%         end
%         leg{a} = allAcc{a}.setname;
%     end
%     
%     subplot(211)
%     xlim([-10, length(sigtmpAcc)+10]);
%     xlabel('samples')
%     title({'Alignment of streams after synchronization';''; 'Beginning'})
%     
%     subplot(2, 1, 2)
%     xlim([length(allAcc{a}.data)-length(sigtmpAcc)-10, length(allAcc{a}.data+10)]);
%     xlabel('samples')
%     title('End')
%     
%     legend(leg,'location', [0.87 0.5 0.12 0.06]);
% end

%% merge
disp('Merging datastreams')
allAcctmp = cell(1,numAcc);
for a = 1:numAcc
    
    % delete TTL chan
    if ~ismember(a,idxEEG)
        evalc('allAcctmp{a} = pop_select( allAcc{a}, "nochannel", allAcc{a}.TTL.chan);');
    else
        allAcctmp{a}= allAcc{a};
    end
    
    % delete not shared event fields
    if ismember(a,idxLA) && isfield(allAcctmp{a}.event, 'channel')
        allAcctmp{a}.event= rmfield( allAcctmp{a}.event,{'channel','bvtime','bvmknum','code','urevent'});
    elseif isfield(allAcctmp{a}.event, 'correlation')
        allAcctmp{a}.event = rmfield( allAcctmp{a}.event,'correlation');
    end
    
    % add channel source
    [allAcctmp{a}.chanlocs(:).source] = deal(allAcctmp{a}.setname);
    
    % get fields in same order (to be able to concatinate them)
%     allAcctmp{a}.chanlocs =orderfields(allAcctmp{a}.chanlocs);
    
    % merge
    if a == 1
        allData = allAcctmp{a};
    else
         % participant 007 and the horzcat error..
%         if size(allAcctmp{a}.chanlocs) ~= [1 3];
%             allAcctmp{a}.chanlocs = allAcctmp{a}.chanlocs';
%             allAcctmp{a}.data = allAcctmp{a}.data';
%         end

        % add channel names
        allData.chanlocs = [allData.chanlocs, allAcctmp{a}.chanlocs];
                
        % append data
        allData.data = [allData.data; allAcctmp{a}.data];
        
        % append events
        allData.event = [allData.event, allAcctmp{a}.event];
    end
end

% put in correct order & create event.urevent
evalc('allData = eeg_checkset(allData,"eventconsistency");');
% create index
urevent = num2cell(1:length(allData.event));
% add
[allData.event(:).urevent] = deal(urevent{:});

% adapt channel number...
allData.nbchan = size(allData.data, 1);

% and setname
allData.setname = 'Merged data';
 disp('done!');

%% chopit: if signal and DC offset should be chopped out
if chopit
    disp('Chop signal');
        % get latencies of connection to syncBox
        latplugIn = geteventlat(allData, 'plug in');
        latplugOut = geteventlat(allData, 'plug out');
        ALLEEG = [];
        for i = 1:numTTL-1
            from = latplugOut(numChopStreams*i)+timeChop*allData.srate;
            to = latplugIn(numChopStreams*i-(numChopStreams-1))-timeChop*allData.srate;
            evalc('EEG = pop_select(allData,"point",[from to]);'); % cut signal
            evalc('[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG);');
        end
        
        % merge
        evalc('allData = pop_mergeset(ALLEEG, [1:length(ALLEEG)]);');
     disp('done!');
end
end