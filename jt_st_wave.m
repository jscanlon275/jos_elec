% jt_st_wave() - wavelet transform 
% 
% USAGE:
% [tfr, foi] = st_wave(data, fsample, foi, qmin, qmax, times)
%
% IN:
%  data     	data vector
%  fsample  	sampling rate, e.g., EEG.srate    
%  foi      	frequencies of interest
%  qmin         minimum wavelet width at min frequency
%  qmax         maximum wavelet width at max frequency.
%				For qmax > qmin, wavelet width is being increased
%				for higher frequencies.	Set min=max for running 
%				standard wavelet transform
%  times	    time stamps vector, e.g., EEG.times
%
% OUT
%  trf        	complex time frequency array
%  foi	     	frequencies of interest
%
% st_wave() needs ms_wlt() markus siegel



function [tfr, fres, tres, foi] = jt_st_wave(data, fsample, foi, qmin, qmax, times) %#ok<INUSD>

w = linspace(qmin,qmax,length(foi));  
data=data(:)';

wl = jt_ms_wlt(foi,fsample,w);  % create wavelet family
tfr = zeros(length(foi), length(data));  % hurry-up 

% wavelet transform
for ifr = 1:length(foi)
    tfr(ifr,:) = conv2(data,wl{ifr},'same');
%       tfr(ifr,1:ceil(length(wl{ifr})/2)) = NaN;    
%       tfr(ifr,end-ceil(length(wl{ifr})/2):end) = NaN; 
    fres(ifr) = foi(ifr)/w(ifr); %#ok<AGROW>
    tres(ifr) = 1/(2*pi*fres(ifr)); %#ok<AGROW>
end

% % magnitude (power) 
% trfb = abs(tfr).^2;
% 
% % baseline correction (percent change)
% bstart = find(round(times.*1000)./1000 == bstart);  
% bend = find(round(times.*1000)./1000 == bend);
% base = mean(tfrb(:,bstart:bend),2);
% base = repmat(base,1,size(tfrb,2));
% tfrb = (tfrb - base)/base;
%     [tfr,fres,tres,foi]= st_wave(EEG.databp(:,j), 256,1:1:35,5,5,EEG.times);

% figure;
% subplot(3,1,1);
%     imagesc(times,foi,tfr);
%     set(gca,'YDir','normal');colorbar;
% subplot(3,1,2);
%     imagesc(EEG.times,foi,abs(tfr).^2);
%     set(gca,'YDir','normal');colorbar;
%     axis tight;
%     set(gca, 'ZLimMode', manual, 'ZLim', [-2 2]);
% subplot(3,1,3);
%     plot(times,data);
%     set(gca,'YDir','reverse'); axis tight; 
% 
% figure; plot(fres,1:length(fres));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OR, for standard wavelets of same width:
% w = ones(1,length(foi)).*q;
%
% scale to power: �V.^2
