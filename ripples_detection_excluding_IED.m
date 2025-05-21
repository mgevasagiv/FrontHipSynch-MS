
function [ripples,ripples_stat, exc]=ripples_detection_excluding_IED(signal,BP,t,th,minDistance,minRippleDuration,maxRippleDuration,noise_ch,IED_ch,IED_onsets, Fs)
% Detecting ripples based on Stark et al. 2014 (Neuron) method.
% input:
% signal = normalized squared bandpassed LFP signal from the hippocampus
%          (assuming sampling rate of 500Hz)
% BP = raw bandpass signal
% t = corresponding time in sec
% th = threshold in stdev [onset/offset peak]
% minDistance = in sec
% minRippleDuration = min ripple duration in Sec
% maxRippleDuration = max ripple duration in Sec
% noise_ch = normalized squared bandpassed LFP signal from channels to exclude global artefacts identifed as rippples
% IED_ch = normalized squared bandpassed LFP signal from channel to exclude IEDs that were identifed as rippples
% Fs = sampling rate
%
% output: 
% ripples = output table with ripple-timing and other features (start,
% peak, end, amplitude, etc.)
% ripples_stat = statistics of excluded events
%
% Author: Itzik Norman 15/12/17

% modified by Kamin Kim & Carlos Carrasco, 2020, Feb
% modified by Maya GS, Dec 2023

if size(signal,2)<size(signal,1),signal=signal'; end
if size(noise_ch,2)<size(noise_ch,1),noise_ch=noise_ch'; end
if size(IED_ch,2)<size(IED_ch,1),IED_ch=IED_ch'; end

[pks,locs] = findpeaks(signal, 'MINPEAKHEIGHT', th(2));
ENV=abs(hilbert(BP)); ENV=ENV./nanmedian(ENV); ENV = 10*log10(ENV); % to calculate ripple amplitude in dB relative to median

% Noise rejection:
if ~isempty(noise_ch)
    [r,p]=corr(signal',noise_ch');
    fprintf('\n --- correlation with noise ch.: %.2f \n',r)
    rej_count=0;
    pre_rej_count=size(locs,2);    
    fprintf('\n => rejecting noise artefacts... \n')
    [~,noise_locs] = findpeaks(noise_ch, 'MINPEAKHEIGHT', th(2));
    % Ignore global electrical artifacts:
    exc.GN = []; 
    for i=1:numel(noise_locs)
        tmp=find(abs(locs-noise_locs(i))< (0.05*Fs)); % kk +/- 50ms
        exc.GN = [exc.GN locs(tmp)]; 
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events based on noise channel correlation \n',rej_count,pre_rej_count)
    ripples_stat.SWR_initial = pre_rej_count; 
    ripples_stat.GN_count = rej_count; %Added by CC
    ripples_stat.GN_percent = (rej_count/pre_rej_count) * 100;
end

% IED rejection (Gelinas et al. 2016 Nat. Med.):
if ~isempty(IED_ch)
    rej_count=0;
    pre_rej_count=size(locs,2);    
    fprintf('\n => rejecting IED events... \n')
    [~,IED_locs] = findpeaks(IED_ch, 'MINPEAKHEIGHT', th(2));
    
    % add IED_onsets on this channel    %MGS
    IED_onset_locs = find(IED_onsets);    
    IED_locs = unique(sort([IED_locs, IED_onset_locs]));
    
    % Ignore IED-ripples events (event that coincide within 50 ms):
    exc.IED = []; 
    for i=1:numel(IED_locs)
        tmp=find(abs(locs-IED_locs(i)) < (0.5*Fs)); % kk +/- 500ms
        exc.IED = [exc.IED locs(tmp)]; 
        if ~isempty(tmp)
            locs(tmp)=[];
            pks(tmp)=[];
            rej_count=rej_count+1;
        end
    end
    fprintf('\n *** rejected %d / %d events based on IED channel correlation \n',rej_count,pre_rej_count)
    ripples_stat.IED_count = rej_count; 
    ripples_stat.IED_percent = (rej_count/pre_rej_count) * 100;
end

counter=1;
ripples=nan(1,4);
ripples=array2table(ripples,'VariableNames',{'str','peak','fin','amplitude'});
ripples(1,:)=[];
for k=locs
    % find the starting point of the peak:
    stop=0;
    str=k;
    while ~stop && ~str==0
        if str==1
            break
        end
        str=str-1;
        if signal(str)<th(1), stop=1; end
    end
    % find the ending point of the peak:
    stop=0;
    fin=k;
    while ~stop && ~(fin==numel(signal))
        fin=fin+1;
        if signal(fin)<th(1), stop=1; end
    end
    % =====================================================================
    % Alternative #1:
    % Detect negative peak position for each ripple (closest to ripple's power peak)
    minIndex = [];
    [~,minpos] = findpeaks(-double(BP(str:fin)));
    [~,maxamp] = max(double(ENV(str:fin)));
    if isempty(minpos)
        [~,minpos] = min(BP(str:fin)); 
        minpos = minpos(1);
    else
         minpos=minpos-1; 
    end % cc & kk

    maxamp=maxamp-1;
    [~,tmp] = min(abs(minpos-maxamp));
    minIndex=minpos(tmp);
    peakPosition = min((str + minIndex),numel(signal));

    try
        ripples(counter,:)=array2table([t(str), t(peakPosition), t(fin), ENV(peakPosition)]);
    catch
        disp(ripples);
        fprintf('\n Error has occured in event # %d \n',counter);
    end
    counter=counter+1;
end
disp(['After detection by thresholding: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end


% Merge ripples if inter-ripple period is less than minDistance:
ripples_edit=ripples;
rej=zeros(size(ripples,1),1);
for k = 2:size(ripples,1)
    if (ripples.peak(k)-ripples.peak(k-1)) < minDistance
        % Merge
        ripples_edit.fin(k-1) = ripples.fin(k);
        rej(k)=1;
    end
end
if any(rej), ripples_edit(find(rej),:)=[]; end
ripples=ripples_edit;
disp(['After ripple merge: ' num2str(size(ripples,1)) ' events.']);
if isempty(ripples),return; end

% duration test:
duration = ripples.fin-ripples.str;
ripples(duration<minRippleDuration,:) = [];
disp(['After min duration test: ' num2str(size(ripples,1)) ' events.']);
duration = ripples.fin-ripples.str;
ripples(duration>maxRippleDuration,:) = [];
disp(['After max duration test: ' num2str(size(ripples,1)) ' events.']);

 ripples_stat.SWR_final = size(ripples,1); 
end
