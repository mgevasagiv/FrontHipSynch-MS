%% modified from 
% detect_SWR_events_in_hippocampus.m [Yitzhak Norman] and 
% find_ripples [Carlos Carrasco] 
% by Kamin Kim, 2020 Feb [Ranganath Lab @ UCD] 

function [ripples,ripples_stat, exc] = detect_SWR_MGS(data, coi, cref)
%         addpath(genpath('/Users/kaminkim/Documents/toolboxes/eeglab13_4_4b/'));
        % SWR detection algorithm.
        % Detecting ripples after exclusion of electrical/muscular artifacts and interictal epileptic discharges (IEDs).
        % The algorithm uses three main signals: 
        % (1) hippocampal (CA1/subiculum) ripple band (70-180 Hz)
        % (2) common average of all iEEG cahnnels (filtered between 70-180 Hz for control detection)
        % (3) hippocampal IED band (25-60 Hz)
        addpath(genpath('C:\Users\mgeva\Documents\GitHub\Norman_et_al_2023\'))
        
        % ripple detection  parameters 
        minDistance=0.030; % in sec
        minRippleDuration=0.020; % in sec 
        maxRippleDuration=0.200; % in sec
        th=[2 4]; % ripple detection thresholds (in std from the mean) [onset/offset, peak]
              
        % data parameters 
        time = double(data.time{1, 1}); 
        srate = data.fsample; 
        traces =  data.trial{1}(coi,:); 
        
        if nargin < 3
            cref = [];
        end
        if ~isempty(cref)
            traces_ref = data.trial{1}(cref,:);
            traces = traces - traces_ref;
        end
        elecs = data.elec.label(coi);
        
        % filter parameters (used previously): 
        ripple_band = [70 180]; % in Hz
        IED_band = [25 60]; % in Hz
        Fs = data.fsample; % in Hz
        wtype = 'hamming'; % window type
        
        % Low-pass filter for smoothing the ripple-band envelope:
        LPcutoff = round(mean(ripple_band)/pi); % in Hz (see Stark et al. 2014 for a similar method)
        % Note: the cutoff frequency is 40 Hz when ripple-band frequency window is 70-180 Hz 
        
        % Set up a lowpass filter for smoothing: FIR kaiserwin lowpass filter 
        lpFilt = designfilt('lowpassfir','PassbandFrequency',LPcutoff, ...
            'StopbandFrequency',LPcutoff+10,'PassbandRipple',0.001, ...
            'StopbandAttenuation',60,'DesignMethod','kaiserwin','SampleRate',Fs);
        % fvtool(lpFilt,'OverlayedAnalysis','phase')

        %% Filtered signal
        SWR_signal = eegfilt(traces, srate,70,180); 
        IED_signal = eegfilt(traces, srate, 25, 60);
        Global_signal = eegfilt(nanmean(data.trial{1}), srate, 70, 180); % kk mean to nanmean; bchans must be nanned out in a previous step

        for elec = 1:length(coi)
            
            % Compute stdev of ripple-band envelope across the entire experiment
            
            % hilbert envelope of 70-180Hz filtered data
            absSignal=double(abs(hilbert(SWR_signal(elec,:))));
            
            % Clipping the signal:
            topLim=nanmean(absSignal)+th(2)*nanstd(absSignal);
            absSignal(absSignal>topLim)=topLim;
            
            % Square & smooth
            squaredSignal = absSignal.^2;
            squaredSignal = filtfilt(lpFilt,double(squaredSignal));
            
            % Compute means and std:
            avg = nanmean(squaredSignal);
            stdev = nanstd(squaredSignal,[],2);
            
            % filtered signals
            signalA = SWR_signal(elec,:); %swr signal
            signalB = Global_signal; %common ave
            signalC = IED_signal(elec,:); %ied signal
            
            % analytic signal:
            squaredSignalA = double(abs(hilbert(signalA))).^2;
            squaredSignalB = double(abs(hilbert(signalB))).^2;
            squaredSignalC = double(abs(hilbert(signalC))).^2;
            
            % smoothing :
            squaredSignalA = filtfilt(lpFilt,double(squaredSignalA));
            squaredSignalB = filtfilt(lpFilt,double(squaredSignalB));
            squaredSignalC = filtfilt(lpFilt,double(squaredSignalC));
            
            % ZSCORE:
            squaredSignalNormA = (squaredSignalA-avg)/stdev; % Zscore hippocampal channel - ripples
            squaredSignalNormB = (squaredSignalB-nanmean(squaredSignalB))/nanstd(squaredSignalB); % Zscore common average
            squaredSignalNormC = (squaredSignalC-nanmean(squaredSignalC))/nanstd(squaredSignalC); % Zscore hippocmapal channel - IED
            
            % MayaGS- Exclude sharp transient in the raw EEG:
            rawSignal = traces(elec,:);
            IED_onsets = zeros(size(rawSignal));  % prepare timeseries for IEDs
            IEDthr = 8;
            dv = [0 diff(rawSignal)];
            [avgD,stdD] = robustMean(dv,[],IEDthr);
            dvnorm = (dv-avgD)./stdD;
            % sharpTransientsInd = dv > 50 * (1000/EEG.srate);  % excluding sharp transients (>50muV/ms)
            sharpTransientsInd = dvnorm > IEDthr;  % excluding sharp transients
            IED_onsets(sharpTransientsInd) = 1;
            
            llw = 0.04; % in Sec
            prc = 99.9;
            [ets,ech]=LLspikedetector(rawSignal,srate,llw,prc,[]);
            IED_onsets(ets(:,1)) = 1;
            IED_onsets(ets(:,2)) = 1;
            
            % Find ripples:
            [ripples.(elecs{elec}), ripples_stat.(elecs{elec}), exc.(elecs{elec})] = ripples_detection_excluding_IED_MGS(squaredSignalNormA,signalA,time,th,minDistance,...
                minRippleDuration,maxRippleDuration,squaredSignalNormB,squaredSignalNormC, IED_onsets, Fs);
            
            %                 figure('units','normalized','outerposition',[0 0 1 0.5],'Color','w');
            %                 hold on;
            %                 title(['Ripple Detection, N = ',num2str(height(ripples.(elecs{elec}))), ' Ripples (ripple band envelope)']);
            %                 h1 = plot(time,squaredSignalNormB,'linewidth',0.5,'color',[1,0.7,0.7]); hold on;
            %                 h2 = plot(time,squaredSignalNormA,'k','Linesmoothing','on','Linewidth',1);
            %                 h0 = plot(time,ones(size(time))*th(2),'b--','Linewidth',1);
            %                 h3=scatter(ripples.(elecs{elec}).peak,ones(size(ripples.(elecs{elec}),1),1)*-2,30,'m','fill'); hold on
            %
            %                 ylabel(sprintf('%d-%dHz Amplitude^2 (Z-score)',ripple_band(1),ripple_band(2)))
            %                 xlim([time(1),time(end)])
            %                 ylim([-2,10])
            %                 legend([h2 h1 h3 h0],{'Ripple Band Amplitude (zscore)','Common avg./noise','Ripple','Detection thr.'});legend boxoff;
            %                 drawnow;
            %                 xlabel ('Time (s)');
            %                 set(gca,'linewidth',2);
            %                 set(gca,'fontsize',14);
        end
end 