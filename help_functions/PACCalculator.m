classdef PACCalculator < handle
    properties
    
        samplingRate = 1000;
        
        %general PAC calculation constants
        minSecToCalculatePacOn = 30; %based on Tort et al
        findPreferredPhase = false;
        
        %frequency band constants
        lowLimitPhase = 0.16;
        highLimitPhase = 1.25;
        lowLimitAmp = 9;
        highLimitAmp = 16;
        slowFreqForOutlierRemoval = 1;
        
        %comodulogram constants - taken from Tort et al
        phaseFrequencyStepSize = 2; %the phase frequency bands will be sampled in steps of phaseFrequencyStepSize, each time using the bandwidth phaseFrequencyBandwidth
        phaseFrequencyBandwidth = 4;
        ampFrequencyStepSize = 5; %the amplitude frequency bands will be sampled in steps of ampFrequencyStepSize, each time using the bandwidth ampFrequencyBandwidth
        ampFrequencyBandwidth = 10;
        minFreqToCheckPhase = 2;
        maxFreqToCheckPhase = 50;
        minFreqToCheckAmp = 5;
        maxFreqToCheckAmp = 50;
        
        
        %filtering constants
        defaultFilterOrder = 1;
        nanWarning = 0.01;
        
        %outliers constansts
        removeOutliers = true;
        nSigmasForFastAmpOutlier = 2; %parameter from Lafon et al
        percentileForSlowAmpOutlier = 40; %the percentile to leave (i.e. if it's set to 40 we set the threshold such that 40% of the waves will pass it)
        thresholdForSlowAmpOutlier = nan;
        thresholdForFastAmpOutlier = nan;
        usePredefinedThresholdForSlowAmpOutlier = false;
        predefinedThresholdSlowAmpOutlier = 50; %parameter from Lafon et al
        
        %IIS removal constants
        windowAroundIIS = 500; %ms
        
        %MI calculation constants
        numBins = 18; %parameter from Tort et al
        plotDist = false;
        
        %coupling index method to use
        couplingIndexType = 'MI'; %options: MI (Tort et al), meanVector (used in Lafon et al)
        useHistEqualizationInMeanVector = true;
    
        %simulation constants
        slowFreq = 1;
        fastFreq = 15;
        segLength = 40; %seconds
        spikeFileName = 'spikeData';
        gradLength = 20; %ms
        nItersPernSpikes = 50;
        
        %p value calculation constants
        nItersForPValue = 1000;
        
        %constants for calculation of SI - from Staresina et al
        timeAroundEvent = 1; %seconds
        
        %constants for spectrum display
        freqRange = [5:30];
        timeBeforeAfterEvent = 1.5; %seconds
        timeForBaseline = 1; %seconds, from Starestina et al
        minNCycles = 5; % from Starestina et al
        minWinSizeSpec = 100; %ms, from Starestina et al
        nanThreshPercEpoch = 0.1; %do not include epoch segments with too many nans
        nanThreshPercBL = 0.3; %do not include baseline segments with too many nans
        
        %constants for artificial data
        maxSpindleAmp = 0.3;
        slowRateTest = 0.5;
        nCyclesTest = 1500;
        percSpindlesTest = 0.5;
        interpTime = 0.05; %seconds
        
        %scoring params
        scoringEpochDuration = 0.001; % How many seconds represented by one individual value in the scoring vector [scalar].
        sleepEpochs = [1]; % all the values in the scoring vector which represent sleep stages for which we want to perform the analysis (like NREM/REM/transitions) [1D vector].
        
    end
    
    methods
        

        function [PAC, prefPhase] = calculatePAC(obj, dataPhase, dataAmp, IIStimes, phaseFreqLimits, ampFreqLimits)
            
            %The methods calculates the phase-amplitude coupling index
            %using the method as specific in the propery couplingIndexType,
            %options:
            %A. MI - a modification of KL distance (from uniform distribution) of the
            %"distribution" of the mean amplitude per phase (based on Tort
            %et al, j neurophysiol 2010)
            %B. meanVector - the absolute value of the mean vector of the
            %time series amp*exp(i*phase), either calculated directly from
            %the data or on the same "amplitude distribution" (calculating
            %it on the amplitude distribution will neutralize the
            %irrelevant contribution of the specific phases distribution in the data
            %to this measure and will normalize the size so it's between
            %0-1).
            %input - 
            % dataPhase - the data from which the phase will be extracted
            % dataAmp - the data from which the amplitude will be extracted
            % if only dataPhase is provided or dataAmp is empty the code
            % will use the same data segment for both purposes
            % optional: 
            % IIStimes - if provided, the code will remove segments around
            % IIStimes (the size of the window to remove is defined by the property
            % windowAroundIIS) after the Hilbert transform. IIStimes can be
            % a list of indices if there is only one data set (i.e.
            % dataPhase and dataAmp are the same dataset), or a cell of 2
            % array, each one corresonding to each data set (IIStimes{1}
            % corresponding to dataPhase, IIStimes{2} corresponding to
            % dataAmp). 
            % phaseFreqLimits - a vector of two values defining the
            % frequency band of the "phase modulating" time series (min and
            % max value). If not provided it will be set by default to the
            % properties [lowLimitPhase highLimitPhase]
            % ampFreqLimits - a vector of two values defining the
            % frequency band of the "amplitude modulated" time series (min and
            % max value). If not provided it will be set by default to the
            % properties [lowLimitAmp highLimitAmp]
            % output - 
            % PAC
            % prefPhase - the preferred phase (phase at which amplitude is
            % maximal), only calculated if the property findPreferredPhase
            % is true (otherwise will return nan)
            
            if nargin < 3 || isempty(dataAmp)
                dataAmp = dataPhase;
            end
            
            if nargin < 5 || isempty(phaseFreqLimits)
                lowLimitPhase = obj.lowLimitPhase;
                highLimitPhase = obj.highLimitPhase;
            else
                lowLimitPhase = phaseFreqLimits(1);
                highLimitPhase = phaseFreqLimits(2);
            end
            
            if nargin < 6 || isempty(ampFreqLimits)
                lowLimitAmp = obj.lowLimitAmp;
                highLimitAmp = obj.highLimitAmp;
            else
                lowLimitAmp = ampFreqLimits(1);
                highLimitAmp = ampFreqLimits(2);
            end
            
            if nargin < 4 || isempty(IIStimes)
                removeIIS = false;
            else
                removeIIS = true;
                if ~iscell(IIStimes)
                    oldTimes = IIStimes;
                    IIStimes = cell(1,2);
                    IIStimes{1} = oldTimes;
                    IIStimes{2} = oldTimes;
                end
            end
            
           
            %filter the data at the phase modulating frequency band
            if lowLimitPhase == 0 || lowLimitAmp == 0; error('butter should get non-zero limits'); end
     
            %remove nan inds before the transform % MGS - change 2025
            nanIndsP = isnan(dataPhase);
            nanIndsA = isnan(dataAmp);
            
            dataPhase(isnan(dataPhase)) = 0;
            dataAmp(isnan(dataAmp)) = 0;
            
            dataFilteredFPhase = obj.bandpass(dataPhase, lowLimitPhase, highLimitPhase);
            dataFilteredFAmp = obj.bandpass(dataAmp, lowLimitAmp, highLimitAmp);
                  
            dataFilteredFPhase(nanIndsP) = 0;
            dataFilteredFAmp(nanIndsA) = 0;
            
            %calculate the time series of the phases of dataFilteredFPhase
            phiFP = angle(hilbert(dataFilteredFPhase));
            %Make original NaN indices NaN again
            phiFP(nanIndsP) = nan;
            
            %calculate the time series of the amplitude of dataFilteredFAmp
            AFP = abs(hilbert(dataFilteredFAmp));
            %Make original NaN indices NaN again
            AFP(nanIndsA) = nan;
            
%             subplot(2,1,1); plot(dataPhase(1:10000))
%             subplot(2,1,2); plot(dataFilteredFAmp(1:10000));
%             hold all; plot(AFP(1:10000),'color','g');
%             xlabel('Time (milisec)')
            
            if obj.removeOutliers
                [phiFP, AFP] = obj.removeOutlierCycles(dataFilteredFPhase, dataFilteredFAmp, phiFP, AFP);
            end
            
            %if IIS times were provided - replace them with NaN
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes{1})
                    pointsBefore = min(IIStimes{1}(iTime),winAroundIIS);
                    pointsAfter = min(length(phiFP)-IIStimes{1}(iTime),winAroundIIS);
                    phiFP(IIStimes{1}(iTime)-pointsBefore+1:IIStimes{1}(iTime)+pointsAfter) = nan;
                end
                
                for iTime = 1:length(IIStimes{2})
                    pointsBefore = min(IIStimes{2}(iTime),winAroundIIS);
                    pointsAfter = min(length(AFP)-IIStimes{2}(iTime),winAroundIIS);
                    AFP(IIStimes{2}(iTime)-pointsBefore+1:IIStimes{2}(iTime)+pointsAfter) = nan;
                end
            end
%             figure(1);
%             hold off;
%             subplot(2,1,1); plot(dataPhase);
%             hold off;
%             subplot(2,1,2); plot(dataFilteredFAmp); hold all; plot(AFP);
%             pause;
            
            if sum(~isnan(phiFP)&~isnan(AFP)) < (obj.minSecToCalculatePacOn * obj.samplingRate)
                disp('Not enough data to calculate PAC on after outlier and IIS removal');
                PAC = nan;
                prefPhase = nan;
                return;
            end

            switch obj.couplingIndexType
                case 'MI'
                    [PAC, prefPhase] = obj.calcMI(phiFP,AFP);
                case 'meanVector'
                    [PAC, prefPhase] = obj.calcMeanVector(phiFP,AFP);
                otherwise
                    PAC = nan;
                    prefPhase = nan;
            end
            
        end        
        
        function [comMat,phaseFreqToIter,ampFreqToIter] = calculateComodulogram(obj, dataPhase, dataAmp, IIStimes, toPlot, dataForThresholdSettingPhase, dataForThresholdSettingAmp)
            
            %calculates the comodulogram matrix of the data:
            %input - 
            % dataPhase - the data from which the phase will be extracted
            % dataAmp - the data from which the amplitude will be extracted
            % if only dataPhase is provided or dataAmp is empty the code
            % will use the same data segment for both purposes
            % IIStimes - if provided, the code will remove segments around
            % IIStimes (the size of the window to remove is defined by the property
            % windowAroundIIS) after the Hilbert transform. IIStimes can be
            % a list of indices if there is only one data set (i.e.
            % dataPhase and dataAmp are the same dataset), or a cell of 2
            % array, each one corresonding to each data set (IIStimes{1}
            % corresponding to dataPhase, IIStimes{2} corresponding to
            % dataAmp). 
            % toPlot - whether to plot imagesc of the comodulogram
            % dataForThresholdSettingPhase - a different set of data to set the
            % threshold for the phase (slow) outliers removal on for each combination of
            % frequencies examined. If not provided, threshold will be set
            % on the input data set (dataPhase & dataAmp)
            % dataForThresholdSettingAmp - a different set of data to set the
            % threshold for the phase (fast) outliers removal on for each combination of
            % frequencies examined. If only dataForThresholdSettingPhase is
            % provided, dataForThresholdSettingAmp will be set to be equal
            % to dataForThresholdSettingPhase
            %comMat - the comodulogram matrix contains the PAC index for each
            %combination of the frequency band of the phase modulating time
            %series and the frequency band of the time modulated time
            %series (based on Tort et al, j neurophysiol 2010).
            %phaseFreqToIter, ampFreqToIter - the center frequencies for the frequency bands which were iterated on.
            %They are set by the properties: phaseFrequencyBandwidth (band width of each frequency band),
            %phaseFrequencyStepSize (step size between centers of frequency bands), 
            %minFreqToCheckPhase (minimal frequency that will be the center of a frequency band)
            %maxFreqToCheckPhase (maximal frequency that will be the center of a frequency band)
            %for the phase modulating time series. Accordingly ampFrequencyBandwidth,
            %ampFrequencyStepSize, maxFreqToCheckAmp, minFreqToCheckAmp set the iterated frequency bands for the amplitude
            %modulated time series.
            
            if nargin < 3 || isempty(dataAmp)
                dataAmp = dataPhase;
            end
            
            if nargin < 4
                IIStimes = [];
            end
            
            if nargin < 5 || isempty(toPlot)
                toPlot = true;
            end
            
            if nargin < 6 || isempty(dataForThresholdSettingPhase)
                setThresholdPerFreq = false;
            else
                setThresholdPerFreq = true;
                if nargin < 7 || isempty(dataForThresholdSettingAmp)
                    dataForThresholdSettingAmp = dataForThresholdSettingPhase;
                end
            end
            
            oldPlotDist = obj.plotDist;
            if obj.plotDist
                obj.plotDist = false;
            end
            
            phaseFreqToIter = [max(obj.phaseFrequencyBandwidth/2,obj.minFreqToCheckPhase):obj.phaseFrequencyStepSize:obj.maxFreqToCheckPhase-obj.phaseFrequencyBandwidth/2];
            ampFreqToIter = [max(obj.ampFrequencyBandwidth/2,obj.minFreqToCheckAmp):obj.ampFrequencyStepSize:obj.maxFreqToCheckAmp-obj.ampFrequencyBandwidth/2];
            
            nPhaseFreqs = length(phaseFreqToIter);
            nAmpFreqs = length(ampFreqToIter);
            
            comMat = zeros(nAmpFreqs,nPhaseFreqs);
            
            oldThresholdSlowOutlier = obj.thresholdForSlowAmpOutlier;
            oldThresholdFastwOutlier = obj.thresholdForFastAmpOutlier;
            
            if ~setThresholdPerFreq
                obj.thresholdForSlowAmpOutlier = nan;
                obj.thresholdForFastAmpOutlier = nan;
            end
            
            %prepare the sets of thresholds for each frequency combination 
            if setThresholdPerFreq
                thresholdsFast = zeros(nAmpFreqs,nPhaseFreqs);
                thresholdsSlow = zeros(1,nPhaseFreqs);
                
                for iPhaseFreq = 1:nPhaseFreqs
                    phaseLims = [phaseFreqToIter(iPhaseFreq)-obj.phaseFrequencyBandwidth/2 phaseFreqToIter(iPhaseFreq)+obj.phaseFrequencyBandwidth/2];
                    obj.setThreshAmpSlow(dataForThresholdSettingPhase, phaseLims);
                    thresholdsSlow(iPhaseFreq) = obj.thresholdForSlowAmpOutlier;
                    for iAmpFreq = 1:nAmpFreqs
                        ampLims = [ampFreqToIter(iAmpFreq)-obj.ampFrequencyBandwidth/2 ampFreqToIter(iAmpFreq)+obj.ampFrequencyBandwidth/2];
                        
                        obj.setThreshAmpFast(dataForThresholdSettingPhase, dataForThresholdSettingAmp, phaseLims, ampLims);
                        thresholdsFast(iAmpFreq,iPhaseFreq) = obj.thresholdForFastAmpOutlier;
                    end
                end
                
            end
            
            for iPhaseFreq = 1:nPhaseFreqs
%                 if mod(iPhaseFreq,10)==0
%                     disp(['completed ', num2str(iPhaseFreq*100/nPhaseFreqs), '%']);
%                 end
                if setThresholdPerFreq
                    obj.thresholdForSlowAmpOutlier = thresholdsSlow(iPhaseFreq);
                end
                
                for iAmpFreq = 1:nAmpFreqs
                    phaseLims = [phaseFreqToIter(iPhaseFreq)-obj.phaseFrequencyBandwidth/2 phaseFreqToIter(iPhaseFreq)+obj.phaseFrequencyBandwidth/2];
                    ampLims = [ampFreqToIter(iAmpFreq)-obj.ampFrequencyBandwidth/2 ampFreqToIter(iAmpFreq)+obj.ampFrequencyBandwidth/2];
                    
                    if setThresholdPerFreq
                        %set the thresholds for each frequency combination
                        obj.thresholdForFastAmpOutlier = thresholdsFast(iAmpFreq,iPhaseFreq);
                    end
                    
                    comMat(iAmpFreq,iPhaseFreq) = obj.calculatePAC(dataPhase, dataAmp, IIStimes, phaseLims, ampLims);
                end
            end
            
            if toPlot
                figure;
                imagesc(ampFreqToIter,phaseFreqToIter,comMat);
                colorbar;
            end
            
            obj.plotDist = oldPlotDist;
            
            obj.thresholdForSlowAmpOutlier = oldThresholdSlowOutlier;
            obj.thresholdForFastAmpOutlier = oldThresholdFastwOutlier;
        end
        
        function meanPAC = calculateMeanPacOnComodulogram(obj, dataPhase, dataAmp, IIStimes)
            
            %Calculates the average of the comodulogram that is created
            %using the class properties
            
            if nargin < 3 || isempty(dataAmp)
                dataAmp = dataPhase;
            end
            
            if nargin < 4
                IIStimes = [];
            end
            
            %calculate comodulogram 
            comMat = obj.calculateComodulogram(dataPhase, dataAmp, IIStimes, false);
            meanPAC = nanmean(comMat(:));
            
        end
        
        function [phiFP, AFP] = removeOutlierCycles(obj, dataFilteredFPhase, dataFilteredFAmp, phiFP, AFP)
            
            %The methods removes outlier cycles of two types:
            %A. cycles for which the mean squared amplitude of the fast
            %frequency band is abnormally large - if obj.thresholdForFastAmpOutlier 
            %is nan then this is defined as a set amount of sigmas (obj.nSigmasForFastAmpOutlier) above the mean 
            %relative either to the current input data. If obj.thresholdForFastAmpOutlier
            %is not nan then it defines the threshold (it can be set for a
            %larger data set using the method setThreshAmpFastPerElectrode)
            %B. cycles for which the peak amplitude of the slow oscillation
            %is too small - as defined by thresholdForSlowAmpOutlier
            %Removing a cycle means replacing it with NaN
            % input - 
            % dataFilteredPhase - the low frequency band data
            % dataFilteredFAmp - the high frequency band data
            % phiFP - the phase of dataFilteredPhase (calculated using
            % hilbert transform)
            % AFP - the amplitude of dataFilteredFAmp (calculated using
            % hilbert transform)
            % output - 
            % updated PhiFP and AFP after the removal (change to NaN) of
            % outlier cycles
            
            %find points where new cycles start (a jump from phase pi to
            %phase -pi)
            diffPhiFP = diff(phiFP);
            newCyclePoints = find(diffPhiFP<-6)+1;
            nCycles = length(newCyclePoints)-1;
            
%             plot(dataFilteredFPhase);
%             hold all;
%             for ii=1:length(newCyclePoints)
%                 plot(newCyclePoints(ii), dataFilteredFPhase(newCyclePoints(ii)), '*', 'color', 'k');
%             end


            ampsFastPerCycle = nan(1, nCycles);
            maxAmpSlowPerCycle = nan(1, nCycles);
            
            %calculate the mean squared amplitude of the HF data per cycle
            %and the peak amplitude of the LF data per cycle
            for iCycle = 1:nCycles
                currCycleSlow = dataFilteredFPhase(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                currCycleFast = dataFilteredFAmp(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                
                ampsFastPerCycle(iCycle) = nanmean(currCycleFast.^2);
                maxAmpSlowPerCycle(iCycle) = nanmax(currCycleSlow);
            end
            
            %if there is no pre-defined threshold for the HF data, use the current data to
            %find a threshold - obj.nSigmasForFastAmpOutlier above the mean
            if isnan(obj.thresholdForFastAmpOutlier)
                cyclesToRemoveFastAmp = find(ampsFastPerCycle>nanmean(ampsFastPerCycle)+obj.nSigmasForFastAmpOutlier*nanstd(ampsFastPerCycle));
            else
                cyclesToRemoveFastAmp = find(ampsFastPerCycle>obj.thresholdForFastAmpOutlier);
            end
            
            %if there is no pre-defined threshold for the LF data, use the current data to
            %find a threshold - at the percentileForSlowAmpOutlier
            %percentile
            if isnan(obj.thresholdForSlowAmpOutlier) && ~obj.usePredefinedThresholdForSlowAmpOutlier
                cyclesToRemoveSlowAmp = find(maxAmpSlowPerCycle<prctile(maxAmpSlowPerCycle,100-obj.percentileForSlowAmpOutlier));
            else
                %if there is a predefined threshold (like in Lafon), use
                %the predefined threshold
                if obj.usePredefinedThresholdForSlowAmpOutlier
                    currThresh = obj.thresholdForSlowAmpOutlier;
                else
                    %if the threshold is not predefined but was calculated
                    %previously (on a larger data set for example), use
                    %this threshold
                    currThresh = obj.thresholdForSlowAmpOutlier;
                end
                cyclesToRemoveSlowAmp = find(maxAmpSlowPerCycle < currThresh);
            end
            
            cyclesToRemove = union(cyclesToRemoveFastAmp,cyclesToRemoveSlowAmp);
            
            %set cycles to remove to NaN
            for iCycle = 1:length(cyclesToRemove)
                phiFP(newCyclePoints(cyclesToRemove(iCycle)):newCyclePoints(cyclesToRemove(iCycle)+1)-1) = nan;
                AFP(newCyclePoints(cyclesToRemove(iCycle)):newCyclePoints(cyclesToRemove(iCycle)+1)-1) = nan;
            end
        end
        
        function [MI, prefPhase] = calcMI(obj, phiFP, AFP)
            %calculate the modulation index of the phase time series and amplitude time series - based on
            %Tort et al, j neurophysiol 2010
            % prefPhase - the preferred phase (phase at which amplitude is
            % maximal), only calculated if the property findPreferredPhase
            % is true (otherwise will return nan)
            
            [ampDist, ~, prefPhase] = obj.buildAmpDistPerPhase(phiFP, AFP);
%             [ampDist, binCenters] = obj.buildAmpDistPerPhase(phiFP, AFP);
%             figure(2);
%             hold off;
%             bar(binCenters,ampDist);
%             pause;
            
            distEnt = -sum(ampDist.*log(ampDist));
            MI = (log(obj.numBins)-distEnt)/log(obj.numBins);
            
        end
        
        function [meanVector, prefPhase] = calcMeanVector(obj,  phiFP, AFP)
            
            %calculate the mean vector of the time series: A*exp(phi)
            % prefPhase - the preferred phase (phase at which amplitude is
            % maximal), only calculated if the property findPreferredPhase
            % is true (otherwise will return nan)
            
            if obj.useHistEqualizationInMeanVector
                %This options can be used in order to:
                %A. Equlize the probability of each phase in order to
                %neutralize bias in the index calcualtion created by the
                %non uniform distribution of phases.
                %B. Normalize the index so that it can't be larger than 1.
                [ampDist, binCenters, prefPhase] = obj.buildAmpDistPerPhase(phiFP, AFP);
                meanVector = abs(mean(ampDist.*exp(i*binCenters)));
            else
                %do not equalize distribution of phases and do not
                %normalize amplitudes
                meanVector = abs(nanmean(AFP.*exp(i*phiFP)));
                prefPhase = nan;
            end
            
        end
        
        function [ampDist, binCenters, prefPhase] = buildAmpDistPerPhase(obj, phiFP, AFP)
            %builds the average amplitude "distribution" per binned phases -
            %based on Tort et al, j neurophysiol 2010
            % prefPhase - the preferred phase (phase at which amplitude is
            % maximal), only calculated if the property findPreferredPhase
            % is true (otherwise will return nan)
            
            %change range to 0-2pi
            phiFP = phiFP + pi;
            
            binSize = 2*pi/obj.numBins;
            binCenters = [0:binSize:2*pi]+(binSize/2)-pi;
            binCenters = binCenters(1:end-1);
            %translate to bin indices
            binIndices = floor(phiFP/binSize)+1;
            ampDist = zeros(1,obj.numBins);
            for iBin = 1:obj.numBins
                ampDist(iBin) = nanmean(AFP(binIndices==iBin));
            end
            
            %normalize mean amplitude
            ampDist = ampDist./sum(ampDist);
            
            binCentersTrip = [binCenters-2*pi binCenters binCenters+2*pi];
            ampDistTrip = [ampDist ampDist ampDist];
            
            if obj.plotDist
                fitRes = obj.createCosFitModel(binCentersTrip', ampDistTrip');
                hold off;
                figure;
                bar(binCenters,ampDist);
                hold all;
                plot(binCenters, fitRes.results.a*cos(binCenters+fitRes.results.b)+fitRes.results.c);
                prefPhase = mod(-fitRes.results.b,2*pi);
                if prefPhase > pi
                    prefPhase = prefPhase - 2*pi;
                end
                plot(prefPhase,0,'*','markersize',14);
                xlabel('Phase (rad)');
                ylabel('Mean amplitude');
                title('Mean amplitude per phase');
            end
            
            if obj.findPreferredPhase
                if ~exist('fitRes','var')
                    %fit a cos model (a*cos(x + b)+c) to the distribution in order to find
                    %the argmax (which is -b)
                    fitRes = obj.createCosFitModel(binCenters', ampDist');
                end
                prefPhase = mod(-fitRes.results.b,2*pi);
                if prefPhase > pi
                    prefPhase = prefPhase - 2*pi;
                end
            else
                prefPhase = nan;
            end
        end
        
        function BP = bandpass(obj, timecourse, lowLimit, highLimit, filterOrder)
            
            %bandpass code - from Maya
            
            if (nargin < 5)
                filterOrder = obj.defaultFilterOrder;
            end
            
            % Maya GS - handle NAN values
            indices = find(isnan(timecourse));
            if length(indices) > obj.nanWarning*length(timecourse)
                warning('many NaN values in filtered signal')
            end
            timecourse(indices) = 0;
            %
            
            [b, a] = butter(filterOrder, [(lowLimit/obj.samplingRate)*2 (highLimit/obj.samplingRate)*2]);
            BP = filtfilt(b, a, timecourse );
            BP(indices) = NaN;
        end
        
        function [signal, spikeLocations] = simulateSignal(obj, multConst, nSpikes, phaseFreq, ampFreq)
            
            %The method creates simulated data which is the summation of
            %two cos signals at different frequncies, where the amplitude of the fast
            %signal (ampPhase) can be affected by the phase of the slow signal (phaseSignal). 
            %The amplitude of ampSignal is multiplied by: (0.5+multConst*phaseSignal). 
            %Thus if multConst is 0 there is no correlation between them,
            %as the larger it is the larger the correlation.
            %The method then optionally adds spikes to the simulated data at random locations.
            %Input:
            %multConst - sets the degree of correlation between the
            %amplitude and the phase of the fast and slow signals
            %nSpikes - the number of spikes that are added to the signal
            %phaseFreq - the frequency of the slow signal, set by default
            %by obj.phaseFreq
            %ampFreq - the frequency of the slow signal, set by default
            %by obj.ampFreq. (Note that the length of the added spike is
            %200 ms, thus a default of ampFreq = 15 had the advantage of
            %being a multiplication of 5, thus the addition of the spikes
            %is "smoother")
            %Output:
            %signal - stimaulted signal
            %spikeLocations - the middle of the added spikes
            
            if nargin < 4 || isempty(slowFreq)
                phaseFreq = obj.slowFreq;
            end
            
            if nargin < 5 || isempty(fastFreq)
                ampFreq = obj.fastFreq;
            end

            domain = [0:1/obj.samplingRate:obj.segLength];
            nPoints = length(domain);
            %creates the two signals and their summation
            phaseSignal = cos(2*pi*phaseFreq*domain);
            ampSignal = cos(2*pi*ampFreq*domain).*(0.5+multConst*phaseSignal);
            signal = phaseSignal+ampSignal;

            if nargin < 3
                nSpikes = 0;
            end
            
            if nSpikes > 0
                
                %the spike shaped is loaded from a file
                spikeData = load(obj.spikeFileName);
                spikeData = spikeData.spikeData;
                lenSpike = length(spikeData);
                %how many points after the spike to gradually return to the
                %value of the signal to which we added the spike
                gradLength = obj.gradLength * 1000 / obj.samplingRate;
                
                currNSpikes = 0;
                spikeLocations = [];
                %the loop chooses random location for the spikes, the
                %distance between two spikes is set to be a minimum of
                %lenSpike
                while currNSpikes < nSpikes
                    randLoc = randi(nPoints-lenSpike+1);
                    if ~isempty(spikeLocations) && any(abs(spikeLocations-randLoc)<lenSpike*2)
                        continue;
                    end
                    spikeLocations(end+1) = randLoc;
                    currNSpikes = currNSpikes+1;
                end
                
                for iLoc = 1:nSpikes
                    %add the spike
                    signal(spikeLocations(iLoc):spikeLocations(iLoc)+lenSpike-1) = signal(spikeLocations(iLoc))+spikeData;
                    
                    %creating a gradient between the end of the spike to a
                    %point in a distance of obj.gradLength ms later
                    endSpikeVal = spikeData(end)+signal(spikeLocations(iLoc));
                    distToNextPoint = min(nPoints-spikeLocations(iLoc)-lenSpike, gradLength);
                    valAtNextPoint = signal(spikeLocations(iLoc)+lenSpike+distToNextPoint);
                    signal(spikeLocations(iLoc)+lenSpike:spikeLocations(iLoc)+lenSpike+distToNextPoint) = [0:distToNextPoint]*(valAtNextPoint-endSpikeVal)/distToNextPoint+endSpikeVal;
                end
                
                %the locations that are returned are the points at the middle of the spikes
                spikeLocations = spikeLocations+round(lenSpike/2);
            end
            
        end
        
        function plotEffectOfNSpikesOnPAC(obj,multConst,maxSpikes,toPlotDist)
            %The method creates signals with the multiplication const mult
            %consts, adds to them [1:maxSpikes] spikes, calculates the PAC
            %for each number of added spikes and plots the results
            
            if nargin < 4 || isempty(toPlotDist)
                toPlotDist = false;
            end
            
           
            oldRemoveOutliers = obj.removeOutliers;
            obj.removeOutliers = false;
            
            %calculate baseline PAC - the "true" PAC of the signal
            tmpSig = obj.simulateSignal(multConst);
            if toPlotDist
                figure;
                plot(tmpSig);
                xlabel('Time');
                ylabel('Amplitude');
                title(['Simulated signal, multConst = ',num2str(multConst)]);
                oldPlotDist = obj.plotDist;
                obj.plotDist = true;
            end
            baselinePac = obj.calculatePAC(tmpSig);
            if toPlotDist
                title(['PAC index type is ',obj.couplingIndexType,' PAC index = ',num2str(baselinePac)]);
                obj.plotDist = oldPlotDist;
            end
            pac = zeros(maxSpikes,obj.nItersPernSpikes);
            pacRem = zeros(maxSpikes,obj.nItersPernSpikes);
            
%             tmpSigsLow = {};
%             tmpSigsHigh = {};
%             pacsLow = [];
%             pacsHigh = [];
%             spLocLow = {};
%             spLocHigh = {};
            for iSpike = 1:maxSpikes
                for iIter = 1:obj.nItersPernSpikes
                    %create signals with spikes
                    [tmpSig, spLoc] = obj.simulateSignal(multConst,iSpike);
                    %calculate the PAC without removing the spikes and with
                    %removing the spikes
                    pac(iSpike,iIter) = obj.calculatePAC(tmpSig);
                    pacRem(iSpike, iIter) = obj.calculatePAC(tmpSig, [], spLoc);
%                     if pacRem(iSpike,iIter)>baselinePac
%                         tmpSigsHigh{end+1} = tmpSig;
%                         pacsHigh(end+1) = pacRem(iSpike, iIter);
%                         spLocHigh{end+1} = spLoc;
%                     end
%                     if pacRem(iSpike,iIter)<baselinePac
%                         tmpSigsLow{end+1} = tmpSig;
%                         pacsLow(end+1) = pacRem(iSpike, iIter);
%                         spLocLow{end+1} = spLoc;
%                     end
                end
            end
            
            %calculates the mean PAC without removing spikes
            meanPac = zeros(1,11);
            meanPac(1) = baselinePac;
            for iSpike=1:maxSpikes
                meanPac(iSpike+1) = mean(pac(iSpike,:));
            end
            
            %calculates the mean PAC with removing spikes
            meanPacRem = zeros(1,11);
            meanPacRem(1) = baselinePac;
            for iSpike=1:maxSpikes
                meanPacRem(iSpike+1) = mean(pacRem(iSpike,:));
            end
            
            figure;
            hold off;
            subplot(1,2,1);
            hold off;
            plot(0,baselinePac,'*','color','k');
            hold all;
            for iSpike = 1:maxSpikes
                plot(iSpike,pac(iSpike,:),'*','color','k');
                hold all;
            end
            
            %plot the mean for each number of spikes
            plot([0:maxSpikes],meanPac,'color','r');
            xlabel('# spikes');
            ylabel('PAC');
            title(['mlutiplication const = ', num2str(multConst), ' without spike removal, segLength = ',num2str(obj.segLength)]);
%             limits1 = axis;
            
            subplot(1,2,2);
            hold off;
            plot(0,baselinePac,'*','color','k');
            hold all;
            for iSpike = 1:maxSpikes
                plot(iSpike,pacRem(iSpike,:),'*','color','k');
                hold all;
            end
            
            %plot the mean for each number of spikes
            plot([0:maxSpikes],meanPacRem,'color','r');
            xlabel('# spikes');
            ylabel('PAC');
            title(['mlutiplication const = ', num2str(multConst), ' with spike removal']);
%             limits2 = axis;
            
%             maxylim = max(limits1(4),limits2(4));
%             subplot(1,2,1); ylim([0 maxylim]);
%             subplot(1,2,2); ylim([0 maxylim]);
            
            obj.removeOutliers = oldRemoveOutliers;
        end
        
        function setThreshAmpFast(obj, dataPhase, dataAmp, phaseFreqLimits, ampFreqLimits)
            
            %The method receives data and sets the property thresholdForFastAmpOutlier
            %to be obj.nSigmasForFastAmpOutlier sigmas above the mean of
            %squared amplitude per cycle of the HF data. The purpose of
            %this method is to find the threshold based on a long stretch
            %of data of the electrode and not separately on each segments
            %for which we calcualte PAC
            %input - 
            %dataPhase - the data from which we extract the LF
            %dataAmp - the data from which we extract the HF (if not
            %provided will be set to be the same set as dataPhase)
            
            if nargin < 3 || isempty(dataAmp)
                dataAmp = dataPhase;
            end
            
            if nargin < 4 || isempty(phaseFreqLimits)
                lowLimitPhase = obj.lowLimitPhase;
                highLimitPhase = obj.highLimitPhase;
            else
                lowLimitPhase = phaseFreqLimits(1);
                highLimitPhase = phaseFreqLimits(2);
            end
            
            if nargin < 5 || isempty(ampFreqLimits)
                lowLimitAmp = obj.lowLimitAmp;
                highLimitAmp = obj.highLimitAmp;
            else
                lowLimitAmp = ampFreqLimits(1);
                highLimitAmp = ampFreqLimits(2);
            end
            
            %filter the data at the phase modulating frequency band
            dataFilteredFPhase = obj.bandpass(dataPhase, lowLimitPhase, highLimitPhase);
            %filter the data at the amplitude frequency band
            dataFilteredFAmp = obj.bandpass(dataAmp, lowLimitAmp, highLimitAmp);
            
            %remove nan inds before the transform
            nanIndsP = isnan(dataPhase);
            nanIndsA = isnan(dataAmp);
            
            dataFilteredFPhase(nanIndsP) = 0;
            dataFilteredFAmp(nanIndsA) = 0;
            
            %calculate the time series of the phases of dataFilteredFPhase
            phiFP = angle(hilbert(dataFilteredFPhase));
            %Make original NaN indices NaN again
            phiFP(nanIndsP) = nan;
            
            %find points where new cycle start (a jump from phase pi to
            %phase -pi)
            diffPhiFP = diff(phiFP);
            newCyclePoints = find(diffPhiFP<-6)+1;
            nCycles = length(newCyclePoints)-1;
            
            ampsFastPerCycle = nan(1, nCycles);
            
            %calculate the mean squared amplitude of the HF data per cycle
            for iCycle = 1:nCycles
                currCycleFast = dataFilteredFAmp(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                ampsFastPerCycle(iCycle) = nanmean(currCycleFast.^2);
            end
           
            %set the threshold to be obj.nSigmasForFastAmpOutlier sigmas
            %above the mean
            obj.thresholdForFastAmpOutlier = nanmean(ampsFastPerCycle)+obj.nSigmasForFastAmpOutlier*nanstd(ampsFastPerCycle);
            
        end
        
        function setThreshAmpSlow(obj, dataPhase, phaseFreqLimits)
            
            %The method receives data and sets the property thresholdForSlowAmpOutlier
            %to be at the 100-obj.percentileForSlowAmpOutlier percentile of
            %the max amplitude per cycle of the LF data. The purpose of
            %this method is to find the threshold based on a long stretch
            %of data of the electrode and not separately on each segments
            %for which we calcualte PAC
            %input - 
            %dataPhase - the data from which we extract the LF
            
            if nargin < 2 || isempty(phaseFreqLimits)
                lowLimitPhase = obj.lowLimitPhase;
                highLimitPhase = obj.highLimitPhase;
            else
                lowLimitPhase = phaseFreqLimits(1);
                highLimitPhase = phaseFreqLimits(2);
            end
            
            %filter the data at the phase modulating frequency band
            dataFilteredFPhase = obj.bandpass(dataPhase, lowLimitPhase, highLimitPhase);
            
            %remove nan inds before the transform
            nanIndsP = isnan(dataPhase);
            
            dataFilteredFPhase(nanIndsP) = 0;
            
            %calculate the time series of the phases of dataFilteredFPhase
            phiFP = angle(hilbert(dataFilteredFPhase));
            %Make original NaN indices NaN again
            phiFP(nanIndsP) = nan;
            
            %find points where new cycle start (a jump from phase pi to
            %phase -pi)
            diffPhiFP = diff(phiFP);
            newCyclePoints = find(diffPhiFP<-6)+1;
            nCycles = length(newCyclePoints)-1;
            
            maxAmpSlowPerCycle = nan(1, nCycles);
            
            %calculate the mean squared amplitude of the HF data per cycle
            for iCycle = 1:nCycles
                currCycleSlow = dataFilteredFPhase(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                maxAmpSlowPerCycle(iCycle) = nanmax(currCycleSlow);
            end
           
            %set the threshold to be at the 100-obj.percentileForSlowAmpOutlier percentile
            obj.thresholdForSlowAmpOutlier = prctile(maxAmpSlowPerCycle,100-obj.percentileForSlowAmpOutlier);
            
        end
        
        function [pvalue, distribution] = calcPValue(obj, PACIndex, dataPhase, dataAmp, IIStimes, toPlot)
            
            %The method receives a PAC index, and returns its pvalue, as
            %calculated by finding the distribution of PAC indices on the
            %input data when the phase is randomly shifted for each cycle
            %(see Lafon et al)
            %input - 
            %PACIndex - the index for which we want to calculate the
            %p-value
            %dataPhase + dataAmp - the dataset on which we want to
            %calculate the PAC distribution. dataPhase - the data from which the phase will be extracted
            % dataAmp - the data from which the amplitude will be extracted
            % if only dataPhase is provided or dataAmp is empty the code
            % will use the same data segment for both purposes
            % IIStimes - if provided, the code will change segments around
            % IIS to nan for ampData
            % output - 
            % pvalue - the calculated pvalue of PACIndex (percentage of
            % calculated PAC indices larger than PACindex when shuffling
            % phases)
            % distribution - the entire distribution of calculated
            % PACindices
            
            if nargin < 4 || isempty(dataAmp)
                dataAmp = dataPhase;
            end
            
            if nargin < 5 || isempty(IIStimes)
                removeIIS = false;
            else
                removeIIS = true;
                if ~iscell(IIStimes)
                    oldTimes = IIStimes;
                    IIStimes = cell(1,2);
                    IIStimes{1} = oldTimes;
                    IIStimes{2} = oldTimes;
                end
            end
            
            if nargin < 6 || isempty(toPlot)
                toPlot = false;
            end
                        
            %filter the data at the phase modulating frequency band
            dataFilteredFPhase = obj.bandpass(dataPhase, obj.lowLimitPhase, obj.highLimitPhase);
            dataFilteredFAmp = obj.bandpass(dataAmp, obj.lowLimitAmp, obj.highLimitAmp);
            
            %remove nan inds before the transform
            nanIndsP = isnan(dataPhase);
            nanIndsA = isnan(dataAmp);
            
            dataFilteredFPhase(nanIndsP) = 0;
            dataFilteredFAmp(nanIndsA) = 0;
            
            %calculate the time series of the phases of dataFilteredFPhase
            phiFP = angle(hilbert(dataFilteredFPhase));
            %Make original NaN indices NaN again
            phiFP(nanIndsP) = nan;
            
            %calculate the time series of the amplitude of dataFilteredFAmp
            AFP = abs(hilbert(dataFilteredFAmp));
            %Make original NaN indices NaN again
            AFP(nanIndsA) = nan;
            
            %if IIS removal is required - remove from AFP
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                for iTime = 1:length(IIStimes{2})
                    pointsBefore = min(IIStimes{2}(iTime),winAroundIIS);
                    pointsAfter = min(length(AFP)-IIStimes{2}(iTime),winAroundIIS);
                    AFP(IIStimes{2}(iTime)-pointsBefore+1:IIStimes{2}(iTime)+pointsAfter) = nan;
                end
            end
            
            %find points where new cycles start (a jump from phase pi to
            %phase -pi)
            diffPhiFP = diff(phiFP);
            newCyclePoints = find(diffPhiFP<-6)+1;
            nCycles = length(newCyclePoints)-1;
            
            %remove the outliers
            if obj.removeOutliers
                [phiFP, AFP] = obj.removeOutlierCycles(dataFilteredFPhase, dataFilteredFAmp, phiFP, AFP);
            end
            
            if removeIIS
                for iTime = 1:length(IIStimes{1})
                    pointsBefore = min(IIStimes{1}(iTime),winAroundIIS);
                    pointsAfter = min(length(phiFP)-IIStimes{1}(iTime),winAroundIIS);
                    phiFP(IIStimes{1}(iTime)-pointsBefore+1:IIStimes{1}(iTime)+pointsAfter) = nan;
                end
            end
            
            if sum(~isnan(phiFP)&~isnan(AFP)) < (obj.minSecToCalculatePacOn * obj.samplingRate)
                disp('Not enough data to calculate PAC on after outlier and IIS removal');
                pvalue = nan;
                distribution = nan(1,obj.nItersForPValue);
                return;
            end
            
            %perform obj.nItersForPValue iterations of randomly changing to
            %phase in order to calculate PAC distribution over the data
            distribution = nan(1,obj.nItersForPValue);
            for iIter = 1:obj.nItersForPValue
                currPhiFP = phiFP;
                %go over the cycles and randomly shift the phase
                for iCycle = 1:nCycles
                    currCyclePhase = phiFP(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1);
                    if all(isnan(currCyclePhase))
                        continue;
                    end
                    currCyclePhase = currCyclePhase + rand*pi;
                    currCyclePhase(currCyclePhase > pi) = currCyclePhase(currCyclePhase > pi)-2*pi;
                    
                    currPhiFP(newCyclePoints(iCycle):newCyclePoints(iCycle+1)-1) = currCyclePhase;
                end
                
                
                switch obj.couplingIndexType
                    case 'MI',
                        distribution(iIter) = obj.calcMI(currPhiFP,AFP);
                    case 'meanVector'
                        distribution(iIter) = obj.calcMeanVector(currPhiFP,AFP);
                end
                
            end
            
            pvalue = sum(distribution>=PACIndex)/obj.nItersForPValue;
            if pvalue==0
                pvalue = eps;
            end
            
            if toPlot
                [F,XI] = ksdensity(distribution);
                figure;
                plot(XI,F);
                title(['Null distribution of PAC index, pvalue = ',num2str(pvalue)]);
                xlabel('PAC index');
                ylabel('Probability Density Function');
                hold all;
                maxlim = ylim;
                maxlim = maxlim(2);
                line([PACIndex PACIndex],[0 maxlim],'color','r');
            end
        end
        
        function fitRes = createCosFitModel(obj, X,Y)
            
            % written by Maya
            % Fit a simple cosine to data (a*cos(x + b)+c)
            % X - x-values
            % Y - y-values
            % Output -
            % fitRes.TC_fit = TC_fit;
            % fitRes.results = fitresult; % returned by fit()
            % fitRes.goodness = goodness;  % returned by fit()
            % fitRes.output_vals = O;  % returned by fit()
            % fitRes.opts = opts;  % fit options when runing fit()
            
            % Complex models, accounting for refractory period, etc.
            %     ft = fittype('(a * (cos(2*pi*f*x)+ 1) +  b) * exp( -abs(x)/t1 ) + c * exp( -x^2/t2^2 )');
            %     ft = fittype('(a * exp( -abs(x)/t1 ) *  (cos(2*pi*f*x)+ 1) +  b) * exp( -abs(x)/t2 ) + c * exp( -x^2/t3^2 ) + d');
            
            % DEBUG - TESTING MODEL
            % X = [0:0.001:10]';
            % Y = 10*cos(X-pi) + 4;
            
            %     create fit model and fit!
            ft = fittype('a*cos(x + b)+c');
            opts = fitoptions( ft );
            m = max(Y);
            % order of params:   a        b        c
            opts.Lower =        [0       -inf      0 ];
            opts.Upper =        [m*1.2    inf      m*1.2];
            %     opts.StartPoint = (opts.Lower + opts.Upper) ./ 2;
            opts.MaxIter = 1000000;
            opts.algorithm = 'Trust-Region';
            [fitresult, goodness, O] = fit( X, Y, ft, opts );
            TC_fit = feval(fitresult , X);
            
            % fitRes.TC_IX_to_fit = TC_IX_to_fit;
            fitRes.TC_fit = TC_fit;
            fitRes.results = fitresult;
            fitRes.goodness = goodness;
            fitRes.output_vals = O; % check what does this variable hold...
            fitRes.opts = opts;
        end
        
        function SIangles = calcPreferredAnglesSI(obj, dataPhase, dataAmp, eventTimes, sleepScoring, IIStimes, toPlot)

            %The method calculates the synchronization index (SI) according
            %to Staresina et al (2015). 
            %Input - 
            %dataPhase, dataAmp - two datasets, for the data from which slow waves
            %will be extracted (low frequency band pass) and the data from which 
            %the fast activity (i.e. spindles) will be extracted (high frequency
            %band pass). If dataAmp is empty or not provided it will be
            %set to be equal to dataPhase (i.e. SI is calculated for slow
            %waves and spindles recorded from the same electrode).
            %[[[This is obsolete currently: eventTimes - there are two options: A. it can be in the format of the output of
            %SlowWaveDetector (an array where each row is a slow wave, the
            %first column is the start index and the second is the end), B.
            %an array where each element is the index of the trough of the
            %slow wave]]]
            %IIStimes - if provided should be a cell where the first
            %element corresponds to dataPhase and the seconds to dataAmp
            %toPlot - whether to plot the polar histogram, if not provided
            %will be true by default
            %output - 
            %SIangles - an array of the SI angles per slow wave
            
            
            if nargin < 7
                toPlot = false;
            end
            
            if isempty(dataAmp)
                dataAmp = dataPhase;
            end
            
            if nargin < 5
                sleepScoring = [];
            end
            
            removeIIS = true;
            if nargin < 6 || isempty(IIStimes) || (isempty(IIStimes{1}) && isempty(IIStimes{2}))
                removeIIS = false;
            end
            
            %filter the data at the low and high frequency bands
            dataFilteredSlow = obj.bandpass(dataPhase, obj.lowLimitPhase, obj.highLimitPhase);
            dataFilteredFast = obj.bandpass(dataAmp, obj.lowLimitAmp, obj.highLimitAmp);
            
            % if sleepScoring is nonempty: leave only the
            % segments in which there was sleep at the desired stage for
            % the analysis
%             if ~isempty(sleepScoring)
%                 segLength = obj.scoringEpochDuration*obj.samplingRate;
%                 isSleep = zeros(1,length(sleepScoring)*segLength);
%                 for iEpoch = 1:length(sleepScoring)
%                     if ismember(sleepScoring(iEpoch),obj.sleepEpochs)
%                         isSleep((iEpoch-1)*segLength+1:iEpoch*segLength) = ones(1,segLength);
%                     end
%                 end
%                 %match the length of data with the length of sleepScoring -
%                 %might get rid of some data points at the end if required, assuming it's
%                 %negligible
%                 if length(isSleep)>length(dataFilteredSlow)
%                     isSleep = isSleep(1:length(dataFilteredSlow));
%                 else if length(isSleep)<length(dataFilteredSlow)
%                         dataFilteredSlow = dataFilteredSlow(1:length(isSleep));
%                         dataFilteredFast = dataFilteredFast(1:length(isSleep));
%                     end
%                 end
%                 %only leave segments of "real" sleep in data
%                 dataFilteredSlow(~isSleep) = nan;
%                 dataFilteredFast(~isSleep) = nan;
%             end
            
            isSleep = sleepScoring==obj.sleepEpochs;
            dataFilteredSlow(~isSleep(1:length(dataFilteredSlow))) = nan;
            dataFilteredFast(~isSleep(1:length(dataFilteredFast))) = nan;
            
            %remove windowAroundIIS ms before and after every IIS as
            %provided as input parameter
            if removeIIS
                winAroundIIS = obj.windowAroundIIS*obj.samplingRate/1000;
                IIStimes{1} = IIStimes{1}(IIStimes{1}<=length(dataFilteredSlow));
                for iTime = 1:length(IIStimes{1})
                    pointsBefore = min(IIStimes{1}(iTime),winAroundIIS);
                    pointsAfter = min(length(dataFilteredSlow)-IIStimes{1}(iTime),winAroundIIS);
                    dataFilteredSlow(IIStimes{1}(iTime)-pointsBefore+1:IIStimes{1}(iTime)+pointsAfter) = nan;
                end
                IIStimes{2} = IIStimes{2}(IIStimes{2}<=length(dataFilteredFast));
                for iTime = 1:length(IIStimes{2})
                    pointsBefore = min(IIStimes{2}(iTime),winAroundIIS);
                    pointsAfter = min(length(dataFilteredFast)-IIStimes{2}(iTime),winAroundIIS);
                    dataFilteredFast(IIStimes{2}(iTime)-pointsBefore+1:IIStimes{2}(iTime)+pointsAfter) = nan;
                end
            end
            
% code specific for slow waves
%check whether slow wave trough indices were provided or the
%             %start and end indices of each slow wave (assumes there is more
%             %than one slow wave in the data...) - this code is specific for
%             %one the events are slow waves
%             if size(eventTimes,1)==1 || size(eventTimes,2)==1
%                 singleIndexSW = true;
%                 nSlowWaves = length(eventTimes);
%             else
%                 singleIndexSW = false;
%                 nSlowWaves = size(eventTimes,1);
%             end
%             

            nEvents = length(eventTimes);
            SIangles = nan(1,nEvents);
            
            %each slow wave epoch contains obj.timeAroundEvent seconds before and
            %after the trough
            timeAroundEvent = obj.timeAroundEvent*obj.samplingRate;
            nPointsInEpoch = timeAroundEvent*2+1;
            
            
            %remove nan inds before the transform
            nanIndsP = isnan(dataFilteredSlow);
            nanIndsA = isnan(dataFilteredFast);
            
            dataFilteredSlow(nanIndsP) = 0;
            dataFilteredFast(nanIndsA) = 0;
            
            %calculate the time series of the phases of dataFilteredSlow
            phiSlow = angle(hilbert(dataFilteredSlow));
            %Make original NaN indices NaN again
            phiSlow(nanIndsP) = nan;
            
            %calculate the power of dataFilteredFast
            phiFast = abs(hilbert(dataFilteredFast));
            %calculate the phase of the power series
            phiFast = angle(hilbert(zscore(phiFast)));
            %Make original NaN indices NaN again
            phiFast(nanIndsA) = nan;
            
            
            %go over epoches of slow waves and calculate SI
            for iEvent = 1:nEvents
                
                pTime = eventTimes(iEvent);
                
% code specific for slow waves                 
%find the index of the slow wave trough
%                 if singleIndexSW
%                     %if it's already provided as input
%                     pTime = eventTimes(iSlowWave);
%                 else
%                     %if the limits of the slow wave are provided find the
%                     %minimum of the filtered data
%                     [~,pTime] = min(dataFilteredSlow(eventTimes(iSlowWave,1):eventTimes(iSlowWave,2)));
%                     pTime = eventTimes(iSlowWave,1)+pTime;
%                 end
                
                %skip slow waves at the edges of the data
                if (pTime-timeAroundEvent > 0) && (pTime+timeAroundEvent < length(dataPhase)) && (pTime+timeAroundEvent < length(dataAmp))
                    %the slow wave is defined as starting at timeAroundEvent data points before
                    %the trough until timeAroundSlow wave data points after the trough
                    currEpochInds = pTime-timeAroundEvent:pTime+timeAroundEvent;
                else
                    continue;
                end
                
               %filtered fast and slow data at the slow wave indices
                currPhiSlow = phiSlow(currEpochInds);
                currPhiFast = phiFast(currEpochInds);
                
                %skip when there are too many NaNs
                if sum(isnan(currPhiSlow))/length(currPhiSlow) >= obj.nanThreshPercEpoch
                    continue;
                end

                currPhiSlow(isnan(currPhiSlow)) = 0;
                currPhiFast(isnan(currPhiFast)) = 0;
                               
                %SI angle calculation
                SI = (1/nPointsInEpoch)*sum(exp(i*(currPhiSlow-currPhiFast)));
                SIangles(iEvent) = angle(SI);
                
            end
            
            if toPlot
                figure;
                polarhistogram(SIangles);
            end
        end
        
        function testArtificialData(obj,spindleData)
            %The method creates artificial data for testing PAC
            %calculation. The format of the data is a slow pure cos wave to
            %which normalized spindle is added at either the: A. troughs,
            %B. peaks, C. random locations. The spindle is provided as
            %input as spindleData and then normalized. 
            %The PAC and SI are then calcualted for each type of data (troughs / peaks / random) 
            
            %normalize spindle
            spindleData = spindleData-mean(spindleData);
            spindleData = spindleData/(max(max(spindleData),-min(spindleData)));
            spindleData = spindleData*obj.maxSpindleAmp;
                       
            nSamplePoints = round(obj.nCyclesTest*obj.samplingRate/obj.slowRateTest);
            
            t = 1:nSamplePoints;
            slowData = sin(2*pi*t*obj.slowRateTest/obj.samplingRate);
            
            %locate spindle at troughs
            troughInds = find(slowData==(-1));
            spinInds = randperm(length(troughInds));
            spinInds = spinInds(1:round(length(troughInds)*obj.percSpindlesTest));
            spinInds = troughInds(spinInds);
            
            interpTime = obj.interpTime*obj.samplingRate;
            
            %add spindles
            troughData = slowData;
            for iSpin=1:length(spinInds)
                startP = spinInds(iSpin)-round(length(spindleData)/2);
                if startP+length(spindleData)+interpTime>length(troughData) || startP-interpTime < 1
                    continue;
                end
                
                troughData(startP:startP+length(spindleData)-1) = troughData(startP:startP+length(spindleData)-1)+spindleData;
                diff1 = (troughData(startP)-troughData(startP-interpTime))/interpTime;
                troughData(startP-interpTime:startP) = [troughData(startP-interpTime):diff1:troughData(startP)];
                
                diff2 = (troughData(startP+length(spindleData)+interpTime)-troughData(startP+length(spindleData)))/interpTime;
                troughData(startP+length(spindleData):startP+length(spindleData)+interpTime) = [troughData(startP+length(spindleData)):diff2:troughData(startP+length(spindleData)+interpTime)];
            end
            
            %locate at peaks
            peakInds = find(slowData==1);
            spinInds = randperm(length(peakInds));
            spinInds = spinInds(1:round(length(peakInds)*obj.percSpindlesTest));
            spinInds = peakInds(spinInds);
                      
            peakData = slowData;
            %add spindles
            for iSpin=1:length(spinInds)
                startP = spinInds(iSpin)-round(length(spindleData)/2);
                if startP+length(spindleData)+interpTime>length(peakData) || startP-interpTime < 1
                    continue;
                end
                
                peakData(startP:startP+length(spindleData)-1) = peakData(startP:startP+length(spindleData)-1)+spindleData;
                diff1 = (peakData(startP)-peakData(startP-interpTime))/interpTime;
                peakData(startP-interpTime:startP) = [peakData(startP-interpTime):diff1:peakData(startP)];
                
                diff2 = (peakData(startP+length(spindleData)+interpTime)-peakData(startP+length(spindleData)))/interpTime;
                peakData(startP+length(spindleData):startP+length(spindleData)+interpTime) = [peakData(startP+length(spindleData)):diff2:peakData(startP+length(spindleData)+interpTime)];
            end
            
            
            %locate randomly
            spinInds = randi(length(slowData),length(spinInds));
            
            randData = slowData;
            %add spindles 
            for iSpin=1:length(spinInds)
                startP = spinInds(iSpin)-round(length(spindleData)/2);
                if startP+length(spindleData)+interpTime>length(randData) || startP-interpTime < 1
                    continue;
                end
                
                randData(startP:startP+length(spindleData)-1) = randData(startP:startP+length(spindleData)-1)+spindleData;
                diff1 = (randData(startP)-randData(startP-interpTime))/interpTime;
                randData(startP-interpTime:startP) = [randData(startP-interpTime):diff1:randData(startP)];
                
                diff2 = (randData(startP+length(spindleData)+interpTime)-randData(startP+length(spindleData)))/interpTime;
                randData(startP+length(spindleData):startP+length(spindleData)+interpTime) = [randData(startP+length(spindleData)):diff2:randData(startP+length(spindleData)+interpTime)];
            end
            
%             oldPlotDist = obj.plotDist;
%             oldPrefPhase = obj.findPreferredPhase;
%             
%             obj.plotDist = true;
%             obj.findPreferredPhase = true;
%             
%             %calculate PAC for three data sets
%             [pactrough, prefPtrough] = obj.calculatePAC(troughData);
%             title(['Trough data, PAC = ',num2str(pactrough), ' pref Phase = ',num2str(prefPtrough)]);
%             
%             %calculate PAC for three data sets
%             [pacpeak, prefPpeak] = obj.calculatePAC(peakData);
%             title(['Peak data, PAC = ',num2str(pacpeak), ' pref Phase = ',num2str(prefPpeak)]);
%            
%             %calculate PAC for three data sets
%             [pacrand, prefPrand] = obj.calculatePAC(randData);
%             title(['Peak data, PAC = ',num2str(pacrand), ' pref Phase = ',num2str(prefPrand)]);
%             
%             obj.plotDist = oldPlotDist;
%             obj.findPreferredPhase = oldPrefPhase;
            
            %calculte SI for each condition
            
            obj.calcPreferredAnglesSI(troughData, [], troughInds);
            title('SI Trough data');
            obj.calcPreferredAnglesSI(peakData, [], troughInds);
            title('SI Peak data');
            obj.calcPreferredAnglesSI(randData, [], troughInds);
            title('SI Rand data');
            
            
        end
        
        function [meanSpec, meanEpochs, stdEpochs, allSpec, nEpochs] = plotAvgSpecDiff(obj, data, epochTimes, toPlot, stimV, specFileName)
            
            % Plots the average spectrogram over all the slow waves in the data, as in Staresina et al 2015. 
            % The method receives the data and the event time indices (either as a vector of peak/trough indices, as an array of 
            % start and end times for each event).
            % The method calculates the spectrogram in the epoch around each event  timeBeforeAfterEvent before and after 
            % the event time, and subtracts from it a baseline spectrogram. The baseline spectrogram is the 
            % average spectrogram in the timeForBaseline seconds before the event. The spectrogram is calculated using the 
            % FieldTrip function ft_specest_mtmconcol (as in Staresina et al).
            
            % stimV - adding support for cases where we do not include timeBefore event in TFR (stimulation times)
            % meanSpec - average of the normalized spectograms
            % meanEpochs - average of the event epochs
            % stdEpochs - std of  of the event epochs
            % allSpec - returns all the single spectograms
            % nEpochs - number of epochs in the calculation
            
            % Sep 2020 - adding support in stimV, std of epochs
            
            if nargin < 4 
                toPlot = false;
                stimV = false;
                toSave = false;
            end
           
            if nargin < 5 || isempty(specFileName)
                toSave = false;
            else
               toSave = true;
            end


            timeBeforeAfterEvent = round(obj.timeBeforeAfterEvent*obj.samplingRate); %seconds
            timeForBaseline = round(obj.timeForBaseline*obj.samplingRate);
            
            
            if isempty(epochTimes) 
                disp('no events to calculate mean TFR for')
                allSpec = nan(1,length(obj.freqRange),timeBeforeAfterEvent*2+1);
                meanEpochs = nan(1,timeBeforeAfterEvent*2+1);
                stdEpochs = nan(1,timeBeforeAfterEvent*2+1);
                meanSpec = nan(length(obj.freqRange),timeBeforeAfterEvent*2+1); % mean over allSpec
                nEpochs = 0;
                return
            end
            
            %this is specific for slow waves - this option supports the
            %case where the input events are defined by start+end times
            %rather than one time-point
            if size(epochTimes,2) == 2 && (size(epochTimes,1) > 2) % make sure this is a matrix of start-end times
                dataFilteredSlow = obj.bandpass(data, obj.lowLimitPhase, obj.highLimitPhase);
                
                newET = zeros(1,size(epochTimes,1));
                for iEpoch = 1:size(epochTimes,1)
                    [~,newET(iEpoch)] = max(dataFilteredSlow(epochTimes(iEpoch,1):epochTimes(iEpoch,2)));
                    newET(iEpoch) = epochTimes(iEpoch,1)+newET(iEpoch);
                end
                epochTimes = newET;
            end
            

            
            timeWin = min((1./obj.freqRange)*obj.minNCycles,ones(size(obj.freqRange))*obj.minWinSizeSpec);
            if ~stimV 
                allSpec = nan(length(epochTimes),length(obj.freqRange),timeBeforeAfterEvent*2+1);
                meanEpochs = nan(length(epochTimes),timeBeforeAfterEvent*2+1);
                stdEpochs = nan(length(epochTimes),timeBeforeAfterEvent*2+1);
            else % % adding support for cases where we do not include timeBefore event in TFR
                allSpec = nan(length(epochTimes),length(obj.freqRange),timeBeforeAfterEvent+1);
                meanEpochs = nan(length(epochTimes),timeBeforeAfterEvent+1);
                stdEpochs = nan(length(epochTimes),timeBeforeAfterEvent+1);
            end
            nEpochs = 0;
            ms = [];
            for iEpoch = 1:length(epochTimes)
                if epochTimes(iEpoch)-timeBeforeAfterEvent-timeForBaseline < 1 || epochTimes(iEpoch)+timeBeforeAfterEvent > length(data)
                    continue;
                end
                
                if ~stimV
                    currInds = epochTimes(iEpoch)-timeBeforeAfterEvent:epochTimes(iEpoch)+timeBeforeAfterEvent;
                    currEpoch = data(currInds);
                    meanEpochs(iEpoch,:) = currEpoch;
                    
                    baselineInds = epochTimes(iEpoch)-timeBeforeAfterEvent-timeForBaseline:epochTimes(iEpoch)-timeBeforeAfterEvent;
                    currBL = data(baselineInds);
                    
                    if sum(isnan(currEpoch))/length(currEpoch) > obj.nanThreshPercEpoch
                        continue;
                    end
                else % adding support for cases where we do not include timeBefore event in TFR
                    
                    currInds = epochTimes(iEpoch):epochTimes(iEpoch)+timeBeforeAfterEvent;
                    currEpoch = data(floor(currInds));
                    meanEpochs(iEpoch,:) = currEpoch;
                    
                    baselineInds = epochTimes(iEpoch)-timeForBaseline:epochTimes(iEpoch);
                    currBL = data(floor(baselineInds));
                    
                    if sum(isnan(currEpoch))/length(currEpoch) > obj.nanThreshPercEpoch
                        continue;
                    end

                end
                
                nEpochs = nEpochs+1;
                currEpoch(isnan(currEpoch)) = 0;
                currSpectrum = ft_specest_mtmconvol(currEpoch, [1:length(currEpoch)]/obj.samplingRate,'taper','hanning','freqoi',obj.freqRange,'timwin',timeWin,'verbose',0);
                currSpectrum = squeeze(currSpectrum(1,1,:,:));
                currSpectrum = abs(currSpectrum);

                if sum(isnan(currEpoch))/length(currEpoch) > obj.nanThreshPercBL
                    continue;
                end
                currBL(isnan(currBL)) = 0;
                currBL = ft_specest_mtmconvol(currBL, [1:length(currBL)]/obj.samplingRate,'taper','hanning','freqoi',obj.freqRange,'timwin',timeWin,'verbose',0);
                currBL = squeeze(currBL(1,1,:,:));
                currBL = abs(currBL);
                meanBL = nanmean(currBL,2);               
                
                blMat = repmat(meanBL,1,size(currSpectrum,2));
                currSpectrum = (currSpectrum-blMat)./blMat;
                tmp = currSpectrum(:);
                tmp(tmp > 50) = nan;
                currSpectrum = reshape(tmp,size(currSpectrum,1),size(currSpectrum,2));                
                    
                allSpec(iEpoch,:,:) = currSpectrum;

            end
            if toSave
                save(specFileName,'allSpec')
            end
            
            stdEpochs = nanstd(meanEpochs,1);
            meanEpochs = nanmean(meanEpochs,1);
            meanSpec = nanmean(allSpec,1);
            meanSpec = squeeze(meanSpec(1,:,:));
            
            if toPlot
                figure;
                imagesc([1:size(meanSpec,2)]/obj.samplingRate, obj.freqRange, meanSpec);
                axis xy
                colorbar;
                title('Change in Freq Power Spectrum around epoch');
                xlabel('Time (sec)');
                ylabel('Frequency (Hz)');
                figure;
                plot(meanEpochs);
                title('Mean epoch event');
            end
            
        end
        
        function psdx = getPS(obj,segment)
            %an help method to calcualte the power spectrum of a segment
            
            segLength = length(segment);
            xdft = fft(segment);
            xdft = xdft(1:segLength/2+1);
            psdx = (1/(obj.samplingRate*segLength)) * abs(xdft).^2;
            psdx(2:end-1) = 2*psdx(2:end-1);
            psdx = 10*log10(psdx);
            
        end
    end
end