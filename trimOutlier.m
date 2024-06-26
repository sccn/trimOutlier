% trimOutlier() - rejects 1.channels that are below or above the specified
%                 SD, and 2.datapoints that are above the specified
%                 threshold. Point spread width [ms] determines the range
%                 for rejection. For the calculation of channel-wise SD,
%                 the data are first divided into 100 bins with equal
%                 length with trimming, calculate SD for each bin, then
%                 calculate SD across the 100 SD values from all the bins.
%                 This way, the bar graph of SDs can show error bars
%                 showing SD of SD. If the error bar is short, it indicates
%                 the high std state is stationary; if the bar is long, it
%                 indicates only a small portion of the continous data has
%                 outlier values. The name was suggested by John Iversen.
%                 
% Usage:
%   >> EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, amplitudeThreshold, pointSpreadWidth);

% History:
% 05/14/2024 Makoto. Mireia reported a bug about the definition of std between pop_trimOutlier() and trimOutlier(). Fixed.
% 02/27/2019 Makoto. line 102 logical(ones(EEG.pnts,1)) fixed.
% 07/05/2018 Makoto. Cleaning mask saved under EEG.etc.trimOutlier
% 06/27/2014 ver 1.4 by Makoto and Clement. Drastic speed up by removing for loop (thanks Christian!) Displays log (requested by Kathleen VanBenthem)
% 04/17/2014 ver 1.3 by Makoto. min(badPntsStart)=1, zero not allowed.
% 04/01/2014 ver 1.2 by Makoto. Check inputs.
% 03/26/2014 ver 1.1 by Makoto. Debug and simplify datapoint rejection
% 03/07/2014 ver 1.0 by Makoto. Former firstpassOutlierTrimmer redesigned and renamed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 09/06/2013 ver 1.5 by Makoto. No datapoint rejection reflects channel rejection.
% 08/05/2013 ver 1.4 by Makoto. Supported 3-D data (except for datapoint rejection)
% 06/27/2013 ver 1.3 by Makoto. Error message for inputting 3-D data.
% 06/13/2013 ver 1.2 by Makoto. Scalp topos added (when valid channel data exist)
% 05/22/2013 ver 1.1 by Makoto. Minor brush up.
% 05/16/2013 ver 1.0 by Makoto. Created.

% Author: Makoto Miyakoshi, SCCN,INC,UCSD 2013; Cincinnati Children's Hospital Medical Center, 2024.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, amplitudeThreshold, pointSpreadWidth)

if ~(nargin==5)
    error('trimOutlier() requires 5 input arguments.')
end

% Remove bad channels with SD. Mireia reported an inconsistency of the
% calculation of the std between pop_trimOutlier() and trimOutlier(). Now
% the definition of std is based on the one defined in the
% pop_trimOutlier(). (05/14/2024 Makoto)
    %stdAllPnts  = std(EEG.data(:,:),0,2);
EEGdata2D = single(EEG.data(:,:));
hundredBinUnitLength = floor(EEG.pnts/100);
truncatedData   = EEGdata2D(:,1:hundredBinUnitLength*100);
truncatedData3D = reshape(truncatedData, [size(truncatedData,1) hundredBinUnitLength 100]);
std3D = squeeze(std(truncatedData3D, 0, 2));
stdAllPnts = mean(std3D,2);

badChanMask = (stdAllPnts < channelSdLowerBound) | (stdAllPnts > channelSdUpperBound);
badChanIdx  = find(badChanMask);
if any(badChanIdx)
    badChanName = {EEG.chanlocs(badChanIdx).labels};
    EEG = pop_select(EEG, 'nochannel', badChanIdx);
    
    % Save the clean channel mask.
    EEG.etc.trimOutlier.cleanChannelMask = ~badChanMask;
    
    % display log
    disp(sprintf('\nThe following channels were removed:'))
    disp(badChanName)
else
    % Save the clean channel mask.
    EEG.etc.trimOutlier.cleanChannelMask = logical(ones(EEG.nbchan,1));
    
    disp(sprintf('\nNo channel removed.'))
end

% return if 3-D
if length(size(EEG.data))==3
    disp('Epoched data detected: datapoint rejection is skipped.')
    return
end

%% remove bad datapoints

% obtain the window size
windowSize = pointSpreadWidth; % millisecond
windowSizeInFrame = round(windowSize/(1000/EEG.srate)); % frame

% compute bad datapoints
absMinMaxAllChan = max([abs(min(EEG.data(:,:))); abs(max(EEG.data(:,:)))],[],1);
badPoints  = absMinMaxAllChan > amplitudeThreshold;

if any(badPoints)
    % expand badPoints
    badPointsExpanded = logical(conv(single(badPoints), ones(1,windowSizeInFrame), 'same'));
    
    % start Christian's impressive code
    rejectDataIntervals = reshape(find(diff([false badPointsExpanded false])),2,[])';
    rejectDataIntervals(:,2) = rejectDataIntervals(:,2)-1;
    
    % reject them
    EEG = pop_select(EEG, 'nopoint', [rejectDataIntervals(:,1) rejectDataIntervals(:,2)]);
    
    % Save the clean data points.
    EEG.etc.trimOutlier.cleanDatapointMask = ~badPointsExpanded;
    
    % display log
    badPointsInSec = length(find(badPointsExpanded))*1000/EEG.srate/1000; %#ok<*NASGU>
    disp(sprintf('\n%2.0fuV threshold with %2.0fms spreading rejected %2.1fsec data, added %1.0f boundaries.', amplitudeThreshold, windowSize, badPointsInSec, size(rejectDataIntervals,1)));
else
    % Save the clean data points.
    EEG.etc.trimOutlier.cleanDatapointMask = logical(ones(EEG.pnts,1));
    
    disp('No datapoint rejected.');
end

disp('trimOutlier done. The masks for clean channels and data points are stored under EEG.etc.trimOutlier.')