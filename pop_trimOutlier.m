% pop_trimOutlier - deploys interactive GUIs to trim outliers in channels
%                   and time-series data. Calls trimOutlier. For the
%                   channel rejection, it divides the time-series data into
%                   100 equal-sized bins, calculates std within each bin,
%                   then finally calculates means and std across the stds
%                   from all the bins. This way, the bar graph of SDs can
%                   show error bars showing SD of SD. If the error bar is
%                   short, it indicates the high std state is stationary;
%                   if the bar is long, it indicates only a small portion
%                   of the continous data has outlier values. In channel
%                   rejection plots, black is a mean value across all
%                   channels; red is a 2SD value across all the channels;
%                   blue, envelope across all the channels.

% History:
% 05/14/2024 Makoto. Mireia Torralba reported a bug about the definition of std between pop_trimOutlier() and trimOutlier(). Fixed. General maintenance for the details.
% 05/15/2019 Makoto. Supported std of std across 100 bins.GUI details fixed.
% 06/27/2014 ver 1.5 by Makoto and Clement. Drastic speed up by removing for loop (thanks Christian!)
% 04/17/2014 ver 1.4 by Makoto. Channel rejection interactive part fixed. 
% 04/02/2014 ver 1.3 by Simon Due Kamronn and Makoto. Fixed com = sprintf(... %s) to %d (thanks Simon again!)
% 04/01/2014 ver 1.2 by Simon Due Kamronn and Makoto. Debugged for epoched data (thanks Simon!)
% 03/20/2014 ver 1.1 by Frank Preston and Makoto. 'windowSize = 0;' added (thanks Frank!)
% 03/07/2014 ver 1.0 by Makoto. Former firstpassOutlierTrimmer redesigned and renamed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [EEG,com] = pop_trimOutlier(EEG)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Does EEG.times exist? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(EEG.times)
    disp('Warning: EEG.times does not exist: Created.')
    EEG.times = (1:EEG.pnts)*(1000/EEG.srate);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Is EEG.chanlocs valid? %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch ~isempty(EEG.chanlocs)
    case 0
        validChanFlag = 0;
    case 1
        switch ~isempty(EEG.chanlocs(1,1).X)
            case 0
                validChanFlag = 0;
            case 1
                validChanFlag = 1;
        end
end
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First, the topograph and bar graph for standard diviation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the main figure.
figureHandle1 = figure('position', [800 300 800 700]);
set(gcf, 'Name', 'Original Data; trimOutlier()', 'NumberTitle', 'off','Color', [0.93 0.96 1])

% Compute standard deviation after dividing the data into 100 bins. 05/14/2019 Makoto.
    %stdAllPnts = std(EEG.data(:,:),0,2);
EEGdata2D = single(EEG.data(:,:));
hundredBinUnitLength = floor(EEG.pnts/100);
truncatedData   = EEGdata2D(:,1:hundredBinUnitLength*100);
truncatedData3D = reshape(truncatedData, [size(truncatedData,1) hundredBinUnitLength 100]);
std3D = squeeze(std(truncatedData3D, 0, 2));
stdAllPnts = mean(std3D,2);

if validChanFlag == 1
    % plot topograph - all channels
    subplot(3,2,1)
    topoplot(stdAllPnts, EEG.chanlocs, 'electrodes', 'on', 'emarker', {'.','k',14,1});
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', '(\muV)', 'FontSize', 11)
    caxis([min(stdAllPnts) max(stdAllPnts)]);
    title('Channel SD','FontSize', 11)

    % plot topograph - 90 percentile
    subplot(3,2,2)
    sortedPnts = sort(stdAllPnts);
    cutoffPoint = round(length(stdAllPnts)*0.90);
    cutoffValue = sortedPnts(cutoffPoint);
    std90percentIdx   = find(stdAllPnts<=cutoffValue);
    std90percentPnts  = stdAllPnts(std90percentIdx);
    chanlocs90percent = EEG.chanlocs(std90percentIdx);
    topoplot(std90percentPnts, chanlocs90percent, 'electrodes', 'on', 'emarker', {'.','k',14,1});
    cbarHandle = colorbar;
    set(get(cbarHandle, 'title'), 'string', '(\muV)', 'FontSize', 11)
    caxis([min(std90percentPnts) max(std90percentPnts)]);
    title('90-%tile trimmed','FontSize', 11)
end

% Plot bargraph with error bars. 05/14/2019 Makoto.
switch validChanFlag
    case 0
        subplot(2,1,1)
    case 1
        subplot(3,1,2)
end
bar(stdAllPnts);
hold on
errorbar(1:size(truncatedData,1), mean(std3D,2), zeros(1, size(truncatedData,1)), std(std3D,0,2), 'linestyle', 'none', 'color', [0 0 0])

% Add annotations.
title('Channel SD with SD of SD across 100 equally-divided data bins' ,'FontSize', 11)

xlim([0.5 size(EEG.data,1)+0.5]);
barXLabel = get(gca,'XLabel');
barYLabel = get(gca,'YLabel');
set(barXLabel, 'String', 'Channel Index', 'FontSize', 11)
set(barYLabel, 'String', 'Amplitude (\muV)', 'FontSize', 11)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second, time series of channels  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute
meanAllChan = mean(EEG.data(:,:));
stdAllChan  = std( EEG.data(:,:),0,1);
minAllChan  = min( EEG.data(:,:),[],1);
maxAllChan  = max( EEG.data(:,:),[],1);
posi2SDChan = meanAllChan + 2*stdAllChan;
nega2SDChan = meanAllChan - 2*stdAllChan;

% plot
switch validChanFlag
    case 0
        subplot(2,1,2)
    case 1
        subplot(3,1,3)
end
hold on;
if length(size(EEG.data))==2
    plot(EEG.times/1000, minAllChan,  'b'); plot(EEG.times/1000, maxAllChan, 'b')
    plot(EEG.times/1000, nega2SDChan, 'r'); plot(EEG.times/1000, posi2SDChan, 'r')
    plot(EEG.times/1000, meanAllChan, 'k'); 
    timeScaleEnds = [EEG.xmin EEG.xmax];
else
    tmpTimes = 1/EEG.srate*(1:length(minAllChan));
    plot(tmpTimes, minAllChan,  'b'); plot(tmpTimes, maxAllChan, 'b')
    plot(tmpTimes, nega2SDChan, 'r'); plot(tmpTimes, posi2SDChan, 'r')
    plot(tmpTimes, meanAllChan, 'k'); 
    timeScaleEnds = [0 (1/EEG.srate)*(length(minAllChan))];
end

% annotations
xlim(timeScaleEnds)
title('Mean (Black), +/-2SD (Red), and envelope (Blue) across channels','FontSize', 11)
lineXLabel = get(gca,'XLabel');
lineYLabel = get(gca,'YLabel');
set(lineXLabel, 'String', 'Latency (s)', 'FontSize', 11)
set(lineYLabel, 'String', 'Amplitude (\muV)', 'FontSize', 11)

%%%%%%%%%%%%%%%%%%%%%%
%%% cleaning data? %%%
%%%%%%%%%%%%%%%%%%%%%%
userInput = questdlg('Proceed to clean data?','title');
if strcmp(userInput, 'No') || strcmp(userInput, 'Cancel')
    error('User canceled operation')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain channels with too large SD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
userInput = inputgui('title', 'trimOutlier()', 'geom', ...
   {{2 0 [0 0] [1 1]} {3 0 [2 0] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' sprintf('Chanel SD upper bound (microV)')} {'style' 'edit' 'string' ''}});

if isempty(userInput{1,1})
    userInput{1,1} = Inf;
end

% open new fiugre and replicate bar graph
figureHandle2 = figure('position', [800 300 800 700]);
set(gcf, 'Name', 'Channel rejection; trimOutlier()', 'NumberTitle', 'off', 'Color', [0.93 0.96 1])
h6 = subplot(3,1,1);
bar(stdAllPnts);
hold on
errorbar(1:size(truncatedData,1), mean(std3D,2), zeros(1, size(truncatedData,1)), std(std3D,0,2), 'linestyle', 'none', 'color', [0 0 0])
%h6 = gca;

close(figureHandle1)

% annotations
title('Before channel rejection (Red, upper bound)','FontSize', 11)
xlim([0.5 size(EEG.data(:,:),1)+0.5]);
barXLabel = get(gca,'XLabel');
barYLabel = get(gca,'YLabel');
set(barXLabel, 'String', 'Channels', 'FontSize', 11)
set(barYLabel, 'String', 'Amplitudes (uV)', 'FontSize', 11)


while ~isempty(userInput{1,1})
    if isnumeric(userInput{1,1})
        threshBar = userInput{1,1};
    else
        threshBar = str2num(userInput{1,1}); %#ok<*ST2NM>
    end
    
    % draw a threshold line
    axes(h6) %#ok<*LAXES>
    set(get(h6, 'title'), 'FontSize', 11)
    hold on
    h7 = line([0 size(EEG.data,1)], [threshBar threshBar], 'color', [1 0 0]); % 05/14/2019 Makoto.
    %h7 = plot(0.5:0.1:size(EEG.data(:,:),1)+0.5, threshBar, 'r');
    
    % compute
    badChanIdx  = find(stdAllPnts > threshBar);
    goodChanIdx = setdiff(1:EEG.nbchan, badChanIdx);
    stdAllPntsPost = stdAllPnts(goodChanIdx);
    disp([num2str(length(badChanIdx)) ' channels will be rejected for upper bound threshold.'])

    % plot badChan-removed on bottom 
    afterRejBarHandle = subplot(3,1,2);
    bar(stdAllPnts(goodChanIdx));
    hold on
    errorbar(1:size(truncatedData(goodChanIdx,:),1), mean(std3D(goodChanIdx,:),2), zeros(1, size(truncatedData(goodChanIdx,:),1)), std(std3D(goodChanIdx,:),0,2), 'linestyle', 'none', 'color', [0 0 0])
    plotYlim = max(mean(std3D(goodChanIdx,:),2) + std(std3D(goodChanIdx,:),0,2))*1.03;
    ylim([0 plotYlim])

    % annotations
    title('After channel rejection','FontSize', 11)
    xlim([0.5 length(stdAllPntsPost)+0.5]);
    barXLabel = get(gca,'XLabel');
    barYLabel = get(gca,'YLabel');
    set(barXLabel, 'String', 'Channel Index', 'FontSize', 11)
    set(barYLabel, 'String', 'Amplitude (\muV)', 'FontSize', 11)
    
    userInput = inputgui('title', 'trimOutlier()', 'geom', ...
       {{1 1 [0 0] [1 1]}, ...
        {2 1 [0 1] [1 1]} {2 1 [1 1] [1 1]}}, ... 
    'uilist',...
       {{'style' 'text' 'string' ['Upper bound at ' num2str(threshBar) 'microV rejects ' num2str(length(badChanIdx)) ' channels. Is it ok?']} ...
        {'style' 'text' 'string' 'If not, enter other values'} {'style' 'edit' 'string' ''}});
    if ~isempty(userInput{1,1})
        delete(h7)
        cla(afterRejBarHandle)
    end
end
channelSdUpperBound = threshBar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain channels with too small SD %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sort bars
[stdAllPntsPostSort, sortingIdx] = sort(stdAllPntsPost,'ascend');
stdAcross100bins = std(std3D(goodChanIdx,:),0,2);
stdAcross100binsSort = stdAcross100bins(sortingIdx);

% plot bars
h8 = subplot(3,1,3);
bar(stdAllPntsPostSort);
hold on
errorbar(1:length(stdAllPntsPostSort), stdAllPntsPostSort, zeros(1, length(stdAllPntsPostSort)), stdAcross100binsSort, 'linestyle', 'none', 'color', [0 0 0]); % show 20 channels from the lowest

% annotations
title('Rejecting abnormally small channels (Red, lower bound)','FontSize', 11)
xlim([0.5 20+0.5]);
barXLabel = get(gca,'XLabel');
barYLabel = get(gca,'YLabel');
set(barXLabel, 'String', 'Sorted channel index', 'FontSize', 11)
set(barYLabel, 'String', 'Amplitude (\muV)', 'FontSize', 11)

userInput = inputgui('title', 'trimOutlier()', 'geom', ...
   {{2 0 [0 0] [1 1]} {3 0 [2 0] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Chanel SD lower bound (microV)'} {'style' 'edit' 'string' ''}});

if isempty(userInput{1,1})
    goodChansPost = goodChanIdx;
    threshBar = 0;
else
    axes(h8)
    set(get(h8, 'title'), 'FontSize', 11)
    while ~isempty(userInput{1,1})
        threshBar = str2num(userInput{1,1});
        
        % draw a threshold line
        hold on
        h9 = line([0 size(EEG.data,1)], [threshBar threshBar], 'color', [1 0 0]); % 05/14/2019 Makoto.
        %h9 = plot(0.5:0.1:size(EEG.data(:,:),1)+0.5, threshBar, 'r');
        
        % compute
        badChanIdx    = find(stdAllPntsPost < threshBar);
        goodChansPost = setdiff(goodChanIdx, badChanIdx);
        disp([num2str(length(badChanIdx)) ' channels will be rejected for lower bound threshold.'])
        
        userInput = inputgui('title', 'trimOutlier()', 'geom', ...
            {{1 1 [0 0] [1 1]}, ...
            {2 1 [0 1] [1 1]} {2 1 [1 1] [1 1]}}, ...
            'uilist',...
            {{'style' 'text' 'string' ['Lower bound at ' num2str(threshBar) 'microV rejects ' num2str(length(badChanIdx)) ' channels. Is it ok?']} ...
            {'style' 'text' 'string' 'If not, enter other values'} {'style' 'edit' 'string' ''}});
        if ~isempty(userInput{1,1})
            delete(h9)
        end
    end
end
channelSdLowerBound = threshBar;
if channelSdLowerBound == 0
    channelSdLowerBound = -Inf;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% if data not continuous, return %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(size(EEG.data))==3
    pointSpreadWidth = -Inf;
    disp(sprintf('\nDapapoint rejection skipped because data is not continuous.'))
    EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, Inf, pointSpreadWidth);
    
    % eegh before terminate
    com = sprintf('EEG = trimOutlier(EEG, %d, %d, %d, %d);', channelSdLowerBound, channelSdUpperBound, Inf, pointSpreadWidth);
    EEG = eegh(com, EEG);
    close(figureHandle2)
    return
end
close(figureHandle2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reject datapoint outliers %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute again
postEEG.data    = EEG.data(goodChansPost,:);

meanAllChanPost = mean(postEEG.data);
stdAllChanPost  = std(postEEG.data,0,1);
posi2SDChanPost = meanAllChanPost + 2*stdAllChanPost;
nega2SDChanPost = meanAllChanPost - 2*stdAllChanPost;
minAllChanPost  = min(postEEG.data,[],1);
maxAllChanPost  = max(postEEG.data,[],1);

% plot time-series data
figureHandle3 = figure('position', [800 300 800 700]);
set(gcf, 'Name', 'Datapoint Rejection; trimOutlier()', 'NumberTitle', 'off', 'Color', [0.93 0.96 1])
subplot(2,1,1)
plot(EEG.times/1000, meanAllChanPost, 'k'); hold on;
plot(EEG.times/1000, minAllChanPost,  'b'); plot(EEG.times/1000, maxAllChanPost, 'b')
plot(EEG.times/1000, nega2SDChanPost, 'r'); plot(EEG.times/1000, posi2SDChanPost, 'r')
gca;

% annotations
xlim([EEG.xmin EEG.xmax])
title('Before datapoint rejection','FontSize', 11)
lineXLabel = get(gca,'XLabel');
lineYLabel = get(gca,'YLabel');
set(lineXLabel, 'String', 'Latency (s)', 'FontSize', 11)
set(lineYLabel, 'String', 'Amplitude (\muV)', 'FontSize', 11)

userInput = inputgui('title', 'trimOutlier()', 'geom', ...
   {{2 1 [0 0] [1 1]} {3 1 [2 0] [1 1]}...
    {2 1 [0 1] [1 1]} {3 1 [2 1] [1 1]}}, ... 
'uilist',...
   {{'style' 'text' 'string' 'Enter absolute threshold [microV] (if no need, press ok)'} {'style' 'edit' 'string' ''}...
    {'style' 'text' 'string' 'Point spread width for rejection [ms]'}           {'style' 'edit' 'string' ''}});

if isempty(userInput{1,1})
    threshPnts = Inf;
    pointSpreadWidth = 0;
    windowSize = 0;
else
    while ~isempty(userInput{1,1})
        threshPnts = str2num(userInput{1,1});
        
        % obtain the window size
        windowSize = str2num(userInput{1,2}); % millisecond
        windowSizeInFrame = round(windowSize/(1000/EEG.srate)); % frame
        if isempty(windowSizeInFrame) % 05/14/2019 Makoto. Updated.
            disp(sprintf('\n'))
            warning(sprintf('Spread width not provided. Default minimum of single frame width is used.'))
            windowSizeInFrame = 1000/EEG.srate;
        end
        
        % compute bad datapoints
        absMinMaxAllChan = max([abs(minAllChanPost); abs(maxAllChanPost)],[],1);
        badPoints  = absMinMaxAllChan > threshPnts;
                
        % expand badPoints
        badPointsExpanded = logical(conv(single(badPoints), ones(1,windowSizeInFrame), 'same'));
        
        % start Christian's impressive code
        rejectDataIntervals = reshape(find(diff([false badPointsExpanded false])),2,[])';
        rejectDataIntervals(:,2) = rejectDataIntervals(:,2)-1;

        % plot how much data will be rejected in sec
        badPointsInSec = length(find(badPointsExpanded))*1000/EEG.srate/1000;
        sprintf('%2.1f sec of data will be rejected.', badPointsInSec);
        sprintf('%1.0f boundary will be made.', size(rejectDataIntervals,1));

        goodPoints = setdiff(1:EEG.pnts, find(badPointsExpanded));
        meanAllChanPostPost = meanAllChanPost(goodPoints);
        posi2SDChanPostPost = posi2SDChanPost(goodPoints);
        nega2SDChanPostPost = nega2SDChanPost(goodPoints);
        minAllChanPostPost  = minAllChanPost(goodPoints);
        maxAllChanPostPost  = maxAllChanPost(goodPoints);
        
        % plot badChan-removed on bottom
        subplot(2,1,2)
        timesPost = EEG.times(1:length(goodPoints))/1000;
        hold on;
        plot(timesPost, minAllChanPostPost,  'b'); plot(timesPost, maxAllChanPostPost, 'b')
        plot(timesPost, nega2SDChanPostPost, 'r'); plot(timesPost, posi2SDChanPostPost, 'r')
        plot(timesPost, meanAllChanPostPost, 'k');
        ylim([-threshPnts threshPnts]);
        h13 = gca;
        
        % annotations
        title('After datapoint rejection','FontSize', 11)
        xlim([0 length(timesPost)*1000/EEG.srate/1000]);
        barXLabel = get(gca,'XLabel');
        barYLabel = get(gca,'YLabel');
        set(barXLabel, 'String', 'Latency (s)', 'FontSize', 11)
        set(barYLabel, 'String', 'Amplitude (\muV)', 'FontSize', 11)
        
        userInput = inputgui('title', 'trimOutlier()', 'geom', ...
            {{1 1 [0 0] [1 1]}, ...
             {3 1 [0 1] [1 1]} {3 1 [1 1] [1 1]},...
             {3 1 [0 2] [1 1]} {3 1 [1 2] [1 1]}},...
            'uilist',...
            {{'style' 'text' 'string' sprintf('Threshold %2.0fmicroV point spread %2.0fms rejects %2.1fsec, creates %1.0f boundaries. Is it ok?',threshPnts,windowSize,badPointsInSec,size(rejectDataIntervals,1))} ...
             {'style' 'text' 'string' 'If not, enter threshold [microV]'} {'style' 'edit' 'string' ''}...
             {'style' 'text' 'string' 'Point spread width [ms]'} {'style' 'edit' 'string' ''}});
        if ~isempty(userInput{1,1})
            cla(h13)
        end
    end
end
close(figureHandle3)

% reject data using user-defined thresholds
EEG = trimOutlier(EEG, channelSdLowerBound, channelSdUpperBound, threshPnts, windowSize);

% eegh before terminate
com = sprintf('EEG = trimOutlier(EEG, %s, %s, %s, %s);', num2str(channelSdLowerBound), num2str(channelSdUpperBound), num2str(threshPnts), num2str(windowSize));
EEG = eegh(com, EEG);