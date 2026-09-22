clear; clc; close all; imtool close all;
mainDir = 'G:\Bitong\Mode 1 cells_RFP cannot be detected 1'; numb = 2;
% mainDir = 'D:\Data\Bitong\20230726'; numb = 2;

% Adjusting value for minimum and maximum gfp/rfp signals, colocalisation
% and mobility
minGFP = 4;
maxGFP = 50;
minRFP = 4;
maxRFP = 50;
distTol = 4;
minBeads = 4;
gfpMobTol = 12; % GFP Mobility Tolerance

% Adjust intensity parameters for defining what is real gfp/rfp signals
gfpTgt = 14; % target gfp intensity group
gfpRng = 20; % range of gfp signals grouping
rfpTgt = 16; % target rfp intensity group
rfpRng = 20; % range of rfp signals grouping

frameIntervals = 1; %30; % number of SMT frames before GPF/RFP
frameTime = 0; %0.606; % time taken to cycle through GFP/RFP in seconds

% Adjust what figures to plot
plotGFP = 0;
plotRFP = 0;
plotNumber = 10; % maximum frame numbers to plot
plotStart = 1; % frames to start plotting
plotGap = 1; % set as 1 to plot everything
plotHeatMap = 0;
plotAllGFPs = 0;

% Specific GFP/RFP to use after run the analysis at least once
GFP2Use = [];
RFP2Use = [];

% Set x, y coordinates of GFP position to start with
manualGFP = [];

bleach_rate = 6.37;
LocalizationError = -6.5; % Localization Error: -6 = 10^-6
EmissionWavelength = 580; % wavelength in nm; consider emission max and filter cutoff
ExposureTime = 500; % in milliseconds
NumDeflationLoops = 0; % Generaly keep this to 0; if you need deflation loops, you are imaging at too high a density;
MaxExpectedD = 0.2; % The maximal expected diffusion constant for tracking in units of um^2/s;
NumGapsAllowed = 1; % the number of gaps allowed in trajectories

%%% DEFINE STRUCTURED ARRAY WITH ALL THE SPECIFIC SETTINGS FOR LOC AND TRACK
% imaging parameters
impars.PixelSize=0.11; % um per pixel
impars.psf_scale=1.35; % PSF scaling
impars.wvlnth= EmissionWavelength/1000; %emission wavelength in um
impars.NA=1.49; % NA of detection objective
impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth / 1)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels
impars.FrameRate= ExposureTime/1000; %secs
impars.FrameSize= ExposureTime/1000; %secs

% localization parameters
locpars.wn=9; %detection box in pixels
locpars.errorRate= LocalizationError; % error rate (10^-)
locpars.dfltnLoops= NumDeflationLoops; % number of deflation loops
locpars.minInt=0; %minimum intensity in counts
locpars.maxOptimIter= 50; % max number of iterations
locpars.termTol= -2; % termination tolerance
locpars.isRadiusTol=false; % use radius tolerance
locpars.radiusTol=50; % radius tolerance in percent
locpars.posTol= 1.5;%max position refinement
locpars.optim = [locpars.maxOptimIter,locpars.termTol,locpars.isRadiusTol,locpars.radiusTol,locpars.posTol];
locpars.isThreshLocPrec = false;
locpars.minLoc = 0;
locpars.maxLoc = inf;
locpars.isThreshSNR = false;
locpars.minSNR = 0;
locpars.maxSNR = inf;
locpars.isThreshDensity = false;

% tracking parameters
trackpars.trackStart=1;
trackpars.trackEnd=inf;
trackpars.Dmax= MaxExpectedD;
trackpars.searchExpFac=1.2;
trackpars.statWin=10;
trackpars.maxComp=3;
trackpars.maxOffTime=NumGapsAllowed;
trackpars.intLawWeight=0.9;
trackpars.diffLawWeight=0.5;

%% Categorise files into state and cells
allFiles = dir(fullfile(mainDir, '**/*.tif'));
allFiles = allFiles(~[allFiles.isdir]);
allCells = repmat(struct('name', [], 'bead', [], 'cell', 0, 'state', '0', 'folder', '0'), 1, 1);
fileName = strsplit(allFiles(1).name, '_');
currentCell = regexp(fileName{numb},'[0-9\.]+','match'); %
currentCell = str2double(currentCell{1}); cellFiles = strings(3, 1); beadFiles = strings(3, 1); cellState = '0'; cellFolder = "0";
for n = 1 : length(allFiles)
  fileName = strsplit(allFiles(n).name, '_');
  cellNumber = regexp(fileName{numb},'[0-9\.]+','match'); %
  cellNumber = str2double(cellNumber{1});
  if cellNumber ~= currentCell
    if convertCharsToStrings(cellFolder) ~= "0"
      allCells(currentCell).name = cellFiles;
      allCells(currentCell).bead = beadFiles;
      allCells(currentCell).cell = currentCell;
      allCells(currentCell).state = cellState;
      allCells(currentCell).folder = cellFolder;
    end
    currentCell = cellNumber;
    cellFiles = strings(3, 1); beadFiles = strings(3, 1); cellState = '0'; cellFolder = "0";
  end
  if contains(allFiles(n).folder, 'ON')
    cellState = 'ON';
  elseif contains(allFiles(n).folder, 'OFF')
    cellState = 'OFF';
  end
  folderName = strsplit(allFiles(n).folder, filesep);
  if convertCharsToStrings(folderName{end}) == convertCharsToStrings(fileName{numb}) %
    cellFolder = allFiles(n).folder;
  end
  if contains(allFiles(n).name, '488 nm.tif')
    if contains(upper(allFiles(n).name), 'BEADS')
      beadFiles(1) = allFiles(n).name;
    else
      cellFiles(1) = allFiles(n).name;
    end
  elseif contains(allFiles(n).name, '561 nm.tif')
    if contains(upper(allFiles(n).name), 'BEADS')
      beadFiles(2) = allFiles(n).name;
    else
      cellFiles(2) = allFiles(n).name;
    end
  elseif contains(allFiles(n).name, '640.tif')
    if contains(upper(allFiles(n).name), 'BEADS')
      beadFiles(3) = allFiles(n).name;
    else
      cellFiles(3) = allFiles(n).name;
    end
  end
end
allCells(currentCell).name = cellFiles;
allCells(currentCell).bead = beadFiles;
allCells(currentCell).cell = currentCell;
allCells(currentCell).state = cellState;
allCells(currentCell).folder = cellFolder;
clear cellFiles cellState cellFolder currentCell fileName folderName cellNumber

%%
addpath(genpath(['.' filesep 'Batch_MTT_code' filesep])); 
finalData = repmat(struct('FileName', [], 'Cell_State_Counts', [], 'SMT_Spots', [], 'Extra', []), length(allCells), 1);
dwellDataCompiled = cell(length(allCells) + 1, 3);
dwellDataCompiled{1, 1} = 'Long Dwell Time (s)';
dwellDataCompiled{1, 2} = 'Fraction';
dwellDataCompiled{1, 3} = 'Short Dwell Time (s)';
for n = 1 :  length(allCells)
  if sum(size(allCells(n).name)) == 0 || sum(size(allCells(n).bead)) == 0 || any(allCells(n).name == "") || any(allCells(n).bead == "")
    continue
  end
  gfpFile = convertStringsToChars(allCells(n).name(1));
  rfpFile = convertStringsToChars(allCells(n).name(2));
  fishFile = convertStringsToChars(allCells(n).name(3));

  gfpAFile = convertStringsToChars(allCells(n).bead(1));
  rfpAFile = convertStringsToChars(allCells(n).bead(2));
  smtAFile = convertStringsToChars(allCells(n).bead(3));

  dirPath = allCells(n).folder;
  mkdir([dirPath, filesep, 'Data']);
  stack = tiffread([dirPath, filesep, smtAFile]);
  smtRaw = stack(1).data;
  
  smtThresh = imquantize(smtRaw, multithresh(smtRaw, 5)) >= 2; % needs to account for 2 beads situation

  [smtMap, smtCounts] = bwlabel(smtThresh);
  smtCenter = zeros(smtCounts, 2);
  smtReal = 0;
  for m = 1 : smtCounts
    [y, x] = find(smtMap == m);
    if length(y) > minBeads
      smtReal = smtReal + 1;
      smtCenter(smtReal, :) = [mean(x), mean(y)];
    end
  end
  if smtReal < smtCounts
    smtCenter(smtReal+1 : end, :) = [];
    smtCounts = smtReal;
  end
  
%   [y, x] = find(smtThresh == 1);
% %   idx = boundary(x, y);
% %   figure; plot(x(idx), y(idx)); hold on;
% %   plot(mean(x), mean(y), '+')
%   smtCenter = [mean(x), mean(y)];

  rfpStack = tiffread([dirPath, filesep, rfpAFile]);
  rfpRaw = rfpStack(1).data;
  rfpThresh = imquantize(rfpRaw, multithresh(rfpRaw, 2)) >= 3; % rfp beads

  [rfpMap, rfpCounts] = bwlabel(rfpThresh);
  rfpCenter = zeros(rfpCounts, 2);
  for m = 1 : smtCounts
    [y, x] = find(rfpMap == m);
    rfpCenter(m, :) = [mean(x), mean(y)];
  end
%   [y, x] = find(rfpThresh == 1);
% %   idx = boundary(x, y);
% %   figure; plot(x(idx), y(idx)); hold on;
% %   plot(mean(x), mean(y), '+')
%   rfpCenter = [mean(x), mean(y)];

  gfpStack = tiffread([dirPath, filesep, gfpAFile]);
  gfpRaw = gfpStack(1).data;
  gfpThresh = imquantize(gfpRaw, multithresh(gfpRaw, 2)) >= 3;

  [gfpMap, gfpCounts] = bwlabel(gfpThresh);
  gfpCenter = zeros(gfpCounts, 2);
  for m = 1 : smtCounts
    [y, x] = find(gfpMap == m);
    gfpCenter(m, :) = [mean(x), mean(y)];
  end
%   [y, x] = find(gfpThresh == 1);
% %   idx = boundary(x, y);
% %   figure; plot(x(idx), y(idx)); hold on;
% %   plot(mean(x), mean(y), '+')
%   gfpCenter = [mean(x), mean(y)];

  % Using gfp as the reference point
  smtCorrection = smtCenter - gfpCenter; smtCorrection = round(smtCorrection);
  rfpCorrection = rfpCenter - gfpCenter; rfpCorrection = round(rfpCorrection);

  % Correcting image based on reference
  smtNew = zeros(size(smtRaw));
  newRow = 1 : size(smtRaw, 1); newCol = 1 : size(smtRaw, 2);
  newRow = newRow - smtCorrection(1); newCol = newCol - smtCorrection(2);
  oldRow = newRow > 0; oldCol = newCol > 0;
  smtNew(newRow(oldRow), newCol(oldCol)) = smtRaw(oldRow, oldCol);

  rfpNew = zeros(size(rfpRaw));
  newRow = 1 : size(rfpNew, 1); newCol = 1 : size(rfpNew, 2);
  newRow = newRow - rfpCorrection(1); newCol = newCol - rfpCorrection(2);
  oldRow = newRow > 0; oldCol = newCol > 0;
  rfpNew(newRow(oldRow), newCol(oldCol)) = rfpRaw(oldRow, oldCol);

  smtCorrection = round(mean(smtCorrection, 1));
  rfpCorrection = round(mean(rfpCorrection, 1));

  stackGFP = tiffread([allCells(n).folder, filesep, gfpFile]);
  stackRFP = tiffread([allCells(n).folder, filesep, rfpFile]);
  stackFISH = tiffread([allCells(n).folder, filesep, fishFile]);

  cellData = zeros(length(stackGFP), 7); % gfp, rfp, colocalise, xgfp, ygfp, gfp signal, rfp signal
  fishData = zeros(size(stackFISH(1).data, 1), size(stackFISH(1).data, 2), length(stackGFP));

  % Huygen correction
%   correction = readtable([allCells(n).folder, filesep, adFiles(n).name], opts);
%   correction = table2array(correction);
%   smtCorrection = [sign(correction(1, 1)) * correction(1, 1), sign(correction(1, 2)) * correction(1, 2)];
%   rfpCorrection = [sign(correction(2, 1)) * correction(2, 1), sign(correction(2, 2)) * correction(2, 2)];

  % Going through time-frames
  gfpUseInd = 1; rfpUseInd = 1;
  % Precompute frames to plot GFP/RFP
  plotFrame = plotStart : plotGap : length(stackGFP);
  for m = 1 : length(stackGFP)
    gfpRaw = stackGFP(m).data;
    rfpRaw = stackRFP(m).data;
    fishRaw = stackFISH(m).data;

    fishData(:, :, m) = double(fishRaw);
    try
%       gfpLoc = imquantize(gfpRaw, multithresh(gfpRaw, 4)) >= 4;
      gfpLoc = imquantize(gfpRaw, multithresh(gfpRaw, gfpRng)) >= gfpTgt; %
    catch
      gfpLoc = imquantize(gfpRaw, multithresh(gfpRaw, 5)) >= 5;
    end
    try
%       rfpLoc = imquantize(rfpRaw, multithresh(rfpRaw, 4)) >= 4;
      rfpLoc = imquantize(rfpRaw, multithresh(rfpRaw, rfpRng)) >= rfpTgt;
    catch
      rfpLoc = imquantize(rfpRaw, multithresh(rfpRaw, 5)) >= 5;
    end
%     if m == 121
%       1; % check gfp/rfp catch point
%     end
    imageSize = size(stackGFP(m).data, 1) * size(stackGFP(m).data, 2);
    gfpTol = 0.01 * imageSize; % was 0.1
    if sum(gfpLoc(:)) > gfpTol
      gfpLoc = 0;
    end
    if sum(rfpLoc(:)) > gfpTol
      rfpLoc = 0;
    end
    [gfpPos, gfpCount] = bwlabel(gfpLoc);
    gfpSize = histcounts(gfpPos, gfpCount + 1);
    gfpData = find(gfpSize >= minGFP & gfpSize <= maxGFP);

    [rfpPos, rfpCount] = bwlabel(rfpLoc);
    rfpSize = histcounts(rfpPos, rfpCount + 1);
    rfpData = find(rfpSize >= minRFP & rfpSize <= maxRFP);

    if ~isempty(gfpData)
%       cellData(m, 1) = 1;
      % unsure to do last gfp found location or last gfp colocalise
      % TODO: current setup is last gfp seen location
%       [yGFP, xGFP] = find(gfpPos == gfpData(1) - 1);
%       centerGFP = [mean(xGFP), mean(yGFP)];
      %
      if plotGFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
        gfpFig = figure;
        pause(0.5);
        [counts, intensity] = histcounts(gfpRaw(:));
        intOI = intensity(counts > 10); % background adjust
        imshow(imadjust(gfpRaw, [intOI(2)/65535, intOI(end)/65535])); hold on;
      end
%       plot(centerGFP(1), centerGFP(2), 'rx');
%       savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_GFP_Frame_', num2str(m), '.fig']);
      %
%       cellData(m, 4:5) = [mean(xGFP), mean(yGFP)];
      %
    end
    if ~isempty(rfpData)
      cellData(m, 2) = 1;
      %
      if plotRFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
        rfpFig = figure;
        pause(0.5);
        [counts, intensity] = histcounts(rfpRaw(:));
        intOI = intensity(counts > 10); % background adjust
        imshow(imadjust(rfpRaw, [intOI(2)/65535, intOI(end)/65535])); hold on;
      end
      %
    end

    if isempty(GFP2Use)
      for k = 1 : length(gfpData)
        [yGFP, xGFP] = find(gfpPos == gfpData(k) - 1);
        centerGFP = [mean(xGFP), mean(yGFP)];
        if m > 1
          gfpPre = find(cellData(:, 1) == 1);
          if length(gfpPre) > 1
            gfpPre = gfpPre(end);
          end
          if sum(pdist([centerGFP; cellData(gfpPre, 4:5)])) > gfpMobTol
            continue
          end
        else
          if sum(pdist([manualGFP; centerGFP])) > gfpMobTol && ~isempty(manualGFP)
            continue
          end
        end
        cellData(m, 1) = 1;
        cellData(m, 4:5) = [mean(xGFP), mean(yGFP)];
        %
        if plotGFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
          figure(gfpFig);
          plot(centerGFP(1), centerGFP(2), 'rx');
          text(centerGFP(1), centerGFP(2), num2str(k), 'Color', 'red');
        end
        gfpInd = sub2ind(size(gfpRaw), yGFP, xGFP);
        cellData(m, 6) = mean(gfpRaw(gfpInd));
        %
        for p = 1 : length(rfpData)
          [yRFP, xRFP] = find(rfpPos == rfpData(p) - 1);
          %         plot(xRFP, yRFP, 'yx');
          centerRFP = [mean(xRFP), mean(yRFP)];
          centerRFP = centerRFP - rfpCorrection;
          %
          if plotRFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
            figure(rfpFig);
            plot(centerRFP(1), centerRFP(2), 'bx');
            text(centerRFP(1), centerRFP(2), num2str(p), 'Color', 'blue');
          end
          rfpInd = sub2ind(size(rfpRaw), yRFP, xRFP);
%           cellData(m, 7) = mean(rfpRaw(rfpInd));
          %
          if pdist([centerGFP; centerRFP]) < distTol
            cellData(m, 3:5) = [1, mean(xGFP), mean(yGFP)];
            cellData(m, 7) = mean(rfpRaw(rfpInd)); % moving rfp int to only record if it colocalise with gfp
            break
          end
        end
      end
      % if no gfp detected, record last seen
      if sum(cellData(m, 4:5)) == 0 && m > 1
        gfpPre = find(cellData(:, 3) == 1);
        if length(gfpPre) > 1
          gfpPre = gfpPre(end);
          cellData(m, 4:5) = cellData(gfpPre, 4:5);
        end
      end
    else
      if ~isempty(gfpData)
        k = GFP2Use(gfpUseInd); gfpUseInd = gfpUseInd + 1;
        [yGFP, xGFP] = find(gfpPos == gfpData(k) - 1);
        centerGFP = [mean(xGFP), mean(yGFP)];
        %
        if plotGFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
          figure(gfpFig);
          plot(centerGFP(1), centerGFP(2), 'rx');
          text(centerGFP(1), centerGFP(2), num2str(k), 'Color', 'red');
        end
        gfpInd = sub2ind(size(gfpRaw), yGFP, xGFP);
        cellData(m, 6) = mean(gfpRaw(gfpInd));
        %
        if ~isempty(rfpData)
          p = RFP2Use(rfpUseInd); rfpUseInd = rfpUseInd + 1;
          [yRFP, xRFP] = find(rfpPos == rfpData(p) - 1);
          %         plot(xRFP, yRFP, 'yx');
          centerRFP = [mean(xRFP), mean(yRFP)];
          centerRFP = centerRFP - rfpCorrection;
          %
          if plotRFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
            figure(rfpFig);
            plot(centerRFP(1), centerRFP(2), 'bx');
            text(centerRFP(1), centerRFP(2), num2str(p), 'Color', 'blue');
          end
          rfpInd = sub2ind(size(rfpRaw), yRFP, xRFP);
%           cellData(m, 7) = mean(rfpRaw(rfpInd));
          %
          if pdist([centerGFP; centerRFP]) < distTol
            cellData(m, 3:5) = [1, mean(xGFP), mean(yGFP)];
            cellData(m, 7) = mean(rfpRaw(rfpInd)); % moving rfp int to only record if it colocalise with gfp
            break
          end
        end
      end
    end
    
    % Can check RFP and GFP points here
    if ~isempty(gfpData) && plotGFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
      figure(gfpFig);
      pause(1);
      savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_GFP_Frame_', num2str(m), '.fig']);
    end
    if ~isempty(rfpData) && plotRFP && plotNumber > 0 && m < plotNumber && m >= plotStart && ismember(m, plotFrame)
      figure(rfpFig);
      pause(1);
      savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_RFP_Frame_', num2str(m), '.fig']);
    end
  end
  finalData(n).FileName = gfpFile;
  finalData(n).Cell_State_Counts = cellData;

  data = localizeParticles_ASH("", impars, locpars, fishData);
  tracks = [data.frame, data.ctrsX, data.ctrsY];
  tracks(:, 2:3) = tracks(:, 2:3) - smtCorrection;

  % Dwell Time
  data=buildTracks2_ASH("", "", data, impars, locpars, trackpars, data.ctrsN, fishData);
  data_cell_array = data.tr;
  trackedPar = struct;
  for i=1:length(data_cell_array)
    %convert to um:
    trackedPar(1,i).xy =  impars.PixelSize .* data_cell_array{i}(:,1:2);
    trackedPar(i).Frame = data_cell_array{i}(:,3);
    trackedPar(i).TimeStamp = impars.FrameRate.* data_cell_array{i}(:,3);
    trackedPar(i).TimeStamp = trackedPar(i).TimeStamp + (floor((trackedPar(i).Frame - 1) / frameIntervals) * frameTime);
  end

  tracklength = zeros(length(trackedPar), 1);
  for i = 1 : length(trackedPar)
%     tracklength(i) = length(trackedPar(i).Frame);
    tracklength(i) = trackedPar(i).TimeStamp(end) - trackedPar(i).TimeStamp(1) + (ExposureTime / 1000);
  end

  % Debug - SMT position check
%   frameCheck = 1;
%   imtool(fishData(:, :, frameCheck));
%   figure; imshow(fishDataExport); hold on;
% %   tracksE = []; % time, x, y, Frame, trajectory number
% %   for i = 1 : length(trackedPar)
% %     tracksE = [tracksE; [trackedPar(i).TimeStamp, trackedPar(i).xy, trackedPar(i).Frame, repelem(i, size(trackedPar(i).xy, 1))']];
% %   end
% %   plot(tracksE(tracksE(:, 4) == frameCheck, 2), tracksE(tracksE(:, 4) == frameCheck, 3), 'rx');
%   plot(tracks(tracks(:, 1) == frameCheck, 2), tracks(tracks(:, 1) == frameCheck, 3), 'rx');
  % End debug - SMT position check

%   tracklength = (ExposureTime / 1000) * tracklength;
  a = tracklength;
  binsize = 1;
  binwindows = 0:binsize:1000;
  [num,~]=hist(a,binwindows);
  cdf_n = [sum(num), sum(num) - cumsum(num)]/sum(num);
  logind=find(cdf_n>0.01);
  binwindows = 0:binsize:(logind(end)+1)*binsize;
  [num,xout]=hist(a,binwindows);
  histogram(a, binwindows);
  %   hist(a,binwindows);
  cdf_n = [sum(num), sum(num) - cumsum(num)]/sum(num);
  cdf_n=cdf_n(1:end-1);

  % xout is the time, n is the count.
  % Please define the bleach rate (unit in second) from previsous measurement
  myfun = @(x) cdf_n(1:end) - x(1) - x(2)*exp(-xout(1:end)/x(3));%1-component fitting
  myfun2 = @(r) cdf_n - r(1)*exp(-xout/r(2))-(1-r(1))*exp(-xout/r(3));%2-component fitting
  x0=[0.05, 1, 5];
  r0=[0.5,8,1.2];
  options = optimset('Algorithm',{'levenberg-marquardt',0.000001});
  x=lsqnonlin(myfun,x0,[],[],options);
  r=lsqnonlin(myfun2,r0,[],[],options);
  figure;
  hold on
  plot(xout, cdf_n,'o');
  xfit=xout(1):0.9:max(xout);
  yfit= x(1) + x(2)*exp(-xfit/x(3));
  plot(xfit,yfit,'-r');
  title('1-Component Fitting', 'interpreter', 'latex');
  hold off
  f1=figure;
  hold on
  plot(xout, cdf_n,'o');
  xfit=0:0.9:max(xout);
  yfit= r(1)*exp(-xfit/r(2))+(1-r(1))*exp(-xfit/r(3));
  plot(xfit,yfit,'-r')
  a = r(1);
  %   b = r(2);
  c = r(3);
  TrueR = (bleach_rate - 1)/(bleach_rate/x(3) -1);
  TrueR2 = (bleach_rate - 1)/(bleach_rate/max(r(2),r(3)) -1);
  title('2-Components Fitting', 'interpreter', 'latex');
%   savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_dwell_time.fig']);
  dwellTimeToSave = cell(3, 2);
  dwellTimeToSave{1, 1} = 'Long Dwell Time (s)'; dwellTimeToSave{1, 2} = TrueR2;
  dwellTimeToSave{2, 1} = 'Fraction'; dwellTimeToSave{2, 2} = a;
  dwellTimeToSave{3, 1} = 'Short Dwell Time (s)'; dwellTimeToSave{3, 2} = c;
  writecell(dwellTimeToSave, [allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_dwell_time.xls'])
  polyfit_str = ['R1= ' num2str(TrueR2) ' Fraction= ' num2str(a)...
    'R2= ' num2str(c)];
  %
  dwellDataCompiled{n + 1, 1} = TrueR2; dwellDataCompiled{n + 1, 2} = a; dwellDataCompiled{n + 1, 3} = c;
  %
  text(1,1,polyfit_str,'FontSize',18);

  fprintf('The residence time is %d before correction of one-component fitting.\n', x(3));
  fprintf('The true residence time is %d of one component fitting after correction.\n', TrueR);
  fprintf('The residence time are %d and %d with fraction %d of two-components fitting.\n', r(2),r(3),r(1));
  fprintf('The true residence time is %d of two components fitting after correction.\n', TrueR2);

  full_row = {TrueR2 c a};

  % Ends Dwell Time Computation
  finalData(n).SMT_Spots = tracks;

  figure;
  plotData = finalData(n).SMT_Spots;
  [~, ~, Counts] = unique(plotData(:, 1));
  histogram(Counts, max(Counts));
  title(gfpFile, 'interpreter', 'latex');
  xlabel('Frame', 'interpreter', 'latex');
  ylabel('Number of Sox2 Detected', 'interpreter', 'latex');
  pause(1);
  savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_sox2_history.fig']);

  % adding last known gfp location to cellData
  finalFrame = size(cellData, 1);
  idx = find(cellData(:, 1) == 1); % TODO: last gfp seen location
  allFrames = 1 : finalFrame; allFrames(idx) = [];
  for m = 1: length(allFrames)
    try
      lastLoc = cellData(allFrames(m) - 1, 4:5);
    catch
      lastLoc = [0, 0];
    end
    cellData(allFrames(m), 4:5) = lastLoc;
  end

  extraData = zeros(1, 20);
  idx = find(cellData(:, 3) == 1); % frames that are on state
  onFrame = [];
  for m = 1 : length(idx) - 1
    if idx(m) + frameIntervals == idx(m + 1)
      onFrame = [onFrame; [idx(m) : idx(m + 1)]'];
    else
      onFrame = [onFrame; idx(m)];
    end
  end
  onFrame = [onFrame; idx(m + 1)];
  onFrame = unique(onFrame);
  idx = onFrame;
  coSox = 0; orSox = 0;
  for m = 1 : length(idx)
    frame = find(tracks(:,1) == idx(m));
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      coSox = coSox + sum(distSox < distTol);
      if sum(distSox < distTol) > 0
        orSox = orSox + 1;
      end
    end
  end
  disp(['Total number of Sox2 detected: ' num2str(length(tracks))])
  disp(['Number of Sox2 on GFP spot during ON state: ' num2str(coSox)])
  disp(['Precentage of time GFP spot is covered by Sox2 during ON state: ' num2str((orSox/length(idx)) * 100), '%'])
  disp(['Precentage of Sox2 on GFP spot during ON state: ' num2str((orSox/length(tracks)) * 100), '%'])
  extraData(1:4) = [length(tracks), coSox, (orSox/length(idx)) * 100, (orSox/length(tracks)) * 100];
  extraData(19) = orSox;

  % Plotting timeline plots
  soxInFrame = [];
  figure;
%   plot(1:length(cellData), cellData(:, 1), 'b'); hold on; % GFP
  plot(1:length(cellData), cellData(:, 1), 'Color', '#76D7C4', 'linewidth', 1.25); hold on; % GFP
%   plot(1:length(cellData), cellData(:, 2), 'r'); % RFP
%   plot(1:length(cellData), cellData(:, 2), 'Color', '	#F08080', 'linewidth', 1.25); % RFP
  sox2Appear = zeros(length(cellData), 1);
  sox2AppearLong = zeros(length(cellData), 1);
  sox2Number = zeros(length(cellData), 1);
  minDist = zeros(length(cellData), 1);
  for m = 1 : length(cellData)
    frame = find(tracks(:,1) == m);
    tracksInFrame = tracks(frame, 2:3);
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(m, 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      minDist(m) = min(distSox);
      if sum(distSox < distTol) > 0
        sox2Appear(m) = 1;
        sox2Number(m) = sum(distSox < distTol);
        if size(soxInFrame, 1) > 0
          for p = 1 : size(soxInFrame, 1)
            distSoxH = squareform(pdist([soxInFrame(p, :); tracksInFrame(distSox < distTol, :)]));
            distSoxH = distSoxH(1, 2:end);
            if sum(distSoxH < 1) > 0
              % a sox last longer than 2 frames!
              sox2AppearLong(m) = sox2AppearLong(m) + 1;
            end
          end
        end
        soxInFrame = tracksInFrame(distSox < distTol, :);
      else
        soxInFrame = [];
      end
    end
  end
  minDist = minDist .* impars.PixelSize;
%   plot(1:length(cellData), sox2Appear, 'k'); % Sox 2 on last seen GFP spot
  plot(1:length(cellData), sox2Appear, 'Color', '#BB8FCE', 'linewidth', 1.25); % Sox 2 on last seen GFP spot
  coFrame = find(cellData(:, 3) > 0);
  if ~isempty(coFrame)
    plot(coFrame, repelem(1.1, length(coFrame)), '+', 'linewidth', 1.25);
  end
  title(gfpFile, 'interpreter', 'latex');
  xlabel('Frame', 'interpreter', 'latex');
  ylabel('Presence of Proteins', 'interpreter', 'latex');
%   legend('GFP', 'RFP', 'Sox 2', 'GFP-RFP colocalise', 'interpreter', 'latex') % legend edit
  legend('GFP', 'Sox 2', 'GFP-RFP colocalise', 'interpreter', 'latex') 
  ylim([0, 1.2]);
  pause(1);
  savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_time_lapse.fig']);
  % Sox plot
  figure;
  bar(1:length(sox2Number), [sox2Number'; sox2AppearLong'], 'BarWidth', 1);
  title(gfpFile, 'interpreter', 'latex');
  xlabel('Frame', 'interpreter', 'latex');
  ylabel('Number of Proteins', 'interpreter', 'latex');
  legend('Sox 2', 'Specific Binding', 'interpreter', 'latex');
  pause(1);
  savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_sox_number_frames.fig']);

  figure;
  plot(cumsum(sox2Appear));
  title(gfpFile, 'interpreter', 'latex');
  xlabel('Frame', 'interpreter', 'latex');
  ylabel('Cumulative Sox2 Occurrence', 'interpreter', 'latex');
  pause(1);
  savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_sox_motif.fig']);
  
%   idx = find(cellData(:, 3) == 1);
  idx = onFrame;
  idx = idx(ismember(idx, tracks(:, 1)));
  extraData(15) = mean(minDist(idx));
  %

  onFrame = idx;
%   idx = find(cellData(:, 3) == 0); % off state, need to compute last seen gfp location onto the data
  allFrame = 1:length(cellData);
  offFrame = allFrame(~ismember(allFrame, onFrame));
  idx = offFrame;
  coSox = 0; orSox = 0; detectedFrame = [];
  for m = 1 : length(idx)
    frame = find(tracks(:,1) == idx(m));
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      coSox = coSox + sum(distSox < distTol);
      if sum(distSox < distTol) > 0
        orSox = orSox + 1;
        detectedFrame = [detectedFrame; frame(distSox < distTol)];
      end
    end
  end
  % onFrame, detectedFrame
  disp(['Number of Sox2 on last detected GFP spot during OFF state: ' num2str(coSox)])
  extraData(5) = coSox;
  extraData(14) = (orSox / length(idx)) * 100;
  extraData(20) = orSox;

%   idx = find(cellData(:, 3) == 0);
  idx = offFrame;
  idx = idx(ismember(idx, tracks(:, 1)));
  extraData(16) = mean(minDist(idx));
  %

  idx = find(cellData(:, 1) == 1 & cellData(:, 2) == 0);
  coSox = 0; orSox = 0;
  for m = 1 : length(idx)
    frame = find(tracks(:,1) == idx(m));
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      coSox = coSox + sum(distSox < distTol);
      if sum(distSox < distTol) > 0
        orSox = orSox + 1;
      end
    end
  end
  disp(['Number of Sox2 on GFP spot when only GFP is detected: ' num2str(coSox)])
  extraData(6) = coSox;

  idx = find(cellData(:, 1) == 0 & cellData(:, 2) == 1);
  coSox = 0; orSox = 0;
  for m = 1 : length(idx)
    frame = find(tracks(:,1) == idx(m));
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      coSox = coSox + sum(distSox < distTol);
      if sum(distSox < distTol) > 0
        orSox = orSox + 1;
      end
    end
  end
  disp(['Number of Sox2 on last detected GFP spot when only RFP is detected: ' num2str(coSox)])
  extraData(7) = coSox;

%   idx = find(cellData(:, 1) == 0 & cellData(:, 2) == 0);
%   coSox = 0; orSox = 0;
%   for m = 1 : length(idx)
%     frame = find(tracks(:,1) == idx(m));
%     if ~isempty(frame)
%       distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
%       distSox = distSox(1, 2:end);
%       coSox = coSox + sum(distSox < distTol);
%       if sum(distSox < distTol) > 0
%         orSox = orSox + 1;
%       end
%     end
%   end
  disp(['Number of Sox2 on last detected GFP spot when neither GFP or RFP is detected: ' num2str(coSox)])
%   extraData(8) = coSox;
  extraData(8) = extraData(5) - extraData(6) - extraData(7);

%   idx = find(cellData(:, 3) == 1);
  idx = onFrame;
  coSox = 0; orSox = 0;
  for m = 1 : length(idx)
    frame = find(tracks(:,1) == idx(m));
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      coSox = coSox + length(distSox);
      if ~isempty(distSox)
        orSox = orSox + 1;
      end
    end
  end
  disp(['Number of Sox2 detected during ON state: ' num2str(coSox)])
  disp(['Number of frames in ON state: ' num2str(length(idx))])
  extraData(9:10) = [coSox, length(idx)];

%   idx = find(cellData(:, 3) == 0);
  idx = offFrame;
  coSox = 0; orSox = 0;
  for m = 1 : length(idx)
    frame = find(tracks(:,1) == idx(m));
    if ~isempty(frame)
      distSox = squareform(pdist([cellData(idx(m), 4:5); tracks(frame, 2:3)]));
      distSox = distSox(1, 2:end);
      coSox = coSox + length(distSox);
      if ~isempty(distSox)
        orSox = orSox + 1;
      end
    end
  end
  disp(['Number of Sox2 detected during OFF state: ' num2str(coSox)])
  disp(['Number of frames in OFF state: ' num2str(length(idx))])
  extraData(11:12) = [coSox, length(idx)];
  extraData(13) = min(minDist);
  
  cellDataSave = cell(1, 5);
  cellDataSave{1} = 'GFP Presence';
  cellDataSave{2} = 'RFP Presence';
  cellDataSave{3} = 'Colocalise';
  cellDataSave{4} = 'GFP x Coordinate';
  cellDataSave{5} = 'GFP y Coordinate';
  cellDataSave = [cellDataSave; arrayfun(@num2str,cellData(:, 1:5),'un',0)];
  writecell(cellDataSave, [allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_Data.xls'])
%   writematrix(cellData, [allCells(n).folder, filesep, gfpFile(1:end-4), '.xls'])

  extraData(17) = (extraData(2) / extraData(10)) * 100;
  extraData(18) = (extraData(5) / extraData(12)) * 100;
  cellSummaryData = cell(12, 1);
  cellSummaryData{1} = 'Total number of Sox2 detected:';
  cellSummaryData{2} = 'Number of Sox2 on GFP spot during ON state:';
  cellSummaryData{3} = 'Precentage of time GFP spot is covered by Sox2 during ON state:';
  cellSummaryData{4} = 'Precentage of Sox2 on GFP spot during ON state:';
  cellSummaryData{5} = 'Number of Sox2 on last detected GFP spot during OFF state:';
  cellSummaryData{6} = 'Number of Sox2 on GFP spot when only GFP is detected:';
  cellSummaryData{7} = 'Number of Sox2 on last detected GFP spot when only RFP is detected:';
  cellSummaryData{8} = 'Number of Sox2 on last detected GFP spot when neither GFP or RFP is detected:';
  cellSummaryData{9} = 'Number of Sox2 detected during ON state:';
  cellSummaryData{10} = 'Number of frames in ON state:';
  cellSummaryData{11} = 'Number of Sox2 detected during OFF state:';
  cellSummaryData{12} = 'Number of frames in OFF state:';
  cellSummaryData{13} = 'Minimum Distance of Sox2 to GFP:';
  cellSummaryData{14} = 'Precentage of time GFP spot is covered by Sox2 during OFF state:';
  cellSummaryData{15} = 'Mean Minimum Distance between Sox2 and GFP During ON state:';
  cellSummaryData{16} = 'Mean Minimum Distance between Sox2 and GFP During OFF state:';
  cellSummaryData{17} = 'Visiting frequency of Sox2 on GFP spot during ON state:';
  cellSummaryData{18} = 'Visiting frequency of Sox2 on GFP spot during OFF state:';
  cellSummaryData{19} = 'Number of Frames Sox2 on GFP spot during ON state:';
  cellSummaryData{20} = 'Number of Frames Sox2 on GFP spot during OFF state:';
  cellSummaryData = [cellSummaryData, arrayfun(@num2str,extraData','un',0)];
  writecell(cellSummaryData, [allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_Summary.xls'])
  finalData(n).Extra = extraData;
%   set(gcf, 'Position', [100,100, 2000, 500])

  figure;
  gfpInt = cellData(1:frameIntervals:length(stackGFP) - 1, 6) .* 3;
  rfpInt = cellData(1:frameIntervals:length(stackGFP) - 1, 7);
  bar(1:frameIntervals:length(stackGFP) - 1, gfpInt, 'g', 'FaceAlpha', 0.3); hold on;
  bar(1:frameIntervals:length(stackGFP) - 1, rfpInt, 'm', 'FaceAlpha', 0.3);
  legend('GFP', 'RFP', 'interpreter', 'latex')
  title('Mean Intensity of Channels', 'interpreter', 'latex');
  xlabel('Frame', 'interpreter', 'latex');
  ylabel('Mean Intensity Detected', 'interpreter', 'latex');
  pause(1);
%   gfpIntL = []; gfpIntx = [];
%   rfpIntL = []; rfpIntx = [];
%   for i = 1 : length(gfpInt)
%     if gfpInt(i) ~= 0
%       gfpIntL(end + 1) = gfpInt(i);
%       gfpIntx(end + 1) = i;
%     end
%     if rfpInt(i) ~= 0
%       rfpIntL(end + 1) = rfpInt(i);
%       rfpIntx(end + 1) = i;
%     end
%   end
%   plot(gfpIntx, gfpIntL, 'g');
%   plot(rfpIntx, rfpIntL, 'm');
%   pause(1);
  savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_signal_intensity.fig']);
  1;

  % plot gfp intensity with sox
  figure;
  gfpInt = cellData(1:frameIntervals:length(stackGFP) - 1, 6) .* 3;
  rfpInt = cellData(1:frameIntervals:length(stackGFP) - 1, 7);
  bar(1:frameIntervals:length(stackGFP) - 1, gfpInt, 'g', 'FaceAlpha', 0.3); hold on;
  bar(1:frameIntervals:length(stackGFP) - 1, rfpInt, 'm', 'FaceAlpha', 0.3);
  title('Mean Intensity of Channels with Sox Occurrence', 'interpreter', 'latex');
  xlabel('Frame', 'interpreter', 'latex');
  ylabel('Mean Intensity Detected', 'interpreter', 'latex');
  plot(1:length(sox2Number), sox2Number' .* max([gfpInt; rfpInt])/2, 'k');
  plot(1:length(sox2Number), sox2AppearLong' .* max([gfpInt; rfpInt]), 'k');
  legend('GFP', 'RFP', 'Sox Occurence', 'interpreter', 'latex')
  pause(0.5);
  savefig([allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_signal_intensity_with_sox.fig']);

  maxGfpInd = find(gfpInt == max(gfpInt));
  soxInd = find(sox2Number > 0);
  writematrix(maxGfpInd, [allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_max_gfp_timing.xls']);
  writematrix(soxInd, [allCells(n).folder, filesep, 'Data', filesep, gfpFile(1:end-4), '_sox_timing.xls']);

  % plotting sox
%   figure;
%   soxPlot = stackFISH(1).data;
%   [counts, intensity] = histcounts(soxPlot(:));
%   intOI = intensity(counts > 10);
%   imshow(imadjust(soxPlot, [intOI(2)/65535, intOI(end)/65535])); hold on;
%   for i = 1 : length(stackGFP)
%     trackNum = tracks(tracks(:, 1) == i, 2);
%     if ~isempty(trackNum)
%       plot(tracks(tracks(:, 1) == i, 2) + smtCorrection(1), tracks(tracks(:, 1) == i, 3) + smtCorrection(2), 'rx');
%     end
%     1;
%   end
%   plot(cellData(1, 4) + smtCorrection(1), cellData(1, 5) + smtCorrection(2), 'gx');
  % plot smt with gfp in final state
%   figure;
%   plot(tracks(:, 2), tracks(:, 3), 'rx'); hold on;
%   plot(cellData(1, 4), cellData(1, 5), 'gx');
  % plotting gfp
%   figure;
%   gfpRaw = stackGFP(1).data;
%   [counts, intensity] = histcounts(gfpRaw(:));
%   intOI = intensity(counts > 10);
%   imshow(imadjust(gfpRaw, [intOI(2)/65535, intOI(end)/65535])); hold on;
%   plot(cellData(1, 4), cellData(1, 5), 'gx');
  %
%   figure; hold on;
%   colour = ['r', "#008000", 'b'];
%   frameNow = 1; colourNow = 1;
%   for i = 1 : length(trackedPar)
%     if trackedPar(i).Frame(1) ~= frameNow
%       if colourNow ~= 3
%         colourNow = colourNow + 1;
%       else
%         colourNow = 1;
%       end
%       frameNow = trackedPar(i).Frame(1);
%     end
% %     plot(trackedPar(i).xy(1, 1), trackedPar(i).xy(1, 2), [colour(colourNow), 'x']);
%     plot(trackedPar(i).xy(1, 1), trackedPar(i).xy(1, 2), 'x', 'MarkerEdgeColor', colour(colourNow));
%   end
%   plot(32.5, 25, 'gpentagram','MarkerSize', 12, 'MarkerFaceColor', 'g')
%   figure;
%   fishRaw = stackFISH(1).data;
%   [counts, intensity] = histcounts(fishRaw(:));
%   intOI = intensity(counts > 10); % background adjust
%   fishRawAd = imadjust(fishRaw, [intOI(2)/65535, intOI(end)/65535]);
  if plotHeatMap
    %   xGrid = 512; yGrid = xGrid; tGrid = zeros(yGrid, xGrid);
    % %   xBoundary = min(tracks(:, 2)) : (max(tracks(:, 2)) - min(tracks(:, 2)))/ xGrid : max(tracks(:, 2));
    % %   yBoundary = min(tracks(:, 3)) : (max(tracks(:, 3)) - min(tracks(:, 3)))/ yGrid : max(tracks(:, 3));
    %   xBoundary = 0 : max(tracks(:, 2))/ xGrid : max(tracks(:, 2));
    %   yBoundary = 0 : max(tracks(:, 3))/ yGrid : max(tracks(:, 3));
    pxPerGrid = 1;
    xBoundary = 0 : pxPerGrid : 512;
    yBoundary = 0 : pxPerGrid : 512;
    xGrid = length(xBoundary) - 1; yGrid = xGrid; tGrid = zeros(yGrid, xGrid);
    for xg = 1 :  xGrid
      for yg = 1 : yGrid
        tGrid(yg, xg) = sum(tracks(:, 2) > xBoundary(xg) & tracks(:, 2) < xBoundary(xg + 1) & tracks(:, 3) > yBoundary(yg) & tracks(:, 3) < yBoundary(yg + 1));
      end
    end
    figure;
    [X, Y] = meshgrid(1:xGrid);
    surface(X, Y, tGrid); view(0, 90);
    %   mesh(X, Y, zeros(10, 10), tGrid); view(0, 90);
    hold on;
    plot3(max(find(cellData(1, 4) > xBoundary)) + 0.5, max(find(cellData(1, 5) > yBoundary)) + 0.5, max(tGrid, [], 'all') + 1, 'rx');
    colorbar;
    set(gca, 'YDir','reverse')
    %   th = 0:pi/50:2*pi;
    %   xunit = distTol * cos(th) + max(find(cellData(1, 4) > xBoundary)) + 0.5;
    %   yunit = distTol * sin(th) + max(find(cellData(1, 5) > yBoundary)) + 0.5;
    %   plot3(xunit, yunit, repelem(max(tGrid, [], 'all') + 1, length(xunit)));
  end
  if plotAllGFPs
    figure;
    gfpRaw = stackGFP(m).data;
    [counts, intensity] = histcounts(gfpRaw(:));
    intOI = intensity(counts > 10); % background adjust
    imshow(imadjust(gfpRaw, [intOI(2)/65535, intOI(end)/65535])); hold on;
    plot(cellData(:,4), cellData(:,5), 'rx');
  end
  close all;
end
writecell(dwellDataCompiled, [mainDir, filesep, 'Summary_dwell_time.xls'])