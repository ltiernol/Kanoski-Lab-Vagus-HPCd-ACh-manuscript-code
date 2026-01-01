%% ================== INPUT: assemble your animals ==================
% For each animal, provide:
%   .name    = animal ID (string or char)
%   .keydown = column vector of timestamps (odd=start, even=end)
%   .fp      = [N x 2], col1 = z-score signal, col2 = time (s)
%
% Example (your A8608):
animals(1).name    = "A8608";
animals(1).keydown = A8608_refeeding_keydown(:);
animals(1).fp      = FP.Refeeding.zscore_finalized.A8608_refeeding;

animals(2).name    = "A8609";
animals(2).keydown = A8609_refeeding_keydown(:);
animals(2).fp      = FP.Refeeding.zscore_finalized.A8609_refeeding;

animals(3).name    = "A8611";
animals(3).keydown = A8611_refeeding_keydown(:);
animals(3).fp      = FP.Refeeding.zscore_finalized.A8611_refeeding;

animals(4).name    = "A8612";
animals(4).keydown = A8612_refeeding_keydown(:);
animals(4).fp      = FP.Refeeding.zscore_finalized.A8612_refeeding;

animals(5).name    = "A9601";
animals(5).keydown = A9601_refeeding_keydown(:);
animals(5).fp      = FP.Refeeding.zscore_finalized.A9601_refeeding;

animals(6).name    = "A9603";
animals(6).keydown = A9603_refeeding_keydown(:);
animals(6).fp      = FP.Refeeding.zscore_finalized.A9603_refeeding;

animals(7).name    = "A9605";
animals(7).keydown = A9605_refeeding_keydown(:);
animals(7).fp      = FP.Refeeding.zscore_finalized.A9605_refeeding;

animals(8).name    = "J2001";
animals(8).keydown = J2001_refeeding_keydown(:);
animals(8).fp      = FP.Refeeding.zscore_finalized.J2001_refeeding;

animals(9).name    = "J2019";
animals(9).keydown = J2019_refeeding_keydown(:);
animals(9).fp      = FP.Refeeding.zscore_finalized.J2019_refeeding;

animals(10).name    = "J2021";
animals(10).keydown = J2021_refeeding_keydown(:);
animals(10).fp      = FP.Refeeding.zscore_finalized.J2021_refeeding;
%% ================== SETTINGS ==================
N = 300;                       % samples before/after each timestamp
winIdx = (-N:N);               % sample offsets relative to the event
T = numel(winIdx);             % 601 samples per bout

% Map samples to seconds for axis labeling:
% 300 samples correspond to 15 seconds  -> effective fs = 20 Hz
fs = N / 15;                   % = 20 when N=300
xAxisSec = winIdx / fs;        % seconds axis: -15 ... +15
% ================== ACCUMULATORS ==================
matStart_all = [];        % rows = bouts across all animals
matEnd_all   = [];
rowAnimalIdx = [];        % animal index per row
rowBoutIdx   = [];        % bout index within animal per row
animalRowBlocks = zeros(numel(animals),1);  % # bouts contributed per animal

% ================== PROCESS EACH ANIMAL ==================
for a = 1:numel(animals)
    % Pull z and t
    fpMat = animals(a).fp;
    zSig  = fpMat(:,1);  zSig = zSig(:);
    tSig  = fpMat(:,2);  tSig = tSig(:);

    % Sanitize time (unique, increasing) and sync zSig accordingly
    [tSigU, zSigU] = sanitizeTime(tSig, zSig);
    if any(diff(tSigU) <= 0)
        error('Animal %s: time vector not strictly increasing after sanitize.', string(animals(a).name));
    end

    % Bout starts/ends
    kd = animals(a).keydown(:);
    starts_raw = kd(1:2:end);
    ends_raw   = kd(2:2:end);
    nPairs     = min(numel(starts_raw), numel(ends_raw));
    starts_raw = starts_raw(1:nPairs);
    ends_raw   = ends_raw(1:nPairs);

    if nPairs == 0
        warning('Animal %s: no complete start/end pairs. Skipping.', string(animals(a).name));
        continue;
    end

    % Build matrices with fixed-sample extraction
    matS = nan(nPairs, T);
    matE = nan(nPairs, T);

    badS = 0; badE = 0;
    for i = 1:nPairs
        idxS = nearestIdx(tSigU, starts_raw(i));
        idxE = nearestIdx(tSigU, ends_raw(i));

        matS(i,:) = getRowPadded(zSigU, idxS, N, T);
        matE(i,:) = getRowPadded(zSigU, idxE, N, T);

        badS = badS + all(isnan(matS(i,:)));
        badE = badE + all(isnan(matE(i,:)));
    end

    if badS > 0 || badE > 0
        warning('Animal %s: %d/%d START rows and %d/%d END rows are all NaN (timestamps likely outside recording).', ...
            string(animals(a).name), badS, nPairs, badE, nPairs);
    end

    % Optional: sort within-animal by chronological start
    if sortWithinAnimal
        [~, sIdx] = sort(starts_raw, 'ascend');
        matS = matS(sIdx,:);
        matE = matE(sIdx,:);
    end

    % Accumulate
    matStart_all   = [matStart_all; matS]; %#ok<AGROW>
    matEnd_all     = [matEnd_all;   matE]; %#ok<AGROW>
    rowAnimalIdx   = [rowAnimalIdx; a * ones(nPairs,1)]; %#ok<AGROW>
    rowBoutIdx     = [rowBoutIdx;   (1:nPairs)'];       %#ok<AGROW>
    animalRowBlocks(a) = nPairs;
end

% Bail if nothing accumulated
if isempty(matStart_all)
    error('No bouts found across animals. Check inputs/timestamps.');
end

% ================== PLOTTING (seconds axis, clean labels) ==================
cmin = nanmin([matStart_all(:); matEnd_all(:)]);
cmax = nanmax([matStart_all(:); matEnd_all(:)]);
cumRows   = cumsum(animalRowBlocks);
totalRows = size(matStart_all,1);

figure('Color','w','Units','normalized','Position',[0.08 0.12 0.84 0.56]);
tl = tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% --- START aligned (combined) ---
ax1 = nexttile(tl);
imagesc(xAxisSec, [], matStart_all);
%set(ax1,'YDir','normal');                       % row 1 at top, increases downward
yticks(ax1, 1:totalRows);                       % bout numbers ascending top->bottom
yticklabels(ax1, 1:totalRows);
xlabel(ax1, 'Time from bout start (seconds)');
ylabel(ax1, 'Bout # (all animals)');
title(ax1, sprintf('HPCd ACh (z) | START-aligned | %d bouts | \x00B1%ds', totalRows, 15));
c = colorbar(ax1); c.Label.String = 'z-score';
caxis(ax1, [cmin cmax]);
xlim(ax1, [-15 15]); xticks(ax1, -15:5:15);

hold(ax1,'on');
plot(ax1, [0 0], [0 totalRows+1], 'k-', 'LineWidth', 1);  % zero line
for k = 1:numel(cumRows)-1
    yline(ax1, cumRows(k) + 0.5, 'k-', 'LineWidth', 0.8); % animal separators
end
hold(ax1,'off');

% --- END aligned (combined) ---
ax2 = nexttile(tl);
imagesc(xAxisSec, [], matEnd_all);
%set(ax2,'YDir','normal');
yticks(ax2, 1:totalRows);
yticklabels(ax2, 1:totalRows);
xlabel(ax2, 'Time from bout end (seconds)');
ylabel(ax2, 'Bout # (all animals)');
title(ax2, sprintf('HPCd ACh (z) | END-aligned | %d bouts | \x00B1%ds', totalRows, 15));
c = colorbar(ax2); c.Label.String = 'z-score';
caxis(ax2, [cmin cmax]);
xlim(ax2, [-15 15]); xticks(ax2, -15:5:15);

hold(ax2,'on');
plot(ax2, [0 0], [0 totalRows+1], 'k-', 'LineWidth', 1);
for k = 1:numel(cumRows)-1
    yline(ax2, cumRows(k) + 0.5, 'k-', 'LineWidth', 0.8);
end
hold(ax2,'off');

% --- Blue -> Yellow colormap (pure gradient) ---
nColors = 256;
blueYellow = [linspace(0,1,nColors)' linspace(0,1,nColors)' linspace(1,0,nColors)'];
colormap(blueYellow);

% === Animal block labels as "Animal 1", "Animal 2", ... ===
midRows = cumsum([1; animalRowBlocks(1:end-1)]) + floor(animalRowBlocks/2);
for a = 1:numel(animalRowBlocks)
    if animalRowBlocks(a) > 0
        txt = sprintf('  Animal %d', a);
        text(ax1, xAxisSec(1), midRows(a), txt, 'VerticalAlignment','middle', ...
             'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',1,'EdgeColor','none');
        text(ax2, xAxisSec(1), midRows(a), txt, 'VerticalAlignment','middle', ...
             'FontSize',9,'FontWeight','bold','BackgroundColor','w','Margin',1,'EdgeColor','none');
    end
end






%% Extracting time of maxima and minima relative to bout start and end. 
% ================== PARAMETERS ==================
% Window & sampling assumptions (20 Hz; 30 s window)
fs               = 20;         % samples/second
halfWin_s        = 15;         % seconds before and after the event
halfN            = halfWin_s * fs;         % 300 samples
winN             = 2 * halfN;              % 600 samples
dt_expected_s    = 1/fs;                   % 0.05 s

% IMPORTANT: Units for keydown timestamps (your earlier message said minutes).
% If your keydowns are ALREADY in seconds, set this to false.
keydown_in_minutes = true;

% Storage (cell/struct arrays that grow per animal)
ResultsByAnimal = struct();
StartSummaryAll = table();
EndSummaryAll   = table();

% Helper: relative time vector for each extracted window (seconds)
relTime_s = ((-halfN):(halfN-1)) * dt_expected_s;

% ================== LOOP OVER ANIMALS ==================
for ai = 1:numel(animals)
    id = string(animals(ai).name);

    % ------------- Load FP data and keydowns -------------
    fpMat = animals(ai).fp;               % [T x 2], col1=z, col2=time
    z     = fpMat(:,1);
    t_raw = fpMat(:,2);                   % FP time (unknown units: minutes or seconds)
    kd    = animals(ai).keydown(:);       % keydown times (odd=start, even=end) in minutes (assumed)

    % ------------- Normalize time units -------------
    % Detect FP timebase units by sampling interval:
    if numel(t_raw) < 3
        error('Animal %s: FP time vector too short.', id);
    end
    dt_est = median(diff(t_raw), 'omitnan');

    % If dt is ~0.05 (seconds), FP is in seconds; if ~0.000833... (~1/1200), FP is in minutes.
    if abs(dt_est - 0.05) < 0.005
        t_fp_s = t_raw;                         % FP already in seconds
    elseif abs(dt_est - (0.05/60)) < 5e-4
        t_fp_s = t_raw * 60;                    % FP was in minutes -> convert to seconds
    else
        % Fallback: compare span to guess
        if max(t_raw) < 1000
            % Likely minutes (e.g., 50 min recording)
            t_fp_s = t_raw * 60;
        else
            t_fp_s = t_raw;
        end
        warning('Animal %s: Unclear FP time step (dt = %.6f). Guessed units.', id, dt_est);
    end

    % Keydown units (minutes -> seconds, unless you flip keydown_in_minutes)
    if keydown_in_minutes
        kd_s = kd * 60;
    else
        kd_s = kd;
    end

    % ------------- Label starts/ends -------------
    isStart = false(size(kd_s)); isStart(1:2:end) = true;
    isEnd   = ~isStart;

    startTimes_s = kd_s(isStart);
    endTimes_s   = kd_s(isEnd);

    % ------------- Preallocate per-animal outputs -------------
    start_keep = true(size(startTimes_s));
    end_keep   = true(size(endTimes_s));

    start_windows = nan(numel(startTimes_s), winN);
    end_windows   = nan(numel(endTimes_s),   winN);

    start_min_latency_s  = nan(numel(startTimes_s),1);
    end_max_latency_s    = nan(numel(endTimes_s),1);
    start_min_idx_in_win = nan(numel(startTimes_s),1);
    end_max_idx_in_win   = nan(numel(endTimes_s),1);

    % ------------- Fast mapping from event time -> nearest FP index -------------
    % Use interp1 over the time vector to get the nearest sample index.
    idx_vec = (1:numel(t_fp_s))';
    getCenterIdx = @(t_event_s) interp1(t_fp_s, idx_vec, t_event_s, 'nearest', 'extrap');

    % ------------- START windows: find local minima & latency -------------
    for k = 1:numel(startTimes_s)
        tevent = startTimes_s(k);
        cIx    = getCenterIdx(tevent);

        i0 = cIx - halfN;          % inclusive
        i1 = cIx + halfN - 1;      % inclusive (600 points total: 300 before, 300 after)

        % Guard against edges
        if i0 < 1 || i1 > numel(z)
            start_keep(k) = false;
            continue
        end

        segZ = z(i0:i1);
        start_windows(k,:) = segZ;

        [~, localMinIdx] = min(segZ);
        start_min_idx_in_win(k) = localMinIdx;

        % Absolute time (seconds) of minima sample
        t_min_s = t_fp_s(i0 + localMinIdx - 1);

        % Latency (seconds): positive if min occurs AFTER bout start
        start_min_latency_s(k) = t_min_s - tevent;
    end

    % ------------- END windows: find local maxima & latency -------------
    for k = 1:numel(endTimes_s)
        tevent = endTimes_s(k);
        cIx    = getCenterIdx(tevent);

        i0 = cIx - halfN;
        i1 = cIx + halfN - 1;

        % Guard against edges
        if i0 < 1 || i1 > numel(z)
            end_keep(k) = false;
            continue
        end

        segZ = z(i0:i1);
        end_windows(k,:) = segZ;

        [~, localMaxIdx] = max(segZ);
        end_max_idx_in_win(k) = localMaxIdx;

        % Absolute time (seconds) of maxima sample
        t_max_s = t_fp_s(i0 + localMaxIdx - 1);

        % Latency (seconds): positive if max occurs AFTER bout end
        end_max_latency_s(k) = t_max_s - tevent;
    end

    % ------------- Keep valid rows (full windows only) -------------
    start_windows          = start_windows(start_keep,:);
    end_windows            = end_windows(end_keep,:);
    start_min_latency_s    = start_min_latency_s(start_keep);
    end_max_latency_s      = end_max_latency_s(end_keep);
    start_min_idx_in_win   = start_min_idx_in_win(start_keep);
    end_max_idx_in_win     = end_max_idx_in_win(end_keep);

    startTimes_s           = startTimes_s(start_keep);
    endTimes_s             = endTimes_s(end_keep);

    % ------------- Package per-animal results -------------
    R = struct();
    R.animal_id          = id;
    R.fs_Hz              = fs;
    R.window_samples     = winN;
    R.half_window_s      = halfWin_s;

    R.start = struct( ...
        'timestamps_s',        startTimes_s, ...
        'win_activity',        start_windows, ...   % [nStarts x 600]
        'win_time_rel_s',      relTime_s, ...
        'min_latency_s',       start_min_latency_s, ...
        'min_idx_in_window',   start_min_idx_in_win ...
    );

    R.end = struct( ...
        'timestamps_s',        endTimes_s, ...
        'win_activity',        end_windows, ...     % [nEnds x 600]
        'win_time_rel_s',      relTime_s, ...
        'max_latency_s',       end_max_latency_s, ...
        'max_idx_in_window',   end_max_idx_in_win ...
    );

    ResultsByAnimal.(id) = R;

end

%% ================== PLOT HISTOGRAMS OF LATENCIES ==================
all_min_latencies = [];
all_max_latencies = [];

ids = fieldnames(ResultsByAnimal);

for i = 1:numel(ids)
    R = ResultsByAnimal.(ids{i});
    all_min_latencies = [all_min_latencies; R.start.min_latency_s(:)];
    all_max_latencies = [all_max_latencies; R.end.max_latency_s(:)];
end

% Remove NaNs just in case
all_min_latencies = all_min_latencies(~isnan(all_min_latencies));
all_max_latencies = all_max_latencies(~isnan(all_max_latencies));

% -------- Plot histograms --------
figure('Name','Latency Histograms','Color','w');

subplot(1,2,1);
histogram(all_min_latencies, 'BinWidth', 1, 'FaceColor',[0.2 0.6 0.8]); % 1s bins
xlabel('Latency to Minimum (s)');
ylabel('Count');
title('Bout Start: Minima Latencies');
xline(0,'k--','LineWidth',1); % reference line at bout start

subplot(1,2,2);
histogram(all_max_latencies, 'BinWidth', 1, 'FaceColor',[0.8 0.4 0.2]);
xlabel('Latency to Maximum (s)');
ylabel('Count');
title('Bout End: Maxima Latencies');
xline(0,'k--','LineWidth',1); % reference line at bout end





%% ================== PRISM EXPORT VARIABLES ==================
% Goal: create variables you can copy/paste into GraphPad Prism

% -------- 1) Raw values (best for Prism histograms / scatter) --------
% Two-column table: each column is a dataset
Prism_Latencies_Raw = table(all_min_latencies, all_max_latencies, ...
    'VariableNames', {'Start_MinLatency_s','End_MaxLatency_s'});

% If Prism prefers separate vectors:
Prism_StartMinLatency_s = all_min_latencies;
Prism_EndMaxLatency_s   = all_max_latencies;

% -------- 2) Long-format per-event table (useful for grouping by animal) --------
% Each row = one event latency with animal ID and type
prismRows = {};
ids = fieldnames(ResultsByAnimal);

for i = 1:numel(ids)
    id = string(ids{i});
    R  = ResultsByAnimal.(ids{i});

    v1 = R.start.min_latency_s(:);
    v1 = v1(~isnan(v1));
    prismRows = [prismRows; [repmat({id}, numel(v1),1), repmat({"Start_Min"}, numel(v1),1), num2cell(v1)]]; %#ok<AGROW>

    v2 = R.end.max_latency_s(:);
    v2 = v2(~isnan(v2));
    prismRows = [prismRows; [repmat({id}, numel(v2),1), repmat({"End_Max"},   numel(v2),1), num2cell(v2)]]; %#ok<AGROW>
end

Prism_Latencies_Long = cell2table(prismRows, ...
    'VariableNames', {'AnimalID','Measure','Latency_s'});

% -------- 3) OPTIONAL: binned histogram tables (to match MATLAB exactly) --------
% Your histogram uses BinWidth = 1 second. Let's build matching bins.
binWidth = 1;
edges = (floor(min([all_min_latencies; all_max_latencies])) : binWidth : ceil(max([all_min_latencies; all_max_latencies])))';
centers = edges(1:end-1) + binWidth/2;

counts_min = histcounts(all_min_latencies, edges);
counts_max = histcounts(all_max_latencies, edges);

Prism_Hist_StartMin = table(centers, counts_min(:), 'VariableNames', {'BinCenter_s','Count'});
Prism_Hist_EndMax   = table(centers, counts_max(:), 'VariableNames', {'BinCenter_s','Count'});

% -------- 4) OPTIONAL: write CSVs for Prism import --------
% (Comment these out if you only want copy/paste from Variable Editor)
%writetable(Prism_Latencies_Raw,  'Prism_Latencies_Raw.csv');
%writetable(Prism_Latencies_Long, 'Prism_Latencies_Long.csv');
%writetable(Prism_Hist_StartMin,  'Prism_Hist_StartMin.csv');
%writetable(Prism_Hist_EndMax,    'Prism_Hist_EndMax.csv');



%% ================== LINEAR REGRESSION: EXTREMA LATENCIES VS BOUT DURATION ==================
all_min_latencies = [];
all_min_durations = [];
all_max_latencies = [];
all_max_durations = [];

ids = fieldnames(ResultsByAnimal);
names_cell = {animals.name};
names_str  = string(names_cell);    % robust to char/string
tolerance_s = 1e-3;                 % tolerance to match event timestamps to keydowns

for i = 1:numel(ids)
    id = string(ids{i});
    R  = ResultsByAnimal.(ids{i});

    % --- find the matching animal by name (single index) ---
    idx = find(names_str == id, 1, 'first');
    if isempty(idx)
        warning('Regression: animal "%s" not found in `animals` struct. Skipping.', id);
        continue
    end

    % --- keydowns and bout durations for this animal ---
    kd = animals(idx).keydown(:);
    if exist('keydown_in_minutes','var') && keydown_in_minutes
        kd_s = kd * 60;
    else
        kd_s = kd;
    end

    % make sure we have complete start/end pairs
    nPairs = floor(numel(kd_s)/2);
    if nPairs < 1
        warning('Regression: animal "%s" has no complete start/end pairs. Skipping.', id);
        continue
    end
    bout_starts = kd_s(1:2:(2*nPairs-1));
    bout_ends   = kd_s(2:2:(2*nPairs));
    bout_durations = bout_ends - bout_starts;

    % --- match kept starts to their durations ---
    for k = 1:numel(R.start.timestamps_s)
        t = R.start.timestamps_s(k);
        % find nearest start within tolerance
        [mindiff, ix] = min(abs(bout_starts - t));
        if ~isempty(ix) && mindiff <= tolerance_s
            all_min_latencies(end+1,1) = R.start.min_latency_s(k);
            all_min_durations(end+1,1) = bout_durations(ix);
        end
    end

    % --- match kept ends to their durations ---
    for k = 1:numel(R.end.timestamps_s)
        t = R.end.timestamps_s(k);
        [mindiff, ix] = min(abs(bout_ends - t));
        if ~isempty(ix) && mindiff <= tolerance_s
            all_max_latencies(end+1,1) = R.end.max_latency_s(k);
            all_max_durations(end+1,1) = bout_durations(ix);
        end
    end
end

% --- clean NaNs ---
valid_min = ~isnan(all_min_latencies) & ~isnan(all_min_durations);
valid_max = ~isnan(all_max_latencies) & ~isnan(all_max_durations);

% --- fit models only if we have >=2 points ---
figure('Name','Regression Analyses','Color','w');

subplot(1,2,1);
if nnz(valid_min) >= 2
    mdl_min = fitlm(all_min_latencies(valid_min), all_min_durations(valid_min));
    scatter(all_min_latencies(valid_min), all_min_durations(valid_min), 36, 'filled'); hold on;
    xgrid = linspace(min(all_min_latencies(valid_min)), max(all_min_latencies(valid_min)), 100)';
    yhat  = predict(mdl_min, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);
    xlabel('Latency to Minimum (s)'); ylabel('Bout Duration (s)');
    title(sprintf('Start: min latency vs duration (R^2=%.3f, p=%.3g)', mdl_min.Rsquared.Ordinary, mdl_min.Coefficients.pValue(2)));
else
    text(0.5,0.5,'Not enough points','HorizontalAlignment','center'); axis off
end

subplot(1,2,2);
if nnz(valid_max) >= 2
    mdl_max = fitlm(all_max_latencies(valid_max), all_max_durations(valid_max));
    scatter(all_max_latencies(valid_max), all_max_durations(valid_max), 36, 'filled'); hold on;
    xgrid = linspace(min(all_max_latencies(valid_max)), max(all_max_latencies(valid_max)), 100)';
    yhat  = predict(mdl_max, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);
    xlabel('Latency to Maximum (s)'); ylabel('Bout Duration (s)');
    title(sprintf('End: max latency vs duration (R^2=%.3f, p=%.3g)', mdl_max.Rsquared.Ordinary, mdl_max.Coefficients.pValue(2)));
else
    text(0.5,0.5,'Not enough points','HorizontalAlignment','center'); axis off
end

%%

%% ================== REGRESSION: Latency vs Bout Order Number ==================
% Requires: `animals`, `ResultsByAnimal`, and the flag `keydown_in_minutes`.

all_start_lat = [];   % latency to min (s)
all_start_ord = [];   % bout index (1,2,...)
all_end_lat   = [];   % latency to max (s)
all_end_ord   = [];   % bout index (1,2,...)

ids = fieldnames(ResultsByAnimal);
names_cell = {animals.name};
names_str  = string(names_cell);
tolerance_s = 1e-3;   % tolerance to match event timestamps to keydowns

for i = 1:numel(ids)
    id = string(ids{i});
    R  = ResultsByAnimal.(ids{i});

    % Find this animal in `animals`
    idx = find(names_str == id, 1, 'first');
    if isempty(idx)
        warning('Latency~Order: animal "%s" not found in `animals`. Skipping.', id);
        continue
    end

    % Build bout start/end sequences (in seconds)
    kd = animals(idx).keydown(:);
    if exist('keydown_in_minutes','var') && keydown_in_minutes
        kd_s = kd * 60;
    else
        kd_s = kd;
    end

    nPairs = floor(numel(kd_s)/2);
    if nPairs < 1, continue; end

    bout_starts = kd_s(1:2:(2*nPairs-1));     % (nPairs x 1)
    bout_ends   = kd_s(2:2:(2*nPairs));       % (nPairs x 1)

    % ----- START side: map kept start times -> bout order -----
    kept_starts = R.start.timestamps_s(:);
    for k = 1:numel(kept_starts)
        t0 = kept_starts(k);
        [mindiff, bix] = min(abs(bout_starts - t0));
        if ~isempty(bix) && mindiff <= tolerance_s && ~isnan(R.start.min_latency_s(k))
            all_start_lat(end+1,1) = R.start.min_latency_s(k);
            all_start_ord(end+1,1) = bix;  % bout order number
        end
    end

    % ----- END side: map kept end times -> bout order -----
    kept_ends = R.end.timestamps_s(:);
    for k = 1:numel(kept_ends)
        t1 = kept_ends(k);
        [mindiff, bix] = min(abs(bout_ends - t1));
        if ~isempty(bix) && mindiff <= tolerance_s && ~isnan(R.end.max_latency_s(k))
            all_end_lat(end+1,1) = R.end.max_latency_s(k);
            all_end_ord(end+1,1) = bix;    % bout order number
        end
    end
end

% Clean
valid_start = ~isnan(all_start_lat) & ~isnan(all_start_ord);
valid_end   = ~isnan(all_end_lat)   & ~isnan(all_end_ord);

% Plot + regress
figure('Name','Latency vs Bout Order','Color','w');

subplot(1,2,1);
if nnz(valid_start) >= 2
    x = all_start_ord(valid_start);
    y = all_start_lat(valid_start);
    mdl_s = fitlm(x, y);
    scatter(x, y, 36, 'filled'); hold on;
    xgrid = linspace(min(x), max(x), 100)';
    yhat  = predict(mdl_s, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);
    xlabel('Bout Order #'); ylabel('Latency to Minimum (s)');
    title(sprintf('Start-side: latency vs order (R^2=%.3f, p=%.3g)', ...
        mdl_s.Rsquared.Ordinary, mdl_s.Coefficients.pValue(2)));
else
    axis off; text(0.5,0.5,'Not enough start data','HorizontalAlignment','center');
end

subplot(1,2,2);
if nnz(valid_end) >= 2
    x = all_end_ord(valid_end);
    y = all_end_lat(valid_end);
    mdl_e = fitlm(x, y);
    scatter(x, y, 36, 'filled'); hold on;
    xgrid = linspace(min(x), max(x), 100)';
    yhat  = predict(mdl_e, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);
    xlabel('Bout Order #'); ylabel('Latency to Maximum (s)');
    title(sprintf('End-side: latency vs order (R^2=%.3f, p=%.3g)', ...
        mdl_e.Rsquared.Ordinary, mdl_e.Coefficients.pValue(2)));
else
    axis off; text(0.5,0.5,'Not enough end data','HorizontalAlignment','center');
end

% Print details to console (optional)
if exist('mdl_s','var'); fprintf('\nStart-side regression:\n'); disp(mdl_s); end
if exist('mdl_e','var'); fprintf('\nEnd-side regression:\n');   disp(mdl_e); end












%% ================== LINEAR REGRESSION: Max latency (from end) vs Next Inter-Bout Interval ==================
% For each kept bout end (with a detected max), pair its latency with the IBI to the *next* bout start:
%   IBI_i = start_{i+1} - end_i

all_max_latencies_IBI = [];   % predictor: latency of max relative to bout end (s)
all_next_IBI_durations = [];  % response: next inter-bout interval duration (s)

ids = fieldnames(ResultsByAnimal);
names_cell = {animals.name};
names_str  = string(names_cell);
tolerance_s = 1e-3;   % matching tolerance for timestamps

for i = 1:numel(ids)
    id = string(ids{i});
    R  = ResultsByAnimal.(ids{i});

    % --- find matching animal in `animals` struct ---
    idx = find(names_str == id, 1, 'first');
    if isempty(idx)
        warning('IBI regression: animal "%s" not found in `animals`. Skipping.', id);
        continue
    end

    % --- keydowns for this animal, to build start/end sequences ---
    kd = animals(idx).keydown(:);
    if exist('keydown_in_minutes','var') && keydown_in_minutes
        kd_s = kd * 60;
    else
        kd_s = kd;
    end

    % ensure complete pairs
    nPairs = floor(numel(kd_s)/2);
    if nPairs < 2
        % need at least two starts to define an IBI after the first end
        warning('IBI regression: animal "%s" has <2 bouts; skipping.', id);
        continue
    end

    bout_starts = kd_s(1:2:(2*nPairs-1));   % length nPairs
    bout_ends   = kd_s(2:2:(2*nPairs));     % length nPairs

    % --- match each kept END time to its index, then grab the *next* start for IBI ---
    kept_ends   = R.end.timestamps_s;       % ends that had full ±15s windows
    lat_max_s   = R.end.max_latency_s;      % latency of max relative to end (s)

    for k = 1:numel(kept_ends)
        t_end = kept_ends(k);

        % find the corresponding end index in the original sequence
        [mindiff, ixEnd] = min(abs(bout_ends - t_end));
        if isempty(ixEnd) || mindiff > tolerance_s
            % couldn't confidently match this end to the original sequence
            continue
        end

        % need a *next* start for IBI: start_{ixEnd+1}
        if ixEnd < numel(bout_starts)
            next_start = bout_starts(ixEnd + 1);
            ibi_next   = next_start - bout_ends(ixEnd);   % next IBI (s)

            % keep only sensible positive IBIs
            if ~isnan(ibi_next) && ibi_next > 0
                all_max_latencies_IBI(end+1,1)  = lat_max_s(k);
                all_next_IBI_durations(end+1,1) = ibi_next;
            end
        end
    end
end

% --- clean and fit model ---
valid = ~isnan(all_max_latencies_IBI) & ~isnan(all_next_IBI_durations);
figure('Name','Regression: Max latency (from end) vs Next IBI','Color','w');

if nnz(valid) >= 2
    mdl_ibi = fitlm(all_max_latencies_IBI(valid), all_next_IBI_durations(valid));
    scatter(all_max_latencies_IBI(valid), all_next_IBI_durations(valid), 36, 'filled'); hold on;

    xgrid = linspace(min(all_max_latencies_IBI(valid)), max(all_max_latencies_IBI(valid)), 100)';
    yhat  = predict(mdl_ibi, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);

    xlabel('Latency to Maximum from Bout End (s)');
    ylabel('Next Inter-Bout Interval (s)');
    title(sprintf('End max latency vs Next IBI (R^2=%.3f, p=%.3g)', ...
        mdl_ibi.Rsquared.Ordinary, mdl_ibi.Coefficients.pValue(2)));

    % Print model to console for details
    fprintf('\nRegression: Max latency (from end) vs Next IBI\n');
    disp(mdl_ibi);
else
    scatter([],[]); axis off
    title('Not enough valid pairs for regression');
end


%% ================== REGRESSION: End Max Latency vs Latency (Max -> Next Min) ==================
% For each animal:
%   t_max_abs   = end_time + latency_to_max
%   t_min_abs   = start_time + latency_to_min
% Pair each end's t_max_abs with the *first* t_min_abs that occurs after it.
% Regress: predictor = latency_to_max (from end), response = latency (max -> next min)

all_max_latencies_from_end   = [];   % x: latency of max relative to bout end (s)
all_latencies_max_to_nextmin = [];   % y: time from that max to the next minimum (s)

ids = fieldnames(ResultsByAnimal);

for i = 1:numel(ids)
    R = ResultsByAnimal.(ids{i});

    % Absolute times (seconds) of extrema already detected
    t_end_s      = R.end.timestamps_s(:);
    lat_max_s    = R.end.max_latency_s(:);
    t_max_abs_s  = t_end_s + lat_max_s;   % absolute time of each max

    t_start_s    = R.start.timestamps_s(:);
    lat_min_s    = R.start.min_latency_s(:);
    t_min_abs_s  = t_start_s + lat_min_s; % absolute time of each min

    % Skip if no mins or no ends for this animal
    if isempty(t_min_abs_s) || isempty(t_max_abs_s)
        continue
    end

    % For each max time, find the next minimum occurring after it
    for j = 1:numel(t_max_abs_s)
        tmax = t_max_abs_s(j);

        % index of first minimum after this max
        ixNextMin = find(t_min_abs_s > tmax, 1, 'first');
        if isempty(ixNextMin)
            continue  % no subsequent minimum in this session
        end

        dt_max_to_min = t_min_abs_s(ixNextMin) - tmax;  % seconds

        % Keep only sensible positive intervals
        if ~isnan(dt_max_to_min) && dt_max_to_min > 0 && ~isnan(lat_max_s(j))
            all_max_latencies_from_end(end+1,1)   = lat_max_s(j);
            all_latencies_max_to_nextmin(end+1,1) = dt_max_to_min;
        end
    end
end

% Clean & check counts
valid = ~isnan(all_max_latencies_from_end) & ~isnan(all_latencies_max_to_nextmin);

figure('Name','Regression: End Max Latency vs (Max → Next Min)','Color','w');
if nnz(valid) >= 2
    % Fit simple linear model: y = beta0 + beta1 * x
    mdl = fitlm(all_max_latencies_from_end(valid), all_latencies_max_to_nextmin(valid));

    % Scatter + fitted line
    scatter(all_max_latencies_from_end(valid), all_latencies_max_to_nextmin(valid), 36, 'filled'); hold on;
    xgrid = linspace(min(all_max_latencies_from_end(valid)), max(all_max_latencies_from_end(valid)), 100)';
    yhat  = predict(mdl, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);

    xlabel('Latency to Maximum from Bout End (s)');
    ylabel('Latency from Maximum to Next Minimum (s)');
    title(sprintf('End Max Latency vs Max→Next Min (R^2=%.3f, p=%.3g)', ...
        mdl.Rsquared.Ordinary, mdl.Coefficients.pValue(2)));

    % Console details
    fprintf('\nRegression: End Max Latency vs (Max → Next Min)\n');
    disp(mdl);
else
    axis off
    text(0.5,0.5,'Not enough valid pairs for regression','HorizontalAlignment','center');
end




%%



%% ================== REGRESSION: (Min at Start_i → Max at End_i) vs Bout Duration_i ==================
% For each bout i:
%   t_min_abs(i) = start_i + latency_to_min(i)     (from R.start)
%   t_max_abs(i) = end_i   + latency_to_max(i)     (from R.end)
%   dt_min_to_max(i) = t_max_abs(i) - t_min_abs(i)
%   bout_duration(i) = end_i - start_i
% Then regress:  dt_min_to_max ~ bout_duration

all_dt_min_to_max = [];   % response (y)
all_bout_dur      = [];   % predictor (x) = bout duration

ids = fieldnames(ResultsByAnimal);
names_cell = {animals.name};
names_str  = string(names_cell);
tolerance_s = 1e-3;   % matching tolerance when mapping to original bouts

for ii = 1:numel(ids)
    id = string(ids{ii});
    R  = ResultsByAnimal.(ids{ii});

    % --- find this animal in `animals` struct ---
    idx = find(names_str == id, 1, 'first');
    if isempty(idx)
        warning('Pairwise Min→Max: animal "%s" not found in `animals`. Skipping.', id);
        continue
    end

    % --- keydowns and bout durations for this animal ---
    kd = animals(idx).keydown(:);
    if exist('keydown_in_minutes','var') && keydown_in_minutes
        kd_s = kd * 60;
    else
        kd_s = kd;
    end

    nPairs = floor(numel(kd_s)/2);
    if nPairs < 1
        continue
    end
    bout_starts = kd_s(1:2:(2*nPairs-1));  % length nPairs
    bout_ends   = kd_s(2:2:(2*nPairs));    % length nPairs
    bout_durs   = bout_ends - bout_starts; % length nPairs

    % --- kept extrema (absolute times) from your 600-sample windows ---
    kept_start_times = R.start.timestamps_s(:);          % starts that had full windows
    kept_end_times   = R.end.timestamps_s(:);            % ends   that had full windows
    t_min_abs        = kept_start_times + R.start.min_latency_s(:);
    t_max_abs        = kept_end_times   + R.end.max_latency_s(:);

    % Build quick maps from original bout start/end times -> indices inside kept arrays
    % (Find nearest kept start/end for each original bout time within tolerance)
    map_start_to_kept = nan(nPairs,1);
    map_end_to_kept   = nan(nPairs,1);

    for b = 1:nPairs
        if ~isempty(kept_start_times)
            [mindiff_s, ix_s] = min(abs(kept_start_times - bout_starts(b)));
            if ~isempty(ix_s) && mindiff_s <= tolerance_s
                map_start_to_kept(b) = ix_s;
            end
        end
        if ~isempty(kept_end_times)
            [mindiff_e, ix_e] = min(abs(kept_end_times - bout_ends(b)));
            if ~isempty(ix_e) && mindiff_e <= tolerance_s
                map_end_to_kept(b) = ix_e;
            end
        end
    end

    % For each bout i, if both start_i and end_i had valid 600-pt windows, compute dt_min_to_max
    for b = 1:nPairs
        ks = map_start_to_kept(b);
        ke = map_end_to_kept(b);
        if ~isnan(ks) && ~isnan(ke)
            dt_mm = t_max_abs(ke) - t_min_abs(ks);    % seconds
            if ~isnan(dt_mm) && dt_mm > 0 && ~isnan(bout_durs(b)) && bout_durs(b) > 0
                all_dt_min_to_max(end+1,1) = dt_mm;
                all_bout_dur(end+1,1)      = bout_durs(b);
            end
        end
    end
end

% --- Fit linear model and plot ---
valid = ~isnan(all_dt_min_to_max) & ~isnan(all_bout_dur);

figure('Name','Regression: (Min@Start → Max@End) vs Bout Duration','Color','w');
if nnz(valid) >= 2
    mdl_mm = fitlm(all_bout_dur(valid), all_dt_min_to_max(valid));  % y ~ x
    scatter(all_bout_dur(valid), all_dt_min_to_max(valid), 36, 'filled'); hold on;

    xgrid = linspace(min(all_bout_dur(valid)), max(all_bout_dur(valid)), 100)';
    yhat  = predict(mdl_mm, xgrid);
    plot(xgrid, yhat, 'LineWidth', 1.5);

    xlabel('Bout Duration (s)');
    ylabel('Latency: Min@Start → Max@End (s)');
    title(sprintf('Min→Max latency vs Bout Duration (R^2=%.3f, p=%.3g)', ...
        mdl_mm.Rsquared.Ordinary, mdl_mm.Coefficients.pValue(2)));

    fprintf('\nRegression: (Min@Start → Max@End) ~ Bout Duration\n');
    disp(mdl_mm);
else
    axis off
    text(0.5,0.5,'Not enough valid pairs for regression','HorizontalAlignment','center');
end
%%
%% ================== SUMMARY TABLE: Min→Max latency & Bout duration ==================
Summary_MinToMax = table();

ids = fieldnames(ResultsByAnimal);
names_cell = {animals.name};
names_str  = string(names_cell);
tolerance_s = 1e-3;

for ii = 1:numel(ids)
    id = string(ids{ii});
    R  = ResultsByAnimal.(ids{ii});

    % find this animal in animals struct
    idx = find(names_str == id, 1, 'first');
    if isempty(idx), continue; end

    % bout starts/ends
    kd = animals(idx).keydown(:);
    if exist('keydown_in_minutes','var') && keydown_in_minutes
        kd_s = kd * 60;
    else
        kd_s = kd;
    end

    nPairs = floor(numel(kd_s)/2);
    if nPairs < 1, continue; end

    bout_starts = kd_s(1:2:(2*nPairs-1));
    bout_ends   = kd_s(2:2:(2*nPairs));
    bout_durs   = bout_ends - bout_starts;

    % extrema absolute times
    kept_start_times = R.start.timestamps_s(:);
    kept_end_times   = R.end.timestamps_s(:);
    t_min_abs = kept_start_times + R.start.min_latency_s(:);
    t_max_abs = kept_end_times   + R.end.max_latency_s(:);

    % map bouts to kept extrema
    map_start_to_kept = nan(nPairs,1);
    map_end_to_kept   = nan(nPairs,1);

    for b = 1:nPairs
        if ~isempty(kept_start_times)
            [mindiff_s, ix_s] = min(abs(kept_start_times - bout_starts(b)));
            if ~isempty(ix_s) && mindiff_s <= tolerance_s
                map_start_to_kept(b) = ix_s;
            end
        end
        if ~isempty(kept_end_times)
            [mindiff_e, ix_e] = min(abs(kept_end_times - bout_ends(b)));
            if ~isempty(ix_e) && mindiff_e <= tolerance_s
                map_end_to_kept(b) = ix_e;
            end
        end
    end

    % collect rows for valid bouts
    for b = 1:nPairs
        ks = map_start_to_kept(b);
        ke = map_end_to_kept(b);
        if ~isnan(ks) && ~isnan(ke)
            dt_mm = t_max_abs(ke) - t_min_abs(ks);
            if ~isnan(dt_mm) && dt_mm > 0 && bout_durs(b) > 0
                newRow = table( ...
                    id, b, bout_durs(b), dt_mm, ...
                    'VariableNames', {'AnimalID','BoutIdx','BoutDuration_s','MinToMaxLatency_s'} ...
                );
                Summary_MinToMax = [Summary_MinToMax; newRow]; %#ok<AGROW>
            end
        end
    end
end

disp('===== Min→Max Latency & Bout Duration Summary =====');
disp(Summary_MinToMax(1:min(10,height(Summary_MinToMax)),:)); % show first 10
% Optional: save
% writetable(Summary_MinToMax,'Summary_MinToMax.csv');





%% ================== HISTOGRAMS: Bout Duration vs Min→Max Latency ==================
if ~exist('Summary_MinToMax','var') || isempty(Summary_MinToMax)
    warning('Summary_MinToMax not found. Run the regression/summary section first.');
else
    % Extract variables
    boutDur = Summary_MinToMax.BoutDuration_s;
    dtMinMax = Summary_MinToMax.MinToMaxLatency_s;

    % Remove NaNs or invalids
    boutDur = boutDur(~isnan(boutDur) & boutDur > 0);
    dtMinMax = dtMinMax(~isnan(dtMinMax) & dtMinMax > 0);

    % Plot side-by-side histograms
    figure('Name','Bout Duration & Min→Max Latency Distributions','Color','w');

    subplot(1,2,1);
    histogram(boutDur, 'BinWidth', 5, 'FaceColor',[0.2 0.6 0.8]); % 5s bins
    xlabel('Bout Duration (s)');
    ylabel('Bout Count');
    title('Distribution of Bout Durations');
    box off;

    subplot(1,2,2);
    histogram(dtMinMax, 'BinWidth', 1, 'FaceColor',[0.8 0.4 0.2]); % 1s bins
    xlabel('Min→Max Latency (s)');
    ylabel('Bout Count');
    title('Distribution of Min→Max Latencies');
    box off;
end


%%
%% ============================================================
% Peak timing within bout (minutes scale) + per-animal count table
% ============================================================

minBoutDur_min = 0;        % set >0 to exclude very short bouts
useInterpolation = true;   % true recommended

allRows = [];

for a = 1:numel(animals)

    name = string(animals(a).name);
    kd = animals(a).keydown(:);     % minutes
    fp = animals(a).fp;

    z = fp(:,1); z = z(:);
    t = fp(:,2); t = t(:);          % minutes

    % Make time unique/increasing
    [tU, ia] = unique(t, 'stable');
    zU = z(ia);

    starts = kd(1:2:end);
    ends   = kd(2:2:end);
    nB = min(numel(starts), numel(ends));
    starts = starts(1:nB);
    ends   = ends(1:nB);

    for b = 1:nB
        t0 = starts(b);
        t1 = ends(b);

        if ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
            continue;
        end

        dur = t1 - t0; % minutes
        if dur < minBoutDur_min
            continue;
        end

        % Skip if bout outside recording
        if t1 < tU(1) || t0 > tU(end)
            continue;
        end

        % Clip to recording bounds for extraction
        t0c = max(t0, tU(1));
        t1c = min(t1, tU(end));

        % Extract within bout and find peak
        if useInterpolation
            inBout = (tU >= t0c) & (tU <= t1c);
            tb = tU(inBout);
            zb = zU(inBout);
            if numel(tb) < 2 || all(isnan(zb)), continue; end
            [peakVal, peakIdx] = max(zb, [], 'omitnan');
            tPeak = tb(peakIdx);
        else
            [~, i0] = min(abs(tU - t0c));
            [~, i1] = min(abs(tU - t1c));
            if i1 < i0, tmp=i0; i0=i1; i1=tmp; end
            tb = tU(i0:i1);
            zb = zU(i0:i1);
            if isempty(tb) || all(isnan(zb)), continue; end
            [peakVal, peakIdx] = max(zb, [], 'omitnan');
            tPeak = tb(peakIdx);
        end

        if ~isfinite(peakVal) || ~isfinite(tPeak)
            continue;
        end

        % Fractional position within the ORIGINAL bout (t0..t1)
        frac = (tPeak - t0) / (t1 - t0);
        frac = max(0, min(1, frac)); % clamp

        if frac <= 0.25
            region = "First25";
        elseif frac < 0.75
            region = "Middle50";
        else
            region = "Last25";
        end

        allRows = [allRows; {a, name, b, region}]; %#ok<AGROW>
    end
end

% Per-bout classification table (kept lightweight)
BoutRegionTable = cell2table(allRows, ...
    'VariableNames', {'AnimalIdx','AnimalID','BoutIdx','PeakRegion'});

% ================== PER-ANIMAL COUNT TABLE ==================
% Initialize counts for all animals (even if 0 bouts in a category)
AnimalCountTable = table((1:numel(animals))', string({animals.name})', ...
    'VariableNames', {'AnimalIdx','AnimalID'});

AnimalCountTable.First25_Count  = zeros(numel(animals),1);
AnimalCountTable.Middle50_Count = zeros(numel(animals),1);
AnimalCountTable.Last25_Count   = zeros(numel(animals),1);

% Fill counts from BoutRegionTable
for a = 1:numel(animals)
    if isempty(BoutRegionTable), continue; end
    maskA = BoutRegionTable.AnimalIdx == a;

    AnimalCountTable.First25_Count(a)  = sum(maskA & BoutRegionTable.PeakRegion == "First25");
    AnimalCountTable.Middle50_Count(a) = sum(maskA & BoutRegionTable.PeakRegion == "Middle50");
    AnimalCountTable.Last25_Count(a)   = sum(maskA & BoutRegionTable.PeakRegion == "Last25");
end

disp(AnimalCountTable);



%%

% ============================================================
% Minimum ACh z-score timing within INTERBOUT intervals (minutes)
% + per-animal count table (First25 / Middle50 / Last25)
% Interbout interval = time between end(i) and start(i+1) for existing bouts
% ============================================================

minIBI_min = 0;            % set >0 to exclude very short IBIs (minutes)
useInterpolation = true;   % true recommended

allRowsIBI = [];

for a = 1:numel(animals)

    name = string(animals(a).name);
    kd = animals(a).keydown(:);     % minutes (odd=start, even=end)
    fp = animals(a).fp;

    z = fp(:,1); z = z(:);
    t = fp(:,2); t = t(:);          % minutes

    % Make time unique/increasing
    [tU, ia] = unique(t, 'stable');
    zU = z(ia);

    % Extract bouts
    starts = kd(1:2:end);
    ends   = kd(2:2:end);
    nB = min(numel(starts), numel(ends));
    starts = starts(1:nB);
    ends   = ends(1:nB);

    % Need at least 2 bouts to define an interbout interval
    if nB < 2
        continue;
    end

    for b = 1:(nB-1)
        tEnd   = ends(b);
        tStart = starts(b+1);

        % Interbout interval must be strictly between two valid bouts
        if ~isfinite(tEnd) || ~isfinite(tStart) || (tStart <= tEnd)
            continue;
        end

        ibiDur = tStart - tEnd;   % minutes
        if ibiDur < minIBI_min
            continue;
        end

        % Skip if IBI outside photometry recording
        if tStart < tU(1) || tEnd > tU(end)
            continue;
        end

        % Clip to recording bounds for extraction
        t0c = max(tEnd,  tU(1));
        t1c = min(tStart, tU(end));

        % Extract within IBI and find MINIMUM
        if useInterpolation
            inIBI = (tU >= t0c) & (tU <= t1c);
            tb = tU(inIBI);
            zb = zU(inIBI);
            if numel(tb) < 2 || all(isnan(zb)), continue; end
            [minVal, minIdx] = min(zb, [], 'omitnan');
            tMin = tb(minIdx);
        else
            [~, i0] = min(abs(tU - t0c));
            [~, i1] = min(abs(tU - t1c));
            if i1 < i0, tmp=i0; i0=i1; i1=tmp; end
            tb = tU(i0:i1);
            zb = zU(i0:i1);
            if isempty(tb) || all(isnan(zb)), continue; end
            [minVal, minIdx] = min(zb, [], 'omitnan');
            tMin = tb(minIdx);
        end

        if ~isfinite(minVal) || ~isfinite(tMin)
            continue;
        end

        % Fractional position within the ORIGINAL IBI (tEnd..tStart)
        frac = (tMin - tEnd) / (tStart - tEnd);
        frac = max(0, min(1, frac)); % clamp

        if frac <= 0.25
            region = "First25";
        elseif frac < 0.75
            region = "Middle50";
        else
            region = "Last25";
        end

        % Store one row per IBI
        % IBIIdx corresponds to the interval after bout b (between b and b+1)
        allRowsIBI = [allRowsIBI; {a, name, b, ibiDur, minVal, tMin, frac, region}]; %#ok<AGROW>
    end
end

% Per-IBI classification table (optional but useful)
IBIRegionTable = cell2table(allRowsIBI, ...
    'VariableNames', {'AnimalIdx','AnimalID','IBIIdx','IBIDur_min','MinZ','MinTime_min','MinFrac','MinRegion'});

% ================== PER-ANIMAL COUNT TABLE ==================
AnimalIBICountTable = table((1:numel(animals))', string({animals.name})', ...
    'VariableNames', {'AnimalIdx','AnimalID'});

AnimalIBICountTable.First25_Count  = zeros(numel(animals),1);
AnimalIBICountTable.Middle50_Count = zeros(numel(animals),1);
AnimalIBICountTable.Last25_Count   = zeros(numel(animals),1);

for a = 1:numel(animals)
    if isempty(IBIRegionTable), continue; end
    maskA = IBIRegionTable.AnimalIdx == a;

    AnimalIBICountTable.First25_Count(a)  = sum(maskA & IBIRegionTable.MinRegion == "First25");
    AnimalIBICountTable.Middle50_Count(a) = sum(maskA & IBIRegionTable.MinRegion == "Middle50");
    AnimalIBICountTable.Last25_Count(a)   = sum(maskA & IBIRegionTable.MinRegion == "Last25");
end

disp(AnimalIBICountTable);




%% Mean ACh z-score within active bouts vs interbout intervals
% Saves per-animal, per-interval means for later inspection
% Keydowns + FP time are in MINUTES (per your reminder)
% ============================================================

BoutInterboutMeansByAnimal = struct();
AllBoutMeansTable = table();
AllIBIMeansTable  = table();

minBoutDur_min = 0;   % optionally exclude short bouts, e.g., 5/60 for 5 sec
minIBIDur_min  = 0;   % optionally exclude short IBIs

for ai = 1:numel(animals)
    id = string(animals(ai).name);

    % --- Load FP data ---
    fp = animals(ai).fp;        % [z, t_min]
    z  = fp(:,1); z = z(:);
    t  = fp(:,2); t = t(:);     % minutes

    % Ensure time is unique/increasing for safe indexing
    [tU, ia] = unique(t, 'stable');
    zU = z(ia);

    % --- Load keydowns ---
    kd = animals(ai).keydown(:);         % minutes, odd=start even=end
    starts = kd(1:2:end);
    ends   = kd(2:2:end);
    nB = min(numel(starts), numel(ends));
    starts = starts(1:nB);
    ends   = ends(1:nB);

    % ---------- Per-bout means ----------
    boutMean = nan(nB,1);
    boutDur  = nan(nB,1);
    boutN    = nan(nB,1);

    for b = 1:nB
        t0 = starts(b);
        t1 = ends(b);

        if ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
            continue;
        end

        dur = t1 - t0;
        if dur < minBoutDur_min
            continue;
        end

        % Only consider overlap with recording
        if t1 < tU(1) || t0 > tU(end)
            continue;
        end

        % Clip to recording range
        t0c = max(t0, tU(1));
        t1c = min(t1, tU(end));

        idx = (tU >= t0c) & (tU <= t1c);
        if ~any(idx)
            continue;
        end

        boutMean(b) = mean(zU(idx), 'omitnan');
        boutDur(b)  = dur;
        boutN(b)    = sum(idx);
    end

    % ---------- Per-interbout means (between two existing bouts) ----------
    nIBI = max(nB-1, 0);
    ibiMean = nan(nIBI,1);
    ibiDur  = nan(nIBI,1);
    ibiN    = nan(nIBI,1);

    for b = 1:(nB-1)
        t0 = ends(b);        % IBI starts at end of bout b
        t1 = starts(b+1);    % IBI ends at start of bout b+1

        if ~isfinite(t0) || ~isfinite(t1) || t1 <= t0
            continue;
        end

        dur = t1 - t0;
        if dur < minIBIDur_min
            continue;
        end

        if t1 < tU(1) || t0 > tU(end)
            continue;
        end

        t0c = max(t0, tU(1));
        t1c = min(t1, tU(end));

        idx = (tU >= t0c) & (tU <= t1c);
        if ~any(idx)
            continue;
        end

        ibiMean(b) = mean(zU(idx), 'omitnan');
        ibiDur(b)  = dur;
        ibiN(b)    = sum(idx);
    end

    % ---------- Save per-animal results to struct ----------
    S = struct();
    S.animal_id = id;

    S.bout = table((1:nB)', starts, ends, boutDur, boutN, boutMean, ...
        'VariableNames', {'BoutIdx','Start_min','End_min','Dur_min','NSamples','MeanZ'});

    S.interbout = table((1:nIBI)', ends(1:nIBI), starts(2:nB), ibiDur, ibiN, ibiMean, ...
        'VariableNames', {'IBIIdx','EndBout_min','NextStart_min','Dur_min','NSamples','MeanZ'});

    BoutInterboutMeansByAnimal.(id) = S;

    % ---------- Append to long-format tables ----------
    tmpBout = S.bout;
    tmpBout.AnimalID = repmat(id, height(tmpBout), 1);
    tmpBout = movevars(tmpBout, 'AnimalID', 'Before', 'BoutIdx');
    AllBoutMeansTable = [AllBoutMeansTable; tmpBout]; %#ok<AGROW>

    tmpIBI = S.interbout;
    tmpIBI.AnimalID = repmat(id, height(tmpIBI), 1);
    tmpIBI = movevars(tmpIBI, 'AnimalID', 'Before', 'IBIIdx');
    AllIBIMeansTable = [AllIBIMeansTable; tmpIBI]; %#ok<AGROW>
end

% Inspect:
% - BoutInterboutMeansByAnimal.A8608.bout
% - BoutInterboutMeansByAnimal.A8608.interbout
% - AllBoutMeansTable / AllIBIMeansTable

disp(AllBoutMeansTable(1:min(10,height(AllBoutMeansTable)),:));
disp(AllIBIMeansTable(1:min(10,height(AllIBIMeansTable)),:));

% Optional: save to CSV for Prism or records
%writetable(AllBoutMeansTable, 'BoutMeanZ_ByAnimal.csv');
%writetable(AllIBIMeansTable,  'InterboutMeanZ_ByAnimal.csv');





%%

% ================== LOCAL FUNCTIONS ==================
function idx = nearestIdx(tSig, ts)
% Return index of tSig closest to timestamp ts
    [~, idx] = min(abs(tSig - ts));
end

function row = getRowPadded(zSig, centerIdx, N, T)
% Return a 1xT row centered at centerIdx with NaN padding at edges
    lo_des = centerIdx - N;
    hi_des = centerIdx + N;
    lo_src = max(lo_des, 1);
    hi_src = min(hi_des, numel(zSig));
    lo_dst = 1 + (lo_src - lo_des);
    hi_dst = T - (hi_des - hi_src);
    row = nan(1, T);
    row(lo_dst:hi_dst) = zSig(lo_src:hi_src);
end

function [tU, zU] = sanitizeTime(t, z)
% Make time unique & increasing; carry z along
    t = t(:); z = z(:);
    [tU, ia] = unique(t, 'stable');
    zU = z(ia);
end

%%


