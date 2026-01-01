%% Photometry and Barnes Spatial Food Seeking (SFS) task:
% This script will be dedicated to the analysis of fiber photometry data
% with bcorrectavioral data from the SFS task.

%% Importing bcorrectavioral data from anymaze:
%Bcorrectavioral data from AnyMaze comes in the form of a csv file. It is a
%table containing state values of regions across the maze in time series
%(seconds).
%HeadDistancetoEscapcorrectole # = Distance from the escape hole specified will
%be provided in meters.
%InEscapcorrectole #, ErrortoEscapcorrectole #, AdjtoEscapcorrectole # = binary state
%value of which hole the animal enters. The data stream does not know which
%is the correct escape hole, this is only tracking all possible escape hole
%and tracking entries to all errors relative to each escape hole.
%Information of each animals' escape hole will be contained in a separate
%file where the bcorrectavioral performance is analyzed.
%Once you know which is the correct escape hole for an animal you are
%analyzing you will import only certain columns from the csv file: for
%example, if escape hole (correct) = 3, you will import, Distancetocorrect3, Incorrect3,
%Adjtocorrect3, and Errorstocorrect3.

%When correct columns are selected change the column name from:
%This will make it easier to access the data under it.
% Time -> time
%HeadDistancetoEscapcorrectole # -> distance
%InEscapcorrectole # -> correct
% ErrortoEscapcorrectole # -> error
% AdjtoEscapcorrectole # -> adj


%%
%Before continuing, the following will convert your corrected data into
%zscore and store in its appropriate spot for analysis:

FP.sfs.D2708_sfs = FP.corrected_data;
FP.sfs.zscore.D2708(:,1) = FP.sfs.D2708_sfs(:,1);
FP.sfs.zscore.D2708(:,2) = zscore(FP.sfs.D2708_sfs(:,2));

%% Plotting distance from correct hole and neural activity on the same graph:
%The following will plot the distance of the animal from the correct hole
%and the trace of neural activity on the same graph:

figure
grid on
yyaxis left
ylabel('Z-score value of activity')
plot(FP.sfs.zscore.D2708(:,1), FP.sfs.zscore.D2708(:,2),'g','LineWidth',0.7)
yyaxis right
ylabel('Distance from the correct correct (meters)')
plot(D2708_sfs_behavior.time, D2708_sfs_behavior.distance,'-k','LineWidth',1);
hold on
title('Z-score value of neural activity and distance from the correct escape hole during SFS');
xlabel('Time (seconds)')




%% After running the script barnes_investigation_filtering and parsing through investigations on a predetermined criteria
% Allocate investigation times to structure that fits the naming 
% schema below

FP.sfs.error.D2708.start = error(:,1);
FP.sfs.correct.D2708.start = correct(:,1)


%% Extracting neural activity data from hole investigations:
%Before extracting neural activity data from investigations, we will need
%to filter out invesitgations that are not at least 5 seconds apart, such
%that the extracted activity trace does not overlap in between
%investigations.
%% ================== EXTRACT activity TRACES (WITH TIME CUTOFF) ==================

% ------------------ USER SETTINGS ------------------
nPre  = 40;    % samples before event
nPost = 60;    % samples after event
winLength = nPre + nPost + 1;

cutoffTime_s = 120;   % <-- ONLY keep events occurring before this time (seconds)
% ---------------------------------------------------

% Load photometry data
t_fp = FP.sfs.zscore.D2708(:,1);   % time (s)
z_fp = FP.sfs.zscore.D2708(:,2);   % ACh z-score

% Event timestamps (s)
t_correct_all = FP.sfs.correct.D2708.start(:);
t_error_all   = FP.sfs.error.D2708.start(:);

% -------- Apply cutoff --------
t_correct = t_correct_all(t_correct_all <= cutoffTime_s);
t_error   = t_error_all(t_error_all <= cutoffTime_s);

% Preallocate output matrices
FP.sfs.correct.D2708_ach = nan(numel(t_correct), winLength);
FP.sfs.error.D2708_ach   = nan(numel(t_error),   winLength);

% Helper: nearest FP index to event time
getIdx = @(t) find(abs(t_fp - t) == min(abs(t_fp - t)), 1, 'first');

% -------- CORRECT INVESTIGATIONS --------
for i = 1:numel(t_correct)
    idx0 = getIdx(t_correct(i));
    i0 = idx0 - nPre;
    i1 = idx0 + nPost;

    % Skip if window exceeds recording bounds
    if i0 < 1 || i1 > numel(z_fp)
        continue
    end

    FP.sfs.correct.D2708_ach(i,:) = z_fp(i0:i1)';
end

% -------- ERROR INVESTIGATIONS --------
for i = 1:numel(t_error)
    idx0 = getIdx(t_error(i));
    i0 = idx0 - nPre;
    i1 = idx0 + nPost;

    % Skip if window exceeds recording bounds
    if i0 < 1 || i1 > numel(z_fp)
        continue
    end

    FP.sfs.error.D2708_ach(i,:) = z_fp(i0:i1)';
end

% -------- Sanity check --------
fprintf('Correct trials (≤ %.1f s): %d\n', cutoffTime_s, ...
    sum(~all(isnan(FP.sfs.correct.D2708_ach),2)));
fprintf('Error trials   (≤ %.1f s): %d\n', cutoffTime_s, ...
    sum(~all(isnan(FP.sfs.error.D2708_ach),2)));










%% Converting the extracted investigation activity traces to its own zscore:
[correct_h correct_l] = size(FP.sfs.correct.D2708_ach);
[error_h error_l] = size(FP.sfs.error.D2708_ach);

for i = 1:correct_h
    FP.sfs.correct.D2708_ach(i,:) = zscore(FP.sfs.correct.D2708_ach(i,:));
end

for i = 1:error_h
    FP.sfs.error.D2708_ach(i,:) = zscore(FP.sfs.error.D2708_ach(i,:));
end
    


%%
%% ==================  individual trial traces (heatmap) ==================

M = FP.sfs.correct.D2708_ach;          % [trials x timepoints]
M = double(M);

nT = size(M,2);
t  = linspace(-2, 3, nT);              % seconds axis inferred from window

figure('Color','w','Name','D2602 correct trials (heatmap)');
imagesc(t, 1:size(M,1), M);
set(gca,'YDir','normal');
xlabel('Time from investigation onset (s)');
ylabel('Trial # (in order)');
title(sprintf('D2602 correct: %d trials', size(M,1)));
colorbar; ylabel(colorbar,'ACh z-score');

hold on;
xline(0,'k-','LineWidth',1);
hold off;






%% ================== MIXED-EFFECTS: ACh z-score ~ Distance (m) * Treatment ==================
% Uses behavior timepoints as the master timebase (one row per behavior sample).
% Each behavior timestamp is paired to the nearest-in-time FP sample.

% ---- Animal list + treatment ----
animalIDs  = ["D2601","D2602","D2708","D2604","D2606","D2607","D2702","D2703","D2704"];
treatment  = ["SDV"  ,"SDV"  ,"SDV"  ,"Sham" ,"Sham" ,"Sham" ,"SDV"  ,"SDV"  ,"Sham"];

% Optional: reject pairings if nearest FP sample is too far away in time
max_dt_s = Inf;   % e.g., 0.2 means require within 200 ms; Inf disables cutoff

% ---- Accumulator table ----
T = table();

for ai = 1:numel(animalIDs)
    id  = animalIDs(ai);
    trt = treatment(ai);

    % ---------- Pull FP ----------
    if ~isfield(FP.sfs.zscore, id)
        warning('%s: FP.sfs.zscore.%s not found. Skipping.', id, id);
        continue
    end
    fpMat = FP.sfs.zscore.(id);     % [N x 2]: col1=time(s), col2=z
    t_fp  = fpMat(:,1);  z_fp = fpMat(:,2);
    t_fp = t_fp(:); z_fp = z_fp(:);

    % ---------- Pull behavior table ----------
    behavVar = id + "_sfs_behavior";  % e.g., "D2602_sfs_behavior"
    if ~evalin('base', sprintf('exist(''%s'',''var'')', behavVar))
        warning('%s: behavior table %s not found. Skipping.', id, behavVar);
        continue
    end
    B = evalin('base', behavVar);

    if ~all(ismember(["time","distance"], string(B.Properties.VariableNames)))
        warning('%s: %s must contain variables .time and .distance. Skipping.', id, behavVar);
        continue
    end
    t_b = B.time(:);
    d_b = B.distance(:);  % meters

    % ---------- Clean ----------
    good_fp = ~isnan(t_fp) & ~isnan(z_fp);
    t_fp = t_fp(good_fp); z_fp = z_fp(good_fp);
    [t_fp, ia] = unique(t_fp, 'stable'); z_fp = z_fp(ia);

    good_b = ~isnan(t_b) & ~isnan(d_b);
    t_b = t_b(good_b); d_b = d_b(good_b);
    [t_b, ib] = unique(t_b, 'stable'); d_b = d_b(ib);

    if numel(t_b) < 3 || numel(t_fp) < 3
        warning('%s: not enough data after cleaning. Skipping.', id);
        continue
    end

    % ---------- Pair each behavior time to nearest FP sample ----------
    % Nearest z at behavior times
    z_at_b = interp1(t_fp, z_fp, t_b, 'nearest', NaN);

    % Get actual nearest indices to compute dt (and optionally threshold)
    if exist('knnsearch','file') == 2
        idx_nn = knnsearch(t_fp, t_b);
    else
        idx_nn = arrayfun(@(x) find(abs(t_fp-x)==min(abs(t_fp-x)),1,'first'), t_b);
    end
    dt_nn = abs(t_fp(idx_nn) - t_b);

    keep = ~isnan(z_at_b) & ~isnan(d_b) & ~isnan(dt_nn);
    if isfinite(max_dt_s)
        keep = keep & (dt_nn <= max_dt_s);
    end

    if nnz(keep) < 3
        warning('%s: not enough paired points after dt filter. Skipping.', id);
        continue
    end

    % ---------- Append rows ----------
    Ti = table();
    Ti.AnimalID  = repmat(categorical(id), nnz(keep), 1);
    Ti.Treatment = repmat(categorical(trt), nnz(keep), 1);
    Ti.Distance_m = d_b(keep);
    Ti.ACh_z      = z_at_b(keep);
    Ti.dt_s       = dt_nn(keep);

    T = [T; Ti]; %#ok<AGROW>
end

% ---- Basic sanity check ----
if isempty(T)
    error('No paired data assembled. Check variable names / presence of behavior tables.');
end

fprintf('Total paired samples: %d\n', height(T));
% ---- Summary: how many paired samples per Animal x Treatment (no groupsummary) ----
[G, gAnimal, gTrt] = findgroups(T.AnimalID, T.Treatment);

Npairs = splitapply(@numel, T.ACh_z, G);
SummaryCounts = table(gAnimal, gTrt, Npairs, ...
    'VariableNames', {'AnimalID','Treatment','Npairs'});

disp(SummaryCounts);

% Optional: center distance to improve interpretability/stability
T.Distance_c = T.Distance_m - mean(T.Distance_m, 'omitnan');

% ================== FIT MIXED-EFFECTS MODELS ==================
% Random intercept model (baseline differs by animal)
lme_int = fitlme(T, 'ACh_z ~ Distance_c * Treatment + (1|AnimalID)');

% Random intercept + random slope model (distance effect can differ by animal)
lme_slope = fitlme(T, 'ACh_z ~ Distance_c * Treatment + (Distance_c|AnimalID)');

fprintf('\n=== Mixed model (random intercept) ===\n');
disp(lme_int);

fprintf('\n=== Mixed model (random intercept + random slope) ===\n');
disp(lme_slope);

% Compare models (likelihood ratio test); if p<0.05, prefer random-slope model
cmp = compare(lme_int, lme_slope);
fprintf('\n=== Model comparison (intercept vs intercept+slope) ===\n');
disp(cmp);

% Choose which model to report/use for plotting
useRandomSlope = (cmp.pValue(2) < 0.05);   % row 2 is the more complex model
if useRandomSlope
    lme = lme_slope;
    fprintf('Using random-slope model (better fit by LRT).\n');
else
    lme = lme_int;
    fprintf('Using random-intercept model.\n');
end

% ================== VISUALIZATION: DATA + FIXED-EFFECT PREDICTIONS ==================
% Make a prediction grid over distance for each treatment using FIXED effects only
dGrid_m = linspace(min(T.Distance_m), max(T.Distance_m), 200)';
dGrid_c = dGrid_m - mean(T.Distance_m, 'omitnan');

Gsham = table(categorical(repmat("Sham", numel(dGrid_m), 1)), dGrid_c, ...
    'VariableNames', {'Treatment','Distance_c'});
Gsdv  = table(categorical(repmat("SDV",  numel(dGrid_m), 1)), dGrid_c, ...
    'VariableNames', {'Treatment','Distance_c'});

% Add AnimalID as a placeholder category (needed for predict); we’ll use fixed effects only
Gsham.AnimalID = categorical(repmat(string(T.AnimalID(1)), numel(dGrid_m), 1));
Gsdv.AnimalID  = categorical(repmat(string(T.AnimalID(1)), numel(dGrid_m), 1));

yhat_sham = predict(lme, Gsham, 'Conditional', false);
yhat_sdv  = predict(lme, Gsdv,  'Conditional', false);

figure('Color','w','Name','ACh vs Distance by Treatment (Mixed Model)');
gscatter(T.Distance_m, T.ACh_z, T.Treatment); hold on;
plot(dGrid_m, yhat_sham, 'LineWidth', 2);
plot(dGrid_m, yhat_sdv,  'LineWidth', 2);
xlabel('Distance from correct hole (m)');
ylabel('ACh z-score');
title('Mixed-effects model: fixed-effect trends (lines) over raw paired samples');
box off;

% ================== OPTIONAL: REPORT KEY FIXED EFFECTS ==================
% This is usually what you interpret:
%  - Distance_c: overall distance effect (reference treatment)
%  - Treatment: group shift at mean distance
%  - Distance_c:Treatment: interaction (does SDV change distance–ACh relationship?)
fprintf('\n=== Fixed effects (what to interpret) ===\n');
disp(lme.Coefficients);

% If you want the interaction p-value explicitly:
ixInt = find(contains(string(lme.Coefficients.Name), "Distance_c:Treatment"));
if ~isempty(ixInt)
    fprintf('Interaction (Distance x Treatment) p = %.3g\n', lme.Coefficients.pValue(ixInt));
end



%%
%% ================== BUILD PAIRED TABLE: Distance (behavior) -> nearest ACh ==================
% Output table: Tdist
% Columns: AnimalID, Treatment, Time_s, Distance_m, ACh_z

% ---- Define animals and treatments (as you provided) ----
animalIDs  = {'D2601','D2602','D2708','D2604','D2606','D2607','D2702','D2703','D2704'};
treatments = {'SDV',  'SDV',  'SDV',  'Sham', 'Sham', 'Sham', 'SDV',  'SDV',  'Sham'};

Tdist = table();  % accumulator

for i = 1:numel(animalIDs)
    id  = animalIDs{i};
    trt = treatments{i};

    % ---------- Pull FP data ----------
    % FP.sfs.zscore.<ID> : col1=time(s), col2=ACh z
    if ~isfield(FP.sfs.zscore, id)
        warning('%s: missing FP.sfs.zscore.%s. Skipping.', id, id);
        continue
    end
    fpMat = FP.sfs.zscore.(id);
    t_fp  = fpMat(:,1);  % seconds
    ach_z = fpMat(:,2);

    % sanitize FP time (strictly increasing unique)
    [t_fp_u, ia] = unique(t_fp(:), 'stable');
    ach_u = ach_z(ia);

    % ---------- Pull behavior table ----------
    % expects a variable in workspace named like "D2602_sfs_behavior"
    behVarName = [id '_sfs_behavior'];
    if ~evalin('base', sprintf('exist(''%s'',''var'')', behVarName))
        warning('%s: missing %s table. Skipping.', id, behVarName);
        continue
    end
    beh = evalin('base', behVarName);

    if ~istable(beh) || ~all(ismember({'time','distance'}, beh.Properties.VariableNames))
        warning('%s: %s exists but is not a table with columns .time and .distance. Skipping.', id, behVarName);
        continue
    end

    t_beh = beh.time(:);       % seconds
    dist  = beh.distance(:);   % meters

    % ---------- Nearest-time pairing: behavior time -> nearest FP sample ----------
    % interp1 with 'nearest' returns the ACh value at closest time
    ach_at_beh = interp1(t_fp_u, ach_u, t_beh, 'nearest', NaN);

    % keep only valid pairs
    good = ~isnan(t_beh) & ~isnan(dist) & ~isnan(ach_at_beh);

    if nnz(good) == 0
        warning('%s: no valid paired samples after NaN removal.', id);
        continue
    end

    % ---------- Add rows to Tdist ----------
    newRows = table( ...
        repmat({id}, nnz(good), 1), ...
        repmat({trt}, nnz(good), 1), ...
        t_beh(good), ...
        dist(good), ...
        ach_at_beh(good), ...
        'VariableNames', {'AnimalID','Treatment','Time_s','Distance_m','ACh_z'} );

    Tdist = [Tdist; newRows]; %#ok<AGROW>
end

% Make types consistent (older MATLAB friendly)
Tdist.AnimalID  = categorical(Tdist.AnimalID);
Tdist.Treatment = categorical(Tdist.Treatment, {'Sham','SDV'}, 'Ordinal', true);

disp('=== Paired distance–ACh table built ===');
disp(Tdist(1:min(10,height(Tdist)),:));
fprintf('Total paired samples: %d\n', height(Tdist));



% ================== SIMPLE REGRESSIONS: ACh (Y) vs Distance (X) by treatment ==================
% Outputs for Prism:
%   XYdist_Sham  (AnimalID, Distance_m, ACh_z)
%   XYdist_SDV

colSham = [0.6 0.6 0.6];      % gray
colSDV  = [0.55 0.35 0.75];   % purple

% Split
idxSham = (Tdist.Treatment == 'Sham');
idxSDV  = (Tdist.Treatment == 'SDV');

% Prism tables
XYdist_Sham = table(string(Tdist.AnimalID(idxSham)), Tdist.Distance_m(idxSham), Tdist.ACh_z(idxSham), ...
    'VariableNames', {'AnimalID','Distance_m','ACh_z'});

XYdist_SDV  = table(string(Tdist.AnimalID(idxSDV)),  Tdist.Distance_m(idxSDV),  Tdist.ACh_z(idxSDV), ...
    'VariableNames', {'AnimalID','Distance_m','ACh_z'});

% Fit models
mdl_sham = fitlm(Tdist.Distance_m(idxSham), Tdist.ACh_z(idxSham));
mdl_sdv  = fitlm(Tdist.Distance_m(idxSDV),  Tdist.ACh_z(idxSDV));

% Correlation (Pearson r)
r_sham = corr(Tdist.Distance_m(idxSham), Tdist.ACh_z(idxSham), 'Rows','complete');
r_sdv  = corr(Tdist.Distance_m(idxSDV),  Tdist.ACh_z(idxSDV),  'Rows','complete');

% Plot
figure('Color','w','Name','Simple linear regressions: ACh vs Distance','Position',[200 200 900 400]);

% ---- Sham ----
subplot(1,2,1); hold on
x = Tdist.Distance_m(idxSham); y = Tdist.ACh_z(idxSham);
scatter(x, y, 18, colSham, 'filled', 'MarkerFaceAlpha', 0.35);

xgrid = linspace(min(x), max(x), 200)';
plot(xgrid, predict(mdl_sham, xgrid), 'Color', colSham, 'LineWidth', 2);

xlabel('Distance from correct hole (m)');
ylabel('ACh z-score');
title(sprintf('Sham: ACh vs Distance\nR = %.2f, p(slope)=%.3g', r_sham, mdl_sham.Coefficients.pValue(2)));
box off

% ---- SDV ----
subplot(1,2,2); hold on
x = Tdist.Distance_m(idxSDV); y = Tdist.ACh_z(idxSDV);
scatter(x, y, 18, colSDV, 'filled', 'MarkerFaceAlpha', 0.35);

xgrid = linspace(min(x), max(x), 200)';
plot(xgrid, predict(mdl_sdv, xgrid), 'Color', colSDV, 'LineWidth', 2);

xlabel('Distance from correct hole (m)');
ylabel('ACh z-score');
title(sprintf('SDV: ACh vs Distance\nR = %.2f, p(slope)=%.3g', r_sdv, mdl_sdv.Coefficients.pValue(2)));
box off

% Optional: export for Prism
% writetable(XYdist_Sham, 'XYdist_Sham_ACh_vs_Distance.csv');
% writetable(XYdist_SDV,  'XYdist_SDV_ACh_vs_Distance.csv');










%% ================== DELTA SPEED vs DELTA ACh REGRESSIONS ==================

figure('Color','w','Position',[200 200 900 400]);

treatments = ["Sham","SDV"];
cols = [0.6 0.6 0.6;      % Sham = gray
        0.55 0.35 0.75]; % SDV = purple

% Initialize containers for Prism export
XYd_Sham = table();
XYd_SDV  = table();

for i = 1:2
    trt = treatments(i);
    col = cols(i,:);

    % --- select valid rows ---
    idx = T.Treatment == trt & ...
          ~isnan(T.Speed_mps) & ...
          ~isnan(T.ACh_z);

    % Extract
    speed = T.Speed_mps(idx);
    ach   = T.ACh_z(idx);
    id    = string(T.AnimalID(idx));

    % --- compute first differences ---
    dSpeed = diff(speed);
    dACh   = diff(ach);
    id_d   = id(2:end);   % align IDs after diff

    % Remove NaNs introduced by diff
    good = ~isnan(dSpeed) & ~isnan(dACh);
    dSpeed = dSpeed(good);
    dACh   = dACh(good);
    id_d   = id_d(good);

    % ---- Save XY pairs (for Prism) ----
    tmpTbl = table(id_d, dSpeed, dACh, ...
        'VariableNames', {'AnimalID','DeltaSpeed_mps','DeltaACh_z'});

    if trt == "Sham"
        XYd_Sham = tmpTbl;
    else
        XYd_SDV = tmpTbl;
    end

    % ---- Plot ----
    subplot(1,2,i); hold on;

    scatter(dSpeed, dACh, 20, col, 'filled', 'MarkerFaceAlpha', 0.4);

    mdl = fitlm(dSpeed, dACh);
    xgrid = linspace(min(dSpeed), max(dSpeed), 200)';
    yhat  = predict(mdl, xgrid);

    plot(xgrid, yhat, 'Color', col, 'LineWidth', 2);

    R = corr(dSpeed, dACh, 'Rows','complete');

    xlabel('\Delta Speed (m/s)');
    ylabel('\Delta ACh z-score');
    title(sprintf('%s: ΔACh vs ΔSpeed\nR = %.2f, p = %.3g', ...
        trt, R, mdl.Coefficients.pValue(2)));

    box off
end

%% ================== AVERAGE RAW SPEED PER ANIMAL ==================

AnimalIDs = categories(T.AnimalID);
AvgSpeed  = nan(numel(AnimalIDs),1);
Treatment = strings(numel(AnimalIDs),1);

for i = 1:numel(AnimalIDs)
    idx = T.AnimalID == AnimalIDs{i};
    AvgSpeed(i)  = mean(T.Speed_mps(idx), 'omitnan');
    Treatment(i) = string(T.Treatment(find(idx,1,'first')));
end

AvgSpeedTable = table( ...
    string(AnimalIDs), ...
    Treatment, ...
    AvgSpeed, ...
    'VariableNames', {'AnimalID','Treatment','AvgSpeed_mps'} );

disp('=== Average raw speed per animal (m/s) ===');
disp(AvgSpeedTable);


%% ================== Mixed Effects Model- delta Speed and delta ACh: BUILD DELTA TABLE (within-animal) ==================
% Output: Td (one row per time step per animal)
% Columns: AnimalID, Treatment, dSpeed, dACh, (optional) dt

% Ensure categorical types
if ~iscategorical(T.AnimalID),  T.AnimalID  = categorical(T.AnimalID); end
if ~iscategorical(T.Treatment), T.Treatment = categorical(T.Treatment); end
T.Treatment = categorical(T.Treatment, {'Sham','SDV'}, 'Ordinal', true);

% Sort within animal by time if available; else keep current order
hasTime = ismember('Time_s', T.Properties.VariableNames);

Td = table();

animalsU = categories(T.AnimalID);
for i = 1:numel(animalsU)
    id  = animalsU{i};
    idx = (T.AnimalID == id);

    Ti = T(idx, :);

    % Sort by time if present
    if hasTime
        Ti = sortrows(Ti, 'Time_s', 'ascend');
    end

    % Pull vectors
    sp = Ti.Speed_mps;
    ac = Ti.ACh_z;

    % Optional dt (seconds between samples) if Time_s exists
    if hasTime
        tt = Ti.Time_s;
        dt = diff(tt);
    else
        dt = nan(height(Ti)-1,1);
    end

    % First differences
    dSpeed = diff(sp);
    dACh   = diff(ac);

    % Align metadata after diff
    AnimalID  = Ti.AnimalID(2:end);
    Treatment = Ti.Treatment(2:end);

    % Keep finite rows only
    good = isfinite(dSpeed) & isfinite(dACh);
    if hasTime
        good = good & isfinite(dt) & dt > 0;
    end

    if nnz(good) < 5
        continue
    end

    newRows = table(AnimalID(good), Treatment(good), dSpeed(good), dACh(good), ...
        'VariableNames', {'AnimalID','Treatment','dSpeed','dACh'});

    if hasTime
        newRows.dt = dt(good);
    end

    Td = [Td; newRows]; %#ok<AGROW>
end

disp('=== Delta table summary ===');
fprintf('Total delta pairs: %d\n',  height(Td));
disp(tabulate(Td.Treatment));

% ================== MIXED-EFFECTS MODELS: dACh ~ dSpeed * Treatment ==================
% Centering dSpeed is optional; for interpretability you can leave it raw.
% If you want, you can add: Td.dSpeed_c = Td.dSpeed - mean(Td.dSpeed,'omitnan');

% Random intercept model
lme_d_int = fitlme(Td, 'dACh ~ 1 + Treatment*dSpeed + (1|AnimalID)', ...
    'FitMethod','ML');

disp('=== Mixed model (random intercept) ===');
disp(lme_d_int);

% Random intercept + random slope model
lme_d_slope = fitlme(Td, 'dACh ~ 1 + Treatment*dSpeed + (1 + dSpeed|AnimalID)', ...
    'FitMethod','ML');

disp('=== Mixed model (random intercept + random slope) ===');
disp(lme_d_slope);

% ================== MODEL COMPARISON (LRT) ==================
cmp = compare(lme_d_int, lme_d_slope);   % likelihood ratio test
disp('=== Model comparison (intercept vs intercept+slope) ===');
disp(cmp);

% Choose model with better fit (lower AIC/BIC or significant LRT)
if cmp.pValue(2) < 0.05
    lme_d = lme_d_slope;
    disp('Using random-slope model (better fit by LRT).');
else
    lme_d = lme_d_int;
    disp('Using random-intercept model (no significant LRT improvement).');
end

disp('=== Fixed effects (interpret these) ===');
disp(lme_d.Coefficients);










%% ================== PANEL B: Example animal time series (Speed + ACh) ==================
% Adds shaded blocks for correct/error investigations.

% --------- USER CHOICE: pick example animal ----------
exampleID = "D2607";     % <-- change as desired

% --------- USER CHOICE: block duration (seconds) ----------
blockDur_s = 1.0;        % width of each shaded block after investigation onset

% --------- Colors ----------
colSham = [0.6 0.6 0.6];      % gray
colSDV  = [0.55 0.35 0.75];   % purple
colCorrect = [0.65 0.85 0.65]; % light green
colError   = [0.90 0.65 0.65]; % light red
alphaShade = 0.18;            % transparency for blocks

% --------- Safety checks ----------
reqVars = {'AnimalID','Treatment','Speed_mps','ACh_z'};
missing = setdiff(reqVars, T.Properties.VariableNames);
if ~isempty(missing)
    error('T is missing required columns: %s', strjoin(missing, ', '));
end

hasTime = ismember('Time_s', T.Properties.VariableNames);

% Extract this animal
idx = string(T.AnimalID) == exampleID;
if nnz(idx) < 5
    error('Not enough samples found for %s in T (found %d).', exampleID, nnz(idx));
end

Ti = T(idx, :);

% Sort by time if available
if hasTime
    Ti = sortrows(Ti, 'Time_s', 'ascend');
    t  = Ti.Time_s;
else
    t  = (1:height(Ti))';  % fallback: sample index
end

speed = Ti.Speed_mps;
ach   = Ti.ACh_z;

% Clean NaNs
good = isfinite(speed) & isfinite(ach) & isfinite(t);
t = t(good); speed = speed(good); ach = ach(good);

% Recompute deltas (aligned to timepoints 2:end)
dSpeed = diff(speed);
dACh   = diff(ach);
t_d    = t(2:end);

% Treatment label for title
trt = string(Ti.Treatment(find(good,1,'first')));

if trt == "Sham"
    colTrt = colSham;
else
    colTrt = colSDV;
end

% --------- Get correct/error investigation onset times (seconds) ----------
% Expected: FP.sfs.correct.<exampleID>.start and FP.sfs.error.<exampleID>.start
tCorr = [];
tErr  = [];

if isfield(FP, 'sfs') && isfield(FP.sfs, 'correct') && isfield(FP.sfs.correct, exampleID)
    if isfield(FP.sfs.correct.(exampleID), 'start')
        tCorr = FP.sfs.correct.(exampleID).start(:);
    end
end

if isfield(FP, 'sfs') && isfield(FP.sfs, 'error') && isfield(FP.sfs.error, exampleID)
    if isfield(FP.sfs.error.(exampleID), 'start')
        tErr = FP.sfs.error.(exampleID).start(:);
    end
end

% --------- Plot ----------
figure('Color','w','Name','Panel B: Example animal traces', ...
       'Units','pixels','Position',[150 150 900 550]);

tl = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% ---- 1) Speed trace ----
ax1 = nexttile; hold on;
plot(t, speed, 'LineWidth', 1.6, 'Color', [0.2 0.2 0.2]);
ylabel('Speed (m/s)');
title(sprintf('Example animal %s (%s): Speed and ACh over time', exampleID, trt));
box off

% ---- 2) ACh trace ----
ax2 = nexttile; hold on;
plot(t, ach, 'LineWidth', 1.6, 'Color', colTrt);
ylabel('ACh (z-score)');
box off

% ---- 3) Deltas ----
ax3 = nexttile; hold on;
plot(t_d, dSpeed, 'LineWidth', 1.2, 'Color', [0.2 0.2 0.2 0.8]);
plot(t_d, dACh,   'LineWidth', 1.2, 'Color', [colTrt 0.8]);
yline(0,'k-','LineWidth',0.8);

if hasTime
    xlabel('Time (s)');
else
    xlabel('Sample # (time column not available)');
end
ylabel('\Delta Speed (m/s), \Delta ACh (z)');
legend({'\Delta Speed','\Delta ACh'}, 'Location','best');
box off

% Link x-axes for easy reading
linkaxes([ax1 ax2 ax3], 'x');

% Limit to first 120 s
tStart = t(1);
tStop  = min(t(1)+120, t(end));
xlim(ax1, [tStart tStop]);

% --------- Add shaded investigation blocks (behind traces) ----------
axesList = [ax1 ax2 ax3];

% Keep only events that overlap the displayed window
keepCorr = tCorr >= (tStart - blockDur_s) & tCorr <= tStop;
keepErr  = tErr  >= (tStart - blockDur_s) & tErr  <= tStop;
tCorr = tCorr(keepCorr);
tErr  = tErr(keepErr);

for ax = axesList
    yl = ylim(ax);

    % Correct blocks (light green)
    for k = 1:numel(tCorr)
        x0 = tCorr(k);
        x1 = x0 + blockDur_s;
        p = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], colCorrect, ...
            'EdgeColor','none', 'FaceAlpha', alphaShade);
        uistack(p, 'bottom');
    end

    % Error blocks (light red)
    for k = 1:numel(tErr)
        x0 = tErr(k);
        x1 = x0 + blockDur_s;
        p = patch(ax, [x0 x1 x1 x0], [yl(1) yl(1) yl(2) yl(2)], colError, ...
            'EdgeColor','none', 'FaceAlpha', alphaShade);
        uistack(p, 'bottom');
    end

    % Re-apply y-limits (patch can sometimes trigger autoscale changes)
    ylim(ax, yl);
end

% Optional: add panel label "B"
annotation('textbox',[0.01 0.92 0.05 0.07], 'String','B', ...
    'FontSize',16,'FontWeight','bold','LineStyle','none','Color','k');

% Optional: warn if events weren't found
if isempty(tCorr) && isempty(tErr)
    warning('No correct/error start times found for %s in FP.sfs.correct/error.(ID).start', exampleID);
end




%% ================== PANEL C: Per-animal regression lines (ΔACh vs ΔSpeed) ==================
% Requires:
%   Td with columns: AnimalID (categorical), Treatment (categorical), dSpeed, dACh
% Optional:
%   lme_d (fitlme object) to overlay fixed-effect lines

% ---- checks ----
reqTd = {'AnimalID','Treatment','dSpeed','dACh'};
missingTd = setdiff(reqTd, Td.Properties.VariableNames);
if ~isempty(missingTd)
    error('Td is missing required columns: %s', strjoin(missingTd, ', '));
end

% Force consistent treatment ordering
if ~iscategorical(Td.Treatment), Td.Treatment = categorical(Td.Treatment); end
Td.Treatment = categorical(Td.Treatment, {'Sham','SDV'}, 'Ordinal', true);

% Colors (your scheme)
colSham = [0.6 0.6 0.6];      % gray
colSDV  = [0.55 0.35 0.75];   % purple

% ---- plotting range for x ----
xMin = prctile(Td.dSpeed, 2);
xMax = prctile(Td.dSpeed, 98);
xgrid = linspace(xMin, xMax, 200)';

% ---- make figure ----
figure('Color','w','Name','Panel C: Per-animal ΔACh vs ΔSpeed lines', ...
       'Units','pixels','Position',[150 150 650 520]);
hold on;

% Optional: light background cloud of all points
scatter(Td.dSpeed, Td.dACh, 8, [0.2 0.2 0.2], 'filled', 'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15);

animalsU = categories(Td.AnimalID);

% Store per-animal slopes (useful for Panel D)
animalIDs = strings(numel(animalsU),1);
animalTrt = strings(numel(animalsU),1);
olsSlope  = nan(numel(animalsU),1);
olsInt    = nan(numel(animalsU),1);
nPairs    = nan(numel(animalsU),1);

for i = 1:numel(animalsU)
    id = animalsU{i};
    idx = (Td.AnimalID == id);

    if nnz(idx) < 20   % minimum points to fit a stable line (tweak as desired)
        continue
    end

    % treatment for this animal
    trt = string(Td.Treatment(find(idx,1,'first')));

    x = Td.dSpeed(idx);
    y = Td.dACh(idx);

    % OLS within animal
    mdl_i = fitlm(x, y);
    b0 = mdl_i.Coefficients.Estimate(1);
    b1 = mdl_i.Coefficients.Estimate(2);

    yhat = b0 + b1*xgrid;

    if trt == "Sham"
        plot(xgrid, yhat, 'Color', [colSham 0.35], 'LineWidth', 1);
    else
        plot(xgrid, yhat, 'Color', [colSDV 0.35], 'LineWidth', 1);
    end

    % Save for later (Panel D)
    animalIDs(i) = string(id);
    animalTrt(i) = trt;
    olsSlope(i)  = b1;
    olsInt(i)    = b0;
    nPairs(i)    = nnz(idx);
end

% ---- overlay fixed-effect (marginal) lines from mixed model (optional) ----

if exist('lme_d','var') && isa(lme_d,'LinearMixedModel')

    % Use the SAME categorical "type" as Td (categories + ordinalness)
    trtCats = categories(Td.Treatment);
    trtOrd  = isordinal(Td.Treatment);

    aniCats = categories(Td.AnimalID);
    aniOrd  = isordinal(Td.AnimalID);

    % Dummy animal ID (must be a valid category)
    a0 = categorical({aniCats{1}}, aniCats, 'Ordinal', aniOrd);

    % Treatment vectors that match Td.Treatment exactly
    trt_sham = categorical(repmat({'Sham'}, numel(xgrid), 1), trtCats, 'Ordinal', trtOrd);
    trt_sdv  = categorical(repmat({'SDV'},  numel(xgrid), 1), trtCats, 'Ordinal', trtOrd);

    % AnimalID vector that matches Td.AnimalID exactly
    ani_vec = categorical(repmat({aniCats{1}}, numel(xgrid), 1), aniCats, 'Ordinal', aniOrd);

    % Prediction tables (variable names must match model)
    G_sham = table(trt_sham, ani_vec, xgrid, 'VariableNames', {'Treatment','AnimalID','dSpeed'});
    G_sdv  = table(trt_sdv,  ani_vec, xgrid, 'VariableNames', {'Treatment','AnimalID','dSpeed'});

    % Predict marginal (fixed effects only)
    y_fix_sham = predict(lme_d, G_sham, 'Conditional', false);
    y_fix_sdv  = predict(lme_d, G_sdv,  'Conditional', false);

    plot(xgrid, y_fix_sham, 'Color', colSham, 'LineWidth', 3);
    plot(xgrid, y_fix_sdv,  'Color', colSDV,  'LineWidth', 3);
end

% ---- final axis labels ----
xlabel('\Delta Speed (m/s)');
ylabel('\Delta ACh (z-score)');
title('Panel C: Within-animal coupling (ΔACh vs ΔSpeed)');
yline(0,'k-','LineWidth',0.7);
xline(0,'k-','LineWidth',0.7);
box off

% Panel label "C"
annotation('textbox',[0.01 0.92 0.05 0.07], 'String','C', ...
    'FontSize',16,'FontWeight','bold','LineStyle','none','Color','k');

% ---- Export slope table for Panel D ----
good = ~isnan(olsSlope) & animalIDs ~= "";
SlopeTable_Delta = table(animalIDs(good), animalTrt(good), nPairs(good), olsSlope(good), olsInt(good), ...
    'VariableNames', {'AnimalID','Treatment','Npairs','OLS_Slope','OLS_Intercept'});

disp('=== Panel C: per-animal OLS slope table (for Panel D) ===');
disp(SlopeTable_Delta);

%%
%% ================== PANEL D: Distribution of per-animal slopes ==================

% Enforce treatment order
SlopeTable_Delta.Treatment = categorical( ...
    SlopeTable_Delta.Treatment, {'Sham','SDV'}, 'Ordinal', false);

% Colors
colSham = [0.6 0.6 0.6];      % gray
colSDV  = [0.55 0.35 0.75];   % purple

% Extract data
slopes_sham = SlopeTable_Delta.OLS_Slope(SlopeTable_Delta.Treatment == 'Sham');
slopes_sdv  = SlopeTable_Delta.OLS_Slope(SlopeTable_Delta.Treatment == 'SDV');

nSham = numel(slopes_sham);
nSDV  = numel(slopes_sdv);

% Stats
[p_ttest,~,~] = ttest2(slopes_sham, slopes_sdv);
[p_rs,~,~]    = ranksum(slopes_sham, slopes_sdv);

% ---- Plot ----
figure('Color','w','Name','Panel D: Per-animal slope distributions', ...
       'Units','pixels','Position',[200 200 520 480]);
hold on;

% Boxplots drawn separately (version-safe)
boxchart( ...
    categorical(repmat("Sham", nSham, 1), {'Sham','SDV'}), ...
    slopes_sham, ...
    'BoxFaceColor', colSham, ...
    'BoxFaceAlpha', 0.25, ...
    'LineWidth', 1.5);

boxchart( ...
    categorical(repmat("SDV", nSDV, 1), {'Sham','SDV'}), ...
    slopes_sdv, ...
    'BoxFaceColor', colSDV, ...
    'BoxFaceAlpha', 0.25, ...
    'LineWidth', 1.5);

% Overlay individual animals
swarmchart( ...
    categorical(repmat("Sham", nSham, 1), {'Sham','SDV'}), ...
    slopes_sham, ...
    70, colSham, 'filled', 'MarkerFaceAlpha', 0.9);

swarmchart( ...
    categorical(repmat("SDV", nSDV, 1), {'Sham','SDV'}), ...
    slopes_sdv, ...
    70, colSDV, 'filled', 'MarkerFaceAlpha', 0.9);

% Reference line
yline(0, 'k--', 'LineWidth', 1);

ylabel('\DeltaACh / \DeltaSpeed slope');
xlabel('Treatment');
title('Panel D: Animal-level coupling strength');

% Sample sizes
text(1, min([slopes_sham;slopes_sdv])*0.9, sprintf('n = %d', nSham), ...
    'HorizontalAlignment','center','FontSize',10);
text(2, min([slopes_sham;slopes_sdv])*0.9, sprintf('n = %d', nSDV), ...
    'HorizontalAlignment','center','FontSize',10);

% Stats annotation
statsStr = sprintf('t-test p = %.3f\nranksum p = %.3f', p_ttest, p_rs);
text(1.5, max([slopes_sham;slopes_sdv])*0.9, statsStr, ...
    'HorizontalAlignment','center','FontSize',10);

box off

% Panel label
annotation('textbox',[0.01 0.92 0.05 0.07], 'String','D', ...
    'FontSize',16,'FontWeight','bold','LineStyle','none','Color','k');

%% ================== SUPPLEMENT TABLES: Mixed-effects model stats (robust across MATLAB versions) ==================
finalModel = lme_d;

% ---- 1) Fixed effects with 95% CI (manual, no coefCI) ----
FE = finalModel.Coefficients;
FE.Term = string(FE.Name);

alpha = 0.05;
tcrit = tinv(1 - alpha/2, FE.DF);
FE.CI_Lower = FE.Estimate - tcrit .* FE.SE;
FE.CI_Upper = FE.Estimate + tcrit .* FE.SE;

Sup_FE = FE(:, {'Term','Estimate','SE','tStat','DF','pValue','CI_Lower','CI_Upper'});

disp('=== Supplement: Fixed effects (with 95% CI) ===');
disp(Sup_FE);

% ---- 2) Model fit summary (robust: no Deviance property, no string(Formula)) ----
deviance_val = -2 * finalModel.LogLikelihood;

% Robust formula-to-text conversion
try
    formula_txt = char(finalModel.Formula);
catch
    formula_txt = 'Formula unavailable in this MATLAB version';
end

Sup_Fit = table( ...
    finalModel.NumObservations, ...
    finalModel.ModelCriterion.AIC, ...
    finalModel.ModelCriterion.BIC, ...
    finalModel.LogLikelihood, ...
    deviance_val, ...
    {formula_txt}, ...  % store as cell so table accepts it cleanly
    'VariableNames', {'Nobs','AIC','BIC','LogLik','Deviance_m2LL','Formula'} );

disp('=== Supplement: Model fit ===');
disp(Sup_Fit);

% ---- 3) Random effects / covariance parameters (as provided) ----
try
    Sup_RE = finalModel.CovarianceParameters;
catch
    Sup_RE = table("CovarianceParameters not available as a table in this MATLAB version.", ...
        'VariableNames', {'Note'});
end

disp('=== Supplement: Random effects / covariance parameters ===');
disp(Sup_RE);

% ---- 4) Model comparison (LRT) ----
Sup_LRT = table();
if exist('cmp','var') && istable(cmp)
    Sup_LRT = cmp;
elseif exist('lme_d_int','var') && exist('lme_d_slope','var')
    Sup_LRT = compare(lme_d_int, lme_d_slope);
end

if ~isempty(Sup_LRT)
    disp('=== Supplement: Model comparison (LRT) ===');
    disp(Sup_LRT);
end

%% ---- 5) Optional: export CSVs ----
% writetable(Sup_FE,  'SupTable_FixedEffects_lme_d.csv');
% writetable(Sup_Fit, 'SupTable_ModelFit_lme_d.csv');
% if istable(Sup_RE),  writetable(Sup_RE, 'SupTable_RandomEffects_lme_d.csv'); end
% if ~isempty(Sup_LRT), writetable(Sup_LRT,'SupTable_ModelCompare_LRT.csv'); end