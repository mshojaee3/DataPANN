% ============================================================
% Plot_Case_vs_Homogeneous.m
%
% SCRIPT 2 of 2.  Plots the Abaqus RUC case data (*_FINAL.csv) against the
% homogeneous MATLAB references stored ONE FILE PER CASE in
%   homog_reference/NC_k####_HOMOG.mat    (partner of NC_k####_FINAL.csv)
% It NEVER solves anything: run Build_Homogeneous_Reference_Cache.m first.
%
% What it draws (per quantity family, once per group of cases)
%   STEP 1  raw case data (solid) with the homogeneous reference (dashed)
%   STEP 2  Case / Homogenized ratio            (breaks down where the reference ~ 0)
%   STEP 3  NORMALISED ERROR  100*(Case - Homog)/scale   <-- robust to small components
%   STEP 4  parity plots (Case vs Homog, all cases pooled)
%   STEP 5  error summary: per-case and per-component, plus case_errors.csv
% A ratio Case/Homog explodes whenever a component of the reference passes
% near zero (off-diagonal P, most Q components), even if the absolute error
% is tiny. Steps 3-5 therefore divide the error by a SCALE that does not go
% to zero (default: the largest |Homog| of the whole family for that case),
% see NORM_SCALE.
% Cases are split into two groups: smooth curves vs. curves with abrupt
% changes (or small vs. large ||G||) -- see GROUP_MODE. By default the
% roughness is measured on the normalised error, not on the ratio.
%
% Also writes case_groups.csv and prints which CSV cases have no cache entry.
% ============================================================
clear; clc; close all

%% ================= CONFIG =================
DATA_DIR     = fileparts(mfilename('fullpath'));
FILE_PATTERN = '*_FINAL.csv';
REF_DIR      = fullfile(DATA_DIR, 'homog_reference');   % where the builder wrote the .mat files
CSV_SUFFIX_REGEX = '_FINAL\.csv$';   % NC_k0007_FINAL.csv -> NC_k0007_HOMOG.mat
HOMOG_SUFFIX     = '_HOMOG.mat';

SELECTED_CASES = [];                % [] = every CSV; or e.g. 0:20, [3 7 12]

% ---- what to plot ----
PLOT_RAW          = true;           % raw case data (not divided)
PLOT_RATIO        = true;           % Case / Homogenized
PLOT_ERR          = true;           % normalised error vs load step (robust to small components)
PLOT_PARITY       = true;           % Case-vs-Homog scatter, all cases pooled
PLOT_SUMMARY      = true;           % per-case / per-component error summary (+ case_errors.csv)

% Scale that the error is divided by (error [%] = 100*(Case-Homog)/scale):
%   'family'          : max |Homog| over ALL components and frames of the family, per case (default)
%   'component'       : max |Homog| of that single component over the load path, per case
%   'family_allcases' : max |Homog| of the family over all cases and frames (one number per family)
% The Energy family always uses 'component' (W_mean and ALLSE_mean have different magnitudes).
NORM_SCALE        = 'family';
OVERLAY_REFERENCE = true;           % raw plots: homogeneous reference dashed
PLOT_FAMILIES     = {'P','Q','Energy'};   % any of 'F','G','P','Q','Energy'
SHOW_LEGEND        = true;
MAX_LEGEND_ENTRIES = 40;
COLORMAP_NAME      = 'parula';
SAVE_FIGS          = false;         % true = also export PNG + FIG next to this script

% ---- two-group split ----
% 'roughness' : group 1 = smooth curves, group 2 = abrupt changes / spikes
% 'Gnorm'     : group 1 = small ||G||,   group 2 = large ||G||
% 'none'      : one group with everything
GROUP_MODE       = 'roughness';
GROUP_METRIC_SOURCE = 'err';        % 'err' = roughness of the normalised error (default, no
                                    % small-denominator artefacts); 'ratio' = old Case/Homog metric
GROUP_THRESHOLD  = [];              % [] = automatic (largest gap of log10 metric)
GROUP_COMPONENTS = {'P11','P12','P21','P22', ...
    'Q111','Q112','Q121','Q122','Q211','Q212','Q221','Q222','W_mean'};

% ---- sanity check: cached H and off-diagonal G (k~=l) vs. the CSV's own F(t), G(t) ----
CHECK_LOADCASE_MATCH = true;
LOADCASE_TOL         = 1e-3;        % warn if max|H_cache-H_csv|, max|G_cache-G_csv| exceed this

%% ================= CHECK REFERENCE FOLDER =================
if ~isfolder(REF_DIR)
    error('Plot:NoRefDir', ['Folder %s not found.\nRun Build_Homogeneous_Reference_Cache.m first ' ...
        '(or edit REF_DIR).'], REF_DIR);
end
nRefFiles = numel(dir(fullfile(REF_DIR, ['*' HOMOG_SUFFIX])));
fprintf('Found %d homogeneous reference file(s) in %s\n', nRefFiles, REF_DIR);

%% ================= READ CSV FILES =================
files = dir(fullfile(DATA_DIR, FILE_PATTERN));
files = files(~[files.isdir]);
[~, o] = sort({files.name});  files = files(o);
if isempty(files)
    error('Plot:NoFiles', 'No files matching "%s" in %s', FILE_PATTERN, DATA_DIR);
end

FAMILIES = {};
FAMILIES{end+1} = struct('tag','F','title','Deformation Gradient', ...
    'comps',{{'F11','F12','F21','F22'}}, 'grid',[1 4]);
FAMILIES{end+1} = struct('tag','G','title','Strain Gradient', ...
    'comps',{{'G111','G112','G121','G122','G211','G212','G221','G222'}}, 'grid',[2 4]);
FAMILIES{end+1} = struct('tag','P','title','First-Order Stress', ...
    'comps',{{'P11','P12','P21','P22'}}, 'grid',[1 4]);
FAMILIES{end+1} = struct('tag','Q','title','Higher-Order (Double) Stress', ...
    'comps',{{'Q111','Q112','Q121','Q122','Q211','Q212','Q221','Q222'}}, 'grid',[2 4]);
FAMILIES{end+1} = struct('tag','Energy','title','Strain Energy', ...
    'comps',{{'W_mean','ALLSE_mean'}}, 'grid',[1 2]);
allComp = {};
for ff = 1:numel(FAMILIES), allComp = [allComp, FAMILIES{ff}.comps]; end %#ok<AGROW>
gComps = FAMILIES{2}.comps;

results = struct([]);
missingIDs = [];  badMatch = [];  geomWarn = [];
for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);
    try
        Tc = readtable(fpath);
    catch ME
        warning('Could not read %s: %s', files(i).name, ME.message);  continue;
    end
    if isempty(Tc) || ~ismember('CaseID', Tc.Properties.VariableNames), continue; end
    CaseID = Tc.CaseID(1);
    if ~isempty(SELECTED_CASES) && ~ismember(CaseID, SELECTED_CASES), continue; end

    % partner file: same base name as the csv, _FINAL.csv -> _HOMOG.mat
    refFile = fullfile(REF_DIR, regexprep(files(i).name, CSV_SUFFIX_REGEX, HOMOG_SUFFIX));
    if ~isfile(refFile) && ismember('LoadCase', Tc.Properties.VariableNames)
        refFile = fullfile(REF_DIR, [char(string(Tc.LoadCase(1))) HOMOG_SUFFIX]);   % fallback by LoadCase name
    end
    if ~isfile(refFile)
        missingIDs(end+1) = CaseID; %#ok<AGROW>
        continue;
    end
    try
        ref = load(refFile, 'T', 'H', 'G', 'meta');
    catch ME
        warning('Could not read %s: %s', refFile, ME.message);
        missingIDs(end+1) = CaseID; %#ok<AGROW>
        continue;
    end
    Tref = ref.T;

    % geometry/mesh the reference was solved with vs. the Abaqus case
    if isfield(ref, 'meta')
        m = ref.meta;
        gd = max(abs([m.Lx_tot - Tc.Lx_tot(1), m.Ly_tot - Tc.Ly_tot(1), ...
            m.Rfrac - Tc.Rfrac(1), m.meshFrac - Tc.meshFrac(1), m.NFRAMES - Tc.NFRAMES(1)]));
        if gd > 1e-9
            geomWarn(end+1) = CaseID; %#ok<AGROW>
            warning('Plot:GeometryMismatch', ['CaseID %d: reference was solved with a different ' ...
                'Lx/Ly/Rfrac/meshFrac/NFRAMES than the Abaqus case (max diff %.3g).'], CaseID, gd);
        end
    end

    if CHECK_LOADCASE_MATCH && isfield(ref, 'H') && isfield(ref, 'G')
        [Hf, Gf] = fitLoadCaseFromCSV(Tc);
        % Only H and the off-diagonal G components (k ~= l: G112,G121,G212,G221)
        % are compared. The k == l components of the CSV's Gbar are the
        % least-squares slope of the actual F field (emergent, differ from the
        % prescribed t*G by the mean-displacement term), so they must NOT match.
        dH = max(abs(ref.H(:) - Hf(:)));
        dG = max(abs([ref.G(1,1,2)-Gf(1,1,2), ref.G(1,2,1)-Gf(1,2,1), ...
                      ref.G(2,1,2)-Gf(2,1,2), ref.G(2,2,1)-Gf(2,2,1)]));
        if dH > LOADCASE_TOL || dG > LOADCASE_TOL
            badMatch(end+1) = CaseID; %#ok<AGROW>
            warning('Plot:LoadCaseMismatch', ['CaseID %d: cached H / off-diagonal G differ from the CSV ' ...
                'load history (max|dH|=%.2e, max|dG(k~=l)|=%.2e) -- wrong case numbering?'], ...
                CaseID, dH, dG);
        end
    end

    n = min(height(Tc), height(Tref));
    if height(Tc) ~= height(Tref)
        warning('CaseID %d: frame count differs (case %d vs ref %d) -- truncated.', ...
            CaseID, height(Tc), height(Tref));
    end
    r = struct();
    r.file = files(i).name;  r.CaseID = CaseID;
    r.LoadStep = (1:n)';  r.StepTime = Tc.StepTime(1:n);
    r.caseVal = struct();  r.refVal = struct();  r.ratio = struct();
    r.err = struct();  r.rmsePct = struct();  r.famScale = struct();
    for c = 1:numel(allComp)
        nm = allComp{c};
        if ~ismember(nm, Tc.Properties.VariableNames) || ~ismember(nm, Tref.Properties.VariableNames)
            continue;
        end
        a = Tc.(nm)(1:n);  b = Tref.(nm)(1:n);
        % Abaqus ALLSE_mean is the EXTENSIVE (total) energy of the domain,
        % W_mean = ALLSE_mean / V0 is the density. The reference stores the
        % density in BOTH columns, so bring it to the same basis (x V0 = Area).
        if strcmp(nm, 'ALLSE_mean') && ismember('Area', Tc.Properties.VariableNames)
            b = b .* Tc.Area(1:n);
        end
        r.caseVal.(nm) = a;  r.refVal.(nm) = b;  r.ratio.(nm) = a ./ b;
    end
    g2 = zeros(height(Tc),1);
    for c = 1:numel(gComps)
        if ismember(gComps{c}, Tc.Properties.VariableNames), g2 = g2 + Tc.(gComps{c}).^2; end
    end
    r.Gnorm = max(sqrt(g2));
    if isempty(results), results = r; else, results(end+1) = r; end %#ok<SAGROW>
end

if ~isempty(missingIDs)
    fprintf(2, ['\nNO *%s file for %d CSV case(s): %s\n' ...
        '  -> they are NOT plotted. Run Build_Homogeneous_Reference_Cache.m for them.\n\n'], ...
        HOMOG_SUFFIX, numel(missingIDs), mat2str(sort(missingIDs)));
end
if ~isempty(geomWarn)
    fprintf(2, 'Geometry/mesh mismatch reference-vs-case for CaseID(s): %s\n', mat2str(sort(geomWarn)));
end
if isempty(results)
    error('Plot:NoResults', 'No CSV case has a reference .mat file -- nothing to plot.');
end
fprintf('Plotting %d case(s).\n', numel(results));

%% ================= NORMALISED ERROR (no division by small numbers) =================
% err.<comp> = 100*(Case - Homog)/scale ; scale never passes through zero.
famScaleAll = struct();
for ff = 1:numel(FAMILIES)
    tag = FAMILIES{ff}.tag;  comps = FAMILIES{ff}.comps;
    sAll = 0;
    for i = 1:numel(results)
        s = 0;
        for c = 1:numel(comps)
            if isfield(results(i).refVal, comps{c})
                s = max(s, max(abs(results(i).refVal.(comps{c}))));
            end
        end
        results(i).famScale.(tag) = s;
        sAll = max(sAll, s);
    end
    famScaleAll.(tag) = sAll;
end
for i = 1:numel(results)
    for ff = 1:numel(FAMILIES)
        tag = FAMILIES{ff}.tag;  comps = FAMILIES{ff}.comps;
        mode_i = NORM_SCALE;
        if strcmp(tag, 'Energy'), mode_i = 'component'; end
        for c = 1:numel(comps)
            nm = comps{c};
            if ~isfield(results(i).refVal, nm), continue; end
            switch mode_i
                case 'family',          sc = results(i).famScale.(tag);
                case 'component',       sc = max(abs(results(i).refVal.(nm)));
                case 'family_allcases', sc = famScaleAll.(tag);
                otherwise, error('Plot:BadNormScale', 'NORM_SCALE must be ''family'', ''component'' or ''family_allcases''.');
            end
            d = results(i).caseVal.(nm) - results(i).refVal.(nm);
            if sc > 0, e = 100 * d / sc; else, e = nan(size(d)); end
            results(i).err.(nm)     = e;
            results(i).rmsePct.(nm) = sqrt(mean(e.^2, 'omitnan'));
        end
    end
end
switch NORM_SCALE
    case 'family',          errLabel = 'error [% of max|Homog| of the family]';
    case 'component',       errLabel = 'error [% of max|Homog| of the component]';
    case 'family_allcases', errLabel = 'error [% of max|Homog| of the family, all cases]';
end

%% ================= TWO GROUPS =================
nR = numel(results);
switch GROUP_MODE
    case 'roughness'
        metric = nan(nR,1);
        for i = 1:nR, metric(i) = curveRoughness(results(i), GROUP_COMPONENTS, GROUP_METRIC_SOURCE); end
        if strcmp(GROUP_METRIC_SOURCE, 'err')
            metricLabel = 'max |2nd difference| of the normalised error [%]';
        else
            metricLabel = 'max |2nd difference| of Case/Homogenized';
        end
        groupNames  = {'Smooth curves', 'Abrupt changes'};
    case 'Gnorm'
        metric = [results.Gnorm]';
        metricLabel = 'max ||G||';
        groupNames  = {'Small strain gradient', 'Large strain gradient'};
    case 'none'
        metric = zeros(nR,1);  metricLabel = '(no grouping)';
        groupNames  = {'All cases', 'none'};
    otherwise
        error('Plot:BadGroupMode', 'GROUP_MODE must be ''roughness'', ''Gnorm'' or ''none''.');
end
if strcmp(GROUP_MODE, 'none')
    thr = Inf;  thrSrc = 'n/a';
elseif isempty(GROUP_THRESHOLD)
    thr = twoMeansThreshold(metric);  thrSrc = 'automatic (2-class split of log10 metric)';
else
    thr = GROUP_THRESHOLD;  thrSrc = 'user-defined';
end
groupOf  = 1 + ~(metric <= thr);
groupIdx = {find(groupOf == 1), find(groupOf == 2)};

fprintf('\nGrouping by %s | threshold = %.4g (%s)\n', metricLabel, thr, thrSrc);
[~, ordM] = sort(metric, 'descend');
for k = 1:min(nR, 60)
    q = ordM(k);
    fprintf('  CaseID %-4d metric %-12.4g group %d\n', results(q).CaseID, metric(q), groupOf(q));
end
if nR > 60, fprintf('  ... (%d more; all are in case_groups.csv)\n', nR-60); end
fprintf('  -> group 1 "%s": %d | group 2 "%s": %d\n\n', groupNames{1}, numel(groupIdx{1}), ...
    groupNames{2}, numel(groupIdx{2}));
try
    Tg = table([results.CaseID]', metric, groupOf, string(groupNames(groupOf))', ...
        'VariableNames', {'CaseID','Metric','Group','GroupName'});
    writetable(Tg, fullfile(DATA_DIR, 'case_groups.csv'));
catch ME
    warning('Could not write case_groups.csv (%s).', ME.message);
end

%% ================= PLOTS =================
modes = {};
if PLOT_RAW,   modes{end+1} = 'raw';   end %#ok<UNRCH>
if PLOT_RATIO, modes{end+1} = 'ratio'; end
if PLOT_ERR,   modes{end+1} = 'err';   end
for mm = 1:numel(modes)
    for ff = 1:numel(FAMILIES)
        fam = FAMILIES{ff};
        if ~ismember(fam.tag, PLOT_FAMILIES), continue; end
        for gg = 1:2
            idx = groupIdx{gg};
            if isempty(idx), continue; end
            lbl = sprintf('%s (%d case%s)', groupNames{gg}, numel(idx), repmat('s',1,double(numel(idx)~=1)));
            plotGroupGrid(results, idx, fam, modes{mm}, lbl, sprintf('group%d', gg), ...
                OVERLAY_REFERENCE, COLORMAP_NAME, SHOW_LEGEND, MAX_LEGEND_ENTRIES, SAVE_FIGS, DATA_DIR, errLabel);
        end
    end
end

if PLOT_PARITY
    for ff = 1:numel(FAMILIES)
        if ~ismember(FAMILIES{ff}.tag, PLOT_FAMILIES), continue; end
        plotParity(results, FAMILIES{ff}, SAVE_FIGS, DATA_DIR);
    end
end

if PLOT_SUMMARY
    Terr = table([results.CaseID]', 'VariableNames', {'CaseID'});
    for ff = 1:numel(FAMILIES)
        fam = FAMILIES{ff};
        if ~ismember(fam.tag, PLOT_FAMILIES), continue; end
        [relErr, maxErr] = plotErrorSummary(results, fam, errLabel, SAVE_FIGS, DATA_DIR);
        Terr.(['relErr_' fam.tag '_pct'])  = relErr;
        Terr.(['maxErr_' fam.tag '_pct'])  = maxErr;
    end
    Terr.Group = groupOf;
    try
        writetable(Terr, fullfile(DATA_DIR, 'case_errors.csv'));
        fprintf('Per-case errors written to %s\n', fullfile(DATA_DIR, 'case_errors.csv'));
    catch ME
        warning('Could not write case_errors.csv (%s).', ME.message);
    end
    fprintf(['\nrelErr_<fam>_pct = 100*||Case-Homog||_2 / ||Homog||_2 over all frames and components\n' ...
        '                   (one number per case, cannot blow up on small components)\n' ...
        'maxErr_<fam>_pct = worst single value of the normalised error (%s)\n'], errLabel);
    for ff = 1:numel(FAMILIES)
        fam = FAMILIES{ff};
        if ~ismember(fam.tag, PLOT_FAMILIES), continue; end
        v = Terr.(['relErr_' fam.tag '_pct']);
        fprintf('  %-7s relErr: median %.3g %% | 90th pct %.3g %% | max %.3g %% (CaseID %d)\n', fam.tag, ...
            median(v,'omitnan'), pctl(v,90), max(v), Terr.CaseID(find(v == max(v), 1)));
    end
end
fprintf('Done.\n');


% ============================================================
% LOCAL FUNCTIONS
% ============================================================
function [H, G] = fitLoadCaseFromCSV(Tc)
    t = Tc.StepTime;  d = sum(t.^2);
    H = [sum(t.*(Tc.F11-1)) sum(t.*Tc.F12); sum(t.*Tc.F21) sum(t.*(Tc.F22-1))] / d;
    G = zeros(2,2,2);
    G(1,1,1)=sum(t.*Tc.G111)/d; G(1,1,2)=sum(t.*Tc.G112)/d;
    G(1,2,1)=sum(t.*Tc.G121)/d; G(1,2,2)=sum(t.*Tc.G122)/d;
    G(2,1,1)=sum(t.*Tc.G211)/d; G(2,1,2)=sum(t.*Tc.G212)/d;
    G(2,2,1)=sum(t.*Tc.G221)/d; G(2,2,2)=sum(t.*Tc.G222)/d;
end

function m = curveRoughness(r, comps, source)
% Largest |second difference| along the load path over the listed components,
% of the normalised error (source 'err') or of the ratio (source 'ratio').
% ~0 for a straight line or gentle curve, large for a spike or kink.
    m = 0;
    for c = 1:numel(comps)
        if ~isfield(r.(source), comps{c}), continue; end
        y = r.(source).(comps{c});  y = y(isfinite(y));
        if numel(y) < 3, m = Inf; continue; end
        m = max(m, max(abs(diff(y, 2))));
    end
end

function thr = twoMeansThreshold(metric)
% 1-D two-means split of log10(metric); Inf if not separable.
    v = sort(log10(max(metric(isfinite(metric)), realmin)));
    n = numel(v);
    if n < 2 || (v(end) - v(1)) < eps, thr = Inf; return; end
    best = Inf; kBest = 1;
    for k = 1:n-1
        a = v(1:k); b = v(k+1:end);
        sse = sum((a-mean(a)).^2) + sum((b-mean(b)).^2);
        if sse < best, best = sse; kBest = k; end
    end
    thr = 10^(0.5*(v(kBest) + v(kBest+1)));
end

function plotGroupGrid(results, idx, fam, mode, groupLabel, gTag, overlayRef, cmapName, ...
        showLegend, maxLegendEntries, saveFigs, outDir, errLabel)
    isRaw = strcmp(mode, 'raw');
    switch mode
        case 'raw',   field = 'caseVal';
        case 'ratio', field = 'ratio';
        case 'err',   field = 'err';
    end
    comps = fam.comps;
    valid = false(size(comps));
    for c = 1:numel(comps)
        for s = idx(:)'
            if isfield(results(s).(field), comps{c}), valid(c) = true; break; end
        end
    end
    comps = comps(valid);
    if isempty(comps), return; end

    nSeries = numel(idx);
    try, cmap = feval(cmapName, max(nSeries, 2)); catch, cmap = parula(max(nSeries, 2)); end %#ok<NOCOM>
    doLegend = showLegend && nSeries <= maxLegendEntries;

    if isRaw
        if overlayRef
            ttl = sprintf('%s - raw data, solid = case, dashed = homogeneous ref  |  %s', fam.title, groupLabel);
        else
            ttl = sprintf('%s - raw case data  |  %s', fam.title, groupLabel);
        end
    elseif strcmp(mode, 'err')
        ttl = sprintf('%s - (Case - Homogenized) / scale, %s  |  %s', fam.title, errLabel, groupLabel);
    else
        ttl = sprintf('Normalized %s (Case / Homogenized)  |  %s', fam.title, groupLabel);
    end
    fig = figure('Name', ttl, 'Color', 'w');
    sgtitle(ttl, 'Interpreter', 'none');

    for c = 1:numel(comps)
        name = comps{c};
        ax = subplot(fam.grid(1), fam.grid(2), c); hold(ax, 'on'); grid(ax, 'on');
        hs = gobjects(0);  le = strings(0);
        for k = 1:nSeries
            r = results(idx(k));
            if ~isfield(r.(field), name), continue; end
            y = r.(field).(name);
            if all(isnan(y)), continue; end
            col = cmap(k, :);
            h = plot(ax, r.LoadStep, y, '-o', 'Color', col, 'MarkerSize', 3, ...
                'MarkerFaceColor', col, 'LineWidth', 1);
            if isRaw && overlayRef && isfield(r.refVal, name)
                plot(ax, r.LoadStep, r.refVal.(name), '--', 'Color', col, 'LineWidth', 1);
            end
            hs(end+1) = h; le(end+1) = sprintf('Case %d', r.CaseID); %#ok<AGROW>
        end
        title(ax, name, 'Interpreter', 'none');
        xlabel(ax, 'Load Step');
        if isRaw
            ylabel(ax, name, 'Interpreter', 'none');
        elseif strcmp(mode, 'err')
            yline(ax, 0, '--k');  ylabel(ax, 'Error [% of scale]');
        else
            yline(ax, 1, '--k');  ylabel(ax, 'Case / Homogenized');
        end
        if doLegend && c == 1 && ~isempty(hs)
            legend(ax, hs, le, 'Location', 'eastoutside', 'FontSize', 6);
        end
    end
    if saveFigs
        base = fullfile(outDir, sprintf('%s_%s_%s', mode, fam.tag, gTag));
        try, exportgraphics(fig, [base '.png'], 'Resolution', 300); catch, end %#ok<NOCOM>
        savefig(fig, [base '.fig']);
    end
end

function plotParity(results, fam, saveFigs, outDir)
% Case (y) vs Homogeneous (x), all cases and frames pooled, one panel per
% component, 1:1 line. Small components sit near the origin and are simply
% part of the cloud -- nothing is divided. The panel title carries the
% pooled relative L2 error of that component w.r.t. the largest |Homog|.
    comps = fam.comps;
    fig = figure('Name', ['Parity - ' fam.title], 'Color', 'w');
    sgtitle(sprintf('%s - parity: Case (y) vs Homogenized (x), all cases pooled', fam.title), 'Interpreter', 'none');
    np = 0;
    for c = 1:numel(comps)
        nm = comps{c};
        x = [];  y = [];
        for i = 1:numel(results)
            if isfield(results(i).refVal, nm)
                x = [x; results(i).refVal.(nm)]; %#ok<AGROW>
                y = [y; results(i).caseVal.(nm)]; %#ok<AGROW>
            end
        end
        if isempty(x), continue; end
        np = np + 1;
        ax = subplot(fam.grid(1), fam.grid(2), np); hold(ax, 'on'); grid(ax, 'on');
        scatter(ax, x, y, 8, 'filled', 'MarkerFaceAlpha', 0.35);
        lo = min([x; y]);  hi = max([x; y]);
        if lo == hi, lo = lo - 1; hi = hi + 1; end
        plot(ax, [lo hi], [lo hi], 'k--');
        xlim(ax, [lo hi]);  ylim(ax, [lo hi]);
        rel = 100 * sqrt(mean((y - x).^2, 'omitnan')) / max(max(abs(x)), realmin);
        title(ax, sprintf('%s  (RMSE = %.2g %% of max|Homog|)', nm, rel), 'Interpreter', 'none', 'FontSize', 8);
        xlabel(ax, 'Homogenized');  ylabel(ax, 'Case');
    end
    if saveFigs
        base = fullfile(outDir, sprintf('parity_%s', fam.tag));
        try, exportgraphics(fig, [base '.png'], 'Resolution', 300); catch, end %#ok<NOCOM>
        savefig(fig, [base '.fig']);
    end
end

function [relErr, maxErr] = plotErrorSummary(results, fam, errLabel, saveFigs, outDir)
% (a) one number per case: relErr = 100*||Case-Homog||_2/||Homog||_2 over all
%     frames and components of the family (Energy: W_mean only).
% (b) per component: distribution over cases of the RMSE of the normalised
%     error along the load path (median bar, 90th percentile and max markers).
    nR = numel(results);
    comps = fam.comps;
    if strcmp(fam.tag, 'Energy'), comps = {'W_mean'}; end
    relErr = nan(nR, 1);  maxErr = nan(nR, 1);
    E = nan(nR, numel(comps));
    for i = 1:nR
        num = 0;  den = 0;  mx = 0;
        for c = 1:numel(comps)
            nm = comps{c};
            if ~isfield(results(i).refVal, nm), continue; end
            d = results(i).caseVal.(nm) - results(i).refVal.(nm);
            num = num + sum(d.^2, 'omitnan');
            den = den + sum(results(i).refVal.(nm).^2, 'omitnan');
            mx  = max(mx, max(abs(results(i).err.(nm)), [], 'omitnan'));
            E(i, c) = results(i).rmsePct.(nm);
        end
        if den > 0, relErr(i) = 100 * sqrt(num / den); end
        maxErr(i) = mx;
    end

    fig = figure('Name', ['Error summary - ' fam.title], 'Color', 'w');
    sgtitle(sprintf('%s - error summary (no division by small components)', fam.title), 'Interpreter', 'none');

    ax1 = subplot(1, 2, 1); hold(ax1, 'on'); grid(ax1, 'on');
    ids = [results.CaseID]';
    semilogy(ax1, ids, max(relErr, eps), 'o', 'MarkerSize', 4, 'MarkerFaceColor', [0.2 0.4 0.8], 'Color', [0.2 0.4 0.8]);
    set(ax1, 'YScale', 'log');
    yline(ax1, median(relErr, 'omitnan'), '--k', sprintf('median %.2g %%', median(relErr, 'omitnan')));
    xlabel(ax1, 'CaseID');  ylabel(ax1, '100 ||Case-Homog|| / ||Homog||  [%]');
    title(ax1, 'Relative L2 error per case (all frames, all components)');

    ax2 = subplot(1, 2, 2); hold(ax2, 'on'); grid(ax2, 'on');
    med = median(E, 1, 'omitnan');  p90 = nan(1, numel(comps));
    for c = 1:numel(comps), p90(c) = pctl(E(:, c), 90); end
     mxc = max(E, [], 1, 'omitnan');
    bar(ax2, 1:numel(comps), med, 0.6);
    plot(ax2, 1:numel(comps), p90, 'r^', 'MarkerFaceColor', 'r');
    plot(ax2, 1:numel(comps), mxc, 'kv', 'MarkerFaceColor', 'k');
    set(ax2, 'XTick', 1:numel(comps), 'XTickLabel', comps, 'TickLabelInterpreter', 'none');
    ylabel(ax2, ['RMSE along load path, ' errLabel]);
    legend(ax2, {'median over cases', '90th percentile', 'max'}, 'Location', 'northwest');
    title(ax2, 'Per component (over cases)');

    if saveFigs
        base = fullfile(outDir, sprintf('errsummary_%s', fam.tag));
        try, exportgraphics(fig, [base '.png'], 'Resolution', 300); catch, end %#ok<NOCOM>
        savefig(fig, [base '.fig']);
    end
end

function q = pctl(x, pc)
% percentile without the Statistics Toolbox (linear interpolation, NaNs ignored)
    x = sort(x(isfinite(x)));
    n = numel(x);
    if n == 0, q = NaN; return; end
    pos = 1 + (n - 1) * pc / 100;
    lo = floor(pos);  hi = ceil(pos);
    q = x(lo) + (pos - lo) * (x(hi) - x(lo));
end