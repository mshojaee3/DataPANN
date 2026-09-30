% ============================================================
%  compare_case_vs_homogeneous.m
%
%  Reads every generated RUC case CSV (*_FINAL.csv, e.g. NC_k0000_FINAL.
%  csv .. NC_k0272_FINAL.csv -- any count) in this script's own folder,
%  computes the matching "homogeneous, 1 RUC, kappa = 1" FE reference for
%  each case's own load direction via RUC_FE_2D_homogeneous.m, divides
%  the case's own F/G/P/Q/energy history by that reference frame-by-
%  frame, and plots the ratios (Case / Homogenized) vs load step -- the
%  same kind of plot as the two reference images (Normalized Strain
%  Gradient, Normalized First-Order Stress).
%
%  Put this file (and RUC_FE_2D_homogeneous.m) in the SAME folder as the
%  *_FINAL.csv files and just run it.
%
%  ------------------------------------------------------------
%  2026-09-29 UPDATE -- RAW PLOTS FIRST, THEN NORMALISED; TWO GROUPS
%
%  For every quantity family (see PLOT_FAMILIES) the script now draws
%    1) the RAW case data with no division (optionally with the
%       homogeneous reference dashed on top), and then
%    2) the same data divided by the homogeneous reference,
%  separately for two groups of cases: smooth curves vs curves with abrupt
%  changes (GROUP_MODE='roughness'), or small vs large strain gradient
%  (GROUP_MODE='Gnorm'). The group assignment is printed and saved to
%  case_groups.csv.
%
%  ------------------------------------------------------------
%  2026-09-23 UPDATE -- REFERENCE MESH NOW MATCHES ABAQUS'S OWN GEOMETRY
%
%  Main2D_hyperelastic.py (the script that generated the *_FINAL.csv case
%  data) ALWAYS partitions the circular inclusion boundary into its
%  Abaqus geometry, regardless of kappaMat -- there is no kappa check in
%  its geometry/meshing block. Every case, including the standard sweep's
%  kappaMat=1 runs, sits on a circle-partitioned mesh.
%
%  RUC_FE_2D_homogeneous.m previously defaulted to MeshType='auto', which
%  silently swapped in a plain REGULAR grid (no inclusion partition) for
%  the reference solve whenever kappaMat==1 -- i.e. for every case in the
%  standard ensemble. That made the "homogeneous reference" a different
%  mesh TOPOLOGY than the paired Abaqus case, which this project's own
%  earlier analysis (a t^2-vs-t^1 error-scaling argument, see project
%  memory) predicts will show up disproportionately in the double stress
%  Qbar -- exactly the large Case/Homogeneous spread observed in Q while
%  P/W stayed tight.
%
%  Fixed here two ways:
%    1) REF_MESH_TYPE below is now explicitly passed as 'inclusion' to
%       every RUC_FE_2D_homogeneous call, so the reference always uses
%       the same circle-partitioned geometry Abaqus does (this is also
%       now RUC_FE_2D_homogeneous.m's own new default, so this is
%       belt-and-suspenders).
%    2) refKey below now carries a 'v2_<meshtype>_' prefix, so every
%       reference computed under the OLD (mesh-mismatched) behavior is
%       automatically treated as a cache miss and recomputed -- you do
%       NOT need to manually delete homog_reference_cache.mat, but the
%       first run after this update WILL re-solve every unique load case
%       (same one-time cost as the original cold-cache run).
%
%  This only matches mesh TOPOLOGY (same geometry, same target element
%  size, same element type/order/integration) -- MATLAB's PDE Toolbox
%  mesher and Abaqus's own free triangular mesher are different
%  algorithms and will not produce a literally identical triangulation
%  even given the same domain and nominal element size. For a genuinely
%  node-for-node/element-for-element identical mesh, export Abaqus's own
%  mesh and pass it through RUC_FE_2D_homogeneous's new 'MeshData'
%  option (see that file's header) instead of letting it re-mesh.
%
%  Also added (opt-in, default OFF): USE_TRUE_TENSORS, which feeds the
%  reference solve the load case's TRUE library-prescribed (H,G) (via
%  `python pann_loadcases.py --tensors`) instead of a least-squares fit
%  recovered from the case's own reported F/G history. The fit is
%  normally accurate to ~KubcErrF/KubcErrG (typically 1e-5..1e-8) and is
%  NOT expected to be the dominant source of the Q mismatch -- the mesh
%  fix above is -- but it removes one more (smaller) source of "close but
%  not identical" between case and reference. Requires pann_loadcases.py
%  and its sample_points_d##.txt next to LIBPY_PATH (or on the MATLAB
%  path); falls back to the least-squares fit with a warning if
%  unavailable, so the script stays runnable without Python.
%  ------------------------------------------------------------
%
%  ------------------------------------------------------------
%  HOW THE HOMOGENEOUS REFERENCE IS OBTAINED FOR EACH CASE
%
%  Each case CSV already reports its own boundary-integrated F11..F22 /
%  G111..G222 at every frame (Method-B boundary integration, this
%  project's established pipeline). Under KUBC, macroscopic Fbar = I +
%  t*H holds EXACTLY by the divergence theorem regardless of
%  microstructure or kappa, and the same boundary integration recovers G
%  to the pipeline's established accuracy -- so H and G for a case's load
%  direction are, by default, read directly off that case's own data (a
%  robust origin-through least-squares slope over its frames) rather than
%  re-derived from the sampling point file; set USE_TRUE_TENSORS=true
%  below to use the library's exact values instead. This means:
%    * the reference solve does not depend on kappaMat, only on the load
%      case (H,G) and domain size -- so it only needs to be run ONCE PER
%      UNIQUE CaseID, and is automatically reused for every kappa batch
%      that shares that CaseID (cached to disk in homog_reference_cache.
%      mat next to this script, so an interrupted run can resume).
%
%  ------------------------------------------------------------
%  COST WARNING: each unique load case needs one real nonlinear FE solve
%  (RUC_FE_2D_homogeneous.m, NX=NY=1, numeric FD tangent). With ~270+
%  unique cases this can take a long time serially. See CONFIG below for
%  a quick-test subset (SELECTED_CASES), a coarser reference mesh
%  (REF_MESHFRAC_OVERRIDE), and USE_PARFOR/CPU_FRACTION_TARGET to
%  parallelize across a fraction of your CPU cores (needs Parallel
%  Computing Toolbox).
% ============================================================

clear; clc;

%% ================= CONFIG =================
DATA_DIR      = fileparts(mfilename('fullpath'));   % folder this script lives in
FILE_PATTERN  = '*_FINAL.csv';                       % matches NC_k####_FINAL.csv
CACHE_FILE    = fullfile(DATA_DIR, 'homog_reference_cache.mat');

% MAX_CASES     = inf;     % (superseded by SELECTED_CASES below)
SELECTED_CASES = [1:4];      % [] = ALL CaseIDs; e.g. [0 5 10] or 0:10 to run a subset first

% ---- Reference mesh topology ----
% 'inclusion' matches Abaqus's own geometry construction (circle
% partitioned into the domain regardless of kappaMat) -- this is the
% correct choice for a case-vs-reference comparison and is baked into
% refKey below so old (mismatched) cache entries auto-invalidate. Only
% change this to 'regular' or 'auto' if you explicitly want the old
% fast-but-mesh-mismatched behavior for a quick smoke test; both print a
% loud warning from RUC_FE_2D_homogeneous.m when used.
REF_MESH_TYPE = 'inclusion';

% ---- True prescribed (H,G) instead of a least-squares fit ----
% Opt-in; requires Python + pann_loadcases.py + its sample_points_d##.txt
% to be reachable. Off by default so this script stays self-contained
% (only needs the *_FINAL.csv files) unless you turn it on.
USE_TRUE_TENSORS = false;
LIBPY_PATH       = fullfile(DATA_DIR, 'pann_loadcases.py');  % adjust if it lives elsewhere
PYTHON_EXE       = 'python';

% ---- Parallelization ----
% Each MISSING homogeneous reference is an independent nonlinear FE
% solve -- parfor hands one to each worker. Set USE_PARFOR = true to
% parallelize; CPU_FRACTION_TARGET controls how much of the machine to
% use (0.80 = 80% of detected CPU cores), so the machine stays usable
% for other work while this runs. The pool is sized/started once, below,
% by setupParallelPool() -- you do not need to call parpool yourself.
USE_PARFOR          = true;
CPU_FRACTION_TARGET = 0.80;   % fraction of feature('numcores') to use as parfor workers

SAVE_EVERY    = 5;       % checkpoint the cache to disk every N new reference solves (serial mode only)
PARFOR_CHECKPOINT_CHUNK = 20;  % parfor mode: process the missing references in chunks of
                                % this size, checkpointing the cache after each chunk. A
                                % single uncaught error used to kill the ENTIRE parfor batch
                                % and discard every other worker's completed-but-uncollected
                                % result; chunking bounds the damage to one chunk's worth of
                                % work, and each iteration is also wrapped in try/catch below
                                % so one bad case (NoConverge, or an inverted-element frame)
                                % is logged and skipped rather than crashing the run.

% Reference FE mesh (kappaMat is ALWAYS forced to 1.0 for the reference,
% regardless of the case's own kappaMat -- that is the whole point of
% this normalization). Rfrac/meshFrac/NFRAMES/Lx_tot/Ly_tot are otherwise
% taken from each case's own CSV metadata columns.
REF_MESHFRAC_OVERRIDE = [];   % [] = use each case's own meshFrac; or set a fixed value, e.g. 0.04, for a faster/coarser reference mesh
REF_VERBOSE   = true;

% Plotting
SHOW_LEGEND        = true;     % the reference images show a full per-case legend
MAX_LEGEND_ENTRIES = 40;       % legend is skipped (colorbar only) above this many series, regardless of SHOW_LEGEND
COLORMAP_NAME       = 'parula';
SAVE_FIGS           = false;   % true = also export PNG/FIG next to this script

% ---- What to plot ----
PLOT_RAW          = true;      % STEP 1: raw case data, NOT divided by the homogeneous reference
PLOT_RATIO        = true;      % STEP 2: Case / Homogenized
OVERLAY_REFERENCE = true;      % raw plots: also draw the homogeneous reference dashed (same colour)
PLOT_FAMILIES     = {'P','Q','Energy'};   % any of 'F','G','P','Q','Energy'

% ---- Two-group split (each figure is drawn once per group) ----
% GROUP_MODE = 'roughness' : group 1 = smooth curves, group 2 = curves with
%                            abrupt changes/spikes along the load path
%              'Gnorm'     : group 1 = small strain gradient ||G||,
%                            group 2 = large strain gradient ||G||
% GROUP_THRESHOLD = []  -> chosen automatically (largest gap in log10 of the
%                          metric); or set a number to fix the split by hand.
%                          The metric of every case is printed and written to
%                          case_groups.csv so you can pick a value.
GROUP_MODE       = 'roughness';
GROUP_THRESHOLD  = [];
GROUP_COMPONENTS = {'P11','P12','P21','P22', ...
    'Q111','Q112','Q121','Q122','Q211','Q212','Q221','Q222','W_mean'};   % used by 'roughness'

fprintf(['Reference mesh: MeshType=''%s'' (matches Abaqus''s own circle-partitioned ' ...
    'geometry). refKey now versioned (''v2_%s_...'') so any cache entries computed ' ...
    'under the old mesh-mismatched behavior are recomputed automatically.\n'], ...
    REF_MESH_TYPE, REF_MESH_TYPE);
if USE_TRUE_TENSORS
    fprintf('True prescribed (H,G) via pann_loadcases.py ENABLED (LIBPY_PATH=%s).\n', LIBPY_PATH);
end

%% ================= DISCOVER FILES =================
files = dir(fullfile(DATA_DIR, FILE_PATTERN));
files = files(~[files.isdir]);
[~, ord] = sort({files.name});
files = files(ord);
if isempty(files)
    error('compare_case_vs_homogeneous:NoFiles', ...
        'No files matching "%s" found in %s', FILE_PATTERN, DATA_DIR);
end

%% ============================================================
% READ ALL FILES FIRST SO WE CAN SELECT BY CaseID
% ============================================================
nAllFiles = numel(files);
allCaseID = nan(nAllFiles,1);
for i = 1:nAllFiles
    fpath = fullfile(files(i).folder, files(i).name);
    try
        Ttmp = readtable(fpath);
        if ismember('CaseID', Ttmp.Properties.VariableNames)
            allCaseID(i) = Ttmp.CaseID(1);
        end
    catch ME
        warning('Could not read %s: %s', files(i).name, ME.message);
    end
end

%% ============================================================
% CASE SELECTION
%
% SELECTED_CASES = []        -> process ALL CaseIDs
% SELECTED_CASES = [0 5 10]  -> process only CaseID 0, 5 and 10
% SELECTED_CASES = 0:10      -> process CaseID 0 through 10
% ============================================================
if isempty(SELECTED_CASES)
    keepMask = true(nAllFiles,1);
    fprintf('Case selection: ALL CASES\n');
else
    keepMask = ismember(allCaseID, SELECTED_CASES);
    fprintf('Case selection: ');
    fprintf('%d ', SELECTED_CASES);
    fprintf('\n');
end

files = files(keepMask);
selectedCaseID = allCaseID(keepMask);
if isempty(files)
    error('compare_case_vs_homogeneous:NoSelectedFiles', ...
        'No files found for the selected CaseID(s).');
end

[~, ord] = sort({files.name});
files = files(ord);
selectedCaseID = selectedCaseID(ord);

fprintf('\nFound %d file(s) matching "%s"\n', numel(files), FILE_PATTERN);
fprintf('Selected CaseIDs:\n');
uniqueSelectedCases = unique(selectedCaseID(~isnan(selectedCaseID)));
fprintf('%d ', uniqueSelectedCases);
fprintf('\n');

%% ================= PASS 1: read metadata + fit (H,G) per file =================
nF = numel(files);
caseMeta(nF) = struct('file','','T',[],'CaseID',NaN,'kappaMat',NaN,'LoadCase','', ...
    'Lx_tot',NaN,'Ly_tot',NaN,'Rfrac',NaN,'meshFrac',NaN,'NFRAMES',NaN, ...
    'H',zeros(2,2),'G',zeros(2,2,2),'refKey','');

nTrueTensorsUsed = 0;
for i = 1:nF
    fpath = fullfile(files(i).folder, files(i).name);
    Tc = readtable(fpath);
    if isempty(Tc)
        warning('Skipping empty file: %s', files(i).name);
        continue;
    end

    kappaMat = Tc.kappaMat(1);
    CaseID   = Tc.CaseID(1);
    LoadCase = string(Tc.LoadCase(1));
    Lx_tot   = Tc.Lx_tot(1);
    Ly_tot   = Tc.Ly_tot(1);
    Rfrac    = Tc.Rfrac(1);
    meshFrac = Tc.meshFrac(1);
    NFRAMES  = Tc.NFRAMES(1);

    gotTrue = false;
    if USE_TRUE_TENSORS
        elongCol = {'Elongation','elongation'};
        elongVal = [];
        for ec = 1:numel(elongCol)
            if ismember(elongCol{ec}, Tc.Properties.VariableNames)
                elongVal = Tc.(elongCol{ec})(1);
                break;
            end
        end
        if ~isempty(elongVal)
            [Htrue, Gtrue, gotTrue] = getTrueTensorsIfAvailable(LIBPY_PATH, ...
                char(LoadCase), elongVal, Lx_tot, Ly_tot, PYTHON_EXE);
        end
        if gotTrue
            H = Htrue; G = Gtrue;
            nTrueTensorsUsed = nTrueTensorsUsed + 1;
        else
            [H, G] = fitLoadCaseFromCSV(Tc);
            warning('compare_case_vs_homogeneous:TrueTensorsUnavailable', ...
                ['Could not get the true prescribed (H,G) for %s from pann_loadcases.py ' ...
                 '-- falling back to the least-squares fit off this case''s own F/G history.'], ...
                files(i).name);
        end
    else
        [H, G] = fitLoadCaseFromCSV(Tc);
    end

    if ~isempty(REF_MESHFRAC_OVERRIDE)
        refMeshFrac = REF_MESHFRAC_OVERRIDE;
    else
        refMeshFrac = meshFrac;
    end

    % 'v2_<meshtype>_' prefix: any refKey computed under a different
    % REF_MESH_TYPE (including the old pre-fix implicit 'auto'->regular
    % behavior, which never had this prefix at all) is a different
    % string, so isKey(refCache, ...) below correctly reports a cache
    % miss and recomputes it -- no manual cache-clearing needed.
    refKey = sprintf('v2_%s_case%d_Lx%.6f_Ly%.6f_Rf%.6f_mf%.6f_NF%d', ...
        REF_MESH_TYPE, CaseID, Lx_tot, Ly_tot, Rfrac, refMeshFrac, NFRAMES);

    caseMeta(i).file     = files(i).name;
    caseMeta(i).T        = Tc;
    caseMeta(i).CaseID   = CaseID;
    caseMeta(i).kappaMat = kappaMat;
    caseMeta(i).LoadCase = LoadCase;
    caseMeta(i).Lx_tot   = Lx_tot;
    caseMeta(i).Ly_tot   = Ly_tot;
    caseMeta(i).Rfrac    = Rfrac;
    caseMeta(i).meshFrac = refMeshFrac;
    caseMeta(i).NFRAMES  = NFRAMES;
    caseMeta(i).H        = H;
    caseMeta(i).G        = G;
    caseMeta(i).refKey   = refKey;
end
if USE_TRUE_TENSORS
    fprintf('True prescribed (H,G) used for %d of %d file(s); the rest fell back to the fit.\n', ...
        nTrueTensorsUsed, nF);
end

%% ================= PASS 2: build/reuse homogeneous references =================
if isfile(CACHE_FILE)
    S = load(CACHE_FILE, 'refCache');
    refCache = S.refCache;
    fprintf('Loaded %d cached homogeneous reference(s) from %s\n', refCache.Count, CACHE_FILE);
else
    refCache = containers.Map('KeyType','char','ValueType','any');
end

% unique jobs still needed
[uniqKeys, firstIdx] = unique({caseMeta.refKey}, 'stable');
needIdx = firstIdx(~isKey(refCache, uniqKeys));
fprintf('%d unique load case(s); %d homogeneous reference(s) already cached, %d to compute.\n', ...
    numel(uniqKeys), numel(uniqKeys)-numel(needIdx), numel(needIdx));

if USE_PARFOR && numel(needIdx) > 1
    setupParallelPool(CPU_FRACTION_TARGET);
end

failedRefs = struct('key', {}, 'file', {}, 'message', {});
FAILED_LOG_FILE = fullfile(DATA_DIR, 'FAILED_HOMOG_REFS.txt');
nNeed = numel(needIdx);

if USE_PARFOR && nNeed > 1
    % Chunked so a hard failure inside one worker's iteration only costs
    % that chunk's worth of unsaved work: previously, ONE uncaught error
    % anywhere in the parfor block terminated the whole batch and
    % discarded every other worker's completed-but-uncollected result
    % (this is exactly what happened in the 273-case run -- a single
    % NoConverge threw away all in-flight progress beyond the last
    % checkpoint). Each iteration is also wrapped in its own try/catch so
    % one bad case (NoConverge after exhausting maxSub, or the new
    % Jmin<=0 inverted-element check) is logged and skipped instead of
    % crashing anything.
    chunkSize = max(1, PARFOR_CHECKPOINT_CHUNK);
    nChunks = ceil(nNeed / chunkSize);
    for ch = 1:nChunks
        cStart = (ch-1)*chunkSize + 1;
        cEnd   = min(nNeed, ch*chunkSize);
        chunkIdx = needIdx(cStart:cEnd);
        nInChunk = numel(chunkIdx);
        chunkResults = cell(nInChunk,1);
        parfor jj = 1:nInChunk
            idx = chunkIdx(jj);
            cm = caseMeta(idx); %#ok<PFBNS>
            try
                [Tref, diagOut] = RUC_FE_2D_homogeneous(cm.H, cm.G, ...
                    'Lx_tot', cm.Lx_tot, 'Ly_tot', cm.Ly_tot, 'Rfrac', cm.Rfrac, ...
                    'meshFrac', cm.meshFrac, 'kappaMat', 1.0, 'NX', 1, 'NY', 1, ...
                    'NFRAMES', cm.NFRAMES, 'Verbose', REF_VERBOSE, ...
                    'MeshType', REF_MESH_TYPE); %#ok<PFBNS>
                chunkResults{jj} = struct('ok', true, 'key', cm.refKey, 'file', cm.file, ...
                    'T', Tref, 'H', cm.H, 'G', cm.G, 'diag', diagOut, 'message', '');
                fprintf('  [parfor] %s  (%d nodes, %d elem, %.1fs)\n', ...
                    cm.refKey, diagOut.Nnodes, diagOut.Nelem, diagOut.solveSeconds);
            catch ME
                chunkResults{jj} = struct('ok', false, 'key', cm.refKey, 'file', cm.file, ...
                    'T', [], 'H', cm.H, 'G', cm.G, 'diag', [], 'message', ME.message);
                fprintf('  [parfor] %s  FAILED: %s\n', cm.refKey, ME.message);
            end
        end
        for jj = 1:nInChunk
            r = chunkResults{jj};
            if r.ok
                refCache(r.key) = struct('T', r.T, 'H', r.H, 'G', r.G); %#ok<NASGU>
            else
                failedRefs(end+1) = struct('key', r.key, 'file', r.file, 'message', r.message); %#ok<AGROW>
            end
        end
        save(CACHE_FILE, 'refCache');
        fprintf('Chunk %d/%d done: %d references cached and saved so far (%d failed so far).\n', ...
            ch, nChunks, refCache.Count, numel(failedRefs));
    end
else
    for jj = 1:nNeed
        idx = needIdx(jj);
        cm = caseMeta(idx);
        tSolve = tic;
        try
            [Tref, diagOut] = RUC_FE_2D_homogeneous(cm.H, cm.G, ...
                'Lx_tot', cm.Lx_tot, 'Ly_tot', cm.Ly_tot, 'Rfrac', cm.Rfrac, ...
                'meshFrac', cm.meshFrac, 'kappaMat', 1.0, 'NX', 1, 'NY', 1, ...
                'NFRAMES', cm.NFRAMES, 'Verbose', REF_VERBOSE, ...
                'MeshType', REF_MESH_TYPE);
            fprintf('[%d/%d] %s  (%s) -> %d nodes, %d elem, %.1fs\n', ...
                jj, nNeed, cm.refKey, cm.file, diagOut.Nnodes, diagOut.Nelem, toc(tSolve));
            refCache(cm.refKey) = struct('T', Tref, 'H', cm.H, 'G', cm.G); %#ok<NASGU>
        catch ME
            fprintf('[%d/%d] %s  (%s) -> FAILED: %s\n', jj, nNeed, cm.refKey, cm.file, ME.message);
            failedRefs(end+1) = struct('key', cm.refKey, 'file', cm.file, 'message', ME.message); %#ok<AGROW>
        end
        if mod(jj, SAVE_EVERY) == 0 || jj == nNeed
            save(CACHE_FILE, 'refCache');
            fprintf('  (checkpoint saved: %d references cached, %d failed so far)\n', ...
                refCache.Count, numel(failedRefs));
        end
    end
end

if ~isempty(failedRefs)
    fid = fopen(FAILED_LOG_FILE, 'w');
    if fid > 0
        fprintf(fid, 'Homogeneous reference solves that FAILED (case skipped in Pass 3):\n\n');
        for k = 1:numel(failedRefs)
            fprintf(fid, '%s  (from %s)\n    %s\n\n', failedRefs(k).key, failedRefs(k).file, failedRefs(k).message);
        end
        fclose(fid);
    end
    fprintf(['\n%d homogeneous reference(s) FAILED and were skipped -- see %s.\n' ...
        'Every OTHER reference was still solved and checkpointed to %s; rerun this script ' ...
        'to retry just the failed ones (they stay uncached, so they are picked up again ' ...
        'next time without recomputing anything that already succeeded).\n\n'], ...
        numel(failedRefs), FAILED_LOG_FILE, CACHE_FILE);
end

%% ================= PASS 3: raw + ratio per case =================
% Every quantity is kept THREE ways per case so the raw data can be
% inspected before any normalisation:
%   r.caseVal.<name>  raw case (Abaqus) series          -- no division
%   r.refVal.<name>   homogeneous MATLAB reference      -- what we divide by
%   r.ratio.<name>    caseVal ./ refVal                 -- "Case / Homogenized"
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
for ff = 1:numel(FAMILIES)
    allComp = [allComp, FAMILIES{ff}.comps]; %#ok<AGROW>
end
gComps = FAMILIES{2}.comps;

results = struct([]);
for i = 1:nF
    cm = caseMeta(i);
    if isempty(cm.T) || ~isKey(refCache, cm.refKey)
        warning('No cached reference for %s (key %s) -- skipping.', cm.file, cm.refKey);
        continue;
    end
    ref = refCache(cm.refKey);
    Tref = ref.T;
    Tc   = cm.T;

    if height(Tc) ~= height(Tref)
        warning('%s: frame count mismatch (case %d vs reference %d) -- truncating to the shorter.', ...
            cm.file, height(Tc), height(Tref));
    end
    n = min(height(Tc), height(Tref));

    r = struct();
    r.file     = cm.file;
    r.CaseID   = cm.CaseID;
    r.kappaMat = cm.kappaMat;
    r.LoadCase = cm.LoadCase;
    r.StepTime = Tc.StepTime(1:n);
    r.LoadStep = (1:n)';
    r.caseVal  = struct();
    r.refVal   = struct();
    r.ratio    = struct();

    for c = 1:numel(allComp)
        name = allComp{c};
        if ~ismember(name, Tc.Properties.VariableNames) || ...
           ~ismember(name, Tref.Properties.VariableNames)
            continue; % e.g. ALLSE_mean may be absent from some datasets
        end
        num = Tc.(name)(1:n);
        den = Tref.(name)(1:n);
        r.caseVal.(name) = num;
        r.refVal.(name)  = den;
        r.ratio.(name)   = num ./ den;
    end

    % strain-gradient magnitude of this load case (G is linear in the
    % load step, so the max over frames is the final-frame value)
    g2 = zeros(height(Tc),1);
    for c = 1:numel(gComps)
        if ismember(gComps{c}, Tc.Properties.VariableNames)
            g2 = g2 + Tc.(gComps{c}).^2;
        end
    end
    r.Gnorm = max(sqrt(g2));

    if isempty(results), results = r; else, results(end+1) = r; end %#ok<AGROW>
end
fprintf('Computed raw + ratio series for %d of %d case file(s).\n', numel(results), nF);
if isempty(results)
    error('compare_case_vs_homogeneous:NoResults', ...
        'No case produced a ratio (check the warnings above) -- nothing to plot.');
end

%% ================= SPLIT INTO TWO GROUPS =================
% Group 1 = "quiet" cases, group 2 = cases with abrupt changes.
%   GROUP_MODE = 'roughness' : how non-smooth a case's Case/Homogenized
%                curves are along the load path. Metric = largest absolute
%                second difference of the ratio over GROUP_COMPONENTS
%                (a straight line or a gentle curve gives ~0; a spike or a
%                kink gives a large value).
%   GROUP_MODE = 'Gnorm'     : magnitude of the imposed strain gradient
%                ||G|| (large-bending cases vs small-bending cases).
nR = numel(results);
switch GROUP_MODE
    case 'roughness'
        metric = nan(nR,1);
        for i = 1:nR
            metric(i) = curveRoughness(results(i), GROUP_COMPONENTS);
        end
        metricLabel = 'max |2nd difference| of Case/Homogenized';
        groupNames  = {'Smooth curves', 'Abrupt changes'};
    case 'Gnorm'
        metric = [results.Gnorm]';
        metricLabel = 'max ||G||';
        groupNames  = {'Small strain gradient', 'Large strain gradient'};
    otherwise
        error('compare_case_vs_homogeneous:BadGroupMode', ...
            'GROUP_MODE must be ''roughness'' or ''Gnorm''.');
end

if isempty(GROUP_THRESHOLD)
    thr = twoMeansThreshold(metric);
    thrSrc = 'automatic (largest gap, 2-class split of log10 metric)';
else
    thr = GROUP_THRESHOLD;
    thrSrc = 'user-defined';
end
isRough  = ~(metric <= thr);          % NaN / Inf metrics land in group 2
groupOf  = 1 + isRough;
groupIdx = {find(groupOf == 1), find(groupOf == 2)};

fprintf('\nGrouping by %s\n  threshold = %.4g  (%s)\n', metricLabel, thr, thrSrc);
[~, ordM] = sort(metric, 'descend');
fprintf('  %-8s %-14s %s\n', 'CaseID', 'metric', 'group');
for k = 1:nR
    q = ordM(k);
    fprintf('  %-8d %-14.4g %d  (%s)\n', results(q).CaseID, metric(q), groupOf(q), groupNames{groupOf(q)});
end
fprintf('  -> group 1 "%s": %d case(s);  group 2 "%s": %d case(s)\n\n', ...
    groupNames{1}, numel(groupIdx{1}), groupNames{2}, numel(groupIdx{2}));

try
    Tgroups = table([results.CaseID]', metric, groupOf, string(groupNames(groupOf))', ...
        'VariableNames', {'CaseID','Metric','Group','GroupName'});
    writetable(Tgroups, fullfile(DATA_DIR, 'case_groups.csv'));
    fprintf('Group assignment written to %s\n', fullfile(DATA_DIR, 'case_groups.csv'));
catch ME
    warning('Could not write case_groups.csv (%s).', ME.message);
end

%% ================= PLOTTING =================
% STEP 1: raw case data (no division); STEP 2: divided by the homogeneous
% reference. Each is drawn once per group so the smooth cases are not
% squashed by the spiky ones.
modes = {};
if PLOT_RAW,   modes{end+1} = 'raw';   end %#ok<UNRCH>
if PLOT_RATIO, modes{end+1} = 'ratio'; end
for mm = 1:numel(modes)
    for ff = 1:numel(FAMILIES)
        fam = FAMILIES{ff};
        if ~ismember(fam.tag, PLOT_FAMILIES), continue; end
        for gg = 1:2
            idx = groupIdx{gg};
            if isempty(idx), continue; end
            groupLabel = sprintf('%s (%d case%s)', groupNames{gg}, numel(idx), ...
                repmat('s', 1, double(numel(idx) ~= 1)));
            plotGroupGrid(results, idx, fam, modes{mm}, groupLabel, sprintf('group%d', gg), ...
                OVERLAY_REFERENCE, COLORMAP_NAME, SHOW_LEGEND, MAX_LEGEND_ENTRIES, SAVE_FIGS, DATA_DIR);
        end
    end
end

fprintf('Done.\n');


% ============================================================
% LOCAL FUNCTIONS
% ============================================================

function setupParallelPool(cpuFraction)
% Starts (or resizes) the local parallel pool to use cpuFraction of the
% machine's detected CPU cores, so the parfor loop above uses only that
% share of the machine and the rest stays usable for other work.
% Called once, right before the parfor block. If the Parallel Computing
% Toolbox isn't available, or the pool fails to start for any reason,
% this warns and returns -- parfor still runs, just serially, so the
% script never errors out over this.
    if isempty(ver('parallel'))
        warning(['Parallel Computing Toolbox not found -- USE_PARFOR will run ' ...
                 'serially (no error, just no speedup). Install/enable it to parallelize.']);
        return;
    end
    try
        nCores = feature('numcores');   % physical cores detected on this machine
        nWorkers = max(1, floor(cpuFraction * nCores));
        pool = gcp('nocreate');
        if ~isempty(pool)
            if pool.NumWorkers == nWorkers
                fprintf('Reusing existing parallel pool (%d workers).\n', pool.NumWorkers);
                return;
            end
            fprintf('Closing existing pool (%d workers) to resize to %d...\n', pool.NumWorkers, nWorkers);
            delete(pool);
        end
        fprintf('Starting parallel pool: %d workers (%.0f%% of %d detected CPU cores)\n', ...
            nWorkers, cpuFraction*100, nCores);
        parpool('local', nWorkers);
    catch ME
        warning('Could not start a parallel pool (%s) -- continuing serially.', ME.message);
    end
end

function [H, G] = fitLoadCaseFromCSV(Tc)
% Recover the KUBC load case (H,G) from a case's own reported boundary
% F11..F22 / G111..G222 history. Fbar(t) = I + t*H and Gbar(t) = t*G hold
% (to this pipeline's established accuracy) at every frame regardless of
% kappaMat -- see the header of this script -- so a simple origin-
% through least-squares slope over all frames is robust to any single
% frame's small numerical noise (KubcErrF/KubcErrG in the CSV bound it).
% Used unless USE_TRUE_TENSORS successfully retrieves the library's own
% exact (H,G) instead (see getTrueTensorsIfAvailable below).
    t = Tc.StepTime;
    denom = sum(t.^2);

    H = zeros(2,2);
    H(1,1) = sum(t.*(Tc.F11-1)) / denom;
    H(1,2) = sum(t.*(Tc.F12))   / denom;
    H(2,1) = sum(t.*(Tc.F21))   / denom;
    H(2,2) = sum(t.*(Tc.F22-1)) / denom;

    G = zeros(2,2,2);
    G(1,1,1) = sum(t.*Tc.G111) / denom;
    G(1,1,2) = sum(t.*Tc.G112) / denom;
    G(1,2,1) = sum(t.*Tc.G121) / denom;
    G(1,2,2) = sum(t.*Tc.G122) / denom;
    G(2,1,1) = sum(t.*Tc.G211) / denom;
    G(2,1,2) = sum(t.*Tc.G212) / denom;
    G(2,2,1) = sum(t.*Tc.G221) / denom;
    G(2,2,2) = sum(t.*Tc.G222) / denom;
end

function [Hout, Gout, ok] = getTrueTensorsIfAvailable(libPy, loadCaseName, elongation, Lx_tot, Ly_tot, pythonExe)
% Prescribed KUBC coefficients for one load case, straight from the
% library (adapted from Run_Hyperelastic_Ensemble*.m's getCaseTensors).
% Calls `pann_loadcases.py --tensors`, i.e. the same buildLoadCase path
% Main2D_hyperelastic.py used to build the case's actual ExpressionField,
% so this cannot drift from what was truly imposed. Returns ok=false (and
% zero tensors) on ANY failure -- missing file, python not found, output
% not parseable -- so the caller can fall back to fitLoadCaseFromCSV
% without this ever aborting the run.
    Hout = zeros(2,2); Gout = zeros(2,2,2); ok = false;

    if isempty(libPy) || ~isfile(libPy)
        return;
    end

    args = sprintf('"%s" --tensors %s %.17g %.17g %.17g', ...
        libPy, loadCaseName, elongation, Lx_tot, Ly_tot);

    [st, out] = system(sprintf('%s %s', pythonExe, args));
    if st ~= 0
        [st, out] = system(sprintf('abaqus python %s', args));
    end
    if st ~= 0
        return;
    end

    hTok = regexp(out, 'HVEC([^\n\r]*)', 'tokens', 'once');
    gTok = regexp(out, 'GVEC([^\n\r]*)', 'tokens', 'once');
    if isempty(hTok) || isempty(gTok)
        return;
    end

    h = sscanf(hTok{1}, '%f').';
    g = sscanf(gTok{1}, '%f').';
    if numel(h) ~= 4 || numel(g) ~= 8
        return;
    end

    Hout = [h(1) h(2); h(3) h(4)];
    Gout = zeros(2,2,2);
    Gout(1,1,1)=g(1); Gout(1,1,2)=g(2); Gout(1,2,1)=g(3); Gout(1,2,2)=g(4);
    Gout(2,1,1)=g(5); Gout(2,1,2)=g(6); Gout(2,2,1)=g(7); Gout(2,2,2)=g(8);
    ok = true;
end

function m = curveRoughness(r, comps)
% Largest absolute second difference of the Case/Homogenized ratio over the
% listed components. ~0 for a straight line or a gentle smooth curve, large
% for a spike or a kink. A series with fewer than 3 finite points counts
% as maximally rough (Inf).
    m = 0;
    for c = 1:numel(comps)
        if ~isfield(r.ratio, comps{c}), continue; end
        y = r.ratio.(comps{c});
        y = y(isfinite(y));
        if numel(y) < 3
            m = Inf;
            continue;
        end
        m = max(m, max(abs(diff(y, 2))));
    end
end

function thr = twoMeansThreshold(metric)
% Automatic 2-class split of a positive metric on a log10 scale (1-D
% two-means / Jenks): pick the split that minimises the within-class sum of
% squares. Returns Inf (everything in group 1) if fewer than 2 distinct
% finite values exist. Inf/NaN metrics are ignored here and fall into the
% "abrupt" group by the caller.
    v = metric(isfinite(metric));
    v = sort(log10(max(v, realmin)));
    n = numel(v);
    if n < 2 || (v(end) - v(1)) < eps
        thr = Inf;
        return;
    end
    best = Inf; kBest = 1;
    for k = 1:n-1
        a = v(1:k); b = v(k+1:end);
        sse = sum((a - mean(a)).^2) + sum((b - mean(b)).^2);
        if sse < best
            best = sse; kBest = k;
        end
    end
    thr = 10^(0.5*(v(kBest) + v(kBest+1)));
end

function plotGroupGrid(results, idx, fam, mode, groupLabel, gTag, overlayRef, cmapName, ...
        showLegend, maxLegendEntries, saveFigs, outDir)
% One figure = one quantity family (P, Q, ...) for ONE group of cases.
%   mode 'raw'   : raw case values (solid, with markers); the homogeneous
%                  reference is overlaid dashed in the same colour when
%                  overlayRef is true.
%   mode 'ratio' : Case / Homogenized, reference line at 1.
    isRaw = strcmp(mode, 'raw');
    if isRaw, field = 'caseVal'; else, field = 'ratio'; end

    comps = fam.comps;
    valid = false(size(comps));
    for c = 1:numel(comps)
        for s = idx(:)'
            if isfield(results(s).(field), comps{c})
                valid(c) = true; break;
            end
        end
    end
    comps = comps(valid);
    if isempty(comps), return; end

    nSeries = numel(idx);
    try
        cmap = feval(cmapName, max(nSeries, 2));
    catch
        cmap = parula(max(nSeries, 2));
    end
    doLegend = showLegend && nSeries <= maxLegendEntries;

    if isRaw
        if overlayRef
            ttl = sprintf('%s - raw data, solid = case, dashed = homogeneous ref  |  %s', fam.title, groupLabel);
        else
            ttl = sprintf('%s - raw case data  |  %s', fam.title, groupLabel);
        end
    else
        ttl = sprintf('Normalized %s (Case / Homogenized)  |  %s', fam.title, groupLabel);
    end

    fig = figure('Name', ttl, 'Color', 'w');
    sgtitle(ttl, 'Interpreter', 'none');

    for c = 1:numel(comps)
        name = comps{c};
        ax = subplot(fam.grid(1), fam.grid(2), c); hold(ax, 'on'); grid(ax, 'on');
        legendHandles = gobjects(0);
        legendEntries = strings(0);
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
            legendHandles(end+1) = h; %#ok<AGROW>
            legendEntries(end+1) = sprintf('Case %d', r.CaseID); %#ok<AGROW>
        end
        title(ax, name, 'Interpreter', 'none');
        xlabel(ax, 'Load Step');
        if isRaw
            ylabel(ax, name, 'Interpreter', 'none');
        else
            yline(ax, 1, '--k');
            ylabel(ax, 'Case / Homogenized');
        end
        if doLegend && c == 1 && ~isempty(legendHandles)
            legend(ax, legendHandles, legendEntries, 'Location', 'eastoutside', 'FontSize', 6);
        end
    end

    if saveFigs
        base = fullfile(outDir, sprintf('%s_%s_%s', mode, fam.tag, gTag));
        try
            exportgraphics(fig, [base '.png'], 'Resolution', 300);
        catch
        end
        savefig(fig, [base '.fig']);
    end
end