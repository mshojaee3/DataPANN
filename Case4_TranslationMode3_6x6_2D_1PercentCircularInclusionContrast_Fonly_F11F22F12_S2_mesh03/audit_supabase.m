% ============================================================
%  audit_supabase.m
%
%  Standalone. Put it next to the driver and run it on its own.
%  It needs no other file.
%
%  WHAT IT DOES
%  ------------
%  Lists every NC_k####_FINAL.csv in the training folder, reads
%  every row of the Supabase case_runs table, and compares them:
%
%    STALE   row in Supabase, no CSV on disk
%            -> the case is claimed but produced nothing (crashed
%               worker, killed MATLAB, deleted output). Its row
%               blocks every machine from retrying it, so the row
%               is deleted and the case becomes runnable again.
%
%    ORPHAN  CSV on disk, no row in Supabase
%            -> the data exists but nothing marks it done, so some
%               machine will eventually redo it. Reported, and
%               optionally recorded as completed.
%
%    OK      row and CSV agree. Left alone.
%
%  READ THIS BEFORE RUNNING
%  ------------------------
%  Each machine holds only ITS OWN CSVs. A case completed on AGAVE
%  has no CSV here, so from this machine it looks STALE. Deleting
%  that row would make this machine redo work that is already
%  finished elsewhere.
%
%  ONLY_THIS_COMPUTER = true (the default) therefore restricts the
%  deletion to rows whose computer_id is this machine. That is
%  always safe. Set it to false ONLY when the training folder holds
%  the CSVs gathered from ALL machines.
%
%  DRY_RUN = true is also the default: it shows what would happen
%  and changes nothing. Read the report, then set DRY_RUN = false.
% ============================================================

clear; clc;

% ============================================================
% CONFIGURATION
% ============================================================
SUPABASE_URL = 'https://mreweybhnsyytsfzxxfj.supabase.co';
SUPABASE_KEY = 'sb_publishable_ZtKwPWDEyxQyBy05Kz08oQ_F_agHBZO';
TABLE        = 'case_runs';

% Folder holding the finished CSVs. Empty = TrainingData next to
% this script. To audit the gathered set from all machines, point
% it at that folder and set ONLY_THIS_COMPUTER = false.
DATA_DIR = '';

% Case-ID range of the sweep (zero based, inclusive).
CASE_ID_MIN = 0;
CASE_ID_MAX = 273;

% Safety switches. Both start in the cautious position.
DRY_RUN            = false;   % false actually deletes rows
ONLY_THIS_COMPUTER = false;   % see the warning above

% A 'running' row younger than this is assumed to be a LIVE job on
% another machine and is never deleted. Without this guard, widening
% the scope while a sweep is in progress would release cases that
% are running right now, and two machines would compute them.
%
% Set it a bit above your longest single case. Yours top out near
% 760 s, so 2 h is generous; raise it if a machine is slow or
% paused. Set to 0 to disable the guard entirely (only when every
% other MATLAB is definitely closed).
MIN_RUNNING_AGE_HOURS = 2;

% Delete stale rows whose status is 'completed' too? A completed row
% with no CSV means the file was deleted or moved after the run.
INCLUDE_COMPLETED = false;   % false = only 'running'/'failed' rows

% A CSV shorter than this many lines (header + data) is treated as
% missing: an empty or header-only file is not a finished case.
MIN_CSV_LINES = 2;

% Write a row for each ORPHAN so other machines stop redoing it.
RECORD_ORPHANS = false;
% ============================================================


rootDir = fileparts(mfilename('fullpath'));
if isempty(DATA_DIR)
    DATA_DIR = fullfile(rootDir);
end
assert(isfolder(DATA_DIR), 'Folder not found: %s', DATA_DIR);

thisComputer = getenv('COMPUTERNAME');
if isempty(thisComputer); thisComputer = getenv('HOSTNAME'); end

fprintf('==================================================\n');
fprintf('Supabase audit\n');
fprintf('  folder    : %s\n', DATA_DIR);
fprintf('  table     : %s\n', TABLE);
fprintf('  computer  : %s\n', thisComputer);
if DRY_RUN
    fprintf('  mode      : DRY RUN (nothing will be changed)\n');
else
    fprintf(2, '  mode      : LIVE (rows WILL be deleted)\n');
end
if ONLY_THIS_COMPUTER
    fprintf('  scope     : rows started by %s only\n', thisComputer);
else
    fprintf(2, '  scope     : ALL rows, every computer\n');
end
fprintf('==================================================\n\n');


% ------------------------------------------------------------
% 1. which case IDs have a real CSV on disk
% ------------------------------------------------------------
files = dir(fullfile(DATA_DIR, 'NC_k*_FINAL.csv'));
onDisk = [];
tooShort = {};

for k = 1:numel(files)
    tok = regexp(files(k).name, '^NC_k(\d+)_FINAL\.csv$', 'tokens', 'once');
    if isempty(tok); continue; end
    id = str2double(tok{1});

    if countLines(fullfile(DATA_DIR, files(k).name)) < MIN_CSV_LINES
        tooShort{end+1} = files(k).name; %#ok<SAGROW>
        continue;    % treated as missing
    end
    onDisk(end+1) = id; %#ok<SAGROW>
end
onDisk = unique(onDisk);

fprintf('CSV files found : %d\n', numel(files));
fprintf('Valid cases     : %d\n', numel(onDisk));
if ~isempty(tooShort)
    fprintf(2, 'Empty/short CSVs (counted as MISSING): %d\n', numel(tooShort));
    for k = 1:min(10, numel(tooShort))
        fprintf(2, '    %s\n', tooShort{k});
    end
end


% ------------------------------------------------------------
% 2. every row in the table
% ------------------------------------------------------------
url = sprintf(['%s/rest/v1/%s' ...
               '?select=case_id,status,computer_id,started_at,created_at' ...
               '&limit=10000'], ...
              SUPABASE_URL, TABLE);
cmd = sprintf('curl -sS -H "apikey: %s" -H "Authorization: Bearer %s" "%s"', ...
              SUPABASE_KEY, SUPABASE_KEY, url);
[st, resp] = system(cmd);
resp = strtrim(resp);

if st ~= 0 || isempty(resp) || resp(1) ~= '['
    error('Could not read %s:\n%s', TABLE, resp);
end

data = jsondecode(resp);
if ~iscell(data)
    if isempty(data); data = {}; else; data = num2cell(data); end
end

rowID   = zeros(numel(data), 1);
rowStat = strings(numel(data), 1);
rowWho  = strings(numel(data), 1);
rowAge  = nan(numel(data), 1);      % hours since started_at
for k = 1:numel(data)
    r = data{k};
    rowID(k)   = r.case_id;
    rowStat(k) = getf(r, 'status',      'unknown');
    rowWho(k)  = getf(r, 'computer_id', '');
    % started_at is null on rows written by older helpers, but
    % created_at carries a database default of now(), so it always
    % has a value and dates the row just as well.
    rowAge(k)  = ageHours(getf(r, 'started_at', ''));
    if isnan(rowAge(k))
        rowAge(k) = ageHours(getf(r, 'created_at', ''));
    end
end

fprintf('Rows in table   : %d\n\n', numel(rowID));


% ------------------------------------------------------------
% 3. compare
% ------------------------------------------------------------
inRange   = rowID >= CASE_ID_MIN & rowID <= CASE_ID_MAX;
hasCsv    = ismember(rowID, onDisk);
isStale   = inRange & ~hasCsv;

if ~INCLUDE_COMPLETED
    isStale = isStale & rowStat ~= "completed";
end
if ONLY_THIS_COMPUTER
    isStale = isStale & startsWith(rowWho, thisComputer);
end

% Protect jobs that may still be running elsewhere. An unknown age
% (NaN started_at) counts as young, i.e. it is protected.
isLive = false(size(isStale));
if MIN_RUNNING_AGE_HOURS > 0
    isLive = rowStat == "running" & ...
             (isnan(rowAge) | rowAge < MIN_RUNNING_AGE_HOURS) & ...
             ~startsWith(rowWho, thisComputer);
    isStale = isStale & ~isLive;
end

% Rows from THIS computer are never treated as live, because the
% usual reason to run the audit is that a local run died. But a row
% only minutes old may belong to a driver running right now in
% another MATLAB session on this same machine, and releasing it
% would let a second machine start the same case.
ownYoung = isStale & startsWith(rowWho, thisComputer) & ...
           ~isnan(rowAge) & rowAge < MIN_RUNNING_AGE_HOURS;
if any(ownYoung)
    fprintf(2, ['WARNING: %d row(s) from %s are less than %g h old.\n' ...
                'If a driver is running on this machine right now, ' ...
                'those cases are IN PROGRESS.\n' ...
                'Close it first, or remove them from the list.\n\n'], ...
            nnz(ownYoung), thisComputer, MIN_RUNNING_AGE_HOURS);
end

orphans = setdiff(onDisk(:), rowID(:));
noRun   = setdiff((CASE_ID_MIN:CASE_ID_MAX)', unique(rowID(inRange)));

fprintf('--- STALE (row but no CSV -> will be released) ---\n');
if ~any(isStale)
    fprintf('  none\n');
else
    idxs = find(isStale);
    for k = 1:numel(idxs)
        i = idxs(k);
        fprintf('  caseID %4d  %-10s %-22s %s\n', ...
                rowID(i), rowStat(i), rowWho(i), ageStr(rowAge(i)));
    end
end

if any(isLive)
    fprintf('\n--- LIVE (running elsewhere, younger than %g h) ---\n', ...
            MIN_RUNNING_AGE_HOURS);
    idxs = find(isLive);
    for k = 1:min(15, numel(idxs))
        i = idxs(k);
        fprintf('  caseID %4d  %-22s %s\n', ...
                rowID(i), rowWho(i), ageStr(rowAge(i)));
    end
    if numel(idxs) > 15
        fprintf('    ... and %d more\n', numel(idxs) - 15);
    end
    fprintf('  Left alone in case that machine is still working.\n');
end

fprintf('\n--- ORPHAN (CSV but no row) ---\n');
if isempty(orphans)
    fprintf('  none\n');
else
    fprintf('  %s\n', strjoin(string(orphans(:).'), ', '));
end

fprintf('\n--- NEVER STARTED (no row, no CSV here) ---\n');
fprintf('  %d of %d cases: %s\n', numel(noRun), ...
        CASE_ID_MAX - CASE_ID_MIN + 1, shortList(noRun));

% Rows skipped only because of the scope switches, so you can see
% what a wider run would catch.
hidden = inRange & ~hasCsv & ~isStale;
if any(hidden)
    fprintf('\n--- NOT TOUCHED (excluded by the switches above) ---\n');
    fprintf('  %d rows have no CSV here but were left alone:\n', nnz(hidden));
    idxs = find(hidden);
    for k = 1:min(15, numel(idxs))
        i = idxs(k);
        fprintf('    caseID %4d  %-10s %-22s %s\n', ...
                rowID(i), rowStat(i), rowWho(i), ageStr(rowAge(i)));
    end
    if numel(idxs) > 15
        fprintf('    ... and %d more\n', numel(idxs) - 15);
    end
    fprintf('  Most are cases finished on another machine. Widen the\n');
    fprintf('  scope only if this folder holds every machine''s CSVs.\n');
end


% ------------------------------------------------------------
% 4. act
% ------------------------------------------------------------
fprintf('\n==================================================\n');

nStale = nnz(isStale);
if nStale == 0
    fprintf('Nothing to release.\n');
elseif DRY_RUN
    fprintf('DRY RUN: %d row(s) would be deleted.\n', nStale);
    fprintf('Set DRY_RUN = false and run again to apply.\n');
else
    fprintf('Deleting %d row(s)...\n', nStale);
    nOK = 0;
    idxs = find(isStale);
    for k = 1:numel(idxs)
        id = rowID(idxs(k));
        [okDel, nDel, why] = deleteCase(SUPABASE_URL, SUPABASE_KEY, TABLE, id);
        if okDel && nDel > 0
            nOK = nOK + 1;
            fprintf('  released caseID %d\n', id);
        elseif okDel && nDel == 0
            % HTTP said fine but no row came back. Row Level Security
            % hides the row from DELETE, so Postgres deleted nothing
            % and still answered 204. Trusting the status code here is
            % exactly how 25 "successful" deletes changed nothing.
            fprintf(2, '  caseID %d: request accepted but NO ROW DELETED\n', id);
        else
            fprintf(2, '  FAILED to release caseID %d (%s)\n', id, why);
        end
    end
    fprintf('Released %d of %d.\n', nOK, nStale);
    if nOK == 0
        fprintf(2, ['\nNothing was actually deleted: Row Level Security is\n' ...
                    'blocking DELETE for the anon key. Either add the policy\n' ...
                    'once,\n\n' ...
                    '  create policy case_runs_delete on public.case_runs\n' ...
                    '    for delete using (true);\n' ...
                    '  notify pgrst, ''reload schema'';\n\n' ...
                    'or just paste this into the SQL editor, which runs as a\n' ...
                    'privileged role and ignores RLS:\n\n']);
        fprintf(2, 'delete from public.case_runs\nwhere case_id in (%s);\n\n', ...
                strjoin(string(rowID(isStale)).', ','));
    elseif nOK < nStale
        fprintf(2, '\n%d row(s) were not removed. See the messages above.\n', ...
                nStale - nOK);
    else
        fprintf('They can now be picked up again.\n');
    end
end

if RECORD_ORPHANS && ~isempty(orphans)
    if DRY_RUN
        fprintf('DRY RUN: %d orphan(s) would be recorded as completed.\n', ...
                numel(orphans));
    else
        fprintf('\nRecording %d orphan(s) as completed...\n', numel(orphans));
        for k = 1:numel(orphans)
            payload = struct('case_id', orphans(k), ...
                             'computer_id', thisComputer, ...
                             'status', 'completed');
            if httpJson('POST', SUPABASE_URL, SUPABASE_KEY, TABLE, '', payload)
                fprintf('  recorded caseID %d\n', orphans(k));
            else
                fprintf(2, '  FAILED to record caseID %d\n', orphans(k));
            end
        end
    end
end
fprintf('==================================================\n');


% ============================================================
% LOCAL FUNCTIONS
% ============================================================
function n = countLines(f)
n = 0;
fid = fopen(f, 'r');
if fid < 0; return; end
while ischar(fgetl(fid)); n = n + 1; end
fclose(fid);
end


function h = ageHours(tsIso)
% Hours since an ISO timestamp; NaN when absent or unparseable.
h = NaN;
tsIso = char(tsIso);
if isempty(tsIso); return; end

% Supabase returns e.g. 2026-08-13T03:28:59.656973+00:00 and also
% plain ...T03:28:59Z. Strip the fraction and the offset, then read
% the rest as UTC.
core = regexprep(tsIso, '\.\d+', '');
core = regexprep(core, '([Zz]|[+-]\d{2}:?\d{2})$', '');
try
    t = datetime(core, 'InputFormat', 'uuuu-MM-dd''T''HH:mm:ss', ...
                 'TimeZone', 'UTC');
catch
    try
        t = datetime(core, 'TimeZone', 'UTC');
    catch
        return;
    end
end
h = hours(datetime('now', 'TimeZone', 'UTC') - t);
end


function s = ageStr(h)
if isnan(h)
    s = 'age unknown';
elseif h < 1
    s = sprintf('%.0f min ago', h * 60);
elseif h < 48
    s = sprintf('%.1f h ago', h);
else
    s = sprintf('%.1f days ago', h / 24);
end
end


function v = getf(s, field, dflt)
if isfield(s, field) && ~isempty(s.(field))
    v = string(s.(field));
else
    v = string(dflt);
end
end


function s = shortList(v)
v = v(:).';
if isempty(v); s = '(none)'; return; end
if numel(v) > 20
    s = sprintf('%s ... %s', strjoin(string(v(1:10)), ', '), ...
                strjoin(string(v(end-4:end)), ', '));
else
    s = strjoin(string(v), ', ');
end
end


function [ok, nDeleted, why] = deleteCase(URL, KEY, TABLE, caseID)
% Delete one case row and REPORT HOW MANY ROWS WERE ACTUALLY REMOVED.
%
% "Prefer: return=representation" makes PostgREST echo the deleted
% rows. Without it a DELETE blocked by Row Level Security returns
% 204 with an empty body -- indistinguishable from a real delete --
% so the script would cheerfully report success while the table
% stayed untouched.
ok = false;  nDeleted = 0;  why = '';

bodyFile = [tempname '.txt'];
cmd = sprintf(['curl -sS -o "%s" -w "%%{http_code}" -X DELETE ' ...
               '-H "apikey: %s" -H "Authorization: Bearer %s" ' ...
               '-H "Prefer: return=representation" ' ...
               '"%s/rest/v1/%s?case_id=eq.%d"'], ...
              bodyFile, KEY, KEY, URL, TABLE, caseID);
[st, out] = system(cmd);

body = '';
if isfile(bodyFile); body = strtrim(fileread(bodyFile)); delete(bodyFile); end

if st ~= 0
    why = sprintf('curl status %d', st);
    return;
end

code = str2double(regexp(strtrim(out), '\d{3}$', 'match', 'once'));
if ~any(code == [200 204])
    why = sprintf('HTTP %g: %s', code, body);
    return;
end

ok = true;

if isempty(body) || body(1) ~= '['
    nDeleted = 0;          % nothing echoed back
    return;
end
try
    d = jsondecode(body);
    nDeleted = numel(d);
catch
    nDeleted = 0;
end
end


function ok = httpJson(verb, URL, KEY, TABLE, query, payloadStruct)
ok = false;
% The JSON goes through a temp file: quoting a brace-and-quote
% payload inside cmd.exe is fragile.
pf = [tempname '.json'];
fid = fopen(pf, 'w');
if fid < 0; return; end
fwrite(fid, jsonencode(payloadStruct), 'char');
fclose(fid);

bodyFile = [tempname '.txt'];
cmd = sprintf(['curl -sS -o "%s" -w "%%{http_code}" -X %s ' ...
               '-H "apikey: %s" -H "Authorization: Bearer %s" ' ...
               '-H "Content-Type: application/json" ' ...
               '-H "Prefer: return=minimal" ' ...
               '--data-binary @"%s" "%s/rest/v1/%s%s"'], ...
              bodyFile, verb, KEY, KEY, pf, URL, TABLE, query);
[st, out] = system(cmd);
if isfile(bodyFile); delete(bodyFile); end
if isfile(pf); delete(pf); end
if st ~= 0; return; end
code = str2double(regexp(strtrim(out), '\d{3}$', 'match', 'once'));
ok = any(code == [200 201 204]);
end