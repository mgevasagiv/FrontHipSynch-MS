% audit_native_sampling_rates.m
% Report each manuscript subject's ORIGINAL acquisition sampling rate.
%
% Primary source (as pointed to by user):
%   E:\MAZE_1.0\data\iEEG\raw_orig\   (Neuralynx .ncs headers and/or EDF/BESA)
% Fallback:
%   E:\MAZE_1.0\data\iEEG\preproc\<subj>\<subj>_maze_wtrig.set  (EEG.srate)
%
% Also reports processed fsample from *_maze_wtrig_pp.mat (should be 512).
%
% Reviewer context: hardware anti-aliasing for native ~512 Hz acquisition is often
% near Fs/4 (~128 Hz), close to high-gamma / ripple bands. Subjects acquired at
% ~1 kHz then software-downsampled are less affected (FieldTrip anti-aliases).

clear; clc;

%% -------- paths (override if needed) --------
pj_home = 'E:\MAZE_1.0\';
if exist('maze_set_path', 'file')
    try
        maze_set_path;
        if exist('dirs', 'var') && isfield(dirs, 'prepfile')
            % keep dirs from maze_set_path when hostname is recognized
        end
    catch
        % hostname not recognized in cloud / other machines — use pj_home below
    end
end
if ~exist('dirs', 'var') || ~isfield(dirs, 'prepfile')
    dirs.raw_orig = fullfile(pj_home, 'data', 'iEEG', 'raw_orig');
    dirs.nldata   = fullfile(dirs.raw_orig, 'neuralynx');
    dirs.rawdata  = fullfile(pj_home, 'data', 'iEEG', 'raw');
    dirs.prepfile = fullfile(pj_home, 'data', 'iEEG', 'preproc');
    dirs.analysis_root = fullfile(pj_home, 'analysis');
else
    if ~isfield(dirs, 'raw_orig')
        dirs.raw_orig = fullfile(pj_home, 'data', 'iEEG', 'raw_orig');
    end
    if ~isfield(dirs, 'nldata')
        dirs.nldata = fullfile(dirs.raw_orig, 'neuralynx');
    end
end

subjects = {'s01','s06','s08','s09','s12','s14','s15','s16','s17','s18'};

fprintf('raw_orig: %s  (exists=%d)\n', dirs.raw_orig, exist(dirs.raw_orig, 'dir')==7);
fprintf('neuralynx: %s  (exists=%d)\n', dirs.nldata, exist(dirs.nldata, 'dir')==7);
fprintf('preproc:  %s  (exists=%d)\n', dirs.prepfile, exist(dirs.prepfile, 'dir')==7);

% Map subject folders under raw_orig (s01, S01, pda017-style aliases if present)
raw_listing = {};
if exist(dirs.raw_orig, 'dir')
    d = dir(dirs.raw_orig);
    raw_listing = {d([d.isdir] & ~ismember({d.name},{'.','..'})).name};
end
if exist(dirs.nldata, 'dir')
    d = dir(dirs.nldata);
    raw_listing = unique([raw_listing, {d([d.isdir] & ~ismember({d.name},{'.','..'})).name}]);
end
fprintf('Subject-like folders under raw_orig/neuralynx:\n');
disp(raw_listing');

use_eeglab = exist('pop_loadset', 'file') == 2;

n = numel(subjects);
subj_id = subjects(:);
native_fs = nan(n, 1);
raw_fs = nan(n, 1);
set_fs = nan(n, 1);
processed_fs = nan(n, 1);
source = strings(n, 1);
raw_folder = strings(n, 1);
impact_flag = strings(n, 1);
notes = strings(n, 1);

for i = 1:n
    subj = subjects{i};
    note_parts = strings(0);

    %% 1) raw_orig Neuralynx / EDF
    cand_dirs = localize_raw_subj_dirs(subj, dirs, raw_listing);
    if ~isempty(cand_dirs)
        raw_folder(i) = strjoin(cand_dirs, ' | ');
        [fs_raw, how] = read_fs_from_raw_dirs(cand_dirs);
        if ~isnan(fs_raw)
            raw_fs(i) = fs_raw;
            note_parts(end+1) = "raw:" + string(how); %#ok<AGROW>
        else
            note_parts(end+1) = "raw folders found but no Fs parsed"; %#ok<AGROW>
        end
    else
        note_parts(end+1) = "no raw_orig folder matched"; %#ok<AGROW>
    end

    %% 2) preproc .set (pre-downsample EEGLAB)
    subjdir = fullfile(dirs.prepfile, subj);
    setfile = fullfile(subjdir, [subj '_maze_wtrig.set']);
    if isfile(setfile)
        try
            if use_eeglab
                EEG = pop_loadset('filename', [subj '_maze_wtrig.set'], 'filepath', [subjdir filesep]);
                set_fs(i) = EEG.srate;
                clear EEG;
            else
                S = load(setfile, '-mat');
                if isfield(S, 'EEG'); set_fs(i) = S.EEG.srate;
                elseif isfield(S, 'srate'); set_fs(i) = S.srate;
                else; note_parts(end+1) = "set has no srate"; %#ok<AGROW>
                end
                clear S;
            end
        catch ME
            note_parts(end+1) = "set load failed: " + string(ME.message); %#ok<AGROW>
        end
    else
        note_parts(end+1) = "missing .set"; %#ok<AGROW>
    end

    %% 3) processed pp.mat
    ppfile = fullfile(subjdir, [subj '_maze_wtrig_pp.mat']);
    if isfile(ppfile)
        try
            mf = matfile(ppfile);
            d = mf.data;
            if isstruct(d) && isfield(d, 'fsample'); processed_fs(i) = d.fsample; end
            try
                ofs = mf.original_fsample;
                note_parts(end+1) = "pp.original_fsample=" + string(ofs); %#ok<AGROW>
            catch
            end
            clear mf d;
        catch ME
            note_parts(end+1) = "pp load failed: " + string(ME.message); %#ok<AGROW>
        end
    end

    %% Choose native Fs: prefer raw header, else .set
    if ~isnan(raw_fs(i))
        native_fs(i) = raw_fs(i);
        source(i) = "raw_orig";
    elseif ~isnan(set_fs(i))
        native_fs(i) = set_fs(i);
        source(i) = "preproc_set";
    else
        source(i) = "none";
    end

    if ~isnan(native_fs(i)) && native_fs(i) <= 520 && native_fs(i) >= 500
        impact_flag(i) = "impacted_native_~512Hz";
    elseif ~isnan(native_fs(i)) && native_fs(i) >= 900
        impact_flag(i) = "not_impacted_native_>=1kHz";
    elseif ~isnan(native_fs(i))
        impact_flag(i) = "check_manually";
    else
        impact_flag(i) = "unknown";
    end
    notes(i) = strjoin(note_parts, "; ");
end

T = table(subj_id, native_fs, raw_fs, set_fs, processed_fs, source, raw_folder, impact_flag, notes, ...
    'VariableNames', {'subj','native_fs_Hz','raw_fs_Hz','set_fs_Hz','processed_fs_Hz', ...
    'fs_source','raw_folder','reviewer_AA_concern','notes'});
disp(T);

n_impacted = sum(impact_flag == "impacted_native_~512Hz");
n_ok = sum(impact_flag == "not_impacted_native_>=1kHz");
n_unk = sum(impact_flag == "unknown" | impact_flag == "check_manually");
fprintf('\n=== Summary (N=%d manuscript subjects) ===\n', n);
fprintf('Native ~512 Hz (reviewer AA concern): %d\n', n_impacted);
fprintf('Native >= ~1 kHz: %d\n', n_ok);
fprintf('Unknown / manual check: %d\n', n_unk);
if n_impacted > 0
    fprintf('Impacted: %s\n', strjoin(subjects(impact_flag == "impacted_native_~512Hz"), ', '));
end
if n_ok > 0
    fprintf('Higher-rate: %s\n', strjoin(subjects(impact_flag == "not_impacted_native_>=1kHz"), ', '));
end

outdir = fullfile(dirs.analysis_root, 'sampling_rate_audit');
if ~exist(outdir, 'dir'); mkdir(outdir); end
outfile = fullfile(outdir, 'sampling_rate_audit.csv');
writetable(T, outfile);
fprintf('\nWrote %s\n', outfile);

%% ================= local helpers =================
function cand = localize_raw_subj_dirs(subj, dirs, raw_listing)
    cand = {};
    patterns = unique({subj, upper(subj), lower(subj)}, 'stable');
    roots = unique({dirs.raw_orig, dirs.nldata}, 'stable');
    for r = 1:numel(roots)
        root = roots{r};
        if exist(root, 'dir') ~= 7; continue; end
        for p = 1:numel(patterns)
            hit = fullfile(root, patterns{p});
            if exist(hit, 'dir') == 7
                cand{end+1} = hit; %#ok<AGROW>
            end
        end
        % fuzzy: folder name contains subject id
        for k = 1:numel(raw_listing)
            name = raw_listing{k};
            if contains(lower(name), lower(subj))
                hit = fullfile(root, name);
                if exist(hit, 'dir') == 7
                    cand{end+1} = hit; %#ok<AGROW>
                end
            end
        end
    end
    cand = unique(cand, 'stable');
end

function [fs, how] = read_fs_from_raw_dirs(cand_dirs)
    fs = nan; how = "";
    rates = [];
    sources = strings(0);
    for d = 1:numel(cand_dirs)
        root = cand_dirs{d};
        % Neuralynx Continuously Sampled Channel files
        ncs = dir(fullfile(root, '**', '*.ncs'));
        for i = 1:min(5, numel(ncs)) % a few files are enough; rates should match
            f = fullfile(ncs(i).folder, ncs(i).name);
            r = read_nlx_sampling_frequency(f);
            if ~isnan(r)
                rates(end+1) = r; %#ok<AGROW>
                sources(end+1) = "ncs:" + string(ncs(i).name); %#ok<AGROW>
            end
        end
        % EDF
        edf = [dir(fullfile(root, '**', '*.edf')); dir(fullfile(root, '**', '*.EDF'))];
        for i = 1:min(3, numel(edf))
            f = fullfile(edf(i).folder, edf(i).name);
            r = read_edf_sampling_frequency(f);
            if ~isnan(r)
                rates(end+1) = r; %#ok<AGROW>
                sources(end+1) = "edf:" + string(edf(i).name); %#ok<AGROW>
            end
        end
        % BESA / generic via FieldTrip if available
        besa = [dir(fullfile(root, '**', '*.besa')); dir(fullfile(root, '**', '*.BESA'))];
        if exist('ft_read_header', 'file') == 2
            for i = 1:min(2, numel(besa))
                f = fullfile(besa(i).folder, besa(i).name);
                try
                    hdr = ft_read_header(f);
                    rates(end+1) = hdr.Fs; %#ok<AGROW>
                    sources(end+1) = "besa:" + string(besa(i).name); %#ok<AGROW>
                catch
                end
            end
        end
    end
    if isempty(rates); return; end
    fs = mode(rates);
    how = sprintf('mode=%.3f from {%s} values [%s]', fs, strjoin(cellstr(sources), ', '), num2str(rates));
end

function fs = read_nlx_sampling_frequency(ncs_path)
    % Neuralynx .ncs: first 16 KB is ASCII header with -SamplingFrequency <Hz>
    fs = nan;
    fid = fopen(ncs_path, 'r');
    if fid < 0; return; end
    cleaner = onCleanup(@() fclose(fid));
    hdr = fread(fid, 16384, '*char')';
    tok = regexp(hdr, '-SamplingFrequency\s+([0-9.]+)', 'tokens', 'once');
    if ~isempty(tok)
        fs = str2double(tok{1});
    end
end

function fs = read_edf_sampling_frequency(edf_path)
    % Minimal EDF header parse: Fs = (samples/record) / (record duration)
    fs = nan;
    fid = fopen(edf_path, 'r');
    if fid < 0; return; end
    cleaner = onCleanup(@() fclose(fid));
    hdr = fread(fid, 256, '*char')';
    if numel(hdr) < 256; return; end
    nchan = str2double(strtrim(hdr(253:256)));
    if ~(nchan >= 1 && nchan < 10000); return; end
    dur = str2double(strtrim(hdr(245:252)));
    % After 256-byte fixed header: nchan * 16 label, ... eventually nchan * 8 samples/record
    % Seek to samples-per-record field: 256 + nchan*(16+80+8+8+8+8+8+8+8) = 256 + nchan*256 - nchan*8?
    % EDF per-channel header is 256 bytes * nchan after main header; last 8 bytes of each
    % channel header are nr of samples / record. Spec: fields accumulate as:
    % labels 16*n, transducer 80*n, dim 8*n, physmin 8*n, physmax 8*n, digmin 8*n,
    % digmax 8*n, prefilt 80*n, samples 8*n, reserved 32*n
    skip = 256 + nchan*(16+80+8+8+8+8+8+80);
    fseek(fid, skip, 'bof');
    samp_field = fread(fid, 8*nchan, '*char')';
    if numel(samp_field) < 8; return; end
    ns = str2double(strtrim(samp_field(1:8))); % first channel
    if dur > 0 && ns > 0
        fs = ns / dur;
    end
end
