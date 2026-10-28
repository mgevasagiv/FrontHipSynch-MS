% audit_native_sampling_rates.m
% Report each subject's ORIGINAL acquisition sampling rate (before MAZE3 downsample).
%
% Reads:  <dirs.prepfile>/<subj>/<subj>_maze_wtrig.set   (pre-downsample EEGLAB)
% Also checks processed: <subj>_maze_wtrig_pp.mat         (should be 512 Hz)
% Writes: analysis/sampling_rate_audit.csv  (+ prints a summary table)
%
% Reviewer context: hardware anti-aliasing for native ~512 Hz acquisition is often
% near Fs/4 (~128 Hz), which can encroach on high-gamma / ripple bands analyzed here.
% Subjects acquired at ~1 kHz (then software-downsampled to 512) are less affected,
% because FieldTrip ft_resampledata anti-aliases before decimation.

clear; clc;
maze_set_path;

subjects = {'s01','s06','s08','s09','s12','s14','s15','s16','s17','s18'};

% Prefer EEGLAB if available; otherwise try loading .set as a MATLAB struct.
use_eeglab = exist('pop_loadset', 'file') == 2;
if ~use_eeglab
    warning(['pop_loadset not on path. Will try load()-ing .set as a MAT struct. ', ...
        'If that fails, add EEGLAB to the path and re-run.']);
end

n = numel(subjects);
subj_id = subjects(:);
native_fs = nan(n, 1);
processed_fs = nan(n, 1);
set_exists = false(n, 1);
pp_exists = false(n, 1);
impact_flag = strings(n, 1); % "impacted" if native Fs is near 512
notes = strings(n, 1);

for i = 1:n
    subj = subjects{i};
    subjdir = fullfile(dirs.prepfile, subj);
    setfile = fullfile(subjdir, [subj '_maze_wtrig.set']);
    ppfile  = fullfile(subjdir, [subj '_maze_wtrig_pp.mat']);

    set_exists(i) = isfile(setfile);
    pp_exists(i)  = isfile(ppfile);

    if set_exists(i)
        try
            if use_eeglab
                EEG = pop_loadset('filename', [subj '_maze_wtrig.set'], 'filepath', [subjdir filesep]);
                native_fs(i) = EEG.srate;
                clear EEG;
            else
                S = load(setfile, '-mat');
                if isfield(S, 'EEG')
                    native_fs(i) = S.EEG.srate;
                elseif isfield(S, 'srate')
                    native_fs(i) = S.srate;
                else
                    notes(i) = "set loaded but no srate field";
                end
                clear S;
            end
        catch ME
            notes(i) = "set load failed: " + string(ME.message);
        end
    else
        notes(i) = "missing .set file";
    end

    if pp_exists(i)
        try
            % Prefer matfile to avoid loading huge continuous data into RAM
            mf = matfile(ppfile);
            d = mf.data;
            if isstruct(d) && isfield(d, 'fsample')
                processed_fs(i) = d.fsample;
            end
            % Also surface any original_fsample we may have saved from an updated MAZE3
            try
                ofs = mf.original_fsample; %#ok<*NASGU>
                if isempty(notes(i)) || notes(i) == ""
                    notes(i) = "pp.original_fsample=" + string(ofs);
                else
                    notes(i) = notes(i) + "; pp.original_fsample=" + string(ofs);
                end
            catch
            end
            clear mf d;
        catch ME
            if notes(i) == ""
                notes(i) = "pp load failed: " + string(ME.message);
            else
                notes(i) = notes(i) + "; pp load failed: " + string(ME.message);
            end
        end
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
end

T = table(subj_id, native_fs, processed_fs, set_exists, pp_exists, impact_flag, notes, ...
    'VariableNames', {'subj','native_fs_Hz','processed_fs_Hz','has_set','has_pp','reviewer_AA_concern','notes'});

disp(T);

n_impacted = sum(impact_flag == "impacted_native_~512Hz");
n_ok = sum(impact_flag == "not_impacted_native_>=1kHz");
n_unk = sum(impact_flag == "unknown" | impact_flag == "check_manually");

fprintf('\n=== Summary (N=%d manuscript subjects) ===\n', n);
fprintf('Native ~512 Hz (reviewer AA concern): %d\n', n_impacted);
fprintf('Native >= ~1 kHz (software downsample to 512): %d\n', n_ok);
fprintf('Unknown / needs manual check: %d\n', n_unk);

if n_impacted > 0
    fprintf('Impacted subjects: %s\n', strjoin(subjects(impact_flag == "impacted_native_~512Hz"), ', '));
end
if n_ok > 0
    fprintf('Higher-rate subjects: %s\n', strjoin(subjects(impact_flag == "not_impacted_native_>=1kHz"), ', '));
end

outdir = fullfile(dirs.analysis_root, 'sampling_rate_audit');
if ~exist(outdir, 'dir'); mkdir(outdir); end
outfile = fullfile(outdir, 'sampling_rate_audit.csv');
writetable(T, outfile);
fprintf('\nWrote %s\n', outfile);
