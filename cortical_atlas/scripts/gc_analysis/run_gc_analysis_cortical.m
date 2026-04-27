%% run_gc_analysis_cortical.m
%
% State-space Granger Causality analysis for DKTatlas (62-region)
% source-reconstructed resting-state MEG data from the NIMH Healthy
% Research Volunteer Dataset (OpenNeuro ds005752).
%
% Source reconstruction performed via LCMV beamformer on individual
% anatomy with DKTatlas parcellation — see
% cortical_atlas/scripts/source_recon/run_pipeline_individual.py.
%
% GC toolbox: MVGC2 (Barnett & Seth) — https://github.com/lcbarnett/MVGC2
%
% PIPELINE
% --------
%   For each subject:
%     0. Load source_epochs.mat (converted from .npy via
%        cortical_atlas/scripts/source_recon/npy_to_mat.py)
%     1. Load per-epoch trialized data (label_ts_epochs,
%        n_epochs × n_labels × T). Data are passed to MVGC2 as
%        (n_labels × T × n_epochs); statistics are pooled across epochs
%        without crossing epoch boundaries.
%        Demean and normalise (per-channel std = 1).
%     2. Estimate VAR model order via BIC → set pf = 2 * varmo
%     3. Estimate SS model order via SVC (singular value criterion)
%     4. Fit innovations-form SS model via CCA subspace algorithm
%     5. Joint GC:     all others → each target          (ss_to_iogc)
%     6. Unconditional pairwise GC: source → target       (ss_to_ucpwgc)
%     7. Conditional pairwise GC:   source → target|rest  (ss_to_pwcgc)
%     8. Derived measures:
%          R̄(i) = Σ_j F_uc(j→i)   / F_joint(i)   (redundancy)
%          V̄(i) = Σ_j F_cond(j→i) / F_joint(i)   (vulnerability)
%
% USAGE
% -----
%   1. Ensure source_epochs.mat exists for each subject:
%        conda activate nimh-meg-atlas
%        python scripts/source_recon/npy_to_mat.py
%   2. In MATLAB:
%        cd /path/to/code_release/cortical_atlas
%        run('scripts/gc_analysis/run_gc_analysis_cortical.m')
%
% REQUIRES
% --------
%   - MVGC2 toolbox (auto-initialised from mvgc2_dir)
%   - Custom function: ss_to_ucpwgc.m (in cortical_atlas/scripts/gc_analysis/)

%% ========================================================================
%  PARAMETERS
%  ========================================================================

% --- Subject filter (set externally or leave empty for all) ---------------
% e.g.:  gc_subjects = {'sub-ON08710'};
%        gc_subjects = {};              % process all available
if ~exist('gc_subjects', 'var')
    gc_subjects = {};
end

% --- Test mode -----------------------------------------------------------
% 'off'      = full run (all available subjects)
% 'quick'    = 1 subject only (for debugging)
test_mode = 'off';

% --- Paths ---------------------------------------------------------------
project_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));

% Expected layout: <data_dir>/sub-*/source_epochs.mat
% Both Zenodo users (extract archive here) and full-pipeline users
% (make source + make mat) use this same directory.
data_dir    = fullfile(project_root, 'derivatives', 'source_timeseries');
results_dir = fullfile(project_root, 'results');

% MVGC2 toolbox path — edit if MVGC2 is not at <project_root>/MVGC2/
% Clone from: git clone https://github.com/lcbarnett/MVGC2.git
if ~exist('mvgc2_dir', 'var') || isempty(mvgc2_dir)
    mvgc2_dir = fullfile(project_root, 'MVGC2');
end

% --- Subject selection ---------------------------------------------------
% If true, restrict analysis to IDs in subjects.txt (recommended).
use_subject_list = true;
subject_list_file = fullfile(project_root, 'subjects.txt');

% --- Output mode ---------------------------------------------------------
% false (default): always create a fresh gc_dkt62_<timestamp> folder.
% true: resume/reuse most recent gc_dkt62_* folder.
resume_existing = false;

% --- Model order selection -----------------------------------------------
var_momax   = 30;           % max VAR order to test (BIC typically picks 10-25 at 300 Hz)
var_ic      = 'BIC';        % information criterion: 'AIC', 'BIC', 'HQC'
var_regmode = 'LWR';        % LWR = Levinson-Wiggins-Robinson (much faster than OLS for large n)

% --- SS model order (empty = auto via SVC) --------------------------------
ss_mosel = [];

% --- Preprocessing (data is already bandpass + notch + resampled at 300 Hz)
do_demean    = true;    % subtract per-channel temporal mean
do_normalise = true;    % divide by per-channel std

% --- Output --------------------------------------------------------------
gc_dirs = dir(fullfile(results_dir, 'gc_analysis', 'gc_dkt62_*'));
if resume_existing && ~isempty(gc_dirs)
    [~, idx] = sort({gc_dirs.name});
    output_folder = fullfile(gc_dirs(idx(end)).folder, gc_dirs(idx(end)).name);
    fprintf('Resuming existing output folder: %s\n', output_folder);
else
    output_tag    = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
    output_folder = fullfile(results_dir, 'gc_analysis', ['gc_dkt62_' output_tag]);
    fprintf('Creating fresh output folder: %s\n', output_folder);
end

%% ========================================================================
%  DISCOVER SUBJECTS
%  ========================================================================

% Auto-discover subjects that have source_epochs.mat in the data directory
d = dir(fullfile(data_dir, 'sub-*'));
candidate_ids = {d([d.isdir]).name}';

% Filter to subjects that actually have source_epochs.mat
subject_ids = {};
for i = 1:length(candidate_ids)
    mat_path = fullfile(data_dir, candidate_ids{i}, 'source_epochs.mat');
    if exist(mat_path, 'file')
        subject_ids{end+1} = candidate_ids{i}; %#ok<SAGROW>
    end
end
subject_ids = subject_ids(:);

if isempty(subject_ids)
    error('No subjects with source_epochs.mat found in %s.\nRun: python scripts/source_recon/npy_to_mat.py', data_dir);
end

% Restrict to curated analysis subject list (default: subjects.txt)
if use_subject_list
    listed_ids = read_subject_list(subject_list_file);
    if isempty(listed_ids)
        error('Subject list is empty: %s', subject_list_file);
    end

    present = ismember(listed_ids, subject_ids);
    missing_ids = listed_ids(~present);
    subject_ids = listed_ids(present);

    if isempty(subject_ids)
        error('No listed subjects have source_epochs.mat. List file: %s', subject_list_file);
    end

    if ~isempty(missing_ids)
        fprintf('WARNING: %d listed subjects missing source_epochs.mat and will be skipped.\n', ...
            length(missing_ids));
        fprintf('  Missing: %s\n', strjoin(missing_ids, ', '));
    end
end

% Apply subject filter if specified
if ~isempty(gc_subjects)
    subject_ids = subject_ids(ismember(subject_ids, gc_subjects));
    if isempty(subject_ids)
        error('None of the requested subjects have source_epochs.mat.');
    end
    fprintf('Filtered to %d subject(s): %s\n', length(subject_ids), strjoin(subject_ids, ', '));
end

%% ========================================================================
%  INITIALISATION
%  ========================================================================

% --- MVGC2 ---------------------------------------------------------------
if ~exist('tsdata_to_ss', 'file')
    fprintf('Initialising MVGC2 toolbox...\n');
    old_dir = pwd;
    cd(mvgc2_dir);
    startup;
    cd(old_dir);
end

% --- Custom GC functions (this scripts/ directory) -----------------------
scripts_dir = fileparts(mfilename('fullpath'));
addpath(scripts_dir);

% --- Output directory ----------------------------------------------------
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% --- Test mode overrides -------------------------------------------------
if strcmp(test_mode, 'quick')
    subject_ids = subject_ids(1);
    fprintf('\n*** TEST MODE: QUICK — 1 subject only ***\n\n');
end

n_subjects = length(subject_ids);

% --- Store parameters ----------------------------------------------------
params = struct();
params.var_momax    = var_momax;
params.var_ic       = var_ic;
params.var_regmode  = var_regmode;
params.ss_mosel     = ss_mosel;
params.do_demean    = do_demean;
params.do_normalise = do_normalise;
params.test_mode    = test_mode;
params.data_dir     = fullfile('derivatives', 'source_timeseries');
params.mvgc2_dir    = 'MVGC2';

%% ========================================================================
%  PRINT HEADER
%  ========================================================================

fprintf('\n=========================================================\n');
fprintf('  GC State-Space Analysis — DKTatlas-62 MEG Resting State\n');
fprintf('  NIMH Healthy Research Volunteer Dataset\n');
fprintf('=========================================================\n');
fprintf('  Subjects:        %d\n', n_subjects);
fprintf('  Parcellation:    DKTatlas (aparc.DKTatlas, %d regions)\n', 62);
fprintf('  Data:            source_epochs.mat (pre-resampled 300 Hz)\n');
if use_subject_list
    fprintf('  Subject list:    %s\n', subject_list_file);
else
    fprintf('  Subject list:    auto-discovered from data_dir\n');
end
fprintf('  VAR max order:   %d (%s)\n', var_momax, var_ic);
fprintf('  SS pf:           2 × varmo (Bauer heuristic)\n');
fprintf('  Demean:          %d,  Normalise: %d\n', do_demean, do_normalise);
fprintf('  Resume mode:     %d\n', resume_existing);
fprintf('  Output:          %s\n', output_folder);
fprintf('=========================================================\n\n');

%% ========================================================================
%  MAIN LOOP
%  ========================================================================

all_results = cell(n_subjects, 1);
total_timer = tic;

% --- Load existing results (if any) to skip already-processed subjects ---
existing_mat = fullfile(output_folder, 'gc_results.mat');
existing_sids = {};
if exist(existing_mat, 'file')
    prev = load(existing_mat, 'all_results', 'subject_ids');
    if isfield(prev, 'all_results') && isfield(prev, 'subject_ids')
        for k = 1:length(prev.subject_ids)
            if ~isempty(prev.all_results{k})
                existing_sids{end+1} = prev.subject_ids{k}; %#ok<SAGROW>
            end
        end
        % Pre-fill all_results for subjects carried forward
        for k = 1:length(prev.subject_ids)
            idx = find(strcmp(subject_ids, prev.subject_ids{k}));
            if ~isempty(idx) && ~isempty(prev.all_results{k})
                all_results{idx} = prev.all_results{k};
            end
        end
    end
end

for s = 1:n_subjects
    sid = subject_ids{s};

    % Skip subjects already in existing results
    if ismember(sid, existing_sids)
        fprintf('\n  SKIP %s — already in %s\n', sid, output_folder);
        continue;
    end

    fprintf('\n============ Subject %d/%d: %s ============\n', s, n_subjects, sid);
    subj_timer = tic;

    % --- Load data -------------------------------------------------------
    mat_file = fullfile(data_dir, sid, 'source_epochs.mat');
    t0 = tic;
    S = load(mat_file);

    % Prefer trialized epoch data to avoid cross-epoch boundary lags.
    if isfield(S, 'label_ts_epochs') && ~isempty(S.label_ts_epochs) && ndims(S.label_ts_epochs) == 3
        % source_epochs.mat stores epochs as (n_epochs × n_labels × n_times)
        % MVGC expects (n_labels × n_times × n_epochs)
        X3d = permute(double(S.label_ts_epochs), [2 3 1]);
        input_mode = 'epochs';
    elseif isfield(S, 'label_ts') && ~isempty(S.label_ts)
        % Fallback for legacy files: concatenated single-trial data
        X = double(S.label_ts);          % (n_labels × n_times)
        X3d = reshape(X, size(X, 1), size(X, 2), 1);
        input_mode = 'concatenated_fallback';
    else
        error('No usable input in %s (expected label_ts_epochs or label_ts).', mat_file);
    end

    [n, m, N] = size(X3d);         % n = channels, m = samples/trial, N = trials
    sfreq = double(S.sfreq);

    % Read labels — handle both cell array and char array from .mat
    if iscell(S.labels)
        chan_labels = cellfun(@(x) strtrim(char(x)), S.labels, 'UniformOutput', false);
    else
        chan_labels = cellstr(S.labels);
    end
    chan_labels = chan_labels(:)';

    fprintf('  Loaded in %.1f s\n', toc(t0));
    fprintf('  Channels:  %d\n', n);
    fprintf('  Input:     %s\n', input_mode);
    fprintf('  Trials:    %d\n', N);
    fprintf('  Samples:   %d per trial, %d total (%.1f s total at %.0f Hz)\n', ...
        m, m * N, (m * N) / sfreq, sfreq);

    % --- Preprocessing ---------------------------------------------------
    X3d = preprocess_tsdata(X3d, do_demean, do_normalise);

    % --- Fit SS model ----------------------------------------------------
    res = fit_ss_model(X3d, var_momax, var_ic, var_regmode, ss_mosel);

    if res.info.error
        fprintf(2, '  *** SS model has issues (error=%d). Continuing anyway. ***\n', res.info.error);
    end

    % --- Step 1: Joint GC (all → each target) ----------------------------
    fprintf('     Step 1: Joint GC (all → each target)...\n');
    t0 = tic;
    F_joint = ss_to_iogc(res.A, res.C, res.K, res.V, 'in');
    fprintf('     Step 1 done: mean=%.6f  [%.1f s]\n', mean(F_joint), toc(t0));
    res.F_joint = F_joint;

    % --- Step 2: Unconditional pairwise GC --------------------------------
    fprintf('     Step 2: Unconditional pairwise GC (%d×%d DARE)...\n', n, n);
    t0 = tic;
    F_uc = ss_to_ucpwgc(res.A, res.C, res.K, res.V, 1:n, 1:n);
    fprintf('     Step 2 done: mean=%.6f  [%.1f s]\n', nanmean(F_uc(:)), toc(t0));
    res.F_uc = F_uc;

    % --- Step 3: Conditional pairwise GC ----------------------------------
    fprintf('     Step 3: Conditional pairwise GC (%d×%d)...\n', n, n);
    t0 = tic;
    F_cond = ss_to_pwcgc(res.A, res.C, res.K, res.V);
    fprintf('     Step 3 done: mean=%.6f  [%.1f s]\n', nanmean(F_cond(:)), toc(t0));
    res.F_cond = F_cond;

    % --- Derived measures: redundancy & vulnerability ---------------------
    redundancy    = nan(n, 1);
    vulnerability = nan(n, 1);
    for i = 1:n
        uc_sum   = nansum(F_uc(i, :));       % Σ_j F_uc(j→i)
        cond_sum = nansum(F_cond(i, :));      % Σ_j F_cond(j→i)
        if F_joint(i) > 0
            redundancy(i)    = uc_sum   / F_joint(i);
            vulnerability(i) = cond_sum / F_joint(i);
        end
    end
    res.redundancy    = redundancy;
    res.vulnerability = vulnerability;

    fprintf('     Redundancy:     mean=%.4f  [%.4f, %.4f]\n', ...
        nanmean(redundancy), nanmin(redundancy), nanmax(redundancy));
    fprintf('     Vulnerability:  mean=%.4f  [%.4f, %.4f]\n', ...
        nanmean(vulnerability), nanmin(vulnerability), nanmax(vulnerability));

    % --- Store results ---------------------------------------------------
    subj_res = struct();
    subj_res.subject_id     = sid;
    subj_res.chan_labels     = {chan_labels};
    subj_res.sfreq           = sfreq;
    subj_res.input_mode      = input_mode;
    subj_res.n_channels     = n;
    subj_res.n_trials       = N;
    subj_res.n_samples_per_trial = m;
    subj_res.n_samples      = m * N;       % total samples across trials
    subj_res.duration_s     = (m * N) / sfreq;
    subj_res.duration_per_trial_s = m / sfreq;
    subj_res.results        = res;

    all_results{s} = subj_res;

    % --- Print summary ---------------------------------------------------
    fprintf('\n  Subject %s done (%.1f s):\n', sid, toc(subj_timer));
    fprintf('    VAR=%d, SS=%d, pf=%d\n', res.var_mo, res.ss_mo, res.pf);
    fprintf('    rho(A)=%.6f, rho(A-KC)=%.6f, error=%d\n', ...
        res.info.rhoA, res.info.rhoB, res.info.error);
    fprintf('    mean joint=%.6f, mean uc=%.6f, mean cond=%.6f\n', ...
        mean(F_joint), nanmean(F_uc(:)), nanmean(F_cond(:)));
    fprintf('    mean R=%.4f, mean V=%.4f\n', ...
        nanmean(redundancy), nanmean(vulnerability));
end

%% ========================================================================
%  SAVE RESULTS
%  ========================================================================

total_elapsed = toc(total_timer);
fprintf('\n\nTotal elapsed: %.1f s (%.1f min)\n', total_elapsed, total_elapsed/60);
fprintf('Saving to %s\n', output_folder);

% Full .mat
save(fullfile(output_folder, 'gc_results.mat'), ...
    'all_results', 'params', 'subject_ids', '-v7.3');

% Summary CSV
write_summary_csv(all_results, output_folder);

fprintf('\nDone.\n');

%% ========================================================================
%  HELPER: FIT SS MODEL
%  ========================================================================

function res = fit_ss_model(X, var_momax, var_ic, var_regmode, ss_mosel)
% FIT_SS_MODEL  Fit innovations-form SS model via MVGC2.
%
%   1. VAR order estimation (BIC/AIC/HQC)
%   2. pf = 2 * varmo  (Bauer's heuristic)
%   3. SS order estimation (SVC)
%   4. CCA subspace fit → A, C, K, V
%
    %   Input:  X (n × m × N).
%   Output: struct with var_mo, ss_mo, pf, A, C, K, V, info.

    [n, m, ~] = size(X);

    % --- VAR model order --------------------------------------------------
    t0 = tic;
    eff_momax = min(var_momax, floor(m / 3));
    if eff_momax < 1
        warning('Data too short (m=%d). Setting var_mo=1.', m);
        var_mo = 1;
    else
        [moaic, mobic, mohqc, ~] = tsdata_to_varmo(X, eff_momax, var_regmode, ...
            [], false, [], 0);
        switch upper(var_ic)
            case 'AIC', var_mo = moaic;
            case 'BIC', var_mo = mobic;
            case 'HQC', var_mo = mohqc;
            otherwise,  var_mo = moaic;
        end
    end
    fprintf('     VAR order (%s): %d  [%.1f s]\n', var_ic, var_mo, toc(t0));

    if var_mo >= eff_momax
        fprintf(2, '     *** WARNING: VAR order hit ceiling %d ***\n', eff_momax);
    end

    % --- Past/future horizon -----------------------------------------------
    pf = 2 * var_mo;
    fprintf('     pf = 2 × %d = %d\n', var_mo, pf);

    pf_max = floor(m / 2) - 1;
    if pf > pf_max
        fprintf('     ** Capping pf from %d to %d (data length m=%d) **\n', pf, pf_max, m);
        pf = pf_max;
    end

    % --- SS model order (SVC) ---------------------------------------------
    t0 = tic;
    if isempty(ss_mosel)
        [ss_mo, ~] = tsdata_to_ssmo(X, pf, []);
    else
        ss_mo = ss_mosel;
    end
    if ss_mo < 1
        fprintf('     ** SS order = 0, forcing to 1 **\n');
        ss_mo = 1;
    end
    fprintf('     SS order (SVC): %d  [%.1f s]\n', ss_mo, toc(t0));

    % --- Fit SS model ------------------------------------------------------
    t0 = tic;
    [A, C, K, V] = tsdata_to_ss(X, pf, ss_mo);
    info = ss_info(A, C, K, V);

    if info.error
        warning('SS model issues (error=%d).', info.error);
    end
    if info.rhoA > 0.999
        warning('rho(A)=%.6f near unit circle.', info.rhoA);
    end
    if info.rhoA >= 1.0
        warning('rho(A)=%.6f ≥ 1: UNSTABLE.', info.rhoA);
    end
    fprintf('     SS fit: rho(A)=%.6f, rho(A-KC)=%.6f  [%.1f s]\n', ...
        info.rhoA, info.rhoB, toc(t0));

    % --- Pack ---------------------------------------------------------------
    res = struct();
    res.var_mo = var_mo;
    res.ss_mo  = ss_mo;
    res.pf     = pf;
    res.info   = info;
    res.A      = A;
    res.C      = C;
    res.K      = K;
    res.V      = V;
end

%% ========================================================================
%  HELPER: PREPROCESS MULTI-TRIAL DATA
%  ========================================================================

function X = preprocess_tsdata(X, do_demean, do_normalise)
% Pool observations across all trials for channel-wise preprocessing.
% This follows MVGC conventions for multi-trial data.
    [n, m, N] = size(X);
    X2 = reshape(X, n, m * N);
    if do_demean
        X2 = X2 - mean(X2, 2);
    end
    if do_normalise
        ch_std = std(X2, 0, 2);
        ch_std(ch_std < eps) = 1;
        X2 = X2 ./ ch_std;
    end
    X = reshape(X2, n, m, N);
end

%% ========================================================================
%  HELPER: WRITE SUMMARY CSV
%  ========================================================================

function write_summary_csv(all_results, output_folder)
    csv_file = fullfile(output_folder, 'gc_summary.csv');
    fid = fopen(csv_file, 'w');

    % Header
    fprintf(fid, 'subject_id,n_channels,n_trials,n_samples_per_trial,n_samples_total,duration_s,sfreq,');
    fprintf(fid, 'var_mo,ss_mo,pf,rhoA,rhoB,ss_error,');
    fprintf(fid, 'mean_joint_gc,mean_uc_gc,mean_cond_gc,');
    fprintf(fid, 'mean_redundancy,mean_vulnerability\n');

    for s = 1:length(all_results)
        r = all_results{s};
        if isempty(r), continue; end
        res = r.results;

        if isfield(r, 'n_trials'), n_trials = r.n_trials; else, n_trials = 1; end
        if isfield(r, 'n_samples_per_trial'), n_spt = r.n_samples_per_trial; else, n_spt = r.n_samples; end
        if isfield(r, 'n_samples'), n_stot = r.n_samples; else, n_stot = n_spt * n_trials; end
        if isfield(r, 'duration_s'), duration_s = r.duration_s; else, duration_s = n_stot / r.sfreq; end

        fprintf(fid, '%s,%d,%d,%d,%d,%.1f,%.0f,', ...
            r.subject_id, r.n_channels, n_trials, n_spt, n_stot, duration_s, r.sfreq);
        fprintf(fid, '%d,%d,%d,%.6f,%.6f,%d,', ...
            res.var_mo, res.ss_mo, res.pf, res.info.rhoA, res.info.rhoB, res.info.error);
        fprintf(fid, '%.6f,%.6f,%.6f,', ...
            mean(res.F_joint), nanmean(res.F_uc(:)), nanmean(res.F_cond(:)));
        fprintf(fid, '%.6f,%.6f\n', ...
            nanmean(res.redundancy), nanmean(res.vulnerability));
    end

    fclose(fid);
    fprintf('Summary saved to %s\n', csv_file);
end

%% ========================================================================
%  HELPER: READ SUBJECT LIST
%  ========================================================================

function subject_ids = read_subject_list(list_file)
% Read newline-separated subject IDs, skipping blank/comment lines.
    if ~isfile(list_file)
        error('Subject list not found: %s', list_file);
    end

    txt = fileread(list_file);
    lines = regexp(txt, '\r\n|\n|\r', 'split');
    lines = lines(:);

    subject_ids = {};
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if isempty(line) || line(1) == '#'
            continue;
        end
        % Keep first token only if trailing inline comments/whitespace exist.
        toks = regexp(line, '\s+', 'split');
        subject_ids{end+1,1} = toks{1}; %#ok<AGROW>
    end
end
