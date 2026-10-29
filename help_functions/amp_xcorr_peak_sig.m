function out = amp_xcorr_peak_sig(A, B, varargin)
%AMP_XCORR_PEAK_SIG Amplitude cross-correlation peak + shuffle significance
%
% Implements the Adhikari et al. (2010) / manuscript method for ascertaining
% significance of cross-correlation peaks between two amplitude time courses
% (e.g., block-1 vs block-2 OFC mid-gamma power):
%   1) Mean-subtract each vector
%   2) Cross-correlate over +/- maxLag samples (xcorr)
%   3) Take the lag at which the cross-correlation peaks
%   4) Build a null of peak values by randomly permuting one time course
%      nPerm times and re-computing xcorr peaks
%   5) Call the observed peak significant if it exceeds the 95th percentile
%      of the null (also return an empirical p-value)
%
% Optionally also returns corrcoef(A,B) at lag 0 (as reported in the
% manuscript text). That lag-0 parametric p is NOT a valid test of the
% cross-correlation peak; use out.peakSignificant / out.p_peak instead.
%
% Usage:
%   out = amp_xcorr_peak_sig(A, B)
%   out = amp_xcorr_peak_sig(A, B, 'fs', 512, 'maxLagSec', 0.2, 'nPerm', 1000)
%
% Inputs:
%   A, B - equal-length column/row vectors (amplitude envelopes / time courses)
%
% Name-Value options:
%   'fs'          - sampling rate in Hz (default 512). Used only to convert
%                   lag indices to seconds. Pass [] to skip sec conversion.
%   'maxLagSec'   - max |lag| in seconds (default 0.2). Ignored if maxLag
%                   samples is provided.
%   'maxLag'      - max |lag| in samples (overrides maxLagSec when set).
%   'nPerm'       - number of permutations (default 1000)
%   'alpha'       - one-sided percentile threshold (default 0.05 -> 95th pct)
%   'seed'        - rng seed for reproducibility ([] = do not reseed)
%   'normalize'   - xcorr normalization (default 'coeff'); use 'none' to
%                   match a raw-dot-product style if desired
%
% Output struct fields:
%   lags, xcorr           - lag axis (samples) and cross-correlation
%   lagsSec               - lag axis in seconds (empty if fs not given)
%   peakValue, peakLag, peakLagSec
%   r0, p0                - corrcoef(A,B) at lag 0 (descriptive only)
%   nullPeaks             - nPerm x 1 peak values under permutation
%   nullThresh            - (1-alpha) percentile of nullPeaks
%   peakSignificant       - peakValue > nullThresh
%   p_peak                - (1 + #nullPeaks >= peakValue) / (1 + nPerm)

p = inputParser;
p.addRequired('A', @(x) isnumeric(x) && isvector(x));
p.addRequired('B', @(x) isnumeric(x) && isvector(x));
p.addParameter('fs', 512, @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('maxLagSec', 0.2, @(x) isnumeric(x) && isscalar(x) && x >= 0);
p.addParameter('maxLag', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
p.addParameter('nPerm', 1000, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('alpha', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
p.addParameter('normalize', 'coeff', @(x) ischar(x) || isstring(x));
p.parse(A, B, varargin{:});
opt = p.Results;

A = A(:);
B = B(:);
if numel(A) ~= numel(B)
    error('amp_xcorr_peak_sig:LengthMismatch', 'A and B must have the same length.');
end
if numel(A) < 3
    error('amp_xcorr_peak_sig:TooShort', 'A and B must each have at least 3 samples.');
end

if ~isempty(opt.seed)
    rng(opt.seed);
end

if ~isempty(opt.maxLag)
    maxLag = round(opt.maxLag);
else
    if isempty(opt.fs)
        error('amp_xcorr_peak_sig:NeedFs', ...
            'Provide ''fs'' when specifying maxLagSec, or pass ''maxLag'' in samples.');
    end
    maxLag = round(opt.maxLagSec * opt.fs);
end
maxLag = min(maxLag, numel(A) - 1);

% Mean-subtract (Adhikari et al., 2010; manuscript Methods)
A0 = A - nanmean(A);
B0 = B - nanmean(B);

% Replace NaNs after centering so xcorr stays defined
A0(~isfinite(A0)) = 0;
B0(~isfinite(B0)) = 0;

[xc, lags] = local_xcorr(A0, B0, maxLag, opt.normalize);
[peakValue, peakIdx] = max(xc);
peakLag = lags(peakIdx);

% Lag-0 Pearson correlation (manuscript-reported descriptive p; not peak test)
[R0, P0] = corrcoef(A0, B0);
r0 = R0(1, 2);
p0 = P0(1, 2);

% Permutation null of cross-correlation peaks
nPerm = round(opt.nPerm);
nullPeaks = nan(nPerm, 1);
n = numel(B0);
for iPerm = 1:nPerm
    Bp = B0(randperm(n));
    xcPerm = local_xcorr(A0, Bp, maxLag, opt.normalize);
    nullPeaks(iPerm) = max(xcPerm);
end

nullThresh = prctile(nullPeaks, 100 * (1 - opt.alpha));
peakSignificant = peakValue > nullThresh;
p_peak = (1 + sum(nullPeaks >= peakValue)) / (1 + nPerm);

out = struct();
out.lags = lags(:);
out.xcorr = xc(:);
if isempty(opt.fs)
    out.lagsSec = [];
    out.peakLagSec = [];
else
    out.lagsSec = lags(:) / opt.fs;
    out.peakLagSec = peakLag / opt.fs;
end
out.peakValue = peakValue;
out.peakLag = peakLag;
out.r0 = r0;
out.p0 = p0;
out.nullPeaks = nullPeaks;
out.nullThresh = nullThresh;
out.peakSignificant = peakSignificant;
out.p_peak = p_peak;
out.nPerm = nPerm;
out.alpha = opt.alpha;
out.maxLag = maxLag;
out.normalize = char(opt.normalize);
end


function [xc, lags] = local_xcorr(A, B, maxLag, normalize)
%LOCAL_XCORR Prefer MATLAB Signal Toolbox xcorr; fall back to conv-based impl.
normalize = char(normalize);
lags = (-maxLag:maxLag).';
try
    [xc, lags] = xcorr(A, B, maxLag, normalize);
    xc = xc(:);
    lags = lags(:);
    return
catch
    % Fallbacks for environments without Signal Toolbox / Octave differences
end

% Raw cross-correlation via convolution
xcFull = conv(A, flipud(B), 'full');
mid = numel(A);
xc = xcFull(mid - maxLag:mid + maxLag);

switch lower(normalize)
    case 'none'
        % leave as-is
    case 'coeff'
        denom = sqrt(sum(A.^2) * sum(B.^2));
        if denom == 0
            xc = zeros(size(xc));
        else
            xc = xc / denom;
        end
    case 'biased'
        xc = xc / numel(A);
    case 'unbiased'
        lagsTmp = (-maxLag:maxLag).';
        xc = xc ./ (numel(A) - abs(lagsTmp));
    otherwise
        error('amp_xcorr_peak_sig:Normalize', ...
            'Unsupported normalize flag: %s', normalize);
end
xc = xc(:);
end
