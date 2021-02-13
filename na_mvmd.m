function [u, u_hat, omega] = na_mvmd(signal, alpha, tau, K, DC, init, tol, noise_channels, noise_amp)

% -------------------------------------------------------------------------
%   na_mvmd: Noise-Assisted Multivariate Variatonal Mode Decomposition
%
%   args:   - signal: the signal length
%           - alpha: bandwidth parameter  
%           - tau: time-step of the dual ascent
%           - K: the number of modes to be recovered
%           - DC: true if the first mode is put and kept at DC
%           - init: 0 = all omegas start at 0
%                   1 = all omegas start uniformly distributed
%                   2 = all omegas initialized randomly
%           - tol: tolerance value for convergence of ADMM
%           - noise_channels: number of injected noise channels
%           - noise_amp: noise amplitude based on input signal's StD
%
%   returns: - u: the collection of decomposed modes
%            - u_hat: mode spectra
%            - omega: estimated mode center-frequencies
%
%   developers: Charilaos Zisou, Apostolidis Georgios
% -------------------------------------------------------------------------

% Set noise amplitude to 80% of input signal's StD.
if nargin < 9
    noise_amp = 0.8;
end
    
% Reshape signal if needed
[x, y] = size(signal);
if x > y
    T = x;
	signal = signal';
else
    T = y;
end

% Transform every channel to unit variance
s = std(signal, 0, 2);
signal = signal ./ s;

% Inject noise
noise = noise_amp * wgn(noise_channels, T, 0);
signal = cat(1, signal, noise);

% Apply MVMD and remove noise channel
[u, u_hat, omega] = MVMD(signal, alpha, tau, K, DC, init, tol);
u = permute(u, [2 3 1]);
u = u(:, 1:end-noise_channels, :);

% Transform back
for k = 1:K
    u(:, :, k) = u(:, :, k) .* s';   
end