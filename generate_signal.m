function [noised_sig, atom, sig] = generate_signal(N, coord, SNR)

% -------------------------------------------------------------------------
%   generate_signal: generates signals as a summmation of atoms using hann window.
%
%   args:   - N:     the signal length
%           - coord: [t1, f1, T1, A1; t2, f2, T2, A2; ...],
%                    (t1, f1): time-frequency coordinates of atom 1, 
%                     T1: duration of atom 1,
%                     A1: amplitude of atom 1.
%   returns: - noised_sig: the generated signal
%            - atom:       a matrix that contains the all atoms exist in the final signal
%            - sig:        the generated signal without the noise. 
%
%   developer: Apostolidis Georgios
% -------------------------------------------------------------------------

sig = zeros(N, 1);

[Natoms, ccoord] = size(coord);
if (ccoord~=4)
  error('Bad dimension for coord');
end

atom = zeros(N, Natoms);

for k = 1:Natoms
    
    t0 = round(max(min(coord(k, 1), N), 1));
    f0 = max(min(coord(k, 2), pi), 0.0);
    T  = coord(k, 3); 
    A  = coord(k, 4);
    
    win = zeros(2*N, 1);
    win(ceil(N/2)+ceil(t0-T/2+1):ceil(N/2)+ceil(t0+T/2)) = hann(T, 'periodic');
    win = win(ceil(N/2)+1:ceil(3*N/2));
   
    n = (1:N).' - t0;
    fm = cos(f0 * n);
    
    atom(:, k) = A * win .* fm;
    sig = sig + atom(:, k); 
    
end

 noised_sig = addWhiteNoise(sig, SNR);

function [output, noise] = addWhiteNoise(input, snr)

% -------------------------------------------------------------------------
%   addWhiteNoise: add white gaussian noise - awgn function
%   args:   input-> input signal
%           snr-> signal to noise ratio
%   returns:output-> output signal, added with noise
%
%   validated with built-in MATLAB method awgn().
% -------------------------------------------------------------------------

input_power = var(input);
linear_noise = 10 ^ (snr / 10);
noise_power = input_power / linear_noise;
noise_std = sqrt(noise_power);

noise = randn(size(input)) * noise_std;

output = input + noise;

