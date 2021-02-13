function [final_modes, mode_energies] = post_processing(modes, signal, perc)

% -------------------------------------------------------------------------
%   post_processing: This function eliminates low energy modes based on a
%                    relative energy threshold
%
%   args:   - modes: Extracted modes
%           - signal: Input univariate signal
%           - perc: Threshold percentage
%
%   returns: - final_modes: Post processed modes
%            - mode_energies: Mode energies
%
%   developers: Charilaos Zisou, Apostolidis Georgios
% -------------------------------------------------------------------------

% Initialize constants
num_cmps = size(modes, 2);
signal_energy = sum(abs(signal.^2));
counter = 1;

% Processing
for i=1:num_cmps
    comp_energy = sum(abs(modes(:,i).^2));
    mode_energies(i) = comp_energy;
    if comp_energy >= perc*signal_energy
        final_modes(:, counter) = modes(:, i);
        counter = counter + 1;
    end
end


