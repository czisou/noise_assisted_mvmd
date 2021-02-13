function [io] = orthogonality_index(components, initial_signal)

% -------------------------------------------------------------------------
%   orthogonality_index(): computes the index of orthogonality 
%   args:
%       - components: a matrix with the resulted components from a decomposition
%       algorithm. each component is contained in each column. 
% -------------------------------------------------------------------------

[n, m] = size(components);
if n < m
    components = components';
end

no_components = size(components, 2);
initial_energy = sum(initial_signal.^2);
partial_sum = 0;

for i = 1:no_components
    for j = 1:no_components
        if (i ~= j)
            ci = components(:, i);
            cj = components(:, j);
            temp_sum = (ci.' * cj) / initial_energy;
            partial_sum = partial_sum + temp_sum;   
        end     
    end
end

coeff = no_components * (no_components - 1) / 2;
io = abs(partial_sum) / coeff;

