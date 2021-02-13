function [success] = is_successful(u, x, signal, thresh)

% -------------------------------------------------------------------------
%   is_successful: Determines if the decomposition was successful or not
%
%   args:   - u: Collection of extraced modes
%           - x: Actual synthetic modes
%           - signal: Original input signal
%           - threhs: Percentage threshold
%
%   returns: - success: 0 or 1
%
%   developers: Charilaos Zisou, Apostolidis Georgios
% -------------------------------------------------------------------------

% Correct input dimensions
[n_u, m_u] = size(u);
[n_x, m_x] = size(x);
if n_u < m_u
    m_u = n_u;
    u = u';
end
if n_x < m_x
    m_x = n_x;
    x = x';
end

% Eliminate low energy modes
[u, ~] = post_processing(u, signal, 0.01);
[~, m_u] = size(u);

% If they don't have the same dimensions, the decomposition 
% can't be successful
if m_u ~= m_x
    success = 0;
else
    % Compute pair-wise correlation coefficients
    for i=1:m_x
        correlations(i) = xcorr(u(:, i), x(:, i), 0, 'coeff');
    end
    
    % Count how many are above the threshold
    if length(correlations(correlations > thresh)) == m_x
        success = 1;
    else
        success = 0;
    end
end

