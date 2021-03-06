function [ p ] = gpode( U, m )
%GPODE Oversampling empirical interpolation (gappy POD)
% In
%   U       ...     basis matrix
%   m       ...     number of sampling points
% Out
%   p       ...     sampling points
%%%
% Reference:
% Peherstorfer, B., Drmac, Z. & Gugercin, S. Stability of discrete 
% empirical interpolation and gappy proper orthogonal decomposition with 
% randomized and deterministic sampling points. arXiv:1808.10473, 2018.
%%%

assert(m >= size(U, 2) && m <= size(U, 1));
% compute the first n = size(U, 2) points with QDEIM
[~, ~, p] = qr(U', 'vector');
p = p(1:size(U, 2))';

% select points n, ..., m
for i=length(p)+1:m
    [~, S, W] = svd(U(p, :), 0);
    g = S(end-1, end-1).^2 - S(end, end)^2; % eigengap
    Ub = W'*U';
    r = g + sum(Ub.^2, 1);
    r = r - sqrt((g + sum(Ub.^2, 1)).^2 - 4*g*Ub(end, :).^2);
    [~, I] = sort(r, 'descend');
    e = 1;
    while any(I(e) == p)
        e = e + 1;
    end
    p(end + 1) = I(e);
end

end

