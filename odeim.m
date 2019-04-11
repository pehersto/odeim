function [ p ] = odeim( U, m )
%ODEIM Oversampling empirical interpolation
% In
%   U       ...     basis matrix
%   m       ...     number of sampling points
% Out
%   p       ...     sampling points
%%%
% Reference:
% Peherstorfer, B., Drmac, Z. & Gugercin, S. Stabilizing discrete empirical interpolation via 
% randomized and deterministic oversampling. arXiv:1808.10473, 2018.
%%%

assert(m >= size(U, 2));
% compute the first n = size(U, 2) points with QDEIM
[~, ~, p] = qr(U', 'vector');
p = p(1:size(U, 2))';

% select points n, ..., m with ODEIM
for i=length(p)+1:m
    [~, ~, W] = svd(U(p, :), 0);
    r = (W(:, end)'*U').^2;
    [~, I] = sort(r, 'descend');
    g = 1;
    while any(I(g) == p)
        g = g + 1;
    end
    p(end + 1) = I(g);
end

end

