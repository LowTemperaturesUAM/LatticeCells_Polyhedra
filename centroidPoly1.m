function C = centroidPoly1(P)
% Function that calculates the centroid of a polyhedron with vertices P.
%
% Syntax
%   C = centroidPoly1(P)
% Input
%   P - List of vertex coordinates

if any(isnan(P))
    C = nan(1,3);
    return
end
T = delaunayn(P);
n = size(T,1);
W = zeros(n,1);
C=0;
for m = 1:n
    sp = P(T(m,:),:);
    [~,W(m)]=convhull(sp);
    C = C + W(m) * mean(sp);
end
C=C./sum(W);
end