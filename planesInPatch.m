function [planeCoeffs] = planesInPatch(Vertices)
arguments
    Vertices (:,3) double
end
% Given the vertices of a convex polyhedra (check with patch), gives the
% plane equations for all faces, with the normal vectors pointing outward.
% Output is a matrix with 4 columns. One for each parameter of the eq:
% Ax+By+Cz+D=0;

DT = delaunayTriangulation(Vertices);

% Recalculate Convex Polyhedra
[H] = convexHull(DT);

% Analysis of the convex Hull - Triangulate and obtain normals
DT = triangulation(H,DT.Points);
n = faceNormal(DT); % normals to the triangles of hull (with duplicates)
[n,idxU] = uniquetol(n,1e-10,'ByRows',true); % Remove duplicates

% Obtain 1 triangle for each plane
uniqueTri = DT.ConnectivityList(idxU,:);

% Obtain plane equations from normals and a point of the triangle
d = -sum(n.*Vertices(uniqueTri(:,1),:),2);

planeCoeffs = [n d];

end