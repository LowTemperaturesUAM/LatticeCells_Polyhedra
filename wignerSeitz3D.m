function [Vertices,Faces] = wignerSeitz3D(Points,C,opts)
arguments
    Points (:,3) double
    C (1,3) double = [0 0 0]; % center point
    opts.Method {mustBeMember(opts.Method,{'voronoi','planeIntersection'})} = 'planeIntersection'
end
% Given a set of points in 3D space around a center, gives the vertices
% of the Wigner-Seitz Cell. Second output is a Face matrix, indicating 
% which vertex is in which face, used to plot with patch.
% EXAMPLE: [Vertices,Faces] = wignerSeitz3D([0 0 1;0 1 0;1 0 0;
% 0 0 -1;0 -1 0;-1 0 0], [0 0 0])

% Every combination of 3 points
combin = nchoosek(1:size(Points,1),3);
sol = [];
vertexPlanes = [];
% Calculate every coefficient for the plane equation: n(x-m) = 0
n = Points-C; % Normal vector
m = (Points+C)/2; % Midpoint
d = -sum(n.*m,2); % independent factor - position of plane
switch opts.Method
    case 'planeIntersection'
        for i = 1:size(combin,1)
            P = Points(combin(i,:),:);
            N = P-C;
            M = (P+C)/2;
            D = -sum(N.*M,2);

            %MatrixCoefs = n; if determinant is non-zero, there is a single
            %solution
            if abs(det(N)) > 1e-6
                solved = N\(-D);
                sol = [sol; solved'];
                vertexPlanes = [vertexPlanes; combin(i,:)];
            else
            end
        end
        % Leave only those that satisfy ALL inequations
        valCond = any(sol * n' + d' >eps,2);
        sol(valCond,:) = [];
        vertexPlanes(valCond,:) = [];
end
Vertices = sol;

Faces = facesPatch3D(Vertices);

end


