function [Faces, Vol] = facesPatch3D(Vertices)
arguments
    Vertices (:,3) double
end
% Finds the Connecting map for each face of the Convex polyhedra defined by
% the given points, as is used for plotting patches.
% Second output is the enclosed volume. Normal vectors point outward the polyhedra.

% end function if there is no Vertices. Assign 0 volume.
if isempty(Vertices) | all(isnan(Vertices))
    Faces = nan;
    Vol = 0;
    return
end

% Find connecting lines between vertex in same plane-------------------

DT = delaunayTriangulation(Vertices);
[F,Vol] = convexHull(DT);
tt = triangulation(F,DT.Points);

n_tri = faceNormal(tt); % obtain normal vectors from every triangle
[n_tri,~,uIdx] = uniquetol(n_tri,1e-6,'ByRows',true); % remove duplicates


Faces = nan(max(uIdx),...
        numel(unique(F(uIdx==(mode(uIdx,1)),:))) ...
        ); % Initialize variable with nans

for i = 1:max(uIdx) % Loop over all planes
    idV = F(uIdx==i,:);
    idV = unique(idV); 
    if isrow(idV)
        idV = idV';
    end
    V = Vertices(idV,:);% Vertex in that face

    % Convert Vertex to 2d coordinates. Define basis
    if isempty(V)
        warning('Some of the points do not create bisector plane')
    else
    uz = n_tri(i,:)/vecnorm(n_tri(i,:)); % normal vector
    ux = (V(1,:)-V(2,:))/vecnorm(V(1,:)-V(2,:)); % vector inside the plane
    uy = cross(uz,ux); % third vector is the cross product

    V2d =(V-mean(V)) * [ux;uy]';

    % Sort vertices by angle, so they are consecutive
    a = atan2d(V2d(:,2),V2d(:,1));
    [~,id] = sort(a);
    % use sort to
    Faces(i,1:numel(id)) = idV(id,:)';
    end
end

end
