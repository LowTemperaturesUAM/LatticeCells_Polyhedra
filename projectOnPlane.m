function pointsProj = projectOnPlane(points3D,normal,offset)
arguments
    points3D (:,3) double 
    normal (1,3) double
    offset (1,1) double = 0
end
% pointsProj = projectOnPlane(points,n) projects the 3D space coordinates
% points into a plane perpendicular to normal through the origin. Optional
% arguments can be given for the plane to be someewhere else, or through a
% specific point.

% Normalize vector
nNorm = normal/norm(normal);

% Create projection operator (to plane through 0)
P = eye(3) - nNorm'*nNorm;

% Project coordinates into Plane
pointsProj = points3D*P;


end