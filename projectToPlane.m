function [newP3D,newP2D] = projectToPlane(points3D,normal,offset)
arguments
    points3D (:,3) double
    normal (1,3) double
    offset (1,1) double = 0
end
%  [newP3D,newP2D] = projectToPlane(points3D,n) projects 3D points 
%  into a plane perpendicular to n through the origin. It gives the 
%  coordinates of the projected points in the original 3D coordinates and 
% in 2D coordinates of plane. Optional argument offset allows to place the plane
% somewhere else.

% Normalize vector
nNorm = normal/vecnorm(normal);

% Projection operator (to plane through 0)
P = eye(3)-nNorm'*nNorm;

% Project coordinates
vProj = points3D*P + offset;
vProj = uniquetol(vProj,1e-5,'ByRows',true); % remove duplicates

newP3D = vProj;

%Turn to 2D plane coordinates
v1 = vProj(3,:) - vProj(2,:); % Choose in-plane vector
v1 = v1/norm(v1); % normalize
v2 = cross(nNorm,v1); % Use in-plane vector, perpendicular to normal and v1

newP2D = newP3D * [v1;v2]'; % Change coordinates

end