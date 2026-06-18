function [newP3D,newP2D,newBase] = projectToPlane(points3D,normal,offset)
arguments
    points3D (:,3) double
    normal (1,3) double
    offset (1,:) double = 0
end
%  [newP3D,newP2D] = projectToPlane(points3D,n) projects 3D points 
%  into a plane perpendicular to n through the origin. It gives the 
%  coordinates of the projected points in the original 3D coordinates and 
% in 2D coordinates of plane. Third output will be struct with chosen 2D
% vectors that set the coordinates and an offset to displace the plane.
% Optional argument offset allows to place the plane somewhere else. If
% offset is a scales, it will be moved that distance along the normal
% vector. If it is a 3D point, the plane will pas through it

% Normalize vector
nNorm = normal/vecnorm(normal);

% Projection operator (to plane through 0)
P = eye(3)-nNorm'*nNorm;

% Project coordinates
vProj = points3D*P;
vProj = uniquetol(vProj,1e-5,'ByRows',true); % remove duplicates

% Apply offset to plane, if given.
switch size(offset,2)
    case 1 % Move distance offset along normal
        offsetFinal = offset*nNorm;
    case 3 % Move plane to pass through point offsset
        offsetFinal = (nNorm*offset') * nNorm;
         
end
vProj = vProj + offsetFinal;

newP3D = vProj;

if nargout > 1
%Turn to 2D plane coordinates
if size(vProj,1)>2
    v1 = vProj(3,:) - vProj(2,:); % Choose in-plane vector
    v1 = v1/norm(v1); % normalize
    v2 = cross(nNorm,v1); % Use in-plane vector, perpendicular to normal and v1
else % If there are no more than 2 points, choose rand vectors perp to normal
    v1 = null(nNorm);
    v2 = v1(:,2)';
    v1 = v1(:,1)';
end
newBase.vectors = [v1;v2];
newBase.offset = offsetFinal;

newP2D = newP3D * [v1;v2]'; % Change coordinates

end
end