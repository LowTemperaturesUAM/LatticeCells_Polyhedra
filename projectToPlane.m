function [newP_3D,newP_2D] = projectToPlane(points,n,d)
arguments
    points 
    n (1,3)
    d (1,1) = 0
end
% Projects 3D points into a plane. It gives the coordinates of the
% projected points in the original 3D coordinates and in a 2D coordinates
% of plane (to be finished)

% Projection operator
P = n'*n/vecnorm(n).^2;
P = eye(3)-P;

% Project
vProj = points*P + d;
vProj = uniquetol(vProj,1e-8,'ByRows',true); % remove duplicates

newP_3D = vProj;

end