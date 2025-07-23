function [Points] = solvePlaneIntersection(PlaneCoef,opts)
arguments
    PlaneCoef (:,4) double
    opts.InteriorPoints = true % Limit solutions with hemispaces
end
% Given a set of planes, solves the different intersections at points.
% HOWEVER, it assumes that only the hemispace below the plane is of interest
% Input: Matrix where each row are the four coeafficients (A,B,C,D) that 
% define a plane: Ax + By + Cz + D=0

if size(PlaneCoef) < 3
    warning('The function is not prepared for less than 3 planes')
end

% Separate normal vector and offset
n = PlaneCoef(:,1:3);
d = PlaneCoef(:,4);

% Every combination of 3 planes
combin = nchoosek(1:size(PlaneCoef,1),3);
sol = [];
% Loop for each combination
for i = 1:size(combin,1)
    A = n(combin(i,:),:);
    % Solve only when unique solution exists
    if abs(det(A)) > 1e-10
        sol = [sol; (A\(-d(combin(i,:))))'];
    end
end

% Leave only those solutions that satisfy ALL inequations
if opts.InteriorPoints
    valCond = any(sol * n' + d' >1e-13,2);
    sol(valCond,:) = [];
end
% Remove duplicates
Points = uniquetol(sol,1e-10,'ByRows',true);
% if there are no intersection points (empty array) use a nan value
if isempty(Points)
    Points = nan(1,3);
end
end