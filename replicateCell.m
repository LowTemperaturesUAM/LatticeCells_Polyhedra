function [R] = replicateCell(oldR,n,mvs)
arguments
    oldR
    n
    mvs (3,:)
end
% replicateCell(oldR,n,mvs) moves a set of 3D points in space n times in directions set by mvs.
% Ideally, it should extend the boundary of the cells shown, increasin g
% the number of cells viewed.
R = [];
% Loop in 1st direction
for i = -n:n
    R = uniquetol([R; oldR + i*mvs(3,:)],1e-6,'ByRows',true);
end
% 
Rold = R;
for i = -n:n
    R = uniquetol([R; Rold + i*mvs(2,:)],1e-6,'ByRows',true);
end
Rold=R;
for i = -n:n
    R = uniquetol([R; Rold + i*mvs(1,:)],1e-6,'ByRows',true);
end

end






