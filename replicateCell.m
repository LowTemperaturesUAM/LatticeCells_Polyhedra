function [R] = replicateCell(oldR,n,mvs)
arguments
    oldR
    n
    mvs
end
% Replicate a set of 3D points in space n times in directions set by mvs.
R = [];
% Loop in 1st direction
for i = -n:n
    R = uniquetol([R; oldR + i*mvs(3,:)],1e-10,'ByRows',true);
end
% 
Rold = R;
for i = -n:n
    R = uniquetol([R; Rold + i*mvs(2,:)],1e-10,'ByRows',true);
end
Rold=R;
for i = -n:n
    R = uniquetol([R; Rold + i*mvs(1,:)],1e-10,'ByRows',true);
end

end






