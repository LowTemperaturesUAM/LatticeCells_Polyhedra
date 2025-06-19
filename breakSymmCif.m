function [newR] = breakSymmCif(cifData)
% From a cif struct with valid fields, extract atom coordinates and applies
% all space group symmetry operations in order to obtain equivalent atoms.

% Get together ion coordinates and labels
rIons = table(cifData.atom_site_label,cifData.atom_site_fract_x, ...
    cifData.atom_site_fract_y,cifData.atom_site_fract_z, ...
    'VariableNames',{'Label','X','Y','Z'});

% number of ions
nIons = size(cifData.atom_site_fract_x,1);

% all ion positions in a matrix
R = [rIons.X rIons.Y rIons.Z];

% get symmetry operations (strings)
symOps = cifData.space_group_symop_operation_xyz;

newR = zeros(0,3);
% loop over all symm operations
for i = 1:numel(symOps)
    t = symOps{i};
    %find(t==',')

    % Convert operation into affine transformation
    % Get identity or inversion
    refl = (1-2*[contains(t,'-x') contains(t,'-y') contains(t,'-z')]') .* eye(3);
    % Get translation
    trans = str2num(replace(t,{'x','y','z','-x','-y','-z'},'0'))';

    A = [[refl; zeros(1,3)] [trans;1]];
    tform = affinetform3d(A);

    Rtform = R;
    % Obtain ions after symm operations
    [Rtform(:,1), Rtform(:,2), Rtform(:,3)] = transformPointsForward(tform, ...
        rIons.X, rIons.Y, rIons.Z);
    % newR = uniquetol([newR;Rtform],1e-10,'ByRows',true);
    newR = [newR;Rtform];
end

% Separate each ion into a cell
for n = 1:nIons
    R_ion{n} = uniquetol(newR(n:nIons:end,:),1e-10,'ByRows',true);
end
newR = R_ion;
end