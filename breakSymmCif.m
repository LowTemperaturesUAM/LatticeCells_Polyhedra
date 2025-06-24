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
if isfield(cifData,'space_group_symop_operation_xyz')
    symOps = cifData.space_group_symop_operation_xyz;
elseif isfield(cifData,'symmetry_equiv_pos_as_xyz')
    symOps = cifData.symmetry_equiv_pos_as_xyz;
else
    error('There is not known field with symmetry operations')
end

newR = zeros(0,3);
% loop over all symm operations
for i = 1:numel(symOps)
    t = symOps{i}; % take one operation
    t = replace(t,' ',''); % delete empty spaces
    sepDim = find(t==','); % get comma position to separate directions
    sepDim = [0 sepDim length(t)+1];


    % Convert operation into affine transformation
    refl = nan(3); % Initialize matrix
    % obtain matrix rows for reflection
    % Separate string into each dimension for easier parsing
    for n = 1:3
        subt = t(sepDim(n)+1:sepDim(n+1)-1);
        % obtain 1 or -1 for each correspondance, this should give a matrix
        % corresponding to the equivalent positions, outside of traslations
        a = [contains(subt,'x')-2*contains(subt,'-x'),...
            contains(subt,'y')-2*contains(subt,'-y'),...
            contains(subt,'z')-2*contains(subt,'-z')];
        refl(n,:) = a;
    end
    
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
    R_ion{n} = uniquetol(newR(n:nIons:end,:),1e-5,'ByRows',true);
end
newR = R_ion;
end