function [latticeData,params,angles]=createBravaisLattice(bravais,params,angles)
arguments
    bravais
    params (1,3) double {mustBePositive,mustBeFinite}
    angles (1,3) double {mustBeFinite}
end

bravaisLongList = {'Cubic','FCC','BCC','Tetragonal','Body Centered Tetragonal' ...
    'Orthorhombic','Face Centered Orthorhombic', ...
    'Body Centered Orthorhombic','Side Centered Orthorhombic' ...
    'Hexagonal','Rhombohedral','Monoclinic','Side Centered Monoclinic', ...
    'Triclinic'};
bravaisList = bravaisLongList; % default value

bravaisShortList = {'cP'
'cF'
'cI'
'tP'
'tI'
'oP'
'oF'
'oI'
'oS'
'hP'
'hR'
'mP'
'mS'
'aP'}';

switch class(bravais)
    case 'double'
        warning('Only Lattice Systems are known (simple or primitive).')
        if bravais > 194 && bravais <= 230
            bravais = 'Cubic';
        elseif bravais > 167
            bravais = 'Hexagonal';
        elseif bravais > 142
            bravais = 'Trigonal';
        elseif bravais > 74
            bravais = 'Tetragonal';
        elseif bravais > 15
            bravais = 'Orthorhombic';
        elseif bravais> 2
            bravais = 'Monoclinic';
        elseif bravais > 0
            bravais = 'Triclinic';
        end
        % Check lattice type
% clean string. Consider input as abreviations (2 chars)
    case 'string'
        % clean up the string
        if length(char(bravais)) < 3
            bravaisList = bravaisShortList;
        end
        bravais = validatestring(strrep(bravais,'.',''), bravaisList);
    case 'char'
        if length(bravais) < 3
            bravaisList = bravaisShortList;
        end
    bravais = validatestring(strrep(bravais,'.',''), bravaisList);
    otherwise
        error("1st input is not a correct Bravais lattice. Use Group No. or Full Lattice Name");
end

switch bravais
    case bravaisList(1) % cubic
        [alpha,beta,gamma] = deal(90);
        [params(2:3)] = deal(params(1));
        R = eye(3).*params;

    case bravaisList(2) % Face centered Cubic
        [alpha,beta,gamma] = deal(90);
        [params(2:3)] = deal(params(1));

        R = (ones(3)-eye(3))/2.*params;

    case bravaisList(3) % Body Centered Cube
        % Assign default values to not used parameter. (disable them)
        [alpha,beta,gamma] = deal(90);
        [params(2:3)] = deal(params(1));

        R = [-1 1 1;1 -1 1;1 1 -1]/2.*params;

    case bravaisList(4) % Simple tetragonal
        params(2) = params(1);
        [alpha,beta,gamma] = deal(90);

        R = eye(3).*params;

    case bravaisList(5) % body Centered tetragonal
        params(2) = params(1);
        [alpha,beta,gamma] = deal(90);

        R = (ones(3)-2*eye(3))/2.*params;

    case bravaisList(6) % Simple Orthorhombic
        [alpha,beta,gamma] = deal(90);

        R = eye(3).*params;

    case bravaisList(7) % Face Centered orthorhombic
        [alpha,beta,gamma] = deal(90);

        R = (ones(3)-eye(3))/2.*params;

    case bravaisList(8) % Body Centered orthorhombic
        % Assign default values to not used parameter. (disable them)
        [alpha,beta,gamma] = deal(90);

        R = (ones(3)-2*eye(3))/2.*params;

    case bravaisList(9) % Side Centered orthorhombic
        % Assign default values to not used parameter. (disable them)
        [alpha,beta,gamma] = deal(90);

        R = [1 -1 0; 1 1 0; 0 0 2]/2.*params;

    case bravaisList(10) % Hexagonal
        [alpha,beta] = deal(90);
        gamma = 120;
        params(2) = params(1);

        % R = [-cosd(gamma) -sind(gamma) 0;...
        %     -cosd(gamma) sind(gamma) 0;...
        %     0 0 1].*latticeParam;
    warning('Not available')
        %DO%
    case bravaisList(11) % Rombohedral
        [beta,gamma] = deal(angles(1));
        [params(2:3)] = deal(params(1));
        warning('Not available')
        %DO%

    case bravaisList(12) % Monoclinic
        [beta,gamma] = deal(90);
        if beta>90
            warning("Check angle values. All should be 90º or less")
        end
        warning('Not available')
        %DO%
    case bravaisList(13) % Side Centered Monoclinic
        warning('Not available')
        %DO%
    case bravaisList(14) % Triclinic
        warning('Not available')
        %DO%
    otherwise
        warning('bravaisList is not correct')
end

% useful ID for bravais lattice name
bravaisID = contains(bravaisList,bravais);
LongName = bravaisLongList(bravaisID);

% general set of vectors with a//x
if ~exist('alpha','var')
    alpha = angles(1);
end
if ~exist('beta','var')
    beta = angles(2);
end
if ~exist('gamma','var')
    gamma = angles(3);
end

Rconv = [1 0 0;
    cosd(gamma),sind(gamma),0;
    cosd(beta), 1/sind(gamma)*(cosd(alpha)-cosd(beta)*cosd(gamma)),...
    1/sind(gamma)*sqrt(sind(gamma)^2-cosd(alpha)^2-cosd(beta)^2+...
    2*cosd(alpha)*cosd(beta)*cosd(gamma))].*params';

% Try cell in real space
Pr = [0 0 0;Rconv; sum(Rconv([1 2],:)); sum(Rconv([3 2],:)); sum(Rconv([1 3],:)); sum(Rconv)];
Fr = facesPatch3D(Pr);

% If no R was calculated, use convemtional.
if ~exist("R",'var')
    R = Rconv;
end
latticeData.realLatticeVectors = R;
% Calculate reciprocal vectors (primitive). A 2pi could be added
Vr = dot(R(1,:),cross(R(2,:),R(3,:)));
K = 1/Vr .* [cross(R(2,:),R(3,:));...
    cross(R(3,:),R(1,:));
    cross(R(1,:),R(2,:))];

% Calculate several lattice points
fullK = replicateCell(K,2,K);
fullR = replicateCell(R,2,R);

% Brillouin zone
[Pk,Fk] = wignerSeitz3D(fullK,[ 0 0 0],'Method','voronoi');

% Save data into struct
latticeData = struct('Name',LongName,'Real',fullR,'realCellVertices',Pr, ...
    'realCellFaces',Fr,'realLatticeVector',R,'Reciprocal',fullK, ...
    'reciprCellVertices',Pk,'reciprCellFaces',Fk,'reciprLatticeVector',K);

end