function[Type] = ReadingBravais(Name)

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

switch class(Name)
        % Check lattice type
% clean string. Consider input as abreviations (2 chars)
    case 'string'
        % clean up the string
        if length(char(Name)) < 3
            bravaisList = bravaisShortList;
        end
        Name = validatestring(strrep(Name,'.',''), bravaisList);
    case 'char'
        if length(Name) < 3
            bravaisList = bravaisShortList;
        end
    Name = validatestring(strrep(Name,'.',''), bravaisList);
    otherwise
        error("1st input is not a correct Bravais lattice. Use Group No. or Full Lattice Name");
end

switch Name
    case bravaisList(1) % cubic
        Type = bravaisShortList(1);

    case bravaisList(2) % Face centered Cubic
        Type = bravaisShortList(2);

    case bravaisList(3) % Body Centered Cube

        Type = bravaisShortList(3);
    case bravaisList(4) % Simple tetragonal

        Type = bravaisShortList(4);

    case bravaisList(5) % body Centered tetragonal

        Type = bravaisShortList(5);

    case bravaisList(6) % Simple Orthorhombic

        Type = bravaisShortList(6);

    case bravaisList(7) % Face Centered orthorhombic
        
        Type = bravaisShortList(7);

    case bravaisList(8) % Body Centered orthorhombic
        
        Type = bravaisShortList(8);

    case bravaisList(9) % Side Centered orthorhombic
        
        Type = bravaisShortList(9);

    case bravaisList(10) % Hexagonal

        Type = bravaisShortList(10);
        warning('Not available')
        %DO%
    case bravaisList(11) % Rombohedral
        Type = bravaisShortList(11);
        warning('Not available')
        %DO%

    case bravaisList(12) % Monoclinic
        Type = bravaisShortList(12);
        warning('Not available')
        %DO%
    case bravaisList(13) % Side Centered Monoclinic
        warning('Not available')
        Type = bravaisShortList(13);
        %DO%
        % I need to get primitive vectors for a side centered lattice
        %R = [1 -1 0; 1 1 0; 0 0 2]/2.*params;

    case bravaisList(14) % Triclinic
        Type = bravaisShortList(14);
        warning('Not available')
        %DO%
    otherwise
        warning('bravaisList is not correct')
end

% useful ID for bravais lattice name
% No longer used (?)
% bravaisID = contains(bravaisList,bravais);
% LongName = bravaisLongList(bravaisID);
if length(Name) < 3
    LongName = bravaisLongList(contains(bravaisShortList,Name));
else
    LongName = Name;
end

end
