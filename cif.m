classdef cif
    properties (SetAccess = public)
        Data
        % Lattice
    end
    methods
        % constructor function
        function obj = cif(path)
            if nargin > 0
                obj.Data = readcif(path);
            end
        end
        % break space group symmmetry
        function [obj,newR] = breakSymm(obj)
            data = obj.Data;
            newR = breakSymmCif(data);

            % construct new fields with all atoms
            labels = data.atom_site_label;
            % data.atom_site_x = [];
        end
        function l = lengths(obj)
            l = [obj.Data.cell_length_a,...
                obj.Data.cell_length_b,...
                obj.Data.cell_length_c];
        end
        function a = angles(obj)
            a = [obj.Data.cell_angle_alpha,...
                obj.Data.cell_angle_beta,...
                obj.Data.cell_angle_gamma];
        end
        function latName = bravais(obj)
            % If it contains corresp. field, take it directly
            
            % Otherwise, obtain it from space group
            % From number, obtain the crystal system
            groupNum = obj.Data.symmetry_Int_Tables_number;
            % From name first letter, get the structure (primitive, centered...)
            letter = obj.Data.symmetry_space_group_name_H_M;
        aux = [contains(letter,"P"),...
contains(letter,"F"),contains(letter,"I"),contains(letter,"R"),...
contains(letter,"C")||contains(letter,"A")];


            if groupNum > 194 && groupNum <= 230
                system = 'Cubic';
            elseif groupNum > 167
                system = 'Hexagonal';
            elseif groupNum > 142
                system = 'Trigonal';
            elseif groupNum > 74
                system = 'Tetragonal';
            elseif groupNum > 15
                system = 'Orthorhombic';
            elseif groupNum> 2
                system = 'Monoclinic';
            elseif groupNum > 0
                system = 'Triclinic';
            end

            switch find(aux)
                case 1 % simple
                    name = '';
                    if contains(system,'Trigonal')
                        system = 'Hexagonal';
                    end
                case 2 % Face centered
                    name = 'Face Centered';
                    if contains(system,'Cubic')
                        name = 'FCC';
                    end
                case 3 % Body centered
                    name = 'Body Centered';
                    if contains(system,'Cubic')
                        name = 'BCC';
                    end
                case 4 % rhombohedral
                    system = 'Rhombohedral';
                case 5 % side centered
            end
            
latName = [name ' ' system];

        end
        function latticeData = lattice(obj)
            
            latticeData = createBravaisLattice(obj.bravais, ...
                obj.lengths, ...
                obj.angles);

            % obj.Lattice = latticeData;
        end
        function table_atoms = atoms(obj)
            data = obj.Data;
            
            table_atoms = ...
table(data.atom_site_label,data.atom_site_type_symbol, ...
                data.atom_site_symmetry_multiplicity, ...
                data.atom_site_fract_x,data.atom_site_fract_y,data.atom_site_fract_z, ...
                'VariableNames',{'Label','Symbol','Multiplicity', ...
                'X','Y','Z'});
        end
    end

end