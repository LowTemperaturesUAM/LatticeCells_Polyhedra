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
        function [newT,newR] = breakSymm(obj,Ncell)
            arguments
                obj cif
                Ncell double = 1
            end
            % [T,R] = breakSymm(obj) gives all atom positions, with no
            %  implied symmetry in cell coordinates. T is a table 
            % with atom labels and multiplicity, while R is a cell,
            %  with each element in a cell
            data = obj.Data;
            newR = breakSymmCif(data);

            if numel(Ncell) ~= 3
                Ncell = Ncell(1).*ones(1,3);
            end
            [Na,Nb,Nc] = struct('x',num2cell(Ncell)).x;

            % Limit atom positions to Ncell number of cells in each
            % direction
            R = cellfun(@(x) replicateCell(x,max([Na,Nb,Nc]*2),eye(3)), ...
                newR,'UniformOutput',false);
            % Limit to N cells
            cond = @(x) x(:,1) > Na | x(:,2) > Nb | x(:,3) > Nc | any(x<0,2);
            for n = 1:numel(R)
                R{n}(cond(R{n}),:) = [];
            end
            newR = R;

            % construct new table with all atoms-------------------------
            T = obj.atoms; % symmetrized positions in a table
            
            N = height(T); % Different equivalent atoms
            nAtoms = cellfun(@(x)size(x,1), newR); %number af all atoms

            newT = table();
            for n = 1:N
                % Make rows with new positions
                newPos = [repmat(T(n,1:end-3),[nAtoms(n) 1]), num2cell(newR{n})];
                % Preserve variable names
                newPos.Properties.VariableNames=T.Properties.VariableNames;

                newT = unique([newT; T(n,:);newPos],'stable');
            end
            
        end
        function l = lengths(obj)
            % obj.lengths gives the three lattice parameters
            l = [obj.Data.cell_length_a,...
                obj.Data.cell_length_b,...
                obj.Data.cell_length_c];
        end
        function a = angles(obj)
            % obj.angles gives the three lattice angles between vectors a,b
            % and c. alpha=angle(b,c), beta=angle(c,a) and gamma=angle(a,b)
            a = [obj.Data.cell_angle_alpha,...
                obj.Data.cell_angle_beta,...
                obj.Data.cell_angle_gamma];
        end
        function latName = bravais(obj)
            % obj.bravais gives the Bravais lattice of the cif file
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
   
            % Compose lattice full name
            if isempty(name)
                latName = system;
            else
                latName = [name ' ' system];
            end

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