classdef cif
    properties (SetAccess = public)
        Data
        LatticeInfo
        lengths
        angles
    end
    methods
        % constructor function
        function obj = cif(path,lengths,angles)
            if nargin == 1 && contains(path(end-3:end),'.cif')
                disp('creating from cif file')
                obj.Data = readcif(path);

                obj.angles = [obj.Data.cell_angle_alpha,...
                    obj.Data.cell_angle_beta,...
                    obj.Data.cell_angle_gamma];

                obj.lengths = [obj.Data.cell_length_a,...
                    obj.Data.cell_length_b,...
                    obj.Data.cell_length_c];

                obj.LatticeInfo = createBravaisLattice(obj.bravais, ...
                    obj.lengths,obj.angles);

            elseif nargin > 1 % given the lattice data
                bravais = path; % take 1st input as lattice name
                if isempty(angles) % allocate with right angles
                    angles = [90 90 90];
                end
                % create struct with lattice data
                [latticeData,lengths,angles] = createBravaisLattice(bravais,lengths,angles);
                if ~isempty(latticeData) 
                    obj.LatticeInfo = latticeData;
                    obj.lengths = lengths;
                    obj.angles = angles;
                else % in case it fails
                    disp("There are no valid lattice data")
                end
            else
                warning('There are no valid lattice data');
            end
        end
        % break space group symmmetry
        function [newT,newR] = breakSymm(obj,Ncell,opts)
            arguments
                obj cif
                Ncell double = 1
                opts.IncludeBoundary (1,1) logical = true
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

            % Limit to N cells. Include boundary (part of next cell) if
            % indicated
            if opts.IncludeBoundary == false
                cond = @(x) x(:,1) >= Na | x(:,2) >= Nb | x(:,3) >= Nc | any(x<0,2);
            else
                % [Na,Nb,Nc] = struct('x',num2cell(Ncell+1)).x;
                cond = @(x) x(:,1) > Na | x(:,2) > Nb | x(:,3) > Nc | any(x<0,2);
            end

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

        function latName = bravais(obj)
            % obj.bravais gives the Bravais lattice of the cif file
            % If it contains corresp. field, take it directly
            
            % Otherwise, obtain it from space group
            % From number, obtain the crystal system
            if isfield(obj.Data,'symmetry_Int_Tables_number')
                groupNum = obj.Data.symmetry_Int_Tables_number;
            elseif isfield(obj.Data,'space_group_IT_number')
                groupNum = obj.Data.space_group_IT_number;
            else
                error('Group number was not found.')
            end
            % From name first letter, get the structure (primitive, centered...)
            if isfield(obj.Data,'symmetry_space_group_name_H_M')
                letter = obj.Data.symmetry_space_group_name_H_M;
            elseif isfield(obj.Data,'space_group_name_H_M_alt')
                letter = obj.Data.space_group_name_H_M_alt;
            end
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
                otherwise % side centered
                    name = 'Side Centered';
            end
   
            % Compose lattice full name
            if isempty(name)
                latName = system;
            else
                latName = [name ' ' system];
            end

        end

        function table_atoms = atoms(obj)
            data = obj.Data;
            
            table_atoms = ...
table(data.atom_site_label,data.atom_site_type_symbol, ...
                data.atom_site_fract_x,data.atom_site_fract_y,data.atom_site_fract_z, ...
                'VariableNames',{'Label','Symbol', ...
                'X','Y','Z'});
        end
        function [cellPoly] = cellPatch(obj,rkSpace,cellType)
            arguments
                obj
                rkSpace {mustBeMember(rkSpace,{'real','reciprocal'})}= 'real'
                cellType {mustBeMember(cellType,{'conventional','primitive','brillouin'})}= 'conventional'
            end
            % returns a struct with vertices and faces for the cell

            switch rkSpace % choose bewtween cell in r or k - space
                case 'real'
                    if (contains(obj.LatticeInfo.Name,'Centered')|| ...
                            contains(obj.bravais,'Centered')) &&...
                        contains(cellType,'conventional')
                        verts = [0 0 0;
                            1 1 1;
                            perms([1 0 0]);
                            perms([1 1 0])].*obj.lengths;
                    elseif contains(cellType,'brillouin')
                        verts = wignerSeitz3D(obj.LatticeInfo.Real,[0 0 0], ...
                            'Output','struct','Method','voronoi');
                        verts = verts.vertices;

                    else
                        verts = [0 0 0;
                            1 1 1;
                            perms([1 0 0]);
                            perms([1 1 0])]*obj.LatticeInfo.realLatticeVector;
                    end
                case 'reciprocal'
                    if contains(cellType,'primitive')
                    verts = [0 0 0;
                            1 1 1;
                            perms([1 0 0]);
                            perms([1 1 0])]*obj.LatticeInfo.reciprLatticeVector;
                    elseif contains(cellType,'conventional')
                        verts = [0 0 0;
                            1 1 1;
                            perms([1 0 0]);
                            perms([1 1 0])]./obj.lengths;
                        % move cell to center it in 0, as Brillouin cell
                        verts = verts - .5./obj.lengths;
                    elseif contains(cellType,'brillouin')
                        verts = obj.LatticeInfo.reciprCellVertices;
                    end
            end
            verts = uniquetol(verts,1e-7,'ByRows',true);
            faces = facesPatch3D(verts);
            % Save in a struct. Easier to write patch function
            cellPoly = struct('vertices',verts,'faces',faces, ...
                'FaceColor','none');
        end
    end

end