function CutPoly(crystal, Coefplane, axPoly)

%Esta funcion deberia coger los volumenes y los polihedros de la otra
%funcion y calcular el corte y tal. En la original solo quiero que te pinte
%los polihedros.

n = normal(crystal,Coefplane(1:3)); % Convert to real vector
hkl = Coefplane(1:3); % Miller indices
d = -Coefplane(4);
mvId = [0 0 0;Coefplane(1:3);1 0 0;0 1 0];

pPlane = constantplane(axPoly,hkl,-d+0.01);
pPlane.FaceColor = [.2 .4 .7];
pPlane.FaceAlpha = 1;
hold(axPoly,'off');
axPoly.Color = 'none';
axPoly.FontSize = 24;
axis tight;



C = [0.5 0.5 0.5];
numC = numel(crystal.Data.atom_site_fract_x);
if numC > size(C,1)
    C = tab10(numC+1);
    voidColor = C(end,:);
elseif numC < size(C,1)
    voidColor = C(end,:);
else
    voidColor = [.4 .7 .8];
end
C = C(1:numC,:);
clearvars numC




% Create lattice vectors
LatR = latticeVectors(crystal,"real","conventional");
assignin('base','vectors',LatR)
% create convenional cell patch
[cellP] = cellPatch(crystal);

% Load 1st cell
[~,R0] = breakSymm(crystal,1,'IncludeBoundary',true,'NormalizedPosition',false)
assignin('base','R0',R0)
% Decide atoms I want to paint
Rpaint = R0;
Rpaint = cellfun(@(x) x.*[-1 -1 -1]+sum(LatR),R0,'UniformOutput',false);
R0=Rpaint;
% ALL atoms, for Voronoi. Use 3 cells to each side
Rall = cellfun(@(x) replicateCell(x,1,LatR),R0, ...
    'UniformOutput',false);


% Create polyhedra of R0 atoms
Fpoly = cellfun(@(R0) wignerSeitz3D(vertcat(Rall{:}),R0, ...
     "Method","voronoi","Output","struct"),R0,'UniformOutput',false);
for i = 1:numel(R0) % assign colors to the atoms
    tempColor = repmat({C(i,:)},size(Fpoly{i}));
    [Fpoly{i}.faceColor] = tempColor{:};
end

pp = planeModel([n d]);
Fpoly = [Fpoly{:}];

% Replicate polyhedra around using mvId as indices ----------------------
mvts = mvId*LatR;
mvts = vertcat(mvts,-mvts);
% Repeat calculated polyhedra, it moves only one time in each direction
newFPoly = replicatePolyhedron(Fpoly,mvts);
assignin("base",'newpoly',newFPoly)

% 
% save volumes of full polyhedra.....
[~,vols] = cellfun(@(x) facesPatch3D(x),{Fpoly.vertices},'UniformOutput',false);
[fullVols,idFull] = uniquetol([vols{:}],1e-3);
% reshape idFull
idFull = cellfun(@(x)all(x==C,2)', ...
    num2cell(vertcat(Fpoly(idFull).faceColor),2),'UniformOutput',false);
idFull = vertcat(idFull{:});
% order so that first element corresponds to first atom
[~,tempOrder]=sortrows(idFull,1:size(idFull,2),'descend');
idFull = idFull(tempOrder,:);
fullVols = fullVols(tempOrder);
clearvars tempOrder idFull

% calculate cut polyhedra-----------------------------------------------
% Obtain planes for every polyhedron face and add the defined cutPlane
planes = arrayfun(@(x) [planesInPatch(x.vertices); ...
    pp.Parameters],newFPoly,'UniformOutput',false);
% Resolve intersection of planes for new vertices
points = cellfun(@(x) solvePlaneIntersection(x), planes,'UniformOutput',false);

% Create polyhedron from its vertices. Save faces and volumes
[ff,vols] = cellfun(@(x) facesPatch3D(x),points,'UniformOutput',false);

volsUniq = uniquetol([vols{:}],1e-10); % save unique values, for testing
vols = [vols{:}]';
% construct a struct with cut polyhedra
Spoly = newFPoly; % initialize with same fields


[Spoly.vertices] = points{:};
[Spoly.faces] = ff{:};
assignin("base",'spoly',Spoly)

vector = [];
l = 1;
for i = 1:length(Spoly)

    if abs(Spoly(i).UserData.rAtom(1,1)) <= LatR(1,1) && abs(Spoly(i).UserData.rAtom(1,2)) <= LatR(2,2) && Spoly(i).UserData.rAtom(1,1) >= 0 && Spoly(i).UserData.rAtom(1,2) >= 0 && Spoly(i).UserData.rAtom(1,3) >= 0
        %find(isnan(Spoly.faces))

        vector(l,1) = Spoly(i).UserData.rAtom(1,1);       
        vector(l,2) = Spoly(i).UserData.rAtom(1,2); 
        vector(l,3) = Spoly(i).UserData.rAtom(1,3); 
        Spoly_aux(l) = Spoly(i);
        vols_aux(l) = vols(i);
        newFPoly_aux(l) = newFPoly(i);
        %full_vols_aux(l) = fullVols(i)
        l = l+1;

    end


end
vector_u = unique(vector,'rows');
%Spoly_aux = unique( Spoly.UserData.rAtom,'rows');
Spoly = Spoly_aux;
assignin("base",'spoly_lessatoms',Spoly)
assignin("base",'aver',vector)
assignin("base",'aver2',vector_u)
vols = vols_aux.';
newFPoly = newFPoly_aux;
assignin("base",'aver3',newFPoly)
%fullVols = full_vols_aux;

clearvars ff faces planes points Spoly_aux l

idVoid = []; % ids for void polyhedra
idPolyAtom = false(numel(Spoly),size(C,1));
% Paint only below cut Plane

for j = 1:numel(Spoly)
    % make index to detect which atom it is for every polyhedron
    idPolyAtom(j,:) = all(Spoly(j).faceColor == C, 2)';
    % detect those above plane and change their color
    if Spoly(j).UserData.rAtom * n' + d > 0
        idVoid(end+1) = j;
        Spoly(j).UserData.ogColor = Spoly(j).faceColor;
        Spoly(j).faceColor = voidColor;
        Spoly(j).faceAlpha = .7; % void transparency
    else
        Spoly(j).faceAlpha = .4;
    end
end

% Volume calculation------------------------------------
% Total sum of volumes for void polyhedra below plane
voidVols = vols(idVoid); % select volumes of void polyhedra
idZeros = voidVols==0; % filter out the zero volumes (above plane)
voidVols(idZeros)=[];
idPolyAtomVoid = idPolyAtom(idVoid,:); % repeat for atom index
idPolyAtomVoid(idZeros,:) = []; % This identifies which atom cause void

clearvars idZeros
% Work with the rest to obtain partially cut volumes
cutVols = vols;
assignin("base",'cutvols',cutVols)
% remove void volumes
cutVols(idVoid) = [];
% repeat for index carrying the atom type
idPolyAtomCut = idPolyAtom;
idPolyAtomCut(idVoid,:) = [];
% remove full volumes (10 decimal precision)
% check whether is a complete polyhedra (and from what atom)
isFull = round(cutVols-fullVols,10) == 0;
cutVols(any(isFull,2),:) = [];
idPolyAtomCut(any(isFull,2),:) = [];

% Substract to obtain the upper volume
try
    upCutVols = abs(cutVols - idPolyAtomCut*fullVols');
catch
    upCutVols = cutVols*0;
end

% show the different volumes to choose how to add thempla n
disp(["Void Vols: ",round(voidVols',4)])
disp(["Void Vol Total:", sum(voidVols)])
disp(["CutVols: ", round(upCutVols',4)])
disp(["CutVols Total:", sum(upCutVols)])

% Select cut polyhedra and calculate displacement of center of mass-------
SpolyCut = Spoly;
SpolyCut(idVoid) = [];
SpolyCut(any(isFull,2)) = []; % I am sure this could be done in 1 line
% Repeat for full polyhedra for comparison
newFPolyComp = newFPoly;
newFPolyComp(idVoid) = [];
newFPolyComp(any(isFull,2)) = [];

% CM of cut Poly
mvdCM = cellfun(@(x)centroidPoly1(x), ...
    {SpolyCut.vertices},'UniformOutput',false);
mvdCM = vertcat(mvdCM{:});
ogCM = cellfun(@(x)centroidPoly1(x), ...
    {newFPolyComp.vertices},'UniformOutput',false);
ogCM = vertcat(ogCM{:});

if isempty(mvdCM) || isempty(ogCM)
    display("CM array is empty")
    displCM=nan(1,3);
else
    displCM = (mvdCM-ogCM);
end

% plot polyhedra
arrayfun(@(x) patch(axPoly,x,'facealpha',.3),Spoly);
hold(axPoly,"on")
%arrayfun(@(x) patch(axPoly,x,'facealpha',.3),newFPoly);

% Plot CMs
if ~isempty(ogCM)&& ~isempty(mvdCM)
    scatter3(axPoly,ogCM(:,1), ogCM(:,2), ogCM(:,3),180, ...
        [0 .8 0].*ones(size(ogCM,1),1),'filled','^', ...
        'MarkerEdgeColor','k');
    scatter3(axPoly,mvdCM(:,1), mvdCM(:,2), mvdCM(:,3),180, ...
        [0 .8 .8].*ones(size(ogCM,1),1),'filled','v','MarkerEdgeColor','k');
end

end