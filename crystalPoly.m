function [testInfo] = crystalPoly(crystal,planeCoef,C,mvId)
arguments
    crystal
    planeCoef
    C
    mvId (:,3) = [0 0 0;planeCoef(1:3);1 0 0;0 1 0];
end
% crystalPoly(crystal,planeCoef,colors) calculates and plots polyhedra of
% the crystal, cleaved by the plane defined by planeCoef = [A B C D], where
% Ax+By+Cz+D=0 is the plane equation.

n = normal(crystal,planeCoef(1:3)); % Convert to real vector
hkl = planeCoef(1:3); % Miller indices
d = planeCoef(4);

% START ---------------------------------------------------------------
% Check if there is incorrect number of colors given
% If more colors are given, the next one is used for Void polyhedra
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
% create convenional cell patch
[cellP] = cellPatch(crystal);

% Load 1st cell
[~,R0] = breakSymm(crystal,1,'IncludeBoundary',true,'NormalizedPosition',false);
% Decide atoms I want to paint
Rpaint = R0;
Rpaint = cellfun(@(x) x.*[1 -1 -1]+sum(LatR),R0,'UniformOutput',false);
R0=Rpaint;
% ALL atoms, for Voronoi. Use 3 cells to each side
Rall = cellfun(@(x) replicateCell(x,1,LatR),R0, ...
    'UniformOutput',false);

% Define Cut plane
pp = planeModel([n d]);

% Create polyhedra of R0 atoms
Fpoly = cellfun(@(R0) wignerSeitz3D(vertcat(Rall{:}),R0, ...
     "Method","voronoi","Output","struct"),R0,'UniformOutput',false);
for i = 1:numel(R0) % assign colors to the atoms
    tempColor = repmat({C(i,:)},size(Fpoly{i}));
    [Fpoly{i}.faceColor] = tempColor{:};
end

Fpoly = [Fpoly{:}];
% Replicate polyhedra around using mvId as indices ----------------------
mvts = mvId*LatR;
% Repeat calculated polyhedra
newFPoly = replicatePolyhedron(Fpoly,mvts);

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
clearvars ff faces planes points

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
% remove void volumes
cutVols(idVoid) = [];
% repeat for index carrying the atom type
idPolyAtomCut = idPolyAtom;
idPolyAtomCut(idVoid,:) = [];
% remove full volumes (10 decimal precision)
% check whether is a complete polyhedra (and from what atom)
isFull = round(cutVols-fullVols,10) ==0;
cutVols(any(isFull,2),:) = [];
idPolyAtomCut(any(isFull,2),:) = [];

% Substract to obtain the upper volume
try
upCutVols = abs(cutVols - idPolyAtomCut*fullVols');
catch
    upCutVols = cutVols*0;
end

% show the different volumes to choose how to add them
disp(["Void Vols: ",round(voidVols',4)])
disp(["CutVols: ", round(upCutVols',4)])

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
    warning("CM array is empty")
    displCM=nan(1,3);
else
    displCM = (mvdCM-ogCM);
end
%-------------------------------------------------------------------------

% PLOT--------------------------------------------------------------------
figure();
clf;
axPoly = nexttile;
ax = findobj(gcf,'type','axes');
hold(ax,'on');
set(ax,'DataAspectRatio',[1 1 1],'view',[47 17])
xlabel('X (Å)'); ylabel('Y (Å)');zlabel('Z (Å)')
% zlim([0 crystal.lengths(3)*1.1])

% Plot cell countours
arrayfun(@(x) patch(axPoly,x,'facecolor','none','Linewidth',3),cellP);


% plot atoms
for i = 1:length(Rpaint)    
    scatter3(axPoly,Rpaint{i}(:,1), Rpaint{i}(:,2), Rpaint{i}(:,3),100, ...
        C(i,:).*ones(size(Rpaint{i},1),1),'filled','MarkerEdgeColor','k');
end
% plot polyhedra
arrayfun(@(x) patch(axPoly,x,'facealpha',.3),Spoly);
% arrayfun(@(x) patch(axPoly,x,'facealpha',.3),newFPoly);

% Plot CMs
if ~isempty(ogCM)&& ~isempty(mvdCM)
    scatter3(axPoly,ogCM(:,1), ogCM(:,2), ogCM(:,3),180, ...
        [0 .8 0].*ones(size(ogCM,1),1),'filled','^', ...
        'MarkerEdgeColor','k');
    scatter3(axPoly,mvdCM(:,1), mvdCM(:,2), mvdCM(:,3),180, ...
        [0 .8 .8].*ones(size(ogCM,1),1),'filled','v','MarkerEdgeColor','k');
end

pPlane = plot(planeModel([n d+.001]));
% pPlane.FaceColor = [.2 .4 .7];
pPlane.FaceAlpha = 1;
hold(ax,'off');
ax.Color = 'none';
ax.FontSize = 24;
axis tight;

%SAVE INFO AS OUTPUT
% whos;
testInfo.unitPolyhedra = Fpoly;
testInfo.cleavePolyhedra = removeNaNPolyhedra(Spoly);
testInfo.voidPolyhedraID = idVoid;
testInfo.fullVolume = fullVols;
testInfo.excessVolume = upCutVols;
testInfo.voidVolume = voidVols;
testInfo.ogCM = ogCM;
testInfo.movedCM = mvdCM;
%--------------------------------------------------------------------------
    function poly = removeNaNPolyhedra(poly)
        % removeNaNPolyhedra(poly) will remove those entries in the struct
        % that have NaN values in faces

        poly(cellfun(@(x) all(isnan(x),"all"), ...
            {poly.faces})) = [];

    end
end