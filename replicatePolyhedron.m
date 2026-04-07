function newPoly = replicatePolyhedron(Poly,mvts)
% Creates an array from translating the vertices from patches along the
% direction in mvs.
% EXAMPLE: newPoly = replicatePolyhedron(Poly,mvts)
%
%
newPoly = [];
tempPoly = Poly; % to carry the fields and color values

% calculate each set of moved polyhedra, one movement at a time
for iMv = 1:size(mvts,1)
    % moved vertices for all polyhedra
    temp = cellfun(@(x) x+mvts(iMv,:),{Poly.vertices},'UniformOutput',false);
    [tempPoly.vertices] = temp{:}; % assign new vertices to struct
    
    %if there are atom coord saved, move them too
    if isfield(Poly(1),"UserData") && isfield(Poly(1).UserData,"rAtom")
        % Calculate at once the moved positions
        tempUser = arrayfun(@(x)x.UserData ,Poly , ...
            'UniformOutput' ,false);
        mvdAtoms = cellfun(@(x)x.rAtom+mvts(iMv,:),tempUser, ...
            'UniformOutput',false);
        % I use loop to assign them, cause I dont know how not to.
        for iPoly = 1:length(mvdAtoms)
            tempPoly(iPoly).UserData.rAtom = mvdAtoms{iPoly};
        end
    end
    


    newPoly = [newPoly, tempPoly];
end

%this here kills all the duplicate atoms that may appear as you replicate
%the polyhedra
if isfield(newPoly(1),"UserData") && isfield(newPoly(1).UserData,"rAtom")
    nPoly = length(newPoly);

    % Collect one representative atom per polyhedron
    atomPos = zeros(nPoly, size(newPoly(1).UserData.rAtom,2));
 
    for i = 1:nPoly
        atomPos(i,:) = newPoly(i).UserData.rAtom(1,:); 
    end
    % Find unique rows
    [~, idxUnique] = unique(atomPos, 'rows', 'stable');
   
    % Keep only unique polyhedra
    newPoly = newPoly(idxUnique);
% newPoly = deleteDuplicatePoly(newPoly);
end
end