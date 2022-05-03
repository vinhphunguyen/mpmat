function neighbors = getNeighbors1D(elemId, numx)
% find neighbors of a given element "elemId" in a
% structured one dimensional mesh.
% Used in one dimensional GIMP analyses.
% VP Nguyen

col = rem (elemId,numx);

if col == 0
    col = numx;
end

if     col == 1
    cols = [1 2];
elseif col == numx
    cols = [numx-1 numx];
else
    cols = [col-1 col col + 1];
end

neighbors=cols;