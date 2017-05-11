%--------------------------------------------------------------------------
function gridactual = firpmgrid_cg(nfilt,lgrid,ff,neg,nodd)
% firpmgrid

grid = zeros(1,1e5);% Make large initial memory

%    Generate frequency grid
nfcns = fix(nfilt/2);
if nodd == 1 && neg == 0
    nfcns = nfcns + 1;
end
grid(1) = ff(1);
delf = 1/(lgrid*nfcns);
% If value at frequency 0 is constrained, make sure first grid point
% is not too small:
if neg ~= 0 && grid(1) < delf
    % Handle narrow bands
    if ff(1) > sqrt(eps)
       grid(1) = ff(1);
    elseif delf < ff(2)
        grid(1) = delf;
    else
        grid(1) = 0.5*(ff(2) - ff(1));
    end
end
j = 1;
l = 1; gridSize = 1;
while (l+1) <= length(ff)
    fup = ff(l+1);
    newgrid = (grid(j)+delf):delf:(fup+delf);
    if length(newgrid) < 11
        delf1 = ((fup+delf) - (grid(j)+delf)) / 10;
        newgrid = (grid(j)+delf1):delf1:(fup+delf1);
    end
    %grid = [grid newgrid];
    grid(gridSize+1:gridSize+length(newgrid)) = newgrid;
    gridSize = gridSize + length(newgrid);
    
    %jend = length(grid);
    jend = gridSize;%length(grid);
    if jend > 1,
        grid(jend-1) = fup;
        j = jend;
    else
        j = jend + 1;
    end

    l = l + 2;
    if (l+1) <= length(ff)
        grid(j) = ff(l);
    end
end
ngrid = j - 1;
% If value at frequency 1 is constrained, remove that grid point:
if neg == nodd && (grid(ngrid) > 1-delf)
    if ff(end-1) < 1-delf
        ngrid = ngrid - 1;
    else
        grid(ngrid) = ff(end-1);
    end
end
gridactual = grid(1:ngrid);
