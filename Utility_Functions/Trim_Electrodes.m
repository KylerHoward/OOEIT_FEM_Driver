function E_nodes = Trim_Electrodes(E_nodes)    
    % Step 1: Concatenate all coordinates with cell IDs
    allPts = [];
    cellID = [];
    for k = 1:numel(E_nodes)
        Pts    = E_nodes{k};               % n×3 matrix of coordinates
        allPts = [allPts; Pts];
        cellID = [cellID; k*ones(size(Pts,1),1)];
    end
    
    % Step 2: Find duplicates across cells
    [~, ~, ic] = unique(allPts, 'rows');
    counts = accumarray(ic, 1);
    
    % Step 3: Mark coordinates that appear in more than one cell
    sharedMask = counts(ic) > 1;
    
    % Step 4: Remove shared coordinates from each cell
    offset = 0;
    for k = 1:numel(E_nodes)
        n = size(E_nodes{k},1);
        bad = sharedMask(offset + (1:n));   % rows to delete in this cell
        E_nodes{k}(bad,:) = [];
        offset = offset + n;
    end
end