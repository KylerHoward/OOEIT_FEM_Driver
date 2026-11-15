function sigma_GT = Create_GT_Images_Filled(nodes, organ_faces, organ_nodes, lung_split, sigma, E_connect, flags)
    % Ensure surface face is first! This is necessary for logic
    temp = organ_faces.surface;
    organ_faces = rmfield(organ_faces,'surface');
    organ_faces.surface = temp;


    frame = 1;

    % Find the z heights of each electrode and sort them by rows
    E_nodes   = cellfun(@(x) nodes(x,:), E_connect, 'UniformOutput',false);
    E_heights = cellfun(@(x) mean(x(:,3)), E_nodes);
    if flags.E_type == "patch"
        if flags.CP_choice == 1
            E_heights = sort([mean(E_heights([1:4,17:20])), mean(E_heights([5:8,21:24])), mean(E_heights([9:12,25:28])), mean(E_heights([13:16,29:32]))], "descend");
        else
            E_heights = sort([mean(E_heights(1:8)), mean(E_heights(9:16)), mean(E_heights(17:24)), mean(E_heights(25:32))], "descend");
        end
    else
        E_heights = sort([mean(E_heights(1:16)), mean(E_heights(17:32))], "descend");
    end

    % Prep the size of the image
    if flags.E_type == "patch"
        sigma_GT = cell(4,1);
    else
        sigma_GT = cell(2,1);
    end

    % Find the perimeter/polyshape of each organ and the conductivity value per frame
    for row = 1:length(E_heights)
        sigma_GT{row} = PrepOrgans(organ_faces, organ_nodes, lung_split, nodes, sigma, E_heights(row));
    end
    
    % Create the plotted sigma as either absolutes or difference images
    plot_sigma = sigma_GT;
    if flags.plot_diff == 1
        ref_img = cell(length(E_heights), 1);
        % Make the GT to difference images
        for row = 1:length(E_heights)
            ref_img{row} = sigma_GT{row}(:,2);
    
            for c = 2:size(sigma_GT{row}, 2)
                for r = 1:size(sigma_GT{row},1)
                    if ~isempty(sigma_GT{row}{r,c})
                        plot_sigma{row}{r,c} = sigma_GT{row}{r,c} - ref_img{row}{r,1};
                    end
                end
            end
        end
    end

    % % Find electrode heights, and put each of these in Sigma_GT
    % sigma_GT{1} = PrepOrgans(organ_faces, organ_nodes, nodes, sigma, T5);
    % sigma_GT{2} = PrepOrgans(organ_faces, organ_nodes, nodes, sigma, T8);

    %% Create the image
    if flags.plot_GTs == 1
        % Shift each row down. 0 shift for top row, and then increase shift
        shift   = 0;
        for row = 1:size(plot_sigma, 1)
            % Apply shift to second column of each boundary
            for i = 1:size(plot_sigma{row},1)
                P = plot_sigma{row}{i,1};

                if isa(P,'polyshape') && ~isempty(P.Vertices)
                
                    P.Vertices(:,2) = P.Vertices(:,2)-shift;
                    plot_sigma{row}{i,1} = P;
                
                    if i == size(plot_sigma{row},1) && row == 1
                        bot_max_y = max(P.Vertices(:,2));
                    end
                end
            end
    
            % Update shift: add the max of y values plus 33%
            shift = shift + 4/3*max(bot_max_y);
        end

        %%
        if flags.fixed_range
            if flags.plot_diff == 1
                sig_min = -0.1 - 0.1*1i;
                sig_max = 0.04 + 0.04*1i;
            else
                sig_min = 0 + 0*1i;
                sig_max = 0.8 + 0.8*1i;
            end
        else
            allVals = [];
            for i = 1:numel(plot_sigma)
                data = plot_sigma{i}(:,2:end);
                allVals = [allVals, data{~cellfun('isempty',data)}];
            end
            sig_min = min(real(allVals),[],'all') + min(imag(allVals),[],'all')*1i;
            sig_max = max(real(allVals),[],'all') + max(imag(allVals),[],'all')*1i;
        end

        %%
        % Set figure settings
        fig  = uifigure("color","w","Position",[573,337.67,700,420]);
        cmap = jet(256);
        % --- Grid Layout ---
        % Layout:
        %   Row 1: plots (1 or 2)
        %   Row 2: sliders (3 sliders stacked) autofitted
        mainGrid = uigridlayout(fig,[2 1]);
        mainGrid.RowHeight = {"1x", "fit"};
        mainGrid.ColumnWidth = {"1x"};

        % --- Plot Panel ---
        if flags.set_complex
            plotGrid = uigridlayout(mainGrid,[1 2]);
            plotGrid.ColumnWidth = {"1x","1x"};
        else
            plotGrid = uigridlayout(mainGrid,[1 1]);
        end

        % --- Axes Creation ---
        if flags.set_complex
            ax1 = uiaxes(plotGrid);
            ax2 = uiaxes(plotGrid);
            hold(ax1, "on")
            hold(ax2, "on")
        else
            ax1 = uiaxes(plotGrid);
            ax2 = [];
            hold(ax1, "on")
        end

        % Plot 1
        plts1 = cell([size(plot_sigma,1), size(plot_sigma{1},1)]);
        for row = 1:size(plot_sigma, 1)
            for o = size(plot_sigma{row},1):-1:1
                Poly         = plot_sigma{row}{o,1};
                if ~isempty(Poly) && ~isempty(Poly.Vertices)
                    if o == size(plot_sigma{row}, 1)
                        Poly = RemoveInsides(plot_sigma{row}, Poly, o);
                    end

                    color         = FindColor(cmap, plot_sigma, 1, o, frame, sig_min, sig_max);
                    plts1{row, o} = plot(ax1, Poly, 'FaceColor', real(color), 'EdgeColor', 'none','FaceAlpha',1);
                end
            end
        end
        title(ax1, sprintf("Frame %d Ground Truth\nConductivity", frame))
        axis(ax1, "equal", "off")
        set(ax1, "XDir", "reverse") % set in DICOM standard
        clim(ax1, real([sig_min, sig_max]))
        c = colorbar(ax1, "eastoutside");
        c.Label.String = "S/m";
        colormap(ax1, "jet")

        % Plot 2 (if complex)
        plts2 = cell([size(plot_sigma,1), size(plot_sigma{1},1)]);
        if flags.set_complex
            for row = 1:size(plot_sigma, 1)
                for o = size(plot_sigma{row},1):-1:1
                    Poly         = plot_sigma{row}{o,1};
                    if ~isempty(Poly) && ~isempty(Poly.Vertices)
                        if o == size(plot_sigma{row}, 1)
                            Poly = RemoveInsides(plot_sigma{row}, Poly, o);
                        end

                        color         = FindColor(cmap, plot_sigma, 1, o, frame, sig_min, sig_max);
                        plts2{row, o} = plot(ax2, Poly, 'FaceColor', imag(color), 'EdgeColor', 'none','FaceAlpha',1);
                    end
                end
            end
            title(ax2, sprintf("Frame %d Ground Truth\nSusceptivity", frame))
            axis(ax2, "equal", "off")
            set(ax2, "XDir", "reverse") % set in DICOM standard
            clim(ax2, imag([sig_min, sig_max]))
            c = colorbar(ax2, "eastoutside");
            c.Label.String = "S/m";
            colormap(ax2, "jet")
        else
            plts2 = [];
        end

        % --- Frame Slider Block --- 
        frameBlock = uigridlayout(mainGrid,[2 1]); 
        frameBlock.RowHeight = {15,30};
        
        uilabel(frameBlock, ... 
                "Text","Frame Index", ... 
                "HorizontalAlignment","center", ... 
                "FontWeight","bold");
        uislider(frameBlock, ...
                 "Limits",[0.9 size(sigma,2)], ...
                 "Value",1, ...
                 "MajorTicks",1:size(sigma,2), ...
                 "MinorTicks",[],...
                 "ValueChangedFcn",@(src,evt) updateFrame(ax1, ax2, plts1, plts2, src.Value, plot_sigma, sig_min, sig_max, flags));
    end

end

function updateFrame(ax1, ax2, plts1, plts2, framenum, sigma_GT, sig_min, sig_max, flags)
    plotframe = round(framenum);
    cmap = jet(256);

    % Update plot 1
    for row = 1:size(sigma_GT, 1)
        for o = 1:size(sigma_GT{1},1)
            if ~isempty(plts1{row,o})
                color = FindColor(cmap, sigma_GT, row, o, plotframe, sig_min, sig_max);
                plts1{row,o}.FaceColor = real(color);
            end
        end
    end
    title(ax1, sprintf("Frame %d Ground Truth\nConductivity", plotframe))
    clim(ax1, real([sig_min, sig_max]));

    % Update plot 2
    if flags.set_complex
        for row = 1:size(sigma_GT, 1)
            for o = 1:size(sigma_GT{1},1)
                if ~isempty(plts1{row,o})
                    color = FindColor(cmap, sigma_GT, row, o, plotframe, sig_min, sig_max);
                    plts2{row,o}.FaceColor = imag(color);
                end
            end
        end
        title(ax2, sprintf("Frame %d Ground Truth\nSusceptivity", plotframe))
        clim(ax1, imag([sig_min, sig_max]));
    end
end

function color = FindColor(cmap, sigma_GT, row, o, frame, sig_min, sig_max)
    % Normalize desired color value and range to 0-1
    c = real(sigma_GT{row}{o, frame+1} - sig_min)/real(sig_max - sig_min) + ...
        imag(sigma_GT{row}{o, frame+1} - sig_min)/imag(sig_max - sig_min)*1i;
    if isnan(imag(c))
        c = max(0,min(1,real(c)));
    else
        % c = max(0,min(1+1i,c));
        c = max(0,min(1,real(c))) + max(0,min(1,imag(c)))*1i;
    end
    
    % Find the row in the colormap that corresponds to this
    c_i   = round(1 + c*(size(cmap,1)-1));
    try
        color = cmap(real(c_i),:) + cmap(imag(c_i),:)*1i;
    catch
        color = cmap(real(c_i),:);
    end
end

function P = RemoveInsides(organ_polys, P, o)
    % Remove any inclusions that are in the polyshape
    for i = o-1:-1:1
        P0 = organ_polys{i, 1};
        if ~isempty(P0)
            P = subtract(P, P0);
        end
    end
end

function [organ_polys] = PrepOrgans(organ_faces, organ_nodes, lung_split, nodes, sigma, z0)
    warning off

    organ_list = fieldnames(organ_faces);
    organ_polys = cell([length(organ_list)+1, size(sigma,2)+1]);

    % Go backwards to hit surface first
    for o = numel(organ_list):-1:1

        organ = organ_list{o};

        if strcmp(organ, "body")
            continue
        end

        F = organ_faces.(organ);

        if isempty(F)
            continue
        end

        % [organ_polys{o,1}, ~] = SliceOneOrgan(organ, nodes, F, z0, flags, lung_split);
        % P = organ_polys{o,1};
        [P, ~] = SliceOneOrgan(organ, nodes, F, z0, lung_split);

        if strcmp(organ,'lung')
            left_sigma  = sigma(lung_split.left_nodes,:);
            right_sigma = sigma(lung_split.right_nodes,:);
        
            % Assign left to 2 and right to 3
            organ_polys{o,1}   = P{1};
            organ_polys{o+1,1} = P{2};
        
            for i = 1:size(sigma,2)
                % Figure out the conductivity of the lung
                not_body_left  = abs(left_sigma(:,i)  - body_sigma(:,i)) > 1e-10;
                not_body_right = abs(right_sigma(:,i) - body_sigma(:,i)) > 1e-10;
                region_sigma(1,i) = mode(left_sigma(not_body_left,  i),1);
                region_sigma(2,i) = mode(right_sigma(not_body_right,i),1);

                organ_polys{o,i+1}   = region_sigma(1,i);
                organ_polys{o+1,i+1} = region_sigma(2,i);
            end
        
            continue
        else
            if o > 2
                organ_polys{o+1, 1} = P;
            else
                organ_polys{o,1} = P;
            end
        end

        if isempty(P.Vertices)
            continue
        end

        if strcmp(organ,"surface")
            inner_nodes = organ_nodes.boundary;
        elseif strcmp(organ,"esophagus")
            inner_nodes = organ_nodes.esoph_bdry;
        elseif strcmp(organ,"bone")
            inner_nodes = organ_nodes.bone_bdry;
        else
            inner_nodes = organ_nodes.(organ);
        end

        [~,idxInNodes] = ismember(inner_nodes, nodes, 'rows');

        % Find the sigma value for each polyshape
        organ_sigma = sigma(idxInNodes,:);

        if strcmp(organ,"surface")
            body_sigma   = mean(organ_sigma,1);
            region_sigma = mean(organ_sigma,1);
        else
            for i = 1:size(sigma,2)
                not_body = abs(organ_sigma(:,i) - body_sigma(:,i)) > 1e-10;
                region_sigma(i) = mean(organ_sigma(not_body,i),1);
            end
        end

        for i = 1:size(sigma,2)
            if o > 2
                % o+1 to account for the lungs coming later
                organ_polys{o+1,i+1} = region_sigma(i);
            else
                % We're past the lungs, so now we acutally want o
                organ_polys{o,i+1} = region_sigma(i);
            end
        end
    end
    warning off

end

function [P, bdry_nodes] = SliceOneOrgan(organ, nodes, F, z0, lung_split)

    bdry_nodes = unique(F(:));

    %% ---------------- Slice mesh ---------------- %%
    zvals = nodes(:,3);

    triZ = zvals(F);

    crossing = (min(triZ,[],2) <= z0) & (max(triZ,[],2) >= z0);

    Fcross = F(crossing,:);
    
    segments      = [];
    
    leftSegments  = [];
    rightSegments = [];

    for i = 1:size(Fcross,1)
        verts = Fcross(i,:);
        p     = nodes(verts,:);
        edges = [1 2;
                 2 3;
                 3 1];

        pts = [];

        for e = 1:3
            a = p(edges(e,1),:);
            b = p(edges(e,2),:);

            za = a(3);
            zb = b(3);

            % Check if edge crosses slicing plane
            if (za-z0)*(zb-z0) <= 0 && za ~= zb
                t   = (z0-za)/(zb-za);
                pt  = a + t*(b-a);
                pts = [pts; pt];
            end
        end

        if size(pts,1) == 2
            seg = [pts(1,:), pts(2,:)];
        
            if strcmp(organ,'lung')
                sideVals = lung_split.side(verts);
                sideVals = sideVals(sideVals>0);
                triSide = mode(sideVals);
                if triSide == 1
                    leftSegments = [leftSegments; seg];
                elseif triSide == 2
                    rightSegments = [rightSegments; seg];
                end
            else
                segments = [segments; seg];
            end
        end
    end

    %% ---------------- Empty slice ---------------- %%
    if isempty(segments) && isempty(leftSegments) && isempty(rightSegments)
        P = polyshape();
        return
    end

    % ----------------- Reconstruction ----------------%
    if strcmp(organ,'lung')
        % c_left   = mean(leftSegments);
        % c_right  = mean(rightSegments);
        % a_left   = atan2(leftSegments(:,2)  - c_left(2),  leftSegments(:,1)  - c_left(1));
        % a__right = atan2(rightSegments(:,2) - c_right(2), rightSegments(:,1) - c_right(1));
        % [~, order_left]  = sort(a_left, 'descend');
        % [~, order_right] = sort(a__right, 'descend');
        % leftSegments  = leftSegments(order_left,:);
        % rightSegments = rightSegments(order_right,:);


        Pleft  = SegmentsToPolyshape(leftSegments, 'left_lung');
        Pright = SegmentsToPolyshape(rightSegments, 'right_lung');
        
        P = {Pleft,Pright};
    
        return
    else
        P = SegmentsToPolyshape(segments,organ);
    end

end

function P = SegmentsToPolyshape(segments,organ)
    if isempty(segments)
        P = polyshape();
        return
    end

    % Extract endpoints
    pts = [segments(:,1:2); segments(:,4:5)];

    % Merge coincident vertices
    % [verts,~,ic] = uniquetol(pts, 1e-6, 'ByRows',true);
    tol = 1e-8;
    ptsRound = round(pts/tol)*tol;
    [verts,~,ic] = unique(ptsRound,'rows');

    nSeg = size(segments,1);

    edgeStart = ic(1:nSeg);
    edgeEnd   = ic(nSeg+1:end);

    % Build graph
    G   = graph(edgeStart,edgeEnd);
    deg = degree(G);

    % Find path endpoints
    endNodes = find(deg==1);

    % Open path (lungs)
    if length(endNodes)==2
        pathNodes = shortestpath(G, endNodes(1), endNodes(2));

        boundary = verts(pathNodes,:);

        % close the contour across the cut
        boundary(end+1,:) = boundary(1,:);

        P = polyshape(boundary(:,1), boundary(:,2));
        return
    end

    % Closed contour (heart, body, etc.)
    if all(deg==2)
        adj = neighbors(G,1);
    
        prev = 1;
        curr = adj(1);
    
        pathNodes = 1;
    
        while curr ~= 1
            pathNodes(end+1) = curr;
    
            nbrs = neighbors(G,curr);
    
            next = nbrs(nbrs ~= prev);
    
            prev = curr;
            curr = next;
        end
    
        boundary          = verts(pathNodes,:);
        boundary(end+1,:) = boundary(1,:);
    
        P = polyshape(boundary(:,1), boundary(:,2));
        return
    end

    % Fallback
    fprintf(2,'%s: unexpected graph topology\n',organ)
    P = polyshape();
end
