function sigma_GT = Create_GT_Images_Filled(nodes, organ_faces, organ_nodes, sigma, E_connect, flags)

    % T5 = 214;
    % T8 = 149;
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
        sigma_GT{row} = PrepOrgans(organ_faces, organ_nodes, nodes, sigma, E_heights(row));
    end

    % % Find electrode heights, and put each of these in Sigma_GT
    % sigma_GT{1} = PrepOrgans(organ_faces, organ_nodes, nodes, sigma, T5);
    % sigma_GT{2} = PrepOrgans(organ_faces, organ_nodes, nodes, sigma, T8);

    % Create the image
    if flags.plot_GTs == 1
        % Shift each row down. 0 shift for top row, and then increase shift
        shift    = 0;
        for row = 1:size(sigma_GT, 1)
            % Apply shift to second column of each boundary
            for i = 1:size(sigma_GT{row},1)
                P = sigma_GT{row}{i,1};
                if ~isempty(P)
                    P.Vertices(:,2) = P.Vertices(:,2) - shift;
                    sigma_GT{row}{i,1}  = P;
                end

                if i == size(sigma_GT{row}, 1) && row == 1
                    bot_max_y = max(sigma_GT{row}{i,1}.Vertices(:,2));
                end
            end
    
            % Update shift: add the max of y values plus 33%
            shift = shift + 4/3*max(bot_max_y);
        end

        if flags.fixed_range
            sig_min = 0 + 0*1i;
            sig_max = 0.8 + 0.8*1i;
        else
            sig_min = min(real(sigma),[],'all') + min(imag(sigma),[],'all')*1i;
            sig_max = max(real(sigma),[],'all') + max(imag(sigma),[],'all')*1i;
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
        plts1 = cell([size(sigma_GT,1), size(sigma_GT{1},1)]);
        for row = 1:size(sigma_GT, 1)
            for o = size(sigma_GT{row},1):-1:1
                Poly         = sigma_GT{row}{o,1};
                if ~isempty(Poly) && ~isempty(Poly.Vertices)
                    if o == size(sigma_GT{row}, 1)
                        Poly = RemoveInsides(sigma_GT{row}, Poly, o);
                    end

                    color         = FindColor(cmap, sigma_GT, 1, o, frame, sig_min, sig_max);
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
        plts2 = cell([size(sigma_GT,1), size(sigma_GT{1},1)]);
        if flags.set_complex
            for row = 1:size(sigma_GT, 1)
                for o = size(sigma_GT{row},1):-1:1
                    Poly         = sigma_GT{row}{o,1};
                    if ~isempty(Poly) && ~isempty(Poly.Vertices)
                        if o == size(sigma_GT{row}, 1)
                            Poly = RemoveInsides(sigma_GT{row}, Poly, o);
                        end

                        color         = FindColor(cmap, sigma_GT, 1, o, frame, sig_min, sig_max);
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
                 "ValueChangedFcn",@(src,evt) updateFrame(ax1, ax2, plts1, plts2, src.Value, sigma_GT, sig_min, sig_max, flags));
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
    c = max(0,min(1+1i,c));
    
    % Find the row in the colormap that corresponds to this
    c_i   = round(1 + c*(size(cmap,1)-1));
    try
        color = cmap(real(c_i),:) + cmap(imag(c_i),:)*1i;
    catch
        color = cmap(real(c_i),:);
    end
end

function P = RemoveInsides(organ_polys, P, o)
    for i = o-1:-1:1
        P0 = organ_polys{i, 1};
        if ~isempty(P0)
            P = subtract(P, P0);
        end
    end
end

function [organ_polys] = PrepOrgans(organ_faces, organ_nodes, nodes, sigma, z0)
    warning off
    organ_list = fieldnames(organ_faces);

    organ_polys = cell([length(organ_list),  size(sigma,2)+1]);

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

        [organ_polys{o, 1}, bdry_nodes] = SliceOneOrgan(organ, nodes, F, z0);
        P = organ_polys{o, 1};

        if isempty(P.Vertices)
            continue
        end

        if strcmp(organ,"surface")
            inner_nodes = organ_nodes.boundary;
        elseif strcmp(organ, "esophagus")
            inner_nodes = organ_nodes.esoph_bdry;
        elseif strcmp(organ, "bone")
            inner_nodes = organ_nodes.bone_bdry;
        else
            inner_nodes = organ_nodes.(organ);
        end
        [~, idxInNodes] = ismember(inner_nodes, nodes, 'rows');

        organ_sigma  = sigma(idxInNodes,:);
        if strcmp(organ,"surface")
            body_sigma = mean(organ_sigma,1);
            region_sigma = mean(organ_sigma,1);
        else
            for i = 1:size(sigma,2)
                region_sigma(i) = mean(organ_sigma(abs(organ_sigma(:,i) - body_sigma(:,i)) > 1e-10,i),1);
            end
        end
        for i = 1:size(sigma,2)
            organ_polys{o, i+1} = region_sigma(i);
        end
    end
    warning off
end

function [P, bdry_nodes] = SliceOneOrgan(organ, nodes, F, z0)

    bdry_nodes = unique(F(:));

    zvals = nodes(:,3);
    triZ  = zvals(F);

    crossing = (min(triZ,[],2) <= z0) & (max(triZ,[],2) >= z0);
    Fcross   = F(crossing,:);

    segments = [];

    for i = 1:size(Fcross,1)
        verts = Fcross(i,:);
        p = nodes(verts,:);

        edges = [1 2; 2 3; 3 1];
        pts = [];

        for e = 1:3
            a = p(edges(e,1),:);
            b = p(edges(e,2),:);

            za = a(3);
            zb = b(3);

            % Check if one z is above z0 and the other is below
            % Both above or both postive makes first term positive
            if (za - z0)*(zb - z0) <= 0 && za ~= zb
                t  = (z0 - za) / (zb - za);
                pt = a + t*(b - a);
                pts = [pts; pt];
            end
        end

        if size(pts,1) == 2
            segments = [segments; pts(1,:) pts(2,:)];
        end
    end

    if isempty(segments)
        P = polyshape();
        return
    end

    N  = size(segments,1);
    P1 = segments(:,1:2);
    P2 = segments(:,4:5);

    tol  = 1e-6;
    used = false(N,1);
    loops = {};

    for i = 1:N
        if used(i)
            continue
        end

        loop  = P1(i,:);
        endpt = P2(i,:);
        used(i) = true;

        changed = true;
        while changed
            changed = false;

            for j = 1:N
                if used(j), continue; end

                p1 = P1(j,:);
                p2 = P2(j,:);

                if norm(endpt - p1) < tol
                    loop  = [loop; p1];
                    endpt = p2;
                    used(j) = true;
                    changed = true;

                elseif norm(endpt - p2) < tol
                    loop  = [loop; p2];
                    endpt = p1;
                    used(j) = true;
                    changed = true;
                end
            end
        end

        loops{end+1} = loop;
    end

    numLoops = numel(loops);
    isHole   = false(numLoops,1);

    if strcmp(organ, 'body')
        loops = loops(1);
        isHole = false;
        numLoops = 1;
    else
        for i = 1:numLoops
            Ci = mean(loops{i},1);
    
            for j = 1:numLoops
                if i == j, continue; end
    
                Lj = loops{j};
                if inpolygon(Ci(1), Ci(2), Lj(:,1), Lj(:,2))
                    isHole(i) = true;
                end
            end
        end
    end

    P = polyshape();

    for i = 1:numLoops
        if ~isHole(i)
            outer = loops{i};
            P = addboundary(P, outer(:,1), outer(:,2));

            for j = 1:numLoops
                if isHole(j)
                    Cj = mean(loops{j},1);
                    if inpolygon(Cj(1), Cj(2), outer(:,1), outer(:,2))
                        hole = loops{j};
                        P = addboundary(P, hole(:,1), hole(:,2));
                    end
                end
            end
        end
    end
end
