function sigma_GT = Create_GT_Images(thickness, nodes, E_nodes, sigma, flags)
    %{
    Create images of the ground truth conductivites at each row of electrodes
    1/30/26 Kyler Howard

    param: thickness - Amount of nodes to grab in each slice in mm
    param: nodes     - All global nodes in 3D space
    param: E_nodes   - Cell array containing the nodes in each electrode
    param: sigam     - Conductivity value at each node
    param: flags     - Settings used to run the simulation

    return: sigma_GT - Cell array with the nodal indices for each ground truth row
    %}

    % set a starting frame
    frame = 1;

    % Find the z heights of each electrode and sort them by rows
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
    
    % Find what global nodes are in each image
    for row = 1:length(E_heights)
        plane = (nodes(:,3) > E_heights(row) - thickness/2) & (nodes(:,3) < E_heights(row) + thickness/2);
        plane = plane & sigma(:,frame)~=0;
    
        sigma_GT{row} = plane;
    end
       
    % Create the image
    if flags.plot_GTs == 1    
        % Initialize
        shift = 0;
        rowStart = 1;
        plane_nodes = zeros(sum(cellfun(@sum, sigma_GT)),2);
        plane_sigma = zeros(sum(cellfun(@sum, sigma_GT)),size(sigma,2));
        for k = 1:length(E_heights)
            % Current matrices
            nodemat = [nodes(sigma_GT{k},1), nodes(sigma_GT{k},2)];
            nrows   = size(nodemat,1);
            sigmat  = sigma(sigma_GT{k},:);
            PN      = nodemat;
    
            % Apply shift to second column
            PN(:,2) = PN(:,2) - shift;
    
            % Append to full matrices
            plane_nodes(rowStart:rowStart+nrows-1, :) = PN;
            plane_sigma(rowStart:rowStart+nrows-1, :) = sigmat;
    
            % Update shift: add the max of y values plus 33%
            shift = shift + 4/3*max(nodemat(:,2));
            rowStart = rowStart + nrows;
        end
    
        % Set figure settings
        fig = uifigure('color','w','Position',[573,337.67,700,420]);
        % --- Grid Layout ---
        % Layout:
        %   Row 1: plots (1 or 2)
        %   Row 2: sliders (3 sliders stacked) autofitted
        mainGrid = uigridlayout(fig,[2 1]);
        mainGrid.RowHeight = {'1x', 'fit'};
        mainGrid.ColumnWidth = {'1x'};

        % --- Plot Panel ---
        if flags.set_complex
            plotGrid = uigridlayout(mainGrid,[1 2]);
            plotGrid.ColumnWidth = {'1x','1x'};
        else
            plotGrid = uigridlayout(mainGrid,[1 1]);
        end

        % --- Axes Creation ---
        if flags.set_complex
            ax1 = uiaxes(plotGrid);
            ax2 = uiaxes(plotGrid);
        else
            ax1 = uiaxes(plotGrid);
            ax2 = [];
        end
        
        % Plot 1
        s1 = scatter(ax1, plane_nodes(:,1), plane_nodes(:,2), 10, real(plane_sigma(:,frame)), 'filled');
        title(ax1, sprintf("Frame %d Ground Truth\nConductivity", frame))
        axis(ax1, 'equal', 'off')
        set(ax1, 'XDir', 'reverse') % set in DICOM standard
        colormap("jet")
        if flags.fixed_range == 1
            try
                clim(ax1, [0, 0.8])
            catch
                caxis(ax1, [0, 0.8])
            end
        else
            try
                clim(ax1, [min(real(plane_sigma(:,frame))), max(real(plane_sigma(:,frame)))])
            catch
                caxis(ax1, [min(real(plane_sigma(:,frame))), max(real(plane_sigma(:,frame)))])
            end
        end
        c = colorbar(ax1, "eastoutside");
        c.Label.String = "S/m";
        
        % Plot 2 (if complex)
        if flags.set_complex
            s2 = scatter(ax2, plane_nodes(:,1), plane_nodes(:,2), 10, imag(plane_sigma(:,frame)), 'filled');
            title(ax2, sprintf("Frame %d Ground Truth\nSusceptivity", frame))
            axis(ax2, 'equal', 'off')
            set(ax2, 'XDir', 'reverse') % set in DICOM standard
            colormap("jet")
            if flags.fixed_range == 1
                try
                    clim(ax2, [0, 0.8])
                catch
                    caxis(ax2, [0, 0.8])
                end
            else
                try
                    clim(ax2, [min(imag(plane_sigma(:,frame))), max(imag(plane_sigma(:,frame)))])
                catch
                    caxis(ax2, [min(imag(plane_sigma(:,frame))), max(imag(plane_sigma(:,frame)))])
                end
            end
            c = colorbar(ax2, "eastoutside");
            c.Label.String = "S/m";
        else
            s2 = [];
        end
        
        % --- Frame Slider Block --- 
        frameBlock = uigridlayout(mainGrid,[2 1]); 
        frameBlock.RowHeight = {15,30};
        
        uilabel(frameBlock, ... 
                'Text','Frame Index', ... 
                'HorizontalAlignment','center', ... 
                'FontWeight','bold');
        uislider(frameBlock, ...
                 'Limits',[0.9 size(sigma,2)], ...
                 'Value',1, ...
                 'MajorTicks',1:size(sigma,2), ...
                 'MinorTicks',[],...
                 'ValueChangedFcn',@(src,evt) updateFrame(ax1, ax2, s1, s2, src.Value, plane_sigma, flags));
    end % end plotting

% --- Callback Functions ---
function updateFrame(ax1, ax2, s1, s2, framenum, colordata, flags)
        plotframe = round(framenum);

        % Update plot 1
        s1.CData = real(colordata(:,plotframe));
        title(ax1, sprintf("Frame %d Ground Truth\nConductivity", plotframe))
        if flags.fixed_range == 1
            try
                clim(ax1,[0, 0.8])
            catch
                caxis(ax1,[0, 0.8])
            end
        else
            try
                clim(ax1,[min(real(colordata(:,plotframe))), max(real(colordata(:,plotframe)))])
            catch
                caxis(ax1,[min(real(colordata(:,plotframe))), max(real(colordata(:,plotframe)))])
            end
        end

        % Update plot 2
        if flags.set_complex
            s2.CData = imag(colordata(:,plotframe));
            title(ax2, sprintf("Frame %d Ground Truth\nSusceptivity", plotframe))
            if flags.fixed_range == 1
                try
                    clim(ax2,[0, 0.8])
                catch
                    caxis(ax2,[0, 0.8])
                end
            else
                try
                    clim(ax2,[min(imag(colordata(:,plotframe))), max(imag(colordata(:,plotframe)))])
                catch
                    caxis(ax2,[min(imag(colordata(:,plotframe))), max(imag(colordata(:,plotframe)))])
                end
            end
        end

        
    end

end