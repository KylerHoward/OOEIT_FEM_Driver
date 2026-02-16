function sigma_GT = Create_GT_Images(thickness, nodes, E_nodes, sigma, flags)
    %{
    Create images of the ground truth conductivites at each row of electrodes
    1/30/26 Kyler Howard

    param: thickness - Amount of nodes to grab in each slice in mm
    param: nodes     - All global nodes in 3D space
    param: E_nodes   - Cell array containing the nodes in each electrode
    param: sigam     - Conductivity value at each node
    param: flags     - Settings used to run the simulation
    %}

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
        plane = plane & sigma~=0;
    
        sigma_GT{row} = plane;
    end
       
    % Create the image
    if flags.plot_GTs == 1
        % Determine how many columns to plot
        if flags.set_complex == 1
            num_col = 2;
        else
            num_col = 1;
        end
    
        % Initialize
        shift = 0;
        rowStart = 1;
        plane_nodes = zeros(sum(cellfun(@sum, sigma_GT)),2);
        plane_sigma = zeros(sum(cellfun(@sum, sigma_GT)),1);
        for k = 1:length(E_heights)
            % Current matrices
            nodemat = [nodes(sigma_GT{k},1), nodes(sigma_GT{k},2)];
            sigmat  = sigma(sigma_GT{k});
            nrows   = size(nodemat,1);
            PN = nodemat;
    
            % Apply shift to second column
            PN(:,2) = PN(:,2) - shift;
    
            % Append to full matrices
            plane_nodes(rowStart:rowStart+nrows-1, :) = PN;
            plane_sigma(rowStart:rowStart+nrows-1, :) = sigmat;
    
            % Update shift: add the max of y values plus 33%
            shift = shift + 4/3*max(nodemat(:,2));
            rowStart = rowStart + nrows;
        end
    
        figure('color','w','Position',[573,337.67,700,420]);
        for col = 1:num_col
            subplot(1,num_col,col)
            
            if col == 1
                scatter(plane_nodes(:,1), plane_nodes(:,2), 10, real(plane_sigma), 'filled')
            else
                scatter(plane_nodes(:,1), plane_nodes(:,2), 10, imag(plane_sigma), 'filled')
            end
    
            % Flip the x axis horizontally to be in DICOM standard
            % set(gca, 'XDir', 'reverse'); % KH: it already was I think
            axis equal off
    
            colormap("jet")
            if flags.fixed_range == 1
                clim([0, 0.8])
            else
                if col == 1
                    clim([min(real(plane_sigma)), max(real(plane_sigma))])
                else
                    clim([min(imag(plane_sigma)), max(imag(plane_sigma))])
                end
            end
        
            if col == 1
                title(sprintf("Ground Truth\nConductivity"))
            else
                title(sprintf("Ground Truth\nSusceptivity"))
            end
    
            c = colorbar("eastoutside");
            c.Label.String = "S/m";

            % Make sure it is in DICOM format
            set(gca, 'XDir','reverse')
        end
    end
end