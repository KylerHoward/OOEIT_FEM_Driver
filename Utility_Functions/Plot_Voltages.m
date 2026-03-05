function Plot_Voltages(nodes, Uall, flags)
    %{
    Create a dynamic "video" plot of the global voltages
    2/17/26 Kyler Howard

    param: nodes - All global nodes in 3D space
    param: Uall  - Voltage value at each node
    param: flags - Settings used to run the simulation
    %}

    % --- Create UIFigure ---
    fig = uifigure('Position',[573 100 700 700]);
    
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
    
    % --- Initial Plotting ---
    % Frame index starts at 1
    frame = 1;
    
    % Plot 1
    s1 = scatter3(ax1, nodes(:,1), nodes(:,2), nodes(:,3), 10, real(Uall(:,flags.plot_volts,frame)), 'filled');
    title(ax1, sprintf("Frame %d\nReal Voltage: CP %d", frame, flags.plot_volts))
    xlabel(ax1, "x (mm)")
    ylabel(ax1, "y (mm)")
    zlabel(ax1, "z (mm)")
    axis(ax1, 'equal')
    c = colorbar(ax1);
    c.Label.String = "mV";
    if flags.fixed_range == 1
        clim([min(Uall(:,:,:),[],'all'), max(Uall(:,:,:),[],'all')])
    end
    % c.Location     = 'southoutside';
    
    % Plot 2 (if complex)
    if flags.set_complex
        s2 = scatter3(ax2, nodes(:,1), nodes(:,2), nodes(:,3), 10, imag(Uall(:,flags.plot_volts,frame)), 'filled');
        title(ax1, sprintf("Frame %d\nImaginary Voltage: CP %d", frame, flags.plot_volts))
        xlabel(ax1, "x (mm)")
        ylabel(ax1, "y (mm)")
        zlabel(ax1, "z (mm)")
        axis(ax2, 'equal')
        c = colorbar(ax2);
        c.Label.String = "mV";
        % c.Location     = 'southoutside';
    else
        s2 = [];
    end
    
    % --- Slider Panel ---
    sliderGrid = uigridlayout(mainGrid,[3 1]);
    sliderGrid.RowHeight = {75,70,70};
    
    % --- Frame Slider Block --- 
    frameBlock = uigridlayout(sliderGrid,[2 1]); 
    frameBlock.RowHeight = {15,30};
    
    uilabel(frameBlock, ... 
            'Text','Frame Index', ... 
            'HorizontalAlignment','center', ... 
            'FontWeight','bold');
    uislider(frameBlock, ...
             'Limits',[0.9 size(Uall,3)], ...
             'Value',1, ...
             'MajorTicks',1:size(Uall,3), ...
             'MinorTicks',[],...
             'ValueChangedFcn',@(src,evt) updateFrame(ax1, ax2, s1, s2, src.Value, Uall, flags));
    
    % --- X-axis Slider Block --- 
    xaxisBlock = uigridlayout(sliderGrid,[2 1]); 
    xaxisBlock.RowHeight = {15,30};
    
    uilabel(xaxisBlock, ... 
            'Text','X-Axis Zoom', ... 
            'HorizontalAlignment','center', ... 
            'FontWeight','bold');
    try
        uislider(xaxisBlock, 'range',...
                 'Limits',[min(nodes(:,1)), ceil(max(nodes(:,1)))],...
                 'Value', [min(nodes(:,1)), ceil(max(nodes(:,1)))],...
                 'Step', 5,...
                 'ValueChangedFcn',@(src,evt) updateX(ax1, ax2, src.Value, flags));
    catch
        uislider(xaxisBlock,...
                 'Limits',[min(nodes(:,1)), ceil(max(nodes(:,1)))],...
                 'Value', min(nodes(:,1)),...
                 'ValueChangedFcn',@(src,evt) updateX(ax1, ax2, src.Value, flags));
    end
    
    % --- Y-axis Slider Block --- 
    yaxisBlock = uigridlayout(sliderGrid,[2 1]); 
    yaxisBlock.RowHeight = {15,30};
    
    uilabel(yaxisBlock, ... 
            'Text','Y-Axis Zoom', ... 
            'HorizontalAlignment','center', ... 
            'FontWeight','bold');
    try
        uislider(yaxisBlock, 'range',...
                 'Limits',[min(nodes(:,2)), ceil(max(nodes(:,2)))],...
                 'Value', [min(nodes(:,2)), ceil(max(nodes(:,2)))],...
                 'Step', 5,...
                 'ValueChangedFcn',@(src,evt) updateY(ax1, ax2, src.Value, flags));
    catch
        uislider(yaxisBlock,...
                 'Limits',[min(nodes(:,2)), ceil(max(nodes(:,2)))],...
                 'Value', min(nodes(:,2)),...
                 'ValueChangedFcn',@(src,evt) updateY(ax1, ax2, src.Value, flags));
    end

    % --- Callback Functions ---
    function updateFrame(ax1, ax2, s1, s2, framenum, colordata, flags)
        plotframe = round(framenum);

        % Update plot 1
        s1.CData = real(colordata(:,flags.plot_volts,plotframe));
        title(ax1, sprintf("Frame %d\nReal Voltage: CP %d", plotframe, flags.plot_volts))

        % Update plot 2
        if flags.set_complex
            s2.CData = imag(colordata(:,flags.plot_volts,plotframe));
            title(ax2, sprintf("Frame %d\nImaginary Voltage: CP %d", plotframe, flags.plot_volts))
        end
    end

    function updateX(ax1, ax2, rangeVals, flags)
        if numel(rangeVals) == 2
            minVal = rangeVals(1);
            maxVal = rangeVals(2);
        else
            minVal = rangeVals;
            xl     = xlim(ax1);
            maxVal = xl(2);
        end
        ax1.XLim = [minVal, maxVal];
        if flags.set_complex
            ax2.XLim = [minVal, maxVal];
        end
    end

    function updateY(ax1, ax2, rangeVals, flags)
        if numel(rangeVals) == 2
            minVal = rangeVals(1);
            maxVal = rangeVals(2);
        else
            minVal = rangeVals;
            yl     = ylim(ax1);
            maxVal = yl(2);
        end
        ax1.YLim = [minVal, maxVal];
        if flags.set_complex
            ax2.YLim = [minVal, maxVal];
        end
    end
end