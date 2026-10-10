function Plot_Lung_Heart_Stats(nframes, conds, lung_cond, blood_cond, heart_cond, flags, frames)
    figure('WindowState', 'maximized')
        subplot(3,4,1:2)
            plot(1:nframes, flags.breath_curve, 'b')
            xlim([0,nframes])
            if ~isempty(frames)
                xline(frames,"k","LineWidth",2)
            end
            xlabel("Frame")
            ylabel("% of Cycle")
            title("Breathing Cycle")
            del_y = 0.1*range(flags.breath_curve);
            ylim([min(flags.breath_curve) - del_y, max(flags.breath_curve) + del_y])
        subplot(3,4,3:4)
            plot(1:nframes, flags.heart_curve, 'r')
            xlim([0,nframes])
            if ~isempty(frames)
                xline(frames,"k","LineWidth",2)
            end
            xlabel("Frame")
            ylabel("% of Cycle")
            title("Cardiac Cycle")
            del_y = 0.1*range(flags.heart_curve);
            ylim([min(flags.heart_curve) - del_y, max(flags.heart_curve) + del_y])
        sp5 = subplot(3,4,5);
            plot(1:nframes, conds.lung(1,:) , 'b')
            xlim([0,nframes])
            xlabel("Frame")
            ylabel("Conductivity (S/m)")
            title("Lung Conductivity - Ventilation")
            del_y = 0.1*range(conds.lung(1,:));
            ylim([min(conds.lung(1,:)) - del_y, max(conds.lung(1,:)) + del_y])
        subplot(3,4,6)
            plot(1:nframes, blood_cond(1,:) , 'b')
            xlim([0,nframes])
            xlabel("Frame")
            ylabel("Conductivity (S/m)")
            title("Lung Conductivity - Perfusion")
            del_y = 0.1*range(blood_cond(1,:));
            ylim([min(blood_cond(1,:)) - del_y, max(blood_cond(1,:)) + del_y])
        subplot(3,4,7)
            plot(1:nframes, lung_cond(1,:) , 'b')
            xlim([0,nframes])
            xlabel("Frame")
            ylabel("Conductivity (S/m)")
            title("Mixed Lung Conductivity")
            del_y = 0.1*range(lung_cond(1,:));
            ylim([min(lung_cond(1,:)) - del_y, max(lung_cond(1,:)) + del_y])
        subplot(3,4,8)
            plot(1:nframes, heart_cond(1,:) , 'r')
            xlim([0,nframes])
            xlabel("Frame")
            ylabel("Conductivity (S/m)")
            title("Heart Conductivity")
            del_y = 0.1*range(heart_cond(1,:));
            ylim([min(heart_cond(1,:)) - del_y, max(heart_cond(1,:)) + del_y])
        subplot(3,4,9)
            axis off
            stats = sprintf(['Mean: %.3f S/m\n\n' ...
                             'Max:  %.3f S/m\n\n' ...
                             'Min:  %.3f S/m'], ...
                             mean(lung_cond), max(lung_cond), min(lung_cond));
            text(0.3,1, 'Lung Statistics', ...
                'Units','normalized', 'VerticalAlignment','top', ...
                'Interpreter','tex', 'FontWeight','bold');
            text(0.3,0.55, stats, 'Units','normalized');
        subplot(3,4,10:11)
            hold on
            plot(1:nframes, lung_cond(1,:) , 'b')
            plot(1:nframes, heart_cond(1,:) , 'r')
            xlim([0,nframes])
            if ~isempty(frames)
                xline(frames,"k","LineWidth",2)
            end
            xlabel("Frame")
            ylabel("Conductivity (S/m)")
            title("Lung and Heart Conductivities")
            legend("Lung", "Heart")
            del_y = 0.1*(min(lung_cond(1,:)) + max(heart_cond(1,:))/2);
            ylim([min(lung_cond(1,:)) - del_y, max(heart_cond(1,:)) + del_y])
        subplot(3,4,12)
            axis off
            stats = sprintf(['Mean: %.3f S/m\n\n' ...
                             'Max:  %.3f S/m\n\n' ...
                             'Min:  %.3f S/m'], ...
                             mean(heart_cond), max(heart_cond), min(heart_cond));
            text(0.3,1, 'Heart Statistics', ...
                'Units','normalized', 'VerticalAlignment','top', ...
                'Interpreter','tex', 'FontWeight','bold');
            text(0.3,0.55, stats, 'Units','normalized');
    sgtitle(sprintf('\\bf{Lung and Heart Complex Conductvities}'));
end