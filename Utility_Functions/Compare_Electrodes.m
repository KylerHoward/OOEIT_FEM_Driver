function Compare_Electrodes(L, nodes, E_connect, flags)
    %{
    Compare the number of faces and the area of each electrode
    1/30/26 Kyler Howard

    param: L         - How many electrodes are present
    param: nodes     - All nodes in 3D space
    param: E_connect - Cell array containing connectivity matrix of each electrode
    param: flags     - Settings used to create the simulation
    %}

    num_E_faces = zeros(1,L);
    E_areas     = zeros(1,L);
    for i = 1:length(E_connect)
        num_E_faces(i) = size(E_connect{i},1);
        
        for j = 1:num_E_faces(i)
            face_nodes = E_connect{i}(j,:);
            point1 = nodes(face_nodes(1),:);
            point2 = nodes(face_nodes(2),:);
            point3 = nodes(face_nodes(3),:);
    
            E_areas(i) = E_areas(i) + norm(cross(point3-point1, point3-point2)) / 2;
        end
    end
    clear i j
    
    if flags.E_choice < 3 || (flags.E_choice == 5 && flags.E_type == "patch")
        if flags.CP_choice == 1
            front = 1:16;
            back  = 17:32;
        else
            front = [3:6, 11:14, 19:22, 27:30];
            back  = [1:2,  7:10, 15:18, 23:26, 31:32];
        end
    else
        front = [5:12, 21:28];
        back  = [1:4,13:20,29:32];
    end
    
    if flags.plot_electrodes == 1
        figure()
            subplot(1,2,1)
                hold on
                bar(front,num_E_faces(front),'EdgeColor',"#0072BD",'FaceColor',"#0072BD")
                bar(back,num_E_faces(back),  'EdgeColor',"#4DBEEE",'FaceColor',"#4DBEEE")
                plot(1:L, mean(num_E_faces)*ones(1,L),'r')
                xlabel("Electrode")
                title("Number of Faces")
                title(sprintf("Number of Faces\nMean: %.0f, STD: %.0f", mean(num_E_faces), std(num_E_faces)))
                legend("Front", "Back")
            subplot(1,2,2)
                hold on
                bar(front,E_areas(front),'EdgeColor',"#0072BD",'FaceColor',"#0072BD")
                bar(back,E_areas(back),  'EdgeColor',"#4DBEEE",'FaceColor',"#4DBEEE")
                plot(1:L, mean(E_areas)*ones(1,L),'r')
                xlabel("Electrode")
                title(sprintf("Electrode Area (mm²)\nMean: %.2f, STD: %.2f", mean(E_areas), std(E_areas)))
                legend("Front", "Back")
    end
end