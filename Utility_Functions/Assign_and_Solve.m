function [Umeas, Uall, sigma] = Assign_and_Solve(mesh_info, frame_info, heart_BC, solver, fmesh, noise, flags)
    %{
    In order to run either a parfor or a regular for loop, assigning
    conductivites and solving the problem is sent to a separate function
    5/22/26 - Kyler Howard

    

    param: nodes        - All nodes on the boundary of the surface
    param: connectivity - All nodes in the entire body
    param: labels       - Triangular faces that make up the surface of the body
    param: lung_nodes   - Heights of anatomical markers for electrode placement
    param: bframe       - Heights of anatomical markers for electrode placement
    param: bframe       - Heights of anatomical markers for electrode placement
    param: L       - Heights of anatomical markers for electrode placement
    param: K       - Heights of anatomical markers for electrode placement
    param: n_hframes       - Heights of anatomical markers for electrode placement
    param: flags        - Various flags controlling plotting and other parameters

    return: E_nodes       - 1x32 Cell array containing electrode nodes
    return: perim_mm_high - The perimeter around the body in mm
    %}

    % Unpack structures back to local individual variables
    nodes        = mesh_info.nodes; 
    connectivity = mesh_info.connect;
    labels       = mesh_info.labels;
    lung_nodes   = mesh_info.lung_nodes;
    L            = mesh_info.L;
    K            = mesh_info.K;
    cur_pat      = mesh_info.cur_pat;

    n_hframes        = frame_info.n_hframes;
    n_bframes        = frame_info.n_bframes ;
    bframe           = frame_info.bframe;
    i_permutation    = frame_info.i_permutation;
    num_permutations = frame_info.num_permutations;

    heart_BC_indices = heart_BC.indices;
    heart_BC_vals    = heart_BC.vals;

    % Extract which point on the breath curve we want to simulate
    flags.max_inspiration = flags.breath_curve(bframe);
    flags.cardiac_cycle   = flags.heart_curve(bframe);

    if flags.verbose == 1
        fprintf("   Assigning Conductivities\n")
    end

    % Creating the conductivity vector at the nodes (IN SIEMENS PER METER)
    sigma = Assign_Conductivities(nodes, connectivity, labels, lung_nodes, flags);

% ----------------------------------------------------------------------- %
%%                                 Solve                                  %
% ----------------------------------------------------------------------- %                
    % Make a local version for each core
    Umeas = zeros(L, K, n_hframes);
    Uall  = zeros(size(nodes,1), K, n_hframes);

    for hframe = 1:n_hframes
        if flags.solve_problem == 1
            fprintf("   Solving Forward Problem %d of %d\n", (i_permutation-1)*n_hframes*n_bframes + (bframe-1)*n_hframes + hframe, num_permutations*n_hframes*n_bframes)
        
            % Set the current pattern (IN AMPS)
            solver.Iel = cur_pat;
        
            % Solve the voltage
            solve_start          = tic;
            if flags.do_parfor == 1
                [Umeas_frame, Imeas_frame, Uall_frame] = MF_Simulation(fmesh,  [], sigma, solver.zeta, heart_BC_indices(hframe,:), heart_BC_vals, solver, "current", noise);
            else
                [Umeas_frame, Imeas_frame, Uall_frame] = MF_Simulation2(fmesh, [], sigma, solver.zeta, heart_BC_indices(hframe,:), heart_BC_vals, solver, "current", noise);
            end
            solve_time           = toc(solve_start);
            if flags.verbose == 1
                fprintf("      It took %.2f minutes to solve for voltages\n", solve_time/60)
            end
        
            Umeas_frame = Umeas_frame * 1e3; % Convert V to mV
            Uall_frame  = Uall_frame  * 1e3; % Convert V to mV

            % Save the local frame
            Umeas(:,:,hframe) = reshape(Umeas_frame, L, K); 
            Uall(:,:,hframe) = Uall_frame(1:size(nodes,1),:);
            % Umeas(:,:,(bframe-1)*n_hframes + hframe) = reshape(Umeas_frame, L, K);
            % Uall(:,:,(bframe-1)*n_hframes + hframe)  = Uall_frame(1:size(nodes,1),:);
        end% end running a solve
    end % end looping over heart frames
end