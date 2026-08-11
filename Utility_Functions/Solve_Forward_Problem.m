function [Umeas, Uall] = Solve_Forward_Problem(sigma, mesh_info, frame_info, heart_BC, solver, fmesh, noise, flags)
    %{
    In order to run either a parfor or a regular for loop, solving the 
    forward problem is sent to a separate function
    5/22/26 - Kyler Howard
    8/3/26  - Edited by Kyler Howard

    param: sigma      - Assigned conducitvity on every node [nN x 1]. 
    param: mesh_info  - Structure containing variables about the forwrd mesh. 
                        nodes, connectivity, labels, lung_nodes, L, K, and cur_pat
    param: frame_info - Structure containing variables about the frames to be solved on. 
                        n_hframes, n_bframes, bframe, i_permutation, and num_permutations
    param: heart_BC   - Structure containing the Dirichlet BC indicies and values. 
                        indices and vals
    param: solver     - Forward problem solver from OOEIT.
    param: fmesh      - Mesh-object from OOEIT.
    param: noise      - Noise and error parameters.
                        [relative noise, absolute noise, relative system error, absolute system error]
    param: flags      - Various flags controlling plotting and other parameters

    return: Umeas - Measured voltages on the electrdoes for each current pattern and frame [L x K x nF].
    return: Uall  - Global voltages on all nodes for each current pattern and frame [nN x K x nF].
    %}

    % Unpack structures back to local individual variables
    L       = mesh_info.L;
    K       = mesh_info.K;
    cur_pat = mesh_info.cur_pat;

    n_hframes        = frame_info.n_hframes;
    n_bframes        = frame_info.n_bframes ;
    bframe           = frame_info.bframe;
    i_permutation    = frame_info.i_permutation;
    num_permutations = frame_info.num_permutations;

    heart_BC_indices = heart_BC.indices;
    heart_BC_vals    = heart_BC.vals;

% ----------------------------------------------------------------------- %
%%                                 Solve                                  %
% ----------------------------------------------------------------------- %                
    % Make a local version for each core
    Umeas = zeros(L, K, n_hframes);
    Uall  = zeros(fmesh.ng, K, n_hframes);

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
                fprintf("      It took %.2f seconds to solve for voltages\n", solve_time)
            end
        
            Umeas_frame = Umeas_frame * 1e3; % Convert V to mV
            Uall_frame  = Uall_frame  * 1e3; % Convert V to mV

            % Save the local frame
            Umeas(:,:,hframe) = reshape(Umeas_frame, L, K); 
            Uall(:,:,hframe) = Uall_frame(1:fmesh.ng,:);
        end% end running a solve
    end % end looping over heart frames
end