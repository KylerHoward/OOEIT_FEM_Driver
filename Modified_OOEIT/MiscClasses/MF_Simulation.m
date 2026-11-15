function [Umeas, Imeas, Uall, Iall, Umeas_i, Imeas_i, e] = MF_Simulation2(fm, sigma_i, sigma, z, dirichlet_nodes, dirichlet_vals, fsolver, mode, err)
%{
A simple simulation that will return synthetic data. meas_i refer to
initial reference measurements, whereas Umeas and Imeas are the situation
we are actually interested in. if sigma_i = [], reference measurements are
not calculated.

Updated to use Matrix Free methods (gmres in MF_EITFEM.m)
Updated to return the nodal voltages throughout the mesh

Author: Petri Kuusela 5.4.2024
Edited: Kyler Howard  1.22.2025
%}

    %Use default values for arguments that are missing:
    if nargin < 3 || isempty(sigma)
        sigma = ones(fm.nginv,1);
    end
    if nargin < 4 || isempty(z)
        if ~isempty(fsolver)
            z = fsolver.zeta;
        else
            z = 1e-6*ones(fm.nEl,1);
        end
    end
    if nargin < 6 || isempty(mode)
        if nargin > 4 && ~isempty(fsolver)
            mode = fsolver.mode;
        else
            mode = 'current';
        end
    end
    if nargin < 7 || isempty(err)
        err = [1e-6 1e-3 1e-4 1e-2];
    end

    % These are the noise parameters
    noise_rel = err(1);
    noise_abs = err(2);
    e_systematic_rel = err(3);
    e_systematic_abs = err(4);

    if isempty(fsolver)
        fsolver = MF_EITFEM(fm);
        fsolver.zeta = z;
        fsolver.mode = mode;
        fsolver.sigmaMin = 1e-9;
        if strcmp(mode, 'potential')
            fsolver.Uel = eye(length(fm.E));%Injection pattern
            fsolver.Uel= fsolver.Uel(:);
        elseif strcmp(mode, 'current')
            Imeas = eye(length(fm.E));%Create the injection pattern
            Imeas(2:end,1:end-1) = Imeas(2:end,1:end-1) - eye(length(fm.E)-1);
            Imeas(1,end) = -1;
            fsolver.Iel = Imeas(:);
            fsolver.Iel = fsolver.Iel*1e-3;%Injected currents are usually order of mA
        end    
    end

    if ~isempty(sigma_i)
        %reference measurements:
        if length(sigma_i) < 2
            sigma_i = [sigma_i;sigma_i];
        end

        if strcmp(mode, 'potential')
            Imeas_i = fsolver.SolveForwardVec(sigma_i, dirichlet_nodes, dirichlet_vals);%these are the results
            Umeas_i = fsolver.Uel;
            %add noise and error:
            e_sys = randn(length(Imeas_i),1, 'like', Imeas_i)*e_systematic_abs*(max(Imeas_i)-min(Imeas_i));
            e_sysrel = randn(length(Imeas_i),1, 'like', Imeas_i)*e_systematic_rel;
            Imeas_i = Imeas_i.*(1+e_sysrel) + e_sys;
            Imeas_i = Imeas_i.*(1+randn(length(Imeas_i),1, 'like', Imeas_i)*noise_rel) + randn(length(Imeas_i),1, 'like', Imeas_i)*noise_abs*(max(Imeas_i)-min(Imeas_i));

        elseif strcmp(mode, 'current')
            [Umeas_i, ~] = fsolver.SolveForwardVec(sigma_i, dirichlet_nodes, dirichlet_vals);%these are the results
            Imeas_i = fsolver.Iel;
            %add noise and error:
            e_sys = randn(length(Umeas_i),1, 'like', Umeas_i)*e_systematic_abs*(max(Umeas_i)-min(Umeas_i));
            e_sysrel = randn(length(Umeas_i),1, 'like', Umeas_i)*e_systematic_rel;
            Umeas_i = Umeas_i.*(1+e_sysrel) + e_sys;
            Umeas_i = Umeas_i.*(1+randn(length(Umeas_i),1, 'like', Umeas_i)*noise_rel) + randn(length(Umeas_i),1, 'like', Umeas_i)*noise_abs*(max(Umeas_i)-min(Umeas_i));
        end

    end %end reference measurements
        
    %The actual measurements:
    if strcmp(mode, 'potential')
        [Imeas, Iall] = fsolver.SolveForwardVec(sigma, dirichlet_nodes, dirichlet_vals);%these are the results
        Umeas = fsolver.Uel;
        %add noise and error:
        if ~exist('esys', 'var')%there have been no homogeneous measurements where these have already been calculated
            e_sys = randn(length(Imeas),1, 'like', Imeas)*e_systematic_abs*(max(Imeas)-min(Imeas));
            e_sysrel = randn(length(Imeas),1, 'like', Imeas)*e_systematic_rel;
        end%If previous esys and esysrel exist use them
        Imeas = Imeas.*(1+e_sysrel) + e_sys;
        Imeas = Imeas.*(1+randn(length(Imeas),1, 'like', Imeas)*noise_rel) + randn(length(Imeas),1, 'like', Imeas)*noise_abs*(max(Imeas)-min(Imeas));

    elseif strcmp(mode, 'current')
        [Umeas, Uall] = fsolver.SolveForwardVec(sigma, dirichlet_nodes, dirichlet_vals); %these are the results
        Imeas = fsolver.Iel;
        %add noise and error:
        if ~exist('esys', 'var')%there have been no homogeneous measurements where these have already been calculated
            e_sys    = randn(length(Umeas),1, 'like', Umeas)*e_systematic_abs*(max(Umeas)-min(Umeas));
            e_sysrel = randn(length(Umeas),1, 'like', Umeas)*e_systematic_rel;

            e_sys_all    = randn(size(Uall), 'like', Uall)*e_systematic_abs.*(max(Uall)-min(Uall));
            e_sysrel_all = randn(size(Uall), 'like', Uall)*e_systematic_rel;
        end%If previous esys and esysrel exist use them
        Umeas = Umeas.*(1 + e_sysrel) + e_sys;
        Umeas = Umeas.*(1 + randn(length(Umeas),1, 'like', Umeas)*noise_rel) + randn(length(Umeas),1, 'like', Umeas)*noise_abs*(max(Umeas)-min(Umeas));
        Uall  = Uall.*(1+e_sysrel_all) + e_sys_all;
        Uall  = Uall.*(1+randn(size(Uall), 'like', Uall)*noise_rel) + randn(size(Uall), 'like', Uall)*noise_abs.*(max(Uall)-min(Uall));
    end
    
    %return also the systematic error:
    e = [e_sys e_sysrel];
    if ~exist('Umeas_i', 'var')%reference measurements were not done
        Umeas_i = 0;
        Imeas_i = 0;
    end
    
    
end