clear
clc
% close all

load mat_free_test.mat

tol   = 1e-12;
maxit = 10;
tic
[L,U] = ilu(A,struct('type','nofill','droptol',1e-6));
% self.L = ichol(self.A);
MF_solvec = gmres(A, b,[],tol,maxit,L,U);
% MF_solvec = gmres(A, b,[],tol,maxit,L,L');
toc


elval = QC*solVec; % Extract the electrode values from the FEM solution (which contains also the potentials inside the domain)
MF_elval = QC*MF_solvec;

figure;
hold on
plot(real(elval(:,1)),'r')
plot(real(MF_elval(:,1)),'b')
legend("Standard", "gmres")
title(sprintf("%d Iterations", maxit))