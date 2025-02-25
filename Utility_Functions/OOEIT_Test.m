clear
clc
close all

% Add relevant subfolders to path
addpath OOEIT
InitializeOOEIT;

% Load the mesh file and create a 1st order mesh object
load('OOEIT/meshfiles/Mesh2D_dense.mat');
mesh = ForwardMesh1st(g, H, elfaces);
    % g: a matrix containing on each row the coordinates of a single node of the mesh
    % H: a matrix containing on each row the indices of the nodes forming a single element of the mesh
    % E: A cell array, where cell l contains a matrix similar to H, but the rows define the boundary elements forming the electrode l

% run the forward problem solver
solver = EITFEM(mesh);

sigma = GenerateEllipse(g, 1, 10, 0.04, 0.05, 0.04, 0, 0.01);
    % Generates a conductivity distribution on the given mesh, having the background 
    % conductivity of 1, and featuring an ellipse with center at point (0.04, 0),
    % conductivity of 10 units, and minor and major axes of 0.04 and 0.05 units. 
    % The last argument is the width of the transition zone between the inclusion and the background
Vmeas = solver.SolveForwardVec(sigma);

