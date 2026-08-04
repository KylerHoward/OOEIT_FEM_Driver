# OOEIT_FEM_Driver
## About
- A 3D FEM driver for the OOEIT package produced by Petri Kuusela.
- Only modified/needed code from the OOEIT package is included in the OOEIT folder.
- Includes a driver script and functions to build electrodes, assign conductivities, and then solve the forward problem.
- Modified for parallelization and matrix free solving of large systems.

## Requirements
- MATLAB compatible with R2024a
- Parallel Computing Toolbox
- Partial Differential Equation Toolbox

## Instructions
1) Clone the repo on your local machine.
2) Open both "FEM3D_Single_Sbj_Driver" and "FEM3D_Function" scripts.
3) Modify flags in the settings block of the driver script.
   1) Recommended to first test the code with "flags.solve_problem" set to 0 to test electrode placement and conductivity assignment.
4) After running the driver, select a GMSH .mat file containing a structure "tetmesh" which has "cell", "node", and "field" variables.
   1) If solving the problem, you will also be prompted to select a location to save the data to.
   2) The code is setup to create an anatomical atlas set of folders if you select a random location.
   3) To avoid this, create a "Results" folder inside your cloned "OOEIT_FEM_Driver" folder and save the results there.
5) If/ElseIf statements on line 109 of the main function script will need to be modified to point to the correct page of the "CT Data Boundaries" excel sheet.
   1) Each dataset of similar meshes should have its own page in the sheet.
   2) Each subject MUST have information in columns A, G, I, and K. (Subject, CarinaFromST_mm, T5FromST_mm, and T8FromST_mm).
   3) These point the code to which subject you're simulating (must match mesh filename), as well as the heights for a patch/belts of electrodes.
5) With these changes, the code should run from here and generate plots of the electrodes, conductivities, and voltages if you are solving the forward problem. 
