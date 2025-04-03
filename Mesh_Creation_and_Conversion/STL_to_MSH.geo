// Mesh size information
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthMin = 10;
Mesh.CharacteristicLengthMax = 100; 

// Merge STL files
Mesh.Algorithm = 8; // (1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad) 
Mesh.RecombinationAlgorithm = 2; // (0: simple, 1: blossom, 2: simple full-quad, 3: blossom full-quad)

Merge "G:/Lungmap_EIT/Radiology_Images/R1053/R1053_Segmentation_NoBones00001.stl"; // Lungs
Merge "G:/Lungmap_EIT/Radiology_Images/R1053/R1053_Segmentation_NoBones00002.stl"; // Trachea
Merge "G:/Lungmap_EIT/Radiology_Images/R1053/R1053_Segmentation_NoBones00003.stl"; // Soft Tissue
Merge "G:/Lungmap_EIT/Radiology_Images/R1053/R1053_Segmentation_NoBones00005.stl"; // Esophagus
Merge "G:/Lungmap_EIT/Radiology_Images/R1053/R1053_Segmentation_NoBones00006.stl"; // Heart

Coherence; // Remove duplicate vertices
CreateGeometry;

// Define distance-based size field
Field[1] = Distance;
Field[1].SurfacesList = {1,2,3,4,5}; // Apply distance field to all surface loops

// Define a threshold function: finer near the surface, coarser inside
Field[2] = Threshold;
Field[2].InField = 1;  // Uses the distance field
Field[2].SizeMin = 200; // Minimum tetrahedra size
Field[2].SizeMax = 200; // Maximum tetrahedra size
Field[2].DistMin = 10;   // Min mesh size used until 10 units away from surface
Field[2].DistMax = 30;  // Max mesh size used beyond 100 units away from mesh surface

Background Field = 2; // Apply the field to the entire mesh

// Define volumes
Surface Loop(1) = {1};
Surface Loop(2) = {2};
Surface Loop(3) = {3};
Surface Loop(4) = {4};
Surface Loop(5) = {5};
Volume(1) = {1};
Volume(2) = {2};
Volume(3) = {3};
Volume(4) = {4};
Volume(5) = {5};
Physical Volume("Lungs", 1) = {1};
Physical Volume("Trachea", 2) = {2};
Physical Volume("Body", 3) = {3};
Physical Volume("Esophagus", 4) = {4};
Physical Volume("Heart", 5) = {5};

// Mesh Generation
RecombineMesh; // Reduce tetrahedral count
Mesh 3;
Mesh.Optimize = 1; 
Mesh.OptimizeNetgen = 1; 

Coherence Mesh; // Remove duplicate faces
Save "R1053_Mesh_NoBones2.msh";
