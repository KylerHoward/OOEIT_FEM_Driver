import os

def generate_geo(base_dir):

    # define the base NRRD files directory
    nrrd_dir = os.path.join(base_dir, "NRRD files")

    if not os.path.exists(nrrd_dir):
        print(f"Error: '{nrrd_dir}' not found.")
        return

    # loops through each patient folder inside "NRRD files"
    for patient_id in os.listdir(nrrd_dir):
        patient_path = os.path.join(nrrd_dir, patient_id)

        if not os.path.isdir(patient_path):
            continue

        # looks for "STL files" inside the patient folder
        stl_dir = os.path.join(patient_path, "STL files")
        if not os.path.exists(stl_dir):
            print(f"Warning: 'STL files' folder not found for patient {patient_id}. Skipping.")
            continue

        stl_files = {os.path.splitext(f)[0]: os.path.join(stl_dir, f)
                     for f in os.listdir(stl_dir) if f.endswith(".stl")}

        if not stl_files:
            print(f"Warning: No STL files found for patient {patient_id}. Skipping.")
            continue

        # saves .geo file in the patient ID folder (one level up from "STL files")
        geo_filename = os.path.join(patient_path, f"{patient_id}.geo")
        create_geo_file(stl_files, geo_filename)

def create_geo_file(file_paths, output_filename):

    geo_content = """// Global Resolution
Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthMax = 30;

// Merge all files
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2;
Mesh.SubdivisionAlgorithm = 1;
"""

    # merge STL files
    for name, path in file_paths.items():
        geo_content += f'Merge "{path.replace(os.sep, "/")}";\n'

    geo_content += """Coherence;
//CreateGeometry;

// Add Volume
"""

    # generate Surface Loops and Volumes
    for i, name in enumerate(file_paths.keys(), start=1):
        geo_content += f"Surface Loop({i}) = {{{i}}};\n"
        geo_content += f"Volume({i}) = {{{i}}};\n"

    geo_content += "\n"

    # generate Physical Volumes
    for i, name in enumerate(file_paths.keys(), start=1):
        geo_content += f'Physical Volume("{name}", {i}) = {{{i}}};\n'

    geo_content += """\n// Actually mesh and make quadrature
Mesh 3;
//RecombineMesh;
//RefineMesh;
"""

    # save the .geo file
    with open(output_filename, "w") as geo_file:
        geo_file.write(geo_content)

    print(f".geo file created: {output_filename}")

if __name__ == "__main__":
    script_directory = os.path.dirname(os.path.abspath(__file__))
    generate_gmsh_geo_for_patients(script_directory)