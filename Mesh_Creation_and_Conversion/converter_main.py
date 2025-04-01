import os
from geo_creator import generate_geo
from run_gmsh import run_gmsh_geo
from msh_to_mat import msh_to_mat

#------------------------ CREATING .GEO FILES ------------------------    
def geo_creator():
    script_directory = os.path.dirname(os.path.abspath(__file__))
    generate_geo(script_directory)

#--------------------------- RUNNING GMSH ---------------------------    

def process_geo_files():
    """runs created .geo files through GMSH"""
    script_dir = os.path.dirname(os.path.abspath(__file__))  
    parent_folder = os.path.join(script_dir, "NRRD Files")  # path to NRRD Files
    output_folder = os.path.join(script_dir, "input .msh")  # path to input .msh

    if not os.path.exists(parent_folder):
        print(f"Error: The folder '{parent_folder}' does not exist.")
        return

    # loops through each subfolder inside 'NRRD Files'
    for folder_name in os.listdir(parent_folder):
        folder_path = os.path.join(parent_folder, folder_name)

        if os.path.isdir(folder_path):
            # finds all .geo files in the subfolder
            geo_files = [f for f in os.listdir(folder_path) if f.endswith(".geo")]

            if not geo_files:
                print(f"No .geo files found in '{folder_path}', skipping.")
                continue

            # processes each file
            for geo_file in geo_files:
                geo_path = os.path.join(folder_path, geo_file)
                print(f"Processing: {geo_path}")
                run_gmsh_geo(geo_path)

    print("Meshing complete!")


#---------------------------- MSH TO MAT ----------------------------      

def run_batch_conversion():
    """searches for all .msh files in 'input .msh' and converts them to .mat in 'output .mat'."""
    
    script_dir = os.path.dirname(os.path.abspath(__file__)) 
    input_folder = os.path.join(script_dir, "input .msh")  # path to input .msh
    output_folder = os.path.join(script_dir, "output .mat")  # path to output .mat

    # ensures the output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # get list of .msh files in the input directory
    msh_files = [f for f in os.listdir(input_folder) if f.endswith(".msh")]

    if not msh_files:
        print("No .msh files found in 'input .msh' folder!")
        return

    print(f"Found {len(msh_files)} .msh files. Converting...")

    for msh_file in msh_files:
        # define full paths
        msh_path = os.path.join(input_folder, msh_file)
        mat_filename = os.path.splitext(msh_file)[0] + ".mat" 

        msh_to_mat(msh_path, output_folder)
        print(f"Saved: {mat_filename}\n")

    print(".msh batch conversion complete!")



if __name__ == "__main__":
    script_directory = os.path.dirname(os.path.abspath(__file__))
    #geo_creator()
    #process_geo_files()
    run_batch_conversion()
