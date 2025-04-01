import gmsh
import os

def run_gmsh_geo(geo_file):
    gmsh.initialize()
    gmsh.open(geo_file) 

    script_dir = os.path.dirname(os.path.abspath(__file__))
    output_folder = os.path.join(script_dir, "input .msh")
    msh_file = os.path.join(output_folder, os.path.splitext(os.path.basename(geo_file))[0] + ".msh")

    gmsh.write(msh_file)
    print(f"Mesh saved to: {msh_file}")

    gmsh.finalize()

if __name__ == "__main__":
    run_gmsh_geo()