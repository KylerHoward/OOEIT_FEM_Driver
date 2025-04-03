import scipy.io as sp
import numpy as np
import os

def msh_to_mat(filename, output_folder):
    mat_filename = os.path.join(output_folder, os.path.splitext(os.path.basename(filename))[0].replace("_Segmentation_", "_Mesh_") + ".mat")
    print("Converting", filename, "to", mat_filename)
    with open(filename, 'r') as f:
        curr_section = "Header"
        node = []
        cell = []
        field = []
        for line in f:
            line = line.strip()

            if line.startswith("$Nodes"):
                # extract one field from here - node in tetmesh
                curr_section = "Nodes"
                print("Reached Nodes")
                continue

            elif line.startswith("$Elements"):
                # extract cell (tetrahedra) and field (labels)
                curr_section = "Elements"
                print("Reached Elements")
                next(f)
                continue

            if curr_section == "Header":
                continue

            elif curr_section == "Nodes":
                line_sections = line.split(' ')

                if len(line_sections) == 3:
                    node.append([float(x) for x in line_sections])

            elif curr_section == "Elements":
                line_sections = line.split(' ')

                #if we are at a new elements header line
                if len(line_sections) == 4:
                    print("In Elements, label", line_sections[1])
                    field += [int(line_sections[1])] * int(line_sections[3])

                #if we are at a normal line with the tetrahedra
                elif len(line_sections) == 5:
                    cell.append([int(x) for x in line_sections[1:]])

    np_node = np.array(node, dtype=np.float32).transpose()
    np_cell = np.array(cell, dtype=np.int32).transpose()
    np_field = np.array(field, dtype=np.uint8).transpose()

    print(f"Node shape: ({np_node.shape})") # should be (3, N)
    print("Node sample:", node[:5])

    print(f"Cell shape: ({np_cell.shape})") # should be (4, M)
    print("Cell sample:", cell[:5])

    print(f"Field shape: ({np_field.shape})") # should be (1, M)
    print("Field sample:", field[:5])

    tetmesh = {
        "node": np_node,
        "cell": np_cell,
        "field": np_field
    }

    print("\nSaving .mat file...")
    sp.savemat(mat_filename, {"tetmesh": tetmesh})
    print(".mat file saved")
    return mat_filename


if __name__ == "__main__": # test file, no driver 
    msh_to_mat("input .msh\R1053_Segmentation_NoBones.msh")
    print("Done!")