import numpy as np
import fem
import node,element,material

# ------------------------------------------------------------------- #
def input_mesh(mesh_file):
    with open(mesh_file) as f:
        lines = f.readlines()

        nnode,nelem,nmaterial,dof = [int(s) for s in lines[0].split()]

        irec = 1
        nodes = [None] * nnode
        for inode in range(nnode):
            items = lines[inode+irec].split()

            id = int(items[0])
            xyz = np.array([float(s) for s in items[1:4]])
            freedom = np.array([int(s) for s in items[4:]])

            nodes[inode] = node.Node(id,xyz,freedom)

        irec += nnode
        elements = [None] * nelem
        for ielem in range(nelem):
            items = lines[ielem+irec].split()

            id = int(items[0])
            style = items[1]
            material_id = int(items[2])
            inode = np.array([int(s) for s in items[3:]])

            elements[ielem] = element.Element(id,style,material_id,inode)

        irec += nelem
        materials = [None] * nmaterial
        for imaterial in range(nmaterial):
            items = lines[imaterial+irec].split()

            id = int(items[0])
            style = items[1]
            param = np.array([float(s) for s in items[2:]])

            materials[imaterial] = material.Material(id,style,param)

    return fem.Fem(dof,nodes,elements,materials)

# ------------------------------------------------------------------- #
def input_outputs(output_file):
    with open(output_file) as f:
        lines = f.readlines()

        nnode,nelem = [int(s) for s in lines[0].split()]

        irec = 1
        nodes = [None] * nnode
        for inode in range(nnode):
            items = lines[inode+irec].split()

            id = int(items[0])
            nodes[inode] = id

        irec += nnode
        elements = [None] * nelem
        for ielem in range(nelem):
            items = lines[ielem+irec].split()

            id = int(items[0])
            elements[ielem] = id

    return nodes, elements
