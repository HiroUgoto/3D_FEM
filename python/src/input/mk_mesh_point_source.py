import numpy as np
import os

area_x = 1200.0
area_y = 1200.0
area_z = 1200.0

nx = 15
ny = 15
nz = 15
dof = 3

# Perfectly Matched Layer
pml_layers = 5
pml_order = 2
pml_logR = -5   # reflection coeff: R = 10^(-5)

pml_area_x = area_x*(1 + 2*pml_layers/nx)
pml_area_y = area_y*(1 + 2*pml_layers/ny)
pml_area_z = area_z*(1 +   pml_layers/nz)

pml_nx = nx + 2*pml_layers
pml_ny = ny + 2*pml_layers
pml_nz = nz +   pml_layers

xg,dx = np.linspace(-pml_area_x/2,pml_area_x/2,pml_nx+1,endpoint=True,retstep=True)
yg,dy = np.linspace(-pml_area_y/2,pml_area_y/2,pml_ny+1,endpoint=True,retstep=True)
zg,dz = np.linspace(0,pml_area_z,pml_nz+1,endpoint=True,retstep=True)


### Set node ###
node = np.empty([len(xg),len(yg),len(zg)],dtype=np.int32)
node_lines = []

inode = 0
for k in range(len(zg)):
    for j in range(len(yg)):
        for i in range(len(xg)):
            dofx,dofy,dofz = 1,1,1

            node[i,j,k] = inode
            node_lines += [ "{} {} {} {} {} {} {}\n".format(inode,xg[i],yg[j],zg[k],dofx,dofy,dofz)]
            inode += 1

### Set element ###
element_lines = []
pml_lines = []

# 8-Cubic element 
v = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]

ielem = 0
for k in range(pml_nz):
    for j in range(pml_ny):
        for i in range(pml_nx):
            im = 0

            style = "3d8solid"
            param_line = "{} {} {} ".format(ielem,style,im)
            n_list = [node[i+v[l][0],j+v[l][1],k+v[l][2]] for l in range(8)]
            style_line = " ".join(map(str, n_list))

            element_lines += [param_line + style_line + "\n"]

            px,py,pz = 0,0,0
            if i < pml_layers:
                px = (i - pml_layers)*dx
            elif i >= nx + pml_layers:
                px = (i - (nx+pml_layers) + 1)*dx
            if j < pml_layers:
                py = (j - pml_layers)*dy
            elif j >= ny + pml_layers:
                py = (j - (ny+pml_layers) + 1)*dy
            if k >= nz:
                pz = (k - nz + 1)*dz

            if (px != 0) or (py != 0) or (pz != 0):
                pml_lines += ["{} {:.6f} {:.6f} {:.6f}\n".format(ielem,px,py,pz)]

            ielem += 1


# for j in range(ny):
#     for i in range(nx):
#         im = 0
#         style = "2d4visco"

#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[i,j,-1],node[i+1,j,-1],node[i+1,j+1,-1],node[i,j+1,-1])

#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1

# for k in range(nz):
#     for j in range(ny):
#         im = 0
#         style = "2d4visco"

#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[0,j,k],node[0,j,k+1],node[0,j+1,k+1],node[0,j+1,k])

#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1

#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[-1,j,k],node[-1,j+1,k],node[-1,j+1,k+1],node[-1,j,k+1])

#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1

# for k in range(nz):
#     for i in range(nx):
#         im = 0
#         style = "2d4visco"

#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[i,0,k],node[i+1,0,k],node[i+1,0,k+1],node[i,0,k+1])

#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1

#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[i,-1,k],node[i,-1,k+1],node[i+1,-1,k+1],node[i+1,-1,k])

#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1


nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {}\n".format(0,"vs_vp_rho",1000.0,2500.0,2100.0)]

nmaterial = len(material_lines)


### Set PML ###
pml_nelem = len(element_lines)

### Set output ###
output_node_lines = []
output_node_lines += ["{}\n".format(node[len(xg)//2,len(yg)//2,0])]
output_node_lines += ["{}\n".format(node[len(xg)//2,len(yg)-1,0])]

output_element_lines = []
# for i in range(0,nelem-nx-len(zg)):
#     output_element_lines += ["{} \n".format(i)]

output_nnode = len(output_node_lines)
output_nelem = len(output_element_lines)


with open("mesh.in","w") as f:
    f.write("{} {} {} {} \n".format(nnode,nelem,nmaterial,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
    f.writelines(material_lines)

with open("output.in","w") as f:
    f.write("{} {} \n".format(output_nnode,output_nelem))
    f.writelines(output_node_lines)
    f.writelines(output_element_lines)

with open("pml.in","w") as f:
    f.write("{} {} {} {}\n".format(pml_nelem,pml_layers,pml_order,pml_logR)) 
    f.writelines(pml_lines)