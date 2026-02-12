import numpy as np
import os

area_x = 1200.0
area_y = 1200.0
area_z = 600.0

nx = 16
ny = 16
nz = 8
dof = 3

xg = np.linspace(-area_x/2,area_x/2,nx+1,endpoint=True)
yg = np.linspace(-area_y/2,area_y/2,ny+1,endpoint=True)
zg = np.linspace(0,area_z,nz+1,endpoint=True)


### Set node ###
node = np.empty([len(xg),len(yg),len(zg)],dtype=np.int32)
node_lines = []

inode = 0
for k in range(len(zg)):
    for j in range(len(yg)):
        for i in range(len(xg)):
            dofx,dofy,dofz = 1,1,1

            if i == 0:
                dofz = 0
            elif i == len(xg)-1:
                dofz = 0    
            if j == 0:
                dofz = 0
            elif j == len(yg)-1:
                dofz = 0    
            if k == len(zg)-1:
                dofz = 0

            node[i,j,k] = inode
            node_lines += [ "{} {} {} {} {} {} {}\n".format(inode,xg[i],yg[j],zg[k],dofx,dofy,dofz)]
            inode += 1


### Set element ###
element_lines = []

# 8-Cubic element 
v = [[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]]

ielem = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            im = 0
            style = "3d8solid"
            param_line = "{} {} {} ".format(ielem,style,im)
            n_list = [node[i+v[l][0],j+v[l][1],k+v[l][2]] for l in range(8)]
            style_line = " ".join(map(str, n_list))

            # style_line = "{} {} {} {} {} {} {} {}".format(node[i,j,k],node[i+1,j,k],node[i+1,j+1,k],node[i,j+1,k],
            #                                                 node[i,j,k+1],node[i+1,j,k+1],node[i+1,j+1,k+1],node[i,j+1,k+1])
            element_lines += [param_line + style_line + "\n"]
            ielem += 1

# 6-Tetrahedra decomposition
# t0 = [[0,0,0],[2,0,0],[2,2,2],[2,0,2],[1,0,0],[2,1,1],[1,1,1],[1,0,1],[2,1,2],[2,0,1]]
# t1 = [[0,0,0],[0,0,2],[2,0,2],[2,2,2],[0,0,1],[1,0,2],[1,0,1],[1,1,1],[2,1,2],[1,1,2]]
# t2 = [[0,0,0],[0,0,2],[2,2,2],[0,2,2],[0,0,1],[1,1,2],[1,1,1],[0,1,1],[1,2,2],[0,1,2]]
# t3 = [[0,0,0],[0,2,0],[0,2,2],[2,2,2],[0,1,0],[0,2,1],[0,1,1],[1,1,1],[1,2,2],[1,2,1]]
# t4 = [[0,0,0],[0,2,0],[2,2,2],[2,2,0],[0,1,0],[1,2,1],[1,1,1],[1,1,0],[2,2,1],[1,2,0]]
# t5 = [[0,0,0],[2,0,0],[2,2,0],[2,2,2],[1,0,0],[2,1,0],[1,1,0],[1,1,1],[2,2,1],[2,1,1]]
# t = [t0, t1, t2, t3, t4, t5]

# ielem = 0 
# for k in range(0,nz,2):
#     for j in range(0,ny,2):
#         for i in range(0,nx,2):
#             im = 0
#             style = "3d10solid"

#             for it in range(len(t)):
#                 param_line = "{} {} {} ".format(ielem,style,im)
#                 n_list = [node[i+t[it][l][0],j+t[it][l][1],k+t[it][l][2]] for l in range(10)]
#                 style_line = " ".join(map(str, n_list))
#                 element_lines += [param_line + style_line + "\n"]
#                 ielem += 1


for j in range(ny):
    for i in range(nx):
        im = 0
        style = "2d4input"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[i,j,-1],node[i+1,j,-1],node[i+1,j+1,-1],node[i,j+1,-1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

# for k in range(nz):
#     for j in range(ny):
#         im = -1
#         style = "connect"

#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {}".format(node[0,j,k],node[-1,j,k])

#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1

nnode = inode       #number of nodes
nelem = ielem       #number of elements


### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {}\n".format(0,"vs_vp_rho",1000.0,2500.0,2100.0)]
# material_lines += ["{} {} {} {} {} {} {}\n".format(0,"vs_vp_rho_Q",1000.0,2500.0,2100.0,100,200)]

nmaterial = len(material_lines)

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
