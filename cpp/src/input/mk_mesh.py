import numpy as np
import os

ndiv = 500

area_x = 50 * 50
area_y = 50 * ndiv
area_z = 50

nx = 50
ny = ndiv
nz = 1
dof = 3

print("area x[m]:",-area_x/2,area_x/2)
print("area y[m]:",-area_y/2,area_y/2)
print("area z[m]:",0,area_z)

xg = np.linspace(-area_x/2,area_x/2,nx+1,endpoint=True)
yg = np.linspace(-area_y/2,area_y/2,ny+1,endpoint=True)
zg = np.linspace(0,area_z,nz+1,endpoint=True)
z0 = area_z/2

output_space_l = 1000.0
output_space_w = 500.0

output_point_l = np.linspace(0,area_y/2,int(area_y/2/output_space_l)+1)
output_point_w = np.linspace(z0,0,int(z0/output_space_w)+1)


### Set node ###
node = np.empty([len(xg),len(yg),len(zg)],dtype=np.int32)
node_lines = []

inode = 0
for k in range(len(zg)):
    for j in range(len(yg)):
        for i in range(len(xg)):
            dofx,dofy,dofz = 1,1,0

            node[i,j,k] = inode
            node_lines += [ "{} {} {} {} {} {} {}\n".format(inode,xg[i],yg[j],zg[k],dofx,dofy,dofz)]
            inode += 1

# Fault element #
i_fault = nx//2
node_fault = np.empty([1,len(yg),len(zg)],dtype=np.int32)

for k in range(len(zg)):
    for j in range(len(yg)):
        dofx,dofy,dofz = 1,1,0

        node_fault[0,j,k] = inode
        node_lines += [ "{} {} {} {} {} {} {}\n".format(inode,xg[i_fault],yg[j],zg[k],dofx,dofy,dofz)]
        inode += 1


### Set element ###
element_lines = []

ielem = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            im = 0

            style = "3d8solid"

            param_line = "{} {} {} ".format(ielem,style,im)

            if i == i_fault:
                style_line = "{} {} {} {} {} {} {} {}".format(node_fault[0,j,k],node[i+1,j,k],node[i+1,j+1,k],node_fault[0,j+1,k],
                                                                 node_fault[0,j,k+1],node[i+1,j,k+1],node[i+1,j+1,k+1],node_fault[0,j+1,k+1])
            else:
                style_line = "{} {} {} {} {} {} {} {}".format(node[i,j,k],node[i+1,j,k],node[i+1,j+1,k],node[i,j+1,k],
                                                                 node[i,j,k+1],node[i+1,j,k+1],node[i+1,j+1,k+1],node[i,j+1,k+1])
            element_lines += [param_line + style_line + "\n"]
            ielem += 1

# for j in range(ny):
#     for i in range(nx):
#         im = 0
#         style = "2d4visco"
#
#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[i,j,0],node[i,j+1,0],node[i+1,j+1,0],node[i+1,j,0])
#
#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1
#
#         param_line = "{} {} {} ".format(ielem,style,im)
#         style_line = "{} {} {} {}".format(node[i,j,-1],node[i+1,j,-1],node[i+1,j+1,-1],node[i,j+1,-1])
#
#         element_lines += [param_line + style_line + "\n"]
#         ielem += 1

for k in range(nz):
    for j in range(ny):
        im = 0
        style = "2d4visco"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[0,j,k],node[0,j,k+1],node[0,j+1,k+1],node[0,j+1,k])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[-1,j,k],node[-1,j+1,k],node[-1,j+1,k+1],node[-1,j,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

for k in range(nz):
    for i in range(nx):
        im = 0
        style = "2d4visco"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[i,0,k],node[i+1,0,k],node[i+1,0,k+1],node[i,0,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[i,-1,k],node[i,-1,k+1],node[i+1,-1,k+1],node[i+1,-1,k])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

spring_id = np.empty([len(yg),len(zg)],dtype=np.int32)
for k in range(len(zg)):
    for j in range(len(yg)):
        im = 1
        style = "spring"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {}".format(node[i_fault,j,k],node_fault[0,j,k])
        element_lines += [param_line + style_line + "\n"]
        spring_id[j,k] = ielem
        ielem += 1


####################################################
fault_lines = []
fault_lw = []
id_fault = 0
for k in range(nz):
    for j in range(ny):
        im = 0
        style = "2d4fault"

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node[i_fault,j,k],node[i_fault,j+1,k],node[i_fault,j+1,k+1],node[i_fault,j,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem_fault0 = ielem
        ielem += 1

        param_line = "{} {} {} ".format(ielem,style,im)
        style_line = "{} {} {} {}".format(node_fault[0,j,k],node_fault[0,j+1,k],node_fault[0,j+1,k+1],node_fault[0,j,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem_fault1 = ielem
        ielem += 1

        y = (yg[j] + yg[j+1])/2
        z = (zg[k] + zg[k+1])/2
        # if (np.abs(z-area_z/2) <= 1500) and (np.abs(y-area_y/2) <= 1500):
        if np.abs(y) <= 1500:
        # if np.abs(y) <= 500:
            p0 = 81.6e6  # [Pa]
            tp = 81.24e6 # [Pa]
            tr = 63.0e6  # [Pa]
            dc = 0.4     # [m]
        # elif (np.abs(z-area_z/2) <= 7500) and (np.abs(y-area_y/2) <= 7500):
        elif np.abs(y) <= 15000:
        # elif np.abs(y) <= 2000:
            p0 = 70.0e6  # [Pa]
            tp = 81.24e6 # [Pa]
            tr = 63.0e6  # [Pa]
            dc = 0.4     # [m]
        else:
            p0 = 70.0e6  # [Pa]
            tp = 1.e12   # [Pa]
            tr = 1.e12   # [Pa]
            dc = 10.0    # [m]


        l = np.average([yg[j],yg[j+1]])
        w = np.average([zg[k],zg[k+1]])
        fault_lw += [np.array([l,w])]
        # print(id_fault,fault_lw[id_fault])

        fault_lines += ["{} {} {} {} {} {} {} {} {} {} {}\n".format(id_fault,ielem_fault1,ielem_fault0,spring_id[j,k],spring_id[j+1,k],spring_id[j+1,k+1],spring_id[j,k+1],p0,tp,tr,dc)]
        id_fault += 1


nnode = inode       #number of nodes
nelem = ielem       #number of elements
nfault = id_fault   #number of faults

dl = yg[1]-yg[0]
dw = zg[1]-zg[0]

### Set material ###
material_lines = []
material_lines += ["{} {} {} {} {}\n".format(0,"vs_vp_rho",3464.0,6000.0,2670.0)]
material_lines += ["{} {} {} {}\n".format(1,"spring",1.e15,1.e15)]

nmaterial = len(material_lines)

### Set output ###
output_node_lines = []
output_node_lines += ["{}\n".format(node[len(xg)//2,len(yg)//2,0])]
# for i in range(len(xg)):
#     j = len(zg)//2
#     k = 0
#     output_node_lines += ["{}\n".format(node[i,j,k])]


output_element_lines = []
# for i in range(0,nelem-nx-len(zg)):
#     output_element_lines += ["{} \n".format(i)]

# ---- fault output ---- #
dl = yg[1]-yg[0]
dw = zg[1]-zg[0]
print("+++ output fault id (id, y, z)")

output_fault_lines = []
for id_fault in range(nfault):
    l,w = fault_lw[id_fault]
    is_l = np.any(np.logical_and(l-dl/2 <= output_point_l, output_point_l < l+dl/2))
    is_w = np.any(np.logical_and(w-dw/2 <= output_point_w, output_point_w < w+dw/2))

    if (is_l and is_w):
        print(id_fault,l,w)
        output_fault_lines += ["{} \n".format(id_fault)]


output_nnode = len(output_node_lines)
output_nelem = len(output_element_lines)
output_nfault = len(output_fault_lines)

with open("mesh.in","w") as f:
    f.write("{} {} {} {} {} \n".format(nnode,nelem,nfault,nmaterial,dof))
    f.writelines(node_lines)
    f.writelines(element_lines)
    f.writelines(fault_lines)
    f.writelines(material_lines)

with open("output.in","w") as f:
    f.write("{} {} {} \n".format(output_nnode,output_nelem,output_nfault))
    f.writelines(output_node_lines)
    f.writelines(output_element_lines)
    f.writelines(output_fault_lines)
