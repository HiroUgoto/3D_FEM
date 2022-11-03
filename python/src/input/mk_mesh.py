import numpy as np
import os

######################################################################
def uniform_mesh_layout(nx,dx):
    area0 = nx*dx
    xg = np.linspace(-area0/2,area0/2,nx+1,endpoint=True)
    return xg

def uniform_mesh_layout_z(nz,dz,zh):
    area0 = nz*dz

    if -area0/2+zh < 0:
        zg = np.linspace(0,area0,nz+1,endpoint=True)
    else:
        zg = np.linspace(-area0/2+zh,area0/2+zh,nz+1,endpoint=True)

    return zg

def variable_mesh_layout(nx,dx0,x0):
    # nx : number of elements
    # dx0: uniform mesh size
    # x0 : half width of uniform mesh
    xg = np.zeros(nx+1)
    area0 = nx*dx0
    xg0 = np.linspace(-area0/2,area0/2,nx+1,endpoint=True)

    xg1 = (xg0**2 + x0**2) / (2*x0)
    for i in range(nx+1):
        if xg0[i] < -x0:
            xg[i] = -xg1[i]
        elif xg0[i] < x0:
            xg[i] = xg0[i]
        else:
            xg[i] = xg1[i]

    return xg

def variable_mesh_layout_z(nz,dz0,z0,zh):
    # nz : number of elements
    # dz0: uniform mesh size
    # z0 : half width of uniform mesh
    # zh : hypocenter depth
    zg = np.zeros(nz+1)
    area0 = nz*dz0
    zg0 = np.linspace(0,area0,nz+1,endpoint=True)

    zg1 = (zg0**2 + (zh+z0)**2) / (2*(zh+z0))
    for i in range(nz+1):
        if zg0[i] < (zh+z0):
            zg[i] = zg0[i]
        else:
            zg[i] = zg1[i]

    return zg

######################################################################


nx = 4
ny = 20
nz = 10

dx = 100
dy = 100
dz = 100

# zh = 7500
zh = 3000

x0 = 500
xg = variable_mesh_layout(nx,dx,x0)
# xg = uniform_mesh_layout(nx,dx)
yg = uniform_mesh_layout(ny,dy)
zg = uniform_mesh_layout_z(nz,dz,zh)

area_x = xg[-1]-xg[0]
area_y = yg[-1]-yg[0]
area_z = zg[-1]-zg[0]

dof = 3

print("area x[m]:",xg[0],xg[-1])
print("area y[m]:",yg[0],yg[-1])
print("area z[m]:",zg[0],zg[-1])

# output_space_l = 1000.0
output_space_w = 1000.0

# nl = int(area_y/2/output_space_l)
nl = 1
nw = int(zh/output_space_w)
# nw = 1

# output_point_l = np.linspace(0,nl*output_space_l,nl+1)
output_point_l = np.array([0.0])
output_point_w = np.linspace(zh,zh-nw*output_space_w,nw+1)
# output_point_w = np.array([zh])

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

# Fault element #
i_fault = nx//2
node_fault = np.empty([1,len(yg),len(zg)],dtype=np.int32)

for k in range(len(zg)):
    for j in range(len(yg)):
        dofx,dofy,dofz = 1,1,1

        node_fault[0,j,k] = inode
        node_lines += [ "{} {} {} {} {} {} {}\n".format(inode,xg[i_fault],yg[j],zg[k],dofx,dofy,dofz)]
        inode += 1

### Set element ###
element_lines = []

fault_neighbour_element = np.empty([2,ny,nz],dtype=np.int32)

ielem = 0
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            im = 0

            style = "3d8solid"
            param_line = "{} {} {} ".format(ielem,style,im)

            if i == i_fault-1:
                fault_neighbour_element[0,j,k] = ielem
            elif i == i_fault:
                fault_neighbour_element[1,j,k] = ielem

            if i == i_fault:
                style_line = "{} {} {} {} {} {} {} {}".format(node_fault[0,j,k],node[i+1,j,k],
                                                              node[i+1,j+1,k],node_fault[0,j+1,k],
                                                              node_fault[0,j,k+1],node[i+1,j,k+1],
                                                              node[i+1,j+1,k+1],node_fault[0,j+1,k+1])
            else:
                style_line = "{} {} {} {} {} {} {} {}".format(node[i,j,k],node[i+1,j,k],node[i+1,j+1,k],node[i,j+1,k],
                                                                 node[i,j,k+1],node[i+1,j,k+1],node[i+1,j+1,k+1],node[i,j+1,k+1])
            element_lines += [param_line + style_line + "\n"]
            ielem += 1

for j in range(ny):
    for i in range(nx):
        im = 0
        style = "2d4visco"

        # param_line = "{} {} {} ".format(ielem,style,im)
        # style_line = "{} {} {} {}".format(node[i,j,0],node[i,j+1,0],node[i+1,j+1,0],node[i+1,j,0])
        #
        # element_lines += [param_line + style_line + "\n"]
        # ielem += 1

        param_line = "{} {} {} ".format(ielem,style,im)
        if i == i_fault:
            style_line = "{} {} {} {}".format(node_fault[0,j,-1],node[i+1,j,-1],node[i+1,j+1,-1],node_fault[0,j+1,-1])
        else:
            style_line = "{} {} {} {}".format(node[i,j,-1],node[i+1,j,-1],node[i+1,j+1,-1],node[i,j+1,-1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

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
        if i == i_fault:
            style_line = "{} {} {} {}".format(node_fault[0,0,k],node[i+1,0,k],node[i+1,0,k+1],node_fault[0,0,k+1])
        else:
            style_line = "{} {} {} {}".format(node[i,0,k],node[i+1,0,k],node[i+1,0,k+1],node[i,0,k+1])

        element_lines += [param_line + style_line + "\n"]
        ielem += 1

        param_line = "{} {} {} ".format(ielem,style,im)
        if i == i_fault:
            style_line = "{} {} {} {}".format(node_fault[0,-1,k],node_fault[0,-1,k+1],node[i+1,-1,k+1],node[i+1,-1,k])
        else:
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
fault_l,fault_w = [],[]
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
        if (np.abs(z-zh) <= 1500) and (np.abs(y) <= 1500):
        # if np.abs(y) <= 1500:
            p0 = 81.6e6  # [Pa]
            tp = 81.24e6 # [Pa]
            tr = 63.0e6  # [Pa]
            dc = 0.4     # [m]
        # elif (np.abs(z-zh) <= 7500) and (np.abs(y) <= 15000):
        elif (np.abs(z-zh) <= 3000) and (np.abs(y) <= 6000):
        # elif np.abs(y) <= 15000:
            p0 = 70.0e6  # [Pa]
            tp = 81.24e6 # [Pa]
            tr = 63.0e6  # [Pa]
            dc = 0.4     # [m]
        else:
            p0 = 70.0e6  # [Pa]
            tp = 1.e12   # [Pa]
            tr = 1.e12   # [Pa]
            dc = 10.0    # [m]


        fault_l += [np.array([yg[j],yg[j+1]])]
        fault_w += [np.array([zg[k],zg[k+1]])]
        # print(id_fault,fault_lw[id_fault])

        fault_lines += ["{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(id_fault,ielem_fault1,ielem_fault0, \
            fault_neighbour_element[0,j,k],fault_neighbour_element[1,j,k], \
            spring_id[j,k],spring_id[j+1,k],spring_id[j+1,k+1],spring_id[j,k+1], \
            p0,tp,tr,dc)]
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
print("+++ output fault id (id, y, z)")

output_fault_lines = []
for id_fault in range(nfault):
    l0,l1 = fault_l[id_fault]
    w0,w1 = fault_w[id_fault]
    is_l = np.any(np.logical_and(l0 <= output_point_l, output_point_l < l1))
    is_w = np.any(np.logical_and(w0 <= output_point_w, output_point_w < w1))

    if (is_l and is_w):
        print(id_fault,(l0+l1)/2,(w0+w1)/2)
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
