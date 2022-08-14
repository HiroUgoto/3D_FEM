import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d

#--------------------------------------------------------#
def plot_element(ax,nodes,fc):
    poly = [[nodes[0].xyz,nodes[1].xyz,nodes[2].xyz,nodes[3].xyz],
            [nodes[4].xyz,nodes[5].xyz,nodes[6].xyz,nodes[7].xyz],
            [nodes[0].xyz,nodes[1].xyz,nodes[5].xyz,nodes[4].xyz],
            [nodes[1].xyz,nodes[2].xyz,nodes[6].xyz,nodes[5].xyz],
            [nodes[2].xyz,nodes[3].xyz,nodes[7].xyz,nodes[6].xyz],
            [nodes[3].xyz,nodes[0].xyz,nodes[4].xyz,nodes[7].xyz]]

    ax.add_collection3d(art3d.Poly3DCollection(poly,ec="k",fc=fc,alpha=0.2))

def plot_element_disp(ax,nodes,amp,fc):
    poly = [[nodes[0].xyz+nodes[0].u*amp,nodes[1].xyz+nodes[1].u*amp,nodes[2].xyz+nodes[2].u*amp,nodes[3].xyz+nodes[3].u*amp],
            [nodes[4].xyz+nodes[4].u*amp,nodes[5].xyz+nodes[5].u*amp,nodes[6].xyz+nodes[6].u*amp,nodes[7].xyz+nodes[7].u*amp],
            [nodes[0].xyz+nodes[0].u*amp,nodes[1].xyz+nodes[1].u*amp,nodes[5].xyz+nodes[5].u*amp,nodes[4].xyz+nodes[4].u*amp],
            [nodes[1].xyz+nodes[1].u*amp,nodes[2].xyz+nodes[2].u*amp,nodes[6].xyz+nodes[6].u*amp,nodes[5].xyz+nodes[5].u*amp],
            [nodes[2].xyz+nodes[2].u*amp,nodes[3].xyz+nodes[3].u*amp,nodes[7].xyz+nodes[7].u*amp,nodes[6].xyz+nodes[6].u*amp],
            [nodes[3].xyz+nodes[3].u*amp,nodes[0].xyz+nodes[0].u*amp,nodes[4].xyz+nodes[4].u*amp,nodes[7].xyz+nodes[7].u*amp]]

    ax.add_collection3d(art3d.Poly3DCollection(poly,ec="k",fc=fc,alpha=0.2))

#--------------------------------------------------------#
def plot_mesh(fem):
    pc = ["gray","yellow","green","pink"]

    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111,projection='3d')

    x = [node.xyz[0] for node in fem.nodes]
    y = [node.xyz[1] for node in fem.nodes]
    z = [node.xyz[2] for node in fem.nodes]

    area_x = max(x)-min(x)
    area_y = max(y)-min(y)
    area_z = max(z)-min(z)

    ax.set_xlim([min(x)-0.1*area_x,max(x)+0.1*area_x])
    ax.set_ylim([min(y)-0.1*area_y,max(y)+0.1*area_y])
    ax.set_zlim([max(z)+0.1*area_z,min(z)-0.25*area_z])
    ax.set_box_aspect((1,1,1))

    for element in fem.elements:
        if element.dim == 3:
            ic = element.material_id % len(pc)
            plot_element(ax,element.nodes,fc=pc[ic])

    for node in fem.nodes:
        ax.scatter(node.xyz[0],node.xyz[1],node.xyz[2],color="k")

    plt.show()

#--------------------------------------------------------#
def plot_mesh_update_init():
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_axisbelow(True)
    return ax

def plot_mesh_update(ax,fem,amp=1.0,fin=False):
    pc = ["gray","yellow","green","pink"]

    ax.cla()
    ax.grid()

    x = [node.xyz[0] for node in fem.nodes]
    y = [node.xyz[1] for node in fem.nodes]
    z = [node.xyz[2] for node in fem.nodes]

    area_x = max(x)-min(x)
    area_y = max(y)-min(y)
    area_z = max(z)-min(z)

    ax.set_xlim([min(x)-0.1*area_x,max(x)+0.1*area_x])
    ax.set_ylim([min(y)-0.1*area_y,max(y)+0.1*area_y])
    ax.set_zlim([max(z)+0.1*area_z,min(z)-0.25*area_z])
    ax.set_box_aspect((1,1,1))

    for element in fem.elements:
        if element.dim == 3:
            ic = element.material_id % len(pc)
            plot_element_disp(ax,element.nodes,amp,fc=pc[ic])

    for node in fem.nodes:
        ax.scatter(node.xyz[0]+node.u[0]*amp,node.xyz[1]+node.u[1]*amp,node.xyz[2]+node.u[2]*amp,color="k")

    if fin:
        plt.show()
    else:
        plt.pause(0.001)
