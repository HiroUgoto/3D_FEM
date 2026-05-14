import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib.gridspec as gridspec
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

def plot_tetra_element(ax,nodes,fc):
    poly = [[nodes[0].xyz,nodes[1].xyz,nodes[2].xyz],
            [nodes[0].xyz,nodes[1].xyz,nodes[3].xyz],
            [nodes[0].xyz,nodes[2].xyz,nodes[3].xyz],
            [nodes[1].xyz,nodes[2].xyz,nodes[3].xyz]]

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
def plot_mesh(fem,fmt="cube"):
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

            if fmt == "cube":
                plot_element(ax,element.nodes,fc=pc[ic])
            elif fmt == "tetra":
                plot_tetra_element(ax,element.nodes,fc=pc[ic])

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

#--------------------------------------------------------#
def _estimate_tol(fem, axis_idx):
    coords = np.array([n.xyz[axis_idx] for n in fem.nodes])
    unique_coords = np.unique(np.round(coords, decimals=4))
    if len(unique_coords) > 1:
        min_dx = np.min(np.diff(unique_coords))
        return min_dx * 0.55
    return 0.1

def _get_slice_info(fem, axis, target_coord, tol):
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    slice_idx = axis_map[axis.lower()]
    if tol == 'auto': tol = _estimate_tol(fem, slice_idx)

    if slice_idx == 0:   idx_h, idx_v = 1, 2 # X断面 -> Y, Z
    elif slice_idx == 1: idx_h, idx_v = 0, 2 # Y断面 -> X, Z
    else:                idx_h, idx_v = 1, 0 # Z断面 -> Y, X

    slice_h, slice_v, node_indices = [], [], []
    for i, node in enumerate(fem.nodes):
        if abs(node.xyz[slice_idx] - target_coord) < tol:
            slice_h.append(node.xyz[idx_h]); slice_v.append(node.xyz[idx_v])
            node_indices.append(i)

    if len(slice_h) < 3: return None
    slice_h, slice_v = np.array(slice_h), np.array(slice_v)
    return {
        'triang': mtri.Triangulation(slice_h, slice_v),
        'node_indices': node_indices,
        'xlim': [np.min(slice_h), np.max(slice_h)],
        'ylim': [np.min(slice_v), np.max(slice_v)],
        'idx_h': idx_h, 'idx_v': idx_v
    }


def plot_multislice_update_init(fem, targets={'x':0.0, 'y':0.0, 'z':0.0}, tol='auto'):
    all_xyz = np.array([n.xyz for n in fem.nodes])
    mins, maxs = all_xyz.min(axis=0), all_xyz.max(axis=0)

    x_range = maxs[0] - mins[0]
    y_range = maxs[1] - mins[1]
    z_range = maxs[2] - mins[2]

    plt.ion() 
    fig = plt.figure(figsize=(12, 10))

    gs = gridspec.GridSpec(2, 2, height_ratios=[z_range, x_range], hspace=0.3, wspace=0.3)

    axs = {
        'x': fig.add_subplot(gs[0, 0]),
        'y': fig.add_subplot(gs[0, 1]),
        'z': fig.add_subplot(gs[1, 0]),
        '3d': fig.add_subplot(gs[1, 1], projection='3d')
    }

    config = {
        'slice_data': {ax: _get_slice_info(fem, ax, targets[ax], tol) for ax in ['x', 'y', 'z']},
        'targets': targets,
        'hex_colors': {'x': '#FF0000', 'y': '#008000', 'z': '#0000FF'},
        'auto_vmax': 0.0,  
    }

    ax3d = axs['3d']
    ax3d.set_xlim([mins[0], maxs[0]]); ax3d.set_ylim([mins[1], maxs[1]]); ax3d.set_zlim([mins[2], maxs[2]])
    ax3d.invert_zaxis()
    ax3d.invert_yaxis() 
    ax3d.set_box_aspect(maxs - mins)
    ax3d.view_init(elev=35, azim=200)
    ax3d.set_xlabel('X (North)'); ax3d.set_ylabel('Y (East)'); ax3d.set_zlabel('Z (Depth)')

    # 断面のプレビュー表示
    for ax_key, color in config['hex_colors'].items():
        if config['slice_data'][ax_key] is None: continue
        if ax_key == 'x':
            yy, zz = np.meshgrid([mins[1], maxs[1]], [mins[2], maxs[2]])
            ax3d.plot_surface(np.full_like(yy, targets['x']), yy, zz, alpha=0.2, color=color, shade=False)
        elif ax_key == 'y':
            xx, zz = np.meshgrid([mins[0], maxs[0]], [mins[2], maxs[2]])
            ax3d.plot_surface(xx, np.full_like(xx, targets['y']), zz, alpha=0.2, color=color, shade=False)
        elif ax_key == 'z':
            xx, yy = np.meshgrid([mins[0], maxs[0]], [mins[1], maxs[1]])
            ax3d.plot_surface(xx, yy, np.full_like(xx, targets['z']), alpha=0.2, color=color, shade=False)

    fig.show()
    fig.canvas.draw() 

    return fig, axs, config


def plot_multislice_update(fig, axs, fem, config, comp='norm', fin=False):
    axis_names = ['X (North)', 'Y (East)', 'Z (Depth)']

    all_vals = []  
    for ax_key, sd in config['slice_data'].items():
        if sd is None: continue
        if comp.lower() == 'norm':
            vals = np.array([np.linalg.norm(fem.nodes[i].v) for i in sd['node_indices']])
        else:
            idx = {'x':0, 'y':1, 'z':2}.get(comp.lower(), 1)
            vals = np.abs([fem.nodes[i].v[idx] for i in sd['node_indices']])
        if len(vals) > 0: all_vals.append(vals)

    if all_vals:
        step_vmax = np.max(np.concatenate(all_vals))
        config['auto_vmax'] = step_vmax

    vmax_plot = max(config['auto_vmax'], 1e-10)

    for ax_key in ['x', 'y', 'z']:
        ax, sd = axs[ax_key], config['slice_data'][ax_key]
        if sd is None: continue

        ax.cla(); ax.grid(True)
        
        if comp.lower() == 'norm':
            val = np.array([np.linalg.norm(fem.nodes[i].v) for i in sd['node_indices']])
            vmin, cmap = 0.0, "BuPu"
        else:
            idx = {'x':0, 'y':1, 'z':2}.get(comp.lower(), 1)
            val = np.array([fem.nodes[i].v[idx] for i in sd['node_indices']])
            vmin, cmap = -vmax_plot, "PiYG"

        ax.tripcolor(sd['triang'], val, vmin=vmin, vmax=vmax_plot, cmap=cmap, shading='gouraud')
        
        ax.set_xlim(sd['xlim'])
        ax.set_ylim([sd['ylim'][1], sd['ylim'][0]]) if sd['idx_v'] == 2 else ax.set_ylim(sd['ylim'])
        ax.set_aspect('equal')
        
        comp_str = "(vector norm)" if comp.lower() == 'norm' else f"({comp.upper()} comp.)"
        ax.set_title(f"{ax_key.upper()} Slice {comp_str}\n(at {ax_key}={config['targets'][ax_key]})", 
                     color=config['hex_colors'][ax_key], fontweight='bold')
        ax.set_xlabel(axis_names[sd['idx_h']]); ax.set_ylabel(axis_names[sd['idx_v']])

    if fin: 
        plt.ioff()
        plt.show()
    else: 
        # plt.pause(0.001)
        fig.canvas.draw()
        fig.canvas.flush_events()
    