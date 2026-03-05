import os
# os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import matplotlib.pyplot as plt

import time

import io_data
import input_wave
import source
import plot_model
import vtk

import datetime
from scipy import interpolate

start = time.time()

## --- Input FEM Mesh --- ##
fem = io_data.input_mesh("input/mesh.in")
outputs = io_data.input_outputs("input/output.in")
output_dir = "result/"

## --- FEM Set up --- ##
fem.set_init(n_threads=8)
fem.set_output(outputs)
# plot_model.plot_mesh(fem)
# exit()

## --- Input --- ##
input_type = "double_couple"
# input_type = "plane_wave"

input_time_series = "function"
# input_time_series = "datafile"

# --------------------- #
if input_time_series == "function":
    fsamp = 100
    fp = 1.0
    duration = 5.0

    tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
    # input_acc = input_wave.ricker(tim,fp,tp=1.5/fp,amp=1.0)
    input_acc = input_wave.diff_gauss(tim,fp,tp=1.5/fp,amp=1.0)
    ntim = len(tim)

elif input_time_series == "datafile":
    # displacement (slip) definition
    input_wave_file = "input/scaled_input_STF.txt"
    tim,input_disp = np.loadtxt(input_wave_file,skiprows=1,unpack=True)
    dt = tim[1] - tim[0]
    
    input_vel = np.diff(input_disp)/dt
    input_acc = np.diff(input_vel)/dt
    tim = tim[:-2]
    ntim = len(tim)

else:
    print("Error: 'input_time_series' is not defined")
    exit()

# --------------------- #
if input_type == "plane_wave":
    polarity = 90    # [deg] N[XX]E
    wave_accx = input_acc * np.cos(np.deg2rad(polarity))
    wave_accy = input_acc * np.sin(np.deg2rad(polarity))

    # output_line = np.vstack([tim,wave_accx,wave_accy]).T
    # np.savetxt(output_dir+"input.acc",output_line)

    # plt.figure()
    # plt.plot(tim,wave_accx)
    # plt.plot(tim,wave_accy)
    # plt.show()
    # exit()

elif input_type == "double_couple":
    strike = 0.0   # degree
    dip = 45.0      # degree
    rake = 90.0     # degree

    sx = 0.0
    sy = 0.0
    sz = 600.0
    mw = 3.0

    rmu = 1000*1000*2100
    m0 = 10**(1.5*mw+9.1)
    length = np.sqrt(m0/rmu)
    width = length

    sources = source.set_source(fem.elements,strike,dip,rake,length,width,sx,sy,sz,1,1)

    input_slipr = np.cumsum(input_acc)*dt
    input_stf = np.cumsum(input_slipr)*dt

    # output_line = np.vstack([tim,input_stf]).T
    # np.savetxt(output_dir+"input_stf.dat",output_line)

    # plt.figure()
    # plt.plot(tim,input_stf)
    # plt.show()
    # exit()

else:
    print("Error: 'input_type' is not defined")
    exit()


## --- Prepare time solver --- ##
# ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

output_velx = np.zeros((ntim,fem.output_nnode))
output_vely = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))

output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispy = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

vel0 = np.array([0.0,0.0,0.0])

# --- PML preparation --- #
pml_elems,pml_config = io_data.input_pmls("input/pml.in")
fem.set_pml(pml_elems,pml_config,dt)


## --- Time Iteration --- #
for it in range(len(tim)):
    if input_type == "plane_wave":
        acc0 = np.array([wave_accx[it],wave_accy[it],0.0])
        vel0 += acc0*dt
        fem.update_time(acc0,vel0,input_wave=True)

    elif input_type == "double_couple":
        fem.update_time_source(sources,input_stf[it])

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]
    output_vely[it,:] = [node.v[1] for node in fem.output_nodes]
    output_velz[it,:] = [node.v[2] for node in fem.output_nodes]
    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispy[it,:] = [node.u[1] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[2] for node in fem.output_nodes]

    if it%20 == 0:
        # plot_model.plot_mesh_update(ax,fem,20.)
        print(it,"t=",it*dt,output_vely[it,1])

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,10.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim,output_velx.T]).T
np.savetxt(output_dir+"output_x.vel",output_line)
output_line = np.vstack([tim,output_vely.T]).T
np.savetxt(output_dir+"output_y.vel",output_line)
output_line = np.vstack([tim,output_velz.T]).T
np.savetxt(output_dir+"output_z.vel",output_line)

output_line = np.vstack([tim,output_dispx.T]).T
np.savetxt(output_dir+"output_x.disp",output_line)
output_line = np.vstack([tim,output_dispy.T]).T
np.savetxt(output_dir+"output_y.disp",output_line)
output_line = np.vstack([tim,output_dispz.T]).T
np.savetxt(output_dir+"output_z.disp",output_line)

## --- Write vtk file --- ##
# vtk.output(fem,output_dir+"output.vtk")

## Output result ##
plt.figure()
plt.plot(tim,output_vely[:,1],c='k')
plt.show()
