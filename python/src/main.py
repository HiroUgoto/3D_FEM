import matplotlib.pyplot as plt
import numpy as np
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
fem.set_init()
fem.set_output(outputs)
# plot_model.plot_mesh(fem)
# exit()

## --- Define plane wave --- ##
# fsamp = 1000

# fp = 4.0
# duration = 1

# tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# wave_acc = input_wave.ricker(tim,fp,tp=1.0/fp,amp=1.0)
# ntim = len(tim)

# --- Read input wave --- #
fsamp = 1000
duration = 1.0

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)

input_tim,input_disp = np.loadtxt("input/input_wave.txt",skiprows=1,unpack=True)

input_fp = 100
fp = 

scaling = input_fp/fp
scaled_tim = input_tim*scaling

fd = interpolate.interp1d(scaled_tim,input_disp,kind="cubic")
wave_disp = fd(tim)

wave_vel = np.diff(wave_disp)/dt
wave_acc = np.diff(wave_vel)/dt

# plt.figure()
# plt.plot(tim,wave_disp)
# plt.show()
# exit()

tim = tim[:-2]
ntim = len(tim)


polarity = 0    # [deg] N[XX]E
wave_accx = wave_acc * np.cos(np.deg2rad(polarity))
wave_accy = wave_acc * np.sin(np.deg2rad(polarity))

# plt.figure()
# plt.plot(tim,wave_accx)
# plt.plot(tim,wave_accy)
# plt.show()
# exit()

## --- Define EQ source --- ##
# fsamp = 100

# fp = 0.5
# duration = 6

# tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# slip_rate = input_wave.ricker(tim,fp,tp=1.0/fp,amp=1.0)
# ntim = len(tim)

# strike = 270.0   # degree
# dip = 30.0      # degree
# rake = 90.0     # degree

# width = 1000.0      # m
# length = 1000.0     # m
# nl,nw = 2,2
# sources = source.set_source(fem.elements,strike,dip,rake,length,width,2500,2500,2500,nl,nw)

# plt.figure()
# plt.plot(tim,slip_rate)
# plt.show()
# exit()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispy = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

output_velx = np.zeros((ntim,fem.output_nnode))
output_vely = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))

output_accx = np.zeros((ntim,fem.output_nnode))
output_accy = np.zeros((ntim,fem.output_nnode))
output_accz = np.zeros((ntim,fem.output_nnode))

vel0 = np.array([0.0,0.0,0.0])
# slip0 = 0.0
for it in range(len(tim)):
    acc0 = np.array([wave_accx[it],wave_accy[it],0.0])
    vel0 += acc0*dt
    fem.update_time(acc0,vel0,input_wave=True)

    # slip0 += slip_rate[it]*dt
    # fem.update_time_source(sources,slip0)

    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispy[it,:] = [node.u[1] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[2] for node in fem.output_nodes]

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]
    output_vely[it,:] = [node.v[1] for node in fem.output_nodes]
    output_velz[it,:] = [node.v[2] for node in fem.output_nodes]

    output_accx[it,:] = [node.a[0] for node in fem.output_nodes]
    output_accy[it,:] = [node.a[1] for node in fem.output_nodes]
    output_accz[it,:] = [node.a[2] for node in fem.output_nodes]

    if it%40 == 0:
        plot_model.plot_mesh_update(ax,fem,100.)
        print(it,"t=",it*dt,output_dispx[it,0])

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,10.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim,output_dispx.T]).T
np.savetxt(output_dir+"output_x.disp",output_line)
output_line = np.vstack([tim,output_dispy.T]).T
np.savetxt(output_dir+"output_y.disp",output_line)

output_line = np.vstack([tim,output_velx.T]).T
np.savetxt(output_dir+"output_x.vel",output_line)
output_line = np.vstack([tim,output_vely.T]).T
np.savetxt(output_dir+"output_y.vel",output_line)

output_line = np.vstack([tim,output_accx.T]).T
np.savetxt(output_dir+"output_x.acc",output_line)
output_line = np.vstack([tim,output_accy.T]).T
np.savetxt(output_dir+"output_y.acc",output_line)

output_line = np.vstack([tim,wave_accx,wave_accy]).T
np.savetxt(output_dir+"input.acc",output_line)

input_velx = np.cumsum(wave_accx)*dt
input_vely = np.cumsum(wave_accy)*dt
output_line = np.vstack([tim,input_velx,input_vely]).T
np.savetxt(output_dir+"input.vel",output_line)

input_dispx = np.cumsum(input_velx)*dt
input_dispy = np.cumsum(input_vely)*dt
output_line = np.vstack([tim,input_dispx,input_dispy]).T
np.savetxt(output_dir+"input.disp",output_line)


## --- Write vtk file --- ##
# vtk.output(fem,output_dir+"output.vtk")

## Output result ##
plt.figure()
plt.plot(tim,input_dispx,c='k')
plt.plot(tim,output_dispx[:,0],c='r')
plt.show()
