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
# plot_model.plot_mesh(fem,fmt="tetra")
# exit()

###############################
## --- Define plane wave --- ##
###############################
# fsamp = 1200

# fp = 2
# duration = 3.0

# tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# wave_acc = input_wave.ricker(tim,fp,tp=1.0/fp,amp=1.0)
# ntim = len(tim)

# polarity = 45    # [deg] N[XX]E
# wave_accx = wave_acc * np.cos(np.deg2rad(polarity))
# wave_accy = wave_acc * np.sin(np.deg2rad(polarity))

# --- Read input wave --- #
# fsamp = 1000
# duration = 0.5

# tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
# input_tim,input_disp = np.loadtxt("input/input_wave.txt",skiprows=1,unpack=True)

# input_fp = 100
# fp = 10

# scaling = input_fp/fp
# scaled_tim = input_tim*scaling

# fd = interpolate.interp1d(scaled_tim,input_disp,kind="cubic")
# wave_disp = fd(tim)

# wave_vel = np.diff(wave_disp)/dt
# wave_acc = np.diff(wave_vel)/dt

# plt.figure()
# plt.plot(tim[:-1],wave_vel)
# plt.show()
# exit()

# tim = tim[:-2]
# ntim = len(tim)

# polarity = 45    # [deg] N[XX]E
# wave_accx = wave_acc * np.cos(np.deg2rad(polarity))
# wave_accy = wave_acc * np.sin(np.deg2rad(polarity))

# plt.figure()
# plt.plot(tim,wave_accx)
# plt.plot(tim,wave_accy)
# plt.show()
# exit()

###############################
## --- Define plane wave --- ##
###############################
fsamp = 1000

fp = 2.0
duration = 2.5

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
input_stf = input_wave.ricker(tim,fp,tp=1.5/fp,amp=1.0)
ntim = len(tim)

# tim,input_stf = np.loadtxt("input/scaled_input_STF.txt",skiprows=1,unpack=True)
# dt = tim[1] - tim[0]
# ntim = len(tim)

strike = 0.0   # degree
dip = 45.0      # degree
rake = 90.0     # degree

sx = 250.0
sy = 250.0
sz = 250.0
mw = 3.0

rmu = 1000*1000*2100
m0 = 10**(1.5*mw+9.1)
length = np.sqrt(m0/rmu)
width = length

sources = source.set_source(fem.elements,strike,dip,rake,length,width,sx,sy,sz,1,1)

# output_line = np.vstack([tim,input_stf]).T
# np.savetxt(output_dir+"input_stf.dat",output_line)

# plt.figure()
# plt.plot(tim,input_stf)
# plt.show()
# exit()

## --- Prepare time solver --- ##
# ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_velx = np.zeros((ntim,fem.output_nnode))
output_vely = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))

output_dispx = np.zeros((ntim,fem.output_nnode))
output_dispy = np.zeros((ntim,fem.output_nnode))
output_dispz = np.zeros((ntim,fem.output_nnode))

# vel0 = np.array([0.0,0.0,0.0])

for it in range(len(tim)):
    # acc0 = np.array([wave_accx[it],wave_accy[it],0.0])
    # vel0 += acc0*dt
    # fem.update_time(acc0,vel0,input_wave=True)

    fem.update_time_source(sources,input_stf[it])

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]
    output_vely[it,:] = [node.v[1] for node in fem.output_nodes]
    output_velz[it,:] = [node.v[2] for node in fem.output_nodes]
    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispy[it,:] = [node.u[1] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[2] for node in fem.output_nodes]

    if it%20 == 0:
        # plot_model.plot_mesh_update(ax,fem,20.)
        print(it,"t=",it*dt,output_vely[it,0])

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

# output_line = np.vstack([tim,wave_accx,wave_accy]).T
# np.savetxt(output_dir+"input.acc",output_line)

# input_velx = np.cumsum(wave_accx)*dt
# input_vely = np.cumsum(wave_accy)*dt
# output_line = np.vstack([tim,input_velx,input_vely]).T
# np.savetxt(output_dir+"input.vel",output_line)

# input_dispx = np.cumsum(input_velx)*dt
# input_dispy = np.cumsum(input_vely)*dt
# output_line = np.vstack([tim,input_dispx,input_dispy]).T
# np.savetxt(output_dir+"input.disp",output_line)


## --- Write vtk file --- ##
# vtk.output(fem,output_dir+"output.vtk")

## Output result ##
plt.figure()
# plt.plot(tim,input_velx*2,c='gray')
plt.plot(tim,output_velx[:,0],c='k')
plt.show()
