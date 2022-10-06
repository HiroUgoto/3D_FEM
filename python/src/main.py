import matplotlib.pyplot as plt
import numpy as np
import time

import io_data
import input_wave
import source
import plot_model
import datetime


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

## --- Define EQ source --- ##
fsamp = 1000
duration = 2

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
ntim = len(tim)


## --- Fault setup --- ##
fem.set_initial_fault()

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

for it in range(len(tim)):
    fem.update_time_dynamic_fault()

    output_dispx[it,:] = [node.u[0] for node in fem.output_nodes]
    output_dispy[it,:] = [node.u[1] for node in fem.output_nodes]
    output_dispz[it,:] = [node.u[2] for node in fem.output_nodes]

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]
    output_vely[it,:] = [node.v[1] for node in fem.output_nodes]
    output_velz[it,:] = [node.v[2] for node in fem.output_nodes]

    output_accx[it,:] = [node.a[0] for node in fem.output_nodes]
    output_accy[it,:] = [node.a[1] for node in fem.output_nodes]
    output_accz[it,:] = [node.a[2] for node in fem.output_nodes]

    if it%20 == 0:
        plot_model.plot_mesh_update(ax,fem,10.)
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

## Output result ##
# plt.figure()
# plt.plot(tim,output_velx[:,0],c='r')
# plt.show()
