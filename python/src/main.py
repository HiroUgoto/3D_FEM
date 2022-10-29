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
fsamp = 2500
duration = 1.0

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)
ntim = len(tim)


## --- Fault setup --- ##
fem.set_initial_fault()

## --- Prepare time solver --- ##
ax = plot_model.plot_mesh_update_init()
fem.update_init(dt)

## Iteration ##
output_velx = np.zeros((ntim,fem.output_nnode))
output_vely = np.zeros((ntim,fem.output_nnode))
output_velz = np.zeros((ntim,fem.output_nnode))

output_slip = np.zeros((ntim,fem.output_nfault))
output_sliprate = np.zeros((ntim,fem.output_nfault))
output_traction = np.zeros((ntim,fem.output_nfault))

for it in range(len(tim)):
    fem.update_time_dynamic_fault(tim[it])

    output_velx[it,:] = [node.v[0] for node in fem.output_nodes]
    output_vely[it,:] = [node.v[1] for node in fem.output_nodes]
    output_velz[it,:] = [node.v[2] for node in fem.output_nodes]

    output_slip[it,:] = [fault.slip for fault in fem.output_faults]
    output_sliprate[it,:] = [fault.sliprate for fault in fem.output_faults]
    output_traction[it,:] = [fault.traction + fault.p0 for fault in fem.output_faults]

    if it%20 == 0:
        # plot_model.plot_mesh_update(ax,fem,50.)
        print(it,"t=",tim[it],output_sliprate[it,:])
        print("     ",output_slip[it,:])
        print("     ",output_traction[it,:])

elapsed_time = time.time() - start
print ("elapsed_time: {0}".format(elapsed_time) + "[sec]")

# plot_model.plot_mesh_update(ax,fem,10.,fin=True)

## --- Write output file --- ##
output_line = np.vstack([tim,output_velx.T]).T
np.savetxt(output_dir+"output_x.vel",output_line)
output_line = np.vstack([tim,output_vely.T]).T
np.savetxt(output_dir+"output_y.vel",output_line)

output_line = np.vstack([tim,output_slip.T]).T
np.savetxt(output_dir+"output_slip.dat",output_line)
output_line = np.vstack([tim,output_sliprate.T]).T
np.savetxt(output_dir+"output_sliprate.dat",output_line)
output_line = np.vstack([tim,output_traction.T]).T
np.savetxt(output_dir+"output_traction.dat",output_line)

## Output result ##
# plt.figure()
# plt.plot(tim,output_velx[:,0],c='r')
# plt.show()
