import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

# --- Read input wave --- #
fsamp = 3000
duration = 1.0

tim,dt = np.linspace(0,duration,int(fsamp*duration),endpoint=False,retstep=True)

input_tim,input_disp = np.loadtxt("input_wave.txt",skiprows=1,unpack=True)

input_fp = 100
fp = 4

scaling = input_fp/fp
scaled_tim = input_tim*scaling

fd = interpolate.interp1d(scaled_tim,input_disp,kind="cubic")
wave_disp = fd(tim)

wave_vel = np.diff(wave_disp)/dt
wave_acc = np.diff(wave_vel)/dt

plt.figure()
plt.plot(tim,wave_disp)
plt.show()

with open("scaled_input_acc.txt","w") as f:
    f.write("{}\n".format(len(wave_acc)))
    for i in range(len(wave_acc)):
        f.write("{} {}\n".format(tim[i],wave_acc[i]))
