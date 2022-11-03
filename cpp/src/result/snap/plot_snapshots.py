import numpy as np
import matplotlib.pyplot as plt

######################################################################
def uniform_mesh_layout(nx,dx):
    area0 = nx*dx
    xg = np.linspace(-area0/2,area0/2,nx+1,endpoint=True)
    return xg

def uniform_mesh_layout_z(nz,dz,zh):
    area0 = nz*dz
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

dy = 100 # [m]
dz = 100 # [m]

y0 = 3000
z0 = 3000
zh = 7500

L = 15     # [km]
W =  6     # [km]

# ny = 400
# nz = 200
ny = 150
nz = 70

snap_name = "008"

yn = variable_mesh_layout(ny,dy,y0) / 1000.0
# yn = np.linspace(-L/2,L/2,ny+1)
zn = variable_mesh_layout_z(nz,dz,z0,zh) / 1000.0
# zn = np.linspace(0,W,nz+1)

y = np.convolve(yn,[0.5,0.5],mode="valid")
z = np.convolve(zn,[0.5,0.5],mode="valid")
yg,zg = np.meshgrid(y,z)

slip = np.empty([nz,ny])
sliprate = np.empty([nz,ny])
trac = np.empty([nz,ny])
rupture_time = np.empty([nz,ny])

with open(snap_name + ".dat", "r") as f:
    lines = f.readlines()
    i = 0
    for k in range(nz):
        for j in range(ny):
            items = lines[i].split()

            slip[k,j] = float(items[2])
            sliprate[k,j] = float(items[3])
            trac[k,j] = float(items[4])

            i += 1

with open("../rupture_time.dat", "r") as f:
    lines = f.readlines()
    i = 0
    for k in range(nz):
        for j in range(ny):
            items = lines[i].split()

            rupture_time[k,j] = float(items[2])
            i += 1



plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["font.size"] = 18
fig,ax = plt.subplots(figsize=(8,8))

ax.set_xticks([-L/2,-L/3,-L/6, 0, L/6, L/3, L/2])
ax.set_yticks([0, W/2, W])
ax.set_xlim([-L/2,L/2])
ax.set_ylim([0,W])
ax.set_xlabel("km")
ax.set_aspect('equal')
ax.invert_yaxis()

c = ax.pcolormesh(yg,zg,slip,vmin=0,cmap='OrRd')
cb = fig.colorbar(c,ax=ax,ticks=[0,2,4,6,8],aspect=20,shrink=0.5,label="slip (m)")
# cont = plt.contour(yg,zg,rupture_time,levels=[1,2,3,4,5,6],colors='k')
cont = plt.contour(yg,zg,rupture_time,levels=[0.5,1,1.5,2,2.5,3,3.5,4],colors='k')

# c = ax.pcolormesh(yg,zg,sliprate,vmin=0,cmap='OrRd')
# cb = fig.colorbar(c,ax=ax,ticks=[0,2,4,6,8],aspect=20,shrink=0.5,label="sliprate (m/s)")

ax.grid()

plt.show()
