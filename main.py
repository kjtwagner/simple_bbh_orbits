from bbh.eom_func import *
import matplotlib.pyplot as plt

bh1 = blackHole(10, 0)
bh2 = blackHole(1,1000)

m1, m2, ic = orbit_init(bh1,bh2)

time = np.linspace(0,1000000,5000)

x2,y2,z2 = orbit_ode(two_body_eom_pn, ic, time, G, m1, m2)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter3D(0,0,0)
ax.plot3D(x2, y2, z2, 'green',)
plt.show()