import numpy  as np
import matplotlib.pyplot as plt

d_path = '/home/majong/Documents/simfor/build/examples/ODUtest/'

Eul = np.loadtxt ( d_path + 'lorenz_attr_Eul.txt' )
RK4 = np.loadtxt ( d_path + 'lorenz_attr_RK.txt' )
AdBash = np.loadtxt ( d_path + 'lorenz_attr_AdBash.txt')
AdMltn = np.loadtxt ( d_path + 'lorenz_attr_AdMltn.txt')

ax = plt.figure().add_subplot(projection='3d')
ax.plot(*AdBash.T, lw=0.5)
#ax.plot(*Eul.T, lw=0.6)
#ax.plot(*RK4.T, lw=0.6)
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
ax.set_title("аттрактор Лоренца")
plt.show()
