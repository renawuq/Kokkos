import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('text', usetex=True)

size=np.array([0,1,2])

v0=np.array([5.61956,5.9004,4.81598])

v1=np.array([6.12319,7.07085,6.48815])

v2=np.array([7.21056,7.54525,7.09387])

plt.plot(size,v0, label='serial')
plt.plot(size,v1, label='kokkod 0')
plt.plot(size,v2, label='kokkod 1')
plt.grid(True)
plt.title('draft')
plt.xlabel('N - linear size')
plt.ylabel(r'Bandwidth (GBytes/s)')
plt.legend()
plt.show()
